"""
meteokav-ecm9/ecm9/build.py - the data behind meteokav.gr's precipitation map (wxprecip.php), from ECMWF's free 9 km
open data (2026-10-21, user: "option 1 - 9 km, 3 h steps, zero Open-Meteo calls, the GitHub way"; the same day: "add
a wind layer in km/h ... animated vectors showing the wind flow, Windy style").

WHY HERE AND NOT ON THE HOST: the 0.1 deg GRIB2 files are CCSDS-packed (template 5.42), which only ECMWF's eccodes
library decodes - a Python job, so it runs as a GitHub Action (.github/workflows/ecm9.yml) and the host only pulls the
finished pictures (wxprecip/pull.php). Open-Meteo is not used at all.

WHAT: for the newest IFS HRES run (00 / 06 / 12 / 18 UTC, all four go to 144 h at 3-hour steps) it downloads ONLY the
messages it needs of each step 0..HORIZON (PARAMS: rain, 10 m wind, 2 m temperature, convective rain, CAPE, snowfall)
through the step's .index file (byte ranges: ~14 MB of a ~90 MB file, ~450 MB a run), crops the Greece box at the full
0.1 deg and writes wxprecip's
"data images":
  out/p<run>.png    rain: NX wide, NY x FRAMES tall; frame k = the rain of the 3 hours ENDING at run + 3k h (frame 0 =
                    0), each frame north row first, west to east; pixel = mm x 100 in 16 bits (R = high, G = low byte)
  out/w<run>.png    wind: the same layout; frame k = the 10 m wind AT run + 3k h (frame 0 = the analysis); R = u,
                    G = v (m/s, east / north), each byte = WIND_ZERO + value / WIND_SCALE (+-44 m/s = +-160 km/h in
                    0.35 m/s steps - plenty for a picture, a quarter of 16 bits' size)
  out/t<run>.png    2 m temperature (2026-10-21, user: "a 2m Temp layer"): the same layout, frame k = AT run + 3k h;
                    16 bits like the rain (R = high, G = low byte), deg C = TEMP_MIN + value x TEMP_SCALE
  out/x<run>.png    the rain's MARKS (2026-10-21, user: "lightning symbols where a thunderstorm is likely, snowflakes
                    where it snows"; the open data carry no lightning density, so ECMWF's own ingredients are used):
                    the same layout, frame k = the 3 hours ending at run + 3k h; R = THUNDER STRENGTH 0..255
                    (marks v2, 2026-10-22, user: "size + opacity by the actual intensity"): 0 = no thunder, else the
                    period's convective rain cp >= THUNDER_MIN[0] mm AND the larger CAPE (mucape) of its two ends >=
                    THUNDER_MIN[1] J/kg, and the strength = cp x sqrt(CAPE) (the ingredients of ECMWF's own lightning
                    scheme) on a log scale from THUNDER_SPAN[0] (= 1) to THUNDER_SPAN[1] (= 255); G = the model's snow
                    share of that rain (snowfall sf / tp, 0..255, only where tp >= SNOW_MIN_MM - its own phase
                    physics on its own 9 km terrain)
  out/m<run>.png    mean sea level pressure (2026-10-22, user: "an MSLP layer with smooth pressure lines"): the same
                    layout, frame k = AT run + 3k h; 16 bits, hPa = PRES_MIN + value x PRES_SCALE
  out/latest.json   the run, the grid, frames + step (= 3: the page shows 3-hour periods only - the open data have
                    no hourly steps, and splitting them into hours would be an estimate, not model output), wind {},
                    temp {}, marks {}

usage:
  python ecm9/build.py --check --published published.json   stdlib only: is a newer complete run out? (GitHub output)
  python ecm9/build.py --run 2026102012 --out out             build that run (needs requirements.txt)
  --force (with --check): the newest complete run even if it is published (a format change: rebuild now)
Data: ECMWF open data, CC BY 4.0 - the page credits ECMWF.
"""
import argparse, calendar, json, os, sys, time, urllib.error, urllib.request

# ----------------------------- CONFIG --------------------------------
MIRRORS = [                                    # ECMWF's portal, then the AWS copy - each file falls back to the next one
    "https://data.ecmwf.int/forecasts/",       # (2026-10-21: the AWS copy of the new 0p1 folder answered 503 Slow Down
    "https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com/",   # to every other request; ECMWF limits connections)
]
# the map's box (the page takes it from latest.json). 2026-10-22 (user: "more of the Balkans and a bit more of Italy -
# organised systems and the winter's cold shots from the north"): south Italy + Sicily to the Black Sea's west coast,
# Crete to Hungary / Romania (was 17.5-31.5 E, 34.6-42.0 N). 2026-10-23 (user: "Greece in the centre of the frame"): the
# widest view shows the whole box, so its centre sits on Greece - south to the African coast, east to Cyprus
BOX = dict(lonW=11.0, lonE=35.0, latS=31.5, latN=47.5)
RES = 0.1                                      # deg: the open data's own grid, nothing resampled
STEP = 3                                       # h between the open data's steps (to 144 h)
HORIZON = 96                                   # h of forecast on the map (4 days, as before)
LOOKBACK = 6                                   # runs looked back for the newest complete one (= 36 h)
TRIES = 4                                      # attempts per file and mirror (waits 3, 6, 9 s between them)
PARAMS = ("tp", "10u", "10v", "2t", "cp", "mucape", "sf", "msl")   # the messages taken from each step
PRES_SCALE = 0.01                              # hPa per step of the 16-bit pressure picture (0.05 left terraces in flat fields: broken isobars)
PRES_MIN = 900.0                               # hPa of value 0
# sea-level pressure over mountains is a reduction, noisy at 9 km (small closed loops over the Dinarides / Balkan Mts):
# smoothed like on every weather chart - a Gaussian of this many cells (~30 km) before the isobars are drawn
PRES_SMOOTH = 3.0
THUNDER_MIN = (1.0, 400.0)                     # thunder at all: convective rain mm / 3 h AND CAPE J/kg both reached
THUNDER_SPAN = (20.0, 600.0)                   # cp x sqrt(CAPE) of strength 1 (= the minimum) .. 255 (e.g. 15 mm, 1600 J/kg)
SNOW_MIN_MM = 0.3                              # mm / 3 h: below it no snow share (a trace says nothing)
WIND_SCALE = 0.35                              # m/s per byte step of the wind picture
WIND_ZERO = 128                                # the byte of 0 m/s
TEMP_SCALE = 0.1                               # deg C per step of the 16-bit temperature picture (finer = a bigger PNG)
TEMP_MIN = -60.0                               # deg C of value 0
UA = "MeteokavBot/1.0 (+https://www.meteokav.gr)"
FILL = 65535
# ---------------------------------------------------------------------


def http(url, rng=None, head=False):
    """one request with retries; None on 404 (not out yet)"""
    hdr = {"User-Agent": UA}
    if rng: hdr["Range"] = "bytes=%d-%d" % rng
    for t in range(TRIES):
        try:
            req = urllib.request.Request(url, headers=hdr, method="HEAD" if head else "GET")
            with urllib.request.urlopen(req, timeout=120) as r:
                return b"" if head else r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404): return None   # S3 answers 403 for a missing key
            if t == TRIES - 1: raise
        except Exception:
            if t == TRIES - 1: raise
        time.sleep(3 * (t + 1))


def stem(run, step):
    """the folder + file name of a run's surface file at `step` (run = UTC epoch seconds)"""
    d, h = time.strftime("%Y%m%d", time.gmtime(run)), time.gmtime(run).tm_hour
    return "%s/%02dz/ifs/0p1/oper/%s%02d0000-%dh-oper-fc-sfc" % (d, h, d, h, step)


def fetch(path, rng=None, head=False):
    """one file from the first mirror that has it (a mirror that keeps failing hands over to the next); None = nowhere"""
    err = None
    for m in MIRRORS:
        try:
            r = http(m + path, rng, head)
            if r is not None: return r
        except Exception as e:
            err = e
    if err: raise err
    return None


def complete(run):
    """True when the run's LAST step is out (files arrive in step order)"""
    return fetch(stem(run, HORIZON) + ".index", head=True) is not None


def newest(published):
    """the newest complete run newer than the published one, or None"""
    now = int(time.time()) // 21600 * 21600
    for k in range(LOOKBACK):
        run = now - k * 21600
        if run <= published: break
        if complete(run): return run
    return None


def messages(run, step):
    """{param: GRIB message bytes} of one step for PARAMS, via the .index byte ranges"""
    idx = fetch(stem(run, step) + ".index")
    if idx is None: raise RuntimeError("no index for +%dh" % step)
    lines = {}
    for l in idx.decode().splitlines():
        j = json.loads(l)
        if j.get("param") in PARAMS: lines[j["param"]] = j
    out = {}
    for p in PARAMS:
        if p not in lines: raise RuntimeError("no %s in the index of +%dh" % (p, step))
        o, n = lines[p]["_offset"], lines[p]["_length"]
        msg = fetch(stem(run, step) + ".grib2", (o, o + n - 1))
        if not msg or len(msg) != n: raise RuntimeError("short %s message at +%dh" % (p, step))
        out[p] = msg
    return out


def build(run, out):
    import numpy as np, eccodes
    from PIL import Image

    cut = {}                                   # the box's rows + columns in the global grid (found once)

    def decode(msg):
        h = eccodes.codes_new_from_message(msg)
        try:
            ni, nj = eccodes.codes_get(h, "Ni"), eccodes.codes_get(h, "Nj")
            la1 = eccodes.codes_get(h, "latitudeOfFirstGridPointInDegrees")
            lo1 = eccodes.codes_get(h, "longitudeOfFirstGridPointInDegrees")
            d = eccodes.codes_get(h, "iDirectionIncrementInDegrees")
            v = eccodes.codes_get_values(h).reshape(nj, ni)
        finally:
            eccodes.codes_release(h)
        if not cut:                            # rows north -> south, columns west -> east
            la = np.round(np.arange(BOX["latN"], BOX["latS"] - 1e-6, -RES), 4)
            lo = np.round(np.arange(BOX["lonW"], BOX["lonE"] + 1e-6, RES), 4)
            cut["r"], cut["c"] = np.rint((la1 - la) / d).astype(int), np.rint(((lo - lo1) % 360) / d).astype(int)
        return v[np.ix_(cut["r"], cut["c"])]

    prev, rain, us, vs, ts, ps, got = None, [], [], [], [], [], 0
    pcp = psf = pcape = None; thunder, snow = [], []
    for step in range(0, HORIZON + 1, STEP):
        m = messages(run, step); got += sum(len(x) for x in m.values())
        cur = decode(m["tp"])                  # m accumulated since the run's start
        tp3 = np.zeros_like(cur) if prev is None else np.clip(cur - prev, 0, None) * 1000.0
        rain.append(tp3)
        prev = cur
        us.append(decode(m["10u"])); vs.append(decode(m["10v"])); ts.append(decode(m["2t"]) - 273.15)
        ps.append(decode(m["msl"]) / 100.0)    # Pa -> hPa
        # the marks: this period's convective rain + snowfall (accumulated, m) and the larger CAPE of its two ends
        ccp, csf, cape = decode(m["cp"]), decode(m["sf"]), np.nan_to_num(decode(m["mucape"]))
        lvl = np.zeros(cur.shape, np.uint8); share = np.zeros(cur.shape, np.uint8)
        if pcp is not None:
            cp3, sf3, cmax = np.clip(ccp - pcp, 0, None) * 1000.0, np.clip(csf - psf, 0, None) * 1000.0, np.maximum(cape, pcape)
            on = (cp3 >= THUNDER_MIN[0]) & (cmax >= THUNDER_MIN[1])
            idx = cp3[on] * np.sqrt(cmax[on])
            f = np.log(np.maximum(idx, THUNDER_SPAN[0]) / THUNDER_SPAN[0]) / np.log(THUNDER_SPAN[1] / THUNDER_SPAN[0])
            lvl[on] = (1 + np.rint(np.clip(f, 0, 1) * 254)).astype(np.uint8)
            wet = tp3 >= SNOW_MIN_MM
            share[wet] = np.rint(np.clip(sf3[wet] / tp3[wet], 0, 1) * 255).astype(np.uint8)
        thunder.append(lvl); snow.append(share)
        pcp, psf, pcape = ccp, csf, cape

    ny, nx = rain[0].shape
    os.makedirs(out, exist_ok=True)

    def save16(frames, name):                  # 16 bits per cell: R = high byte, G = low byte
        q = np.clip(np.rint(np.stack(frames)), 0, FILL - 1).astype(np.uint16).reshape(len(frames) * ny, nx)
        img = np.zeros((len(frames) * ny, nx, 3), np.uint8)
        img[..., 0], img[..., 1] = q >> 8, q & 255
        Image.fromarray(img, "RGB").save(os.path.join(out, name), optimize=True)

    name = "p%d.png" % run
    save16([r * 100 for r in rain], name)
    tname = "t%d.png" % run
    save16([(t - TEMP_MIN) / TEMP_SCALE for t in ts], tname)
    def smooth(a, sigma):                      # separable Gaussian, edges padded (no numpy-free scipy needed)
        r = int(3 * sigma); k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2); k /= k.sum()
        p = np.pad(a, r, mode="edge")
        p = np.apply_along_axis(lambda v: np.convolve(v, k, "valid"), 1, p)
        return np.apply_along_axis(lambda v: np.convolve(v, k, "valid"), 0, p)

    ps = [smooth(p, PRES_SMOOTH) for p in ps]
    mname = "m%d.png" % run
    save16([(p - PRES_MIN) / PRES_SCALE for p in ps], mname)

    wb = lambda a: np.clip(np.rint(np.stack(a) / WIND_SCALE + WIND_ZERO), 0, 255).astype(np.uint8).reshape(len(a) * ny, nx)
    img = np.zeros((len(us) * ny, nx, 3), np.uint8)
    img[..., 0], img[..., 1] = wb(us), wb(vs)
    wname = "w%d.png" % run
    Image.fromarray(img, "RGB").save(os.path.join(out, wname), optimize=True)

    img = np.zeros((len(thunder) * ny, nx, 3), np.uint8)
    img[..., 0], img[..., 1] = np.stack(thunder).reshape(-1, nx), np.stack(snow).reshape(-1, nx)
    xname = "x%d.png" % run
    Image.fromarray(img, "RGB").save(os.path.join(out, xname), optimize=True)

    meta = dict(model="ecmwf_ifs_open_0p1", run=run, generated=int(time.time()), file=name, nx=nx, ny=ny,
                frames=len(rain), step=STEP, lonW=BOX["lonW"], latN=BOX["latN"], dx=RES, dy=RES, scale=0.01, fill=FILL,
                wind=dict(file=wname, scale=WIND_SCALE, zero=WIND_ZERO),
                temp=dict(file=tname, scale=TEMP_SCALE, offset=TEMP_MIN), marks=dict(file=xname, v=2),
                pres=dict(file=mname, scale=PRES_SCALE, offset=PRES_MIN),
                source="ECMWF open data (CC BY 4.0)")
    with open(os.path.join(out, "latest.json"), "w") as f: json.dump(meta, f)
    spd = np.hypot(np.stack(us), np.stack(vs)) * 3.6
    kb = lambda n: os.path.getsize(os.path.join(out, n)) / 1024
    th, sn = np.stack(thunder), np.stack(snow)
    print("run %s: %d frames of %d x %d, %.1f MB fetched; rain png %.1f KB (max 3 h %.1f mm), wind png %.1f KB (max %.0f km/h), "
          "temp png %.1f KB (%.1f..%.1f C), marks png %.1f KB (thunder cells %d, strength p50 %d / p90 %d / max %d, snow >= 50 %% cells %d)" % (
          time.strftime("%Y-%m-%d %HZ", time.gmtime(run)), len(rain), nx, ny, got / 1e6, kb(name), float(np.max(rain)),
          kb(wname), float(spd.max()), kb(tname), float(np.min(ts)), float(np.max(ts)), kb(xname),
          int((th > 0).sum()), *(int(np.percentile(th[th > 0], q)) if (th > 0).any() else 0 for q in (50, 90, 100)), int((sn >= 128).sum())))
    print("pressure png %.1f KB (%.1f..%.1f hPa)" % (kb(mname), float(np.min(ps)), float(np.max(ps))))


def output(**kv):
    """step outputs for the workflow (GITHUB_OUTPUT), printed too"""
    p = os.environ.get("GITHUB_OUTPUT")
    for k, v in kv.items():
        print("%s=%s" % (k, v))
        if p:
            with open(p, "a") as f: f.write("%s=%s\n" % (k, v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--published", default="")
    ap.add_argument("--run", default="")
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    pub = 0
    if a.published and os.path.exists(a.published) and not a.force:
        try: pub = int(json.load(open(a.published)).get("run", 0))
        except Exception: pub = 0
    if a.run:
        run = calendar.timegm(time.strptime(a.run, "%Y%m%d%H"))
        if not complete(run): sys.exit("run %s is not complete on any mirror" % a.run)
    else:
        run = newest(pub)
        if not run:
            print("up to date: %s published" % (time.strftime("%Y-%m-%d %HZ", time.gmtime(pub)) if pub else "nothing"))
            output(new="false"); return
    if a.check:
        output(new="true", run=time.strftime("%Y%m%d%H", time.gmtime(run))); return
    build(run, a.out)
    output(new="true", run=time.strftime("%Y%m%d%H", time.gmtime(run)))


if __name__ == "__main__":
    main()
