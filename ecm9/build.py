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
  out/a<run>.png    850 hPa temperature (2026-10-26): the same layout, frame k = AT its hour; 16 bits, deg C =
                    T850_MIN + value x T850_SCALE
  out/z<run>.png    500 hPa geopotential height (2026-10-26): the same layout, frame k = AT its hour; 16 bits, m =
                    Z500_MIN + value x Z500_SCALE, lightly smoothed (Z500_SMOOTH)
  out/ea|ez|em<run>.png  EUROPE (2026-10-26, see EU): 850 hPa temperature, 500 hPa height, sea-level pressure on the
                    EU_RES grid, same frames and encodings as a / z / m; latest.json "eu" {nx, ny, lonW, latN, dx, dy, ...}
  out/g<run>.png    10 m gusts (2026-10-26): the same layout, frame k = the strongest gust of the window ENDING at its
                    hour (latest.json gust.win[k] h: 1 to +90 h, 3 to +144 h, 6 past - see GUST); 16 bits, m/s = value x
                    GUST_SCALE; frame 0 = the analysis' (not shown)
  out/latest.json   the run, the grid, frames + step (= 3: the page shows 3-hour periods only - the open data have
                    no hourly steps, and splitting them into hours would be an estimate, not model output), wind {},
                    temp {}, marks {}, pres {}, t850 {}, z500 {}; since 2026-10-26 "hours" (each frame's hour from the
                    run: 0, 3 .. 144, then 150 .. 240 every 6 h - a rain frame past 144 h = the 6 hours ending there)
                    and for a 06Z / 18Z run "tail" {run, from}: from frame "from" on the frames come from that earlier
                    00Z / 12Z run (06Z / 18Z end at 144 h; the map says so). Steps are downloaded WORKERS at a time.

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
HORIZON = 144                                  # h: the last 3-hour step (6 days; was 96 until 2026-10-23); 06Z / 18Z end here
# 2026-10-26 (user: "go ahead with the 240 h extension; when the data come from a 06 / 18Z run, label that 150-240 come
# from the previous run"): the 00Z + 12Z runs go on at 6-hour steps (the open data reach 360 h; the map stops at 10 days -
# a single run's rain areas past ~day 7-10 carry little skill). A 06Z / 18Z run takes its tail from the 00Z / 12Z run
# before it (TAIL_BACK earlier): the frames past HORIZON then come from that run, and latest.json says so ("tail").
LONG_RUNS = (0, 12)                            # UTC hours of the runs that go past HORIZON
LONG_STEP = 6                                  # h between the steps past HORIZON
LONG_HORIZON = 240                             # h: the map's last step (10 days)
TAIL_BACK = 6                                  # h: a 06Z / 18Z run's tail = the run this much earlier
LOOKBACK = 6                                   # runs looked back for the newest complete one (= 36 h)
TRIES = 4                                      # attempts per file and mirror (waits 3, 6, 9 s between them)
WORKERS = 4                                    # steps fetched at once (the downloads were most of the ~10 min build)
PARAMS = ("tp", "10u", "10v", "2t", "cp", "mucape", "sf", "msl")   # the surface messages taken from each step
# the upper air (2026-10-26, user: "850 hPa temperature with lines of equal temperature" + "a general 500 hPa geopotential
# map with pressure lines"): the same step's pressure-level file (-pl); key = param + level
PL_PARAMS = (("t", "850"), ("gh", "500"))
T850_SCALE = 0.02                              # deg C per step of the 16-bit 850 hPa temperature (lines every 2 deg C: 1/100 of a
                                               # spacing, as the isobars; 0.01 after the smoothing = a 1.8 MB picture)
T850_MIN = -60.0                               # deg C of value 0
T850_SMOOTH = 1.5                              # cells: a Gaussian on it (at 850 hPa over the mountains the level meets the
                                               # ground - raw isotherms zig-zag there; every 850 hPa chart is smoothed)
Z500_SCALE = 0.5                               # m per step of the 16-bit 500 hPa height (contours every 40 m; 0.05 = a 1.8 MB picture)
Z500_MIN = 4500.0                              # m of value 0
Z500_SMOOTH = 1.0                              # cells: a light Gaussian on the height (contours without 9 km wiggles)
# 10 m GUSTS (2026-10-26, user: "a gust layer as a sub option of the wind layer"): the open data's gust is a period maximum
# whose period changes along the run - "10fg" to +90 h = the strongest gust of the LAST HOUR (IFS post-processes hourly
# there; the open data keep every third), "10fg3" +93..+144 h = of the last 3 hours, "10fg" past +144 h = of the last 6
# hours. Each frame's window (h) goes into latest.json "gust" "win", the page labels it so
GUST = ("10fg3", "10fg")                       # the step's gust message: the 3-hour one where it exists
# EUROPE (2026-10-26, user: "a Europe-scale map for 850 hPa and 500 hPa geopotential / pressure - 0.25 deg, no particles"):
# the same global messages, a second, larger cut - Iceland / Iberia to the edge of the Urals, North Africa to Scandinavia -
# at EU_RES: smoothed on the 0.1 grid first (each field's own Gaussian, at least EU_AA cells: no aliasing), then sampled
# bilinearly. -> ea / ez / em<run>.png (850 hPa temperature, 500 hPa height, sea-level pressure) + latest.json "eu"
EU = dict(lonW=-25.0, lonE=45.0, latS=30.0, latN=72.0)
EU_RES = 0.25                                  # deg
EU_AA = 1.5                                    # cells of 0.1 deg: the least smoothing before the 0.25 deg sampling
GUST_SCALE = 0.1                               # m/s per step of the 16-bit gust picture
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


def stem(run, step, kind="sfc"):
    """the folder + file name of a run's surface ("sfc") or pressure-level ("pl") file at `step` (run = UTC epoch s)"""
    d, h = time.strftime("%Y%m%d", time.gmtime(run)), time.gmtime(run).tm_hour
    return "%s/%02dz/ifs/0p1/oper/%s%02d0000-%dh-oper-fc-%s" % (d, h, d, h, step, kind)


def long_run(run):
    return time.gmtime(run).tm_hour in LONG_RUNS


def plan(run):
    """the steps to read, in order: (source run, its step, the frame's hour from `run`, or None = only the base the next
    period's totals start from). 00Z / 12Z: 0..144 every 3 h, then 150..240 every 6 h. 06Z / 18Z: 0..144, then the
    previous long run - its step at the same instant as run + 144 h as the base, then run + 150 .. + 234 h"""
    seq = [(run, s, s) for s in range(0, HORIZON + 1, STEP)]
    if long_run(run):
        seq += [(run, s, s) for s in range(HORIZON + LONG_STEP, LONG_HORIZON + 1, LONG_STEP)]
    else:
        prev = run - TAIL_BACK * 3600
        seq.append((prev, HORIZON + TAIL_BACK, None))
        seq += [(prev, s + TAIL_BACK, s) for s in range(HORIZON + LONG_STEP, LONG_HORIZON - TAIL_BACK + 1, LONG_STEP)]
    return seq


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
    """True when the run's LAST step is out (files arrive in step order; a long run's 6-hour steps come in the same go)"""
    return fetch(stem(run, LONG_HORIZON if long_run(run) else HORIZON) + ".index", head=True) is not None


def newest(published):
    """the newest complete run newer than the published one, or None"""
    now = int(time.time()) // 21600 * 21600
    for k in range(LOOKBACK):
        run = now - k * 21600
        if run <= published: break
        if complete(run): return run
    return None


def pick(run, step, kind, want):
    """{key: GRIB message bytes} from one file via its .index byte ranges; want = {key: (param, level or None)} or
    {key: [(param, level), ...]} = alternatives, the first one the index has wins (its param: out[key + ":param"])"""
    idx = fetch(stem(run, step, kind) + ".index")
    if idx is None: raise RuntimeError("no %s index for +%dh" % (kind, step))
    alts = {key: (w if isinstance(w, list) else [w]) for key, w in want.items()}
    lines = {}
    for l in idx.decode().splitlines():
        j = json.loads(l)
        for key, opts in alts.items():
            for rank, (p, lv) in enumerate(opts):
                if j.get("param") == p and (lv is None or j.get("levelist") == lv) and rank < lines.get(key, (99, None))[0]:
                    lines[key] = (rank, j)
    lines = {key: j for key, (rank, j) in lines.items()}
    out = {}
    for key in want:
        if key not in lines: raise RuntimeError("no %s in the %s index of +%dh" % (key, kind, step))
        o, n = lines[key]["_offset"], lines[key]["_length"]
        msg = fetch(stem(run, step, kind) + ".grib2", (o, o + n - 1))
        if not msg or len(msg) != n: raise RuntimeError("short %s message at +%dh" % (key, step))
        out[key] = msg
        out[key + ":param"] = lines[key]["param"]
    return out


def messages(run, step):
    """{key: GRIB message bytes} of one step: PARAMS from the surface file, PL_PARAMS ("t850", "gh500") from the -pl one"""
    want = {p: (p, None) for p in PARAMS}
    want["gust"] = [(g, None) for g in GUST]
    out = pick(run, step, "sfc", want)
    out.update(pick(run, step, "pl", {p + lv: (p, lv) for p, lv in PL_PARAMS}))
    return out


def fetched(seq):
    """messages() of each (run, step, hour) in order, WORKERS steps downloading at once, at most 2 x WORKERS ahead (memory)"""
    from collections import deque
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(WORKERS) as ex:
        ahead, it = deque(), iter(seq)
        for item in it:
            ahead.append((item, ex.submit(messages, item[0], item[1])))
            if len(ahead) >= 2 * WORKERS: break
        while ahead:
            item, fu = ahead.popleft()
            nxt = next(it, None)
            if nxt is not None: ahead.append((nxt, ex.submit(messages, nxt[0], nxt[1])))
            yield item, fu.result()


def build(run, out):
    import numpy as np, eccodes
    from PIL import Image

    cut = {}                                   # the box's rows + columns in the global grid (found once)

    def smooth(a, sigma):                      # separable Gaussian, edges padded (no numpy-free scipy needed)
        r = int(3 * sigma); k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2); k /= k.sum()
        p = np.pad(a, r, mode="edge")
        p = np.apply_along_axis(lambda v: np.convolve(v, k, "valid"), 1, p)
        return np.apply_along_axis(lambda v: np.convolve(v, k, "valid"), 0, p)

    elat = np.round(np.arange(EU["latN"], EU["latS"] - 1e-6, -EU_RES), 4)
    elon = np.round(np.arange(EU["lonW"], EU["lonE"] + 1e-6, EU_RES), 4)
    def europe(fine, sigma):                   # the 0.1 deg Europe cut -> smoothed -> sampled at EU_RES (bilinear)
        a = smooth(fine, max(sigma, EU_AA))
        fr, fc = (EU["latN"] - elat) / RES, (elon - EU["lonW"]) / RES
        r0 = np.clip(np.floor(fr).astype(int), 0, a.shape[0] - 2); c0 = np.clip(np.floor(fc).astype(int), 0, a.shape[1] - 2)
        wr, wc = (fr - r0)[:, None], (fc - c0)[None, :]
        R, C = np.ix_(r0, c0)
        return (a[R, C] * (1 - wr) * (1 - wc) + a[R, C + 1] * (1 - wr) * wc + a[R + 1, C] * wr * (1 - wc) + a[R + 1, C + 1] * wr * wc)

    def decode(msg, eu=False):
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
            la = np.round(np.arange(EU["latN"], EU["latS"] - 1e-6, -RES), 4)
            lo = np.round(np.arange(EU["lonW"], EU["lonE"] + 1e-6, RES), 4)
            cut["er"], cut["ec"] = np.rint((la1 - la) / d).astype(int), np.rint(((lo - lo1) % 360) / d).astype(int)
        box = v[np.ix_(cut["r"], cut["c"])]
        return (box, v[np.ix_(cut["er"], cut["ec"])]) if eu else box

    prev, rain, us, vs, ts, ps, a8, z5, gs, gw, got = None, [], [], [], [], [], [], [], [], [], 0
    ep, ea, ez = [], [], []                    # Europe: pressure, 850 hPa temperature, 500 hPa height (at EU_RES)
    pcp = psf = pcape = None; thunder, snow = [], []
    hours, tail_from, last_h = [], None, 0
    for (srun, step, hour), m in fetched(plan(run)):
        got += sum(len(x) for x in m.values() if isinstance(x, bytes))
        cur = decode(m["tp"])                  # m accumulated since the SOURCE run's start
        ccp, csf, cape = decode(m["cp"]), decode(m["sf"]), np.nan_to_num(decode(m["mucape"]))
        if hour is None:                       # the tail's base: only where its first period's totals start from
            prev, pcp, psf, pcape = cur, ccp, csf, cape
            continue
        if srun != run and tail_from is None: tail_from = len(hours)
        per = hour - last_h if hours else STEP # this frame's period, h (3, past HORIZON 6)
        tp = np.zeros_like(cur) if prev is None else np.clip(cur - prev, 0, None) * 1000.0
        rain.append(tp)
        prev = cur
        us.append(decode(m["10u"])); vs.append(decode(m["10v"])); ts.append(decode(m["2t"]) - 273.15)
        b, e = decode(m["msl"], True); ps.append(b / 100.0); ep.append(europe(e / 100.0, PRES_SMOOTH))   # Pa -> hPa
        b, e = decode(m["t850"], True); a8.append(b - 273.15); ea.append(europe(e - 273.15, T850_SMOOTH))
        b, e = decode(m["gh500"], True); z5.append(b); ez.append(europe(e, Z500_SMOOTH))
        gs.append(np.nan_to_num(decode(m["gust"])))
        gw.append(3 if m["gust:param"] == "10fg3" else LONG_STEP if step > HORIZON else 1)   # its window, h (see GUST)
        # the marks: this period's convective rain + snowfall (accumulated, m) and the larger CAPE of its two ends; the
        # thresholds are per 3 hours, so a 6-hour period's totals count at their 3-hour rate (x STEP / per)
        lvl = np.zeros(cur.shape, np.uint8); share = np.zeros(cur.shape, np.uint8)
        if pcp is not None:
            rate = STEP / float(per)
            cp3, sf3, cmax = np.clip(ccp - pcp, 0, None) * 1000.0 * rate, np.clip(csf - psf, 0, None) * 1000.0, np.maximum(cape, pcape)
            on = (cp3 >= THUNDER_MIN[0]) & (cmax >= THUNDER_MIN[1])
            idx = cp3[on] * np.sqrt(cmax[on])
            f = np.log(np.maximum(idx, THUNDER_SPAN[0]) / THUNDER_SPAN[0]) / np.log(THUNDER_SPAN[1] / THUNDER_SPAN[0])
            lvl[on] = (1 + np.rint(np.clip(f, 0, 1) * 254)).astype(np.uint8)
            wet = tp * rate >= SNOW_MIN_MM
            share[wet] = np.rint(np.clip(sf3[wet] / tp[wet], 0, 1) * 255).astype(np.uint8)
        thunder.append(lvl); snow.append(share)
        pcp, psf, pcape = ccp, csf, cape
        hours.append(hour); last_h = hour

    ny, nx = rain[0].shape
    os.makedirs(out, exist_ok=True)

    def save16(frames, name):                  # 16 bits per cell: R = high byte, G = low byte (any grid: its own size)
        fy, fx = frames[0].shape
        q = np.clip(np.rint(np.stack(frames)), 0, FILL - 1).astype(np.uint16).reshape(len(frames) * fy, fx)
        img = np.zeros((len(frames) * fy, fx, 3), np.uint8)
        img[..., 0], img[..., 1] = q >> 8, q & 255
        Image.fromarray(img, "RGB").save(os.path.join(out, name), optimize=True)

    name = "p%d.png" % run
    save16([r * 100 for r in rain], name)
    tname = "t%d.png" % run
    save16([(t - TEMP_MIN) / TEMP_SCALE for t in ts], tname)
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

    aname = "a%d.png" % run                    # 850 hPa temperature
    a8 = [smooth(t, T850_SMOOTH) for t in a8]
    save16([(t - T850_MIN) / T850_SCALE for t in a8], aname)
    gname = "g%d.png" % run                    # 10 m gusts (m/s, each frame's own window)
    save16([g / GUST_SCALE for g in gs], gname)
    zname = "z%d.png" % run                    # 500 hPa height, lightly smoothed
    z5 = [smooth(z, Z500_SMOOTH) for z in z5]
    save16([(z - Z500_MIN) / Z500_SCALE for z in z5], zname)

    # "hours" = each frame's hour from the run (3-hourly to HORIZON, then 6-hourly); "tail" (06Z / 18Z) = from frame
    # "from" on, the frames come from the run "run" (TAIL_BACK earlier). "step" stays the base step for older pages
    meta = dict(model="ecmwf_ifs_open_0p1", run=run, generated=int(time.time()), file=name, nx=nx, ny=ny,
                frames=len(rain), step=STEP, hours=hours, lonW=BOX["lonW"], latN=BOX["latN"], dx=RES, dy=RES, scale=0.01, fill=FILL,
                wind=dict(file=wname, scale=WIND_SCALE, zero=WIND_ZERO),
                temp=dict(file=tname, scale=TEMP_SCALE, offset=TEMP_MIN), marks=dict(file=xname, v=2),
                pres=dict(file=mname, scale=PRES_SCALE, offset=PRES_MIN),
                t850=dict(file=aname, scale=T850_SCALE, offset=T850_MIN),
                z500=dict(file=zname, scale=Z500_SCALE, offset=Z500_MIN),
                gust=dict(file=gname, scale=GUST_SCALE, offset=0, win=gw),
                source="ECMWF open data (CC BY 4.0)")
    if tail_from is not None: meta["tail"] = {"run": run - TAIL_BACK * 3600, "from": tail_from}
    # Europe: its own grid, the same frames / hours
    eaname, ezname, emname = "ea%d.png" % run, "ez%d.png" % run, "em%d.png" % run
    save16([(t - T850_MIN) / T850_SCALE for t in ea], eaname)
    save16([(z - Z500_MIN) / Z500_SCALE for z in ez], ezname)
    save16([(p - PRES_MIN) / PRES_SCALE for p in ep], emname)
    meta["eu"] = dict(nx=len(elon), ny=len(elat), lonW=EU["lonW"], latN=EU["latN"], dx=EU_RES, dy=EU_RES,
                      t850=dict(file=eaname, scale=T850_SCALE, offset=T850_MIN),
                      z500=dict(file=ezname, scale=Z500_SCALE, offset=Z500_MIN),
                      pres=dict(file=emname, scale=PRES_SCALE, offset=PRES_MIN))
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
    print("Europe %d x %d: 850 hPa png %.1f KB, 500 hPa png %.1f KB, pressure png %.1f KB (%.1f..%.1f hPa)" % (
          len(elon), len(elat), kb(eaname), kb(ezname), kb(emname), float(np.min(ep)), float(np.max(ep))))
    print("gust png %.1f KB (max %.0f km/h; windows %s)" % (kb(gname), float(np.max(gs)) * 3.6, ",".join(str(w) for w in gw)))
    print("850 hPa temp png %.1f KB (%.1f..%.1f C), 500 hPa height png %.1f KB (%.0f..%.0f m); hours %s; tail %s" % (
          kb(aname), float(np.min(a8)), float(np.max(a8)), kb(zname), float(np.min(z5)), float(np.max(z5)),
          ",".join(str(h) for h in hours), meta.get("tail")))


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
