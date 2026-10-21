"""
meteokav-ecm9/ecm9/build.py - the data behind meteokav.gr's precipitation map (wxprecip.php), from ECMWF's free 9 km
open data (2026-10-21, user: "option 1 - 9 km, 3 h steps, zero Open-Meteo calls, the GitHub way").

WHY HERE AND NOT ON THE HOST: the 0.1 deg GRIB2 files are CCSDS-packed (template 5.42), which only ECMWF's eccodes
library decodes - a Python job, so it runs as a GitHub Action (.github/workflows/ecm9.yml) and the host only pulls the
finished picture (wxprecip/pull.php). Open-Meteo is not used at all.

WHAT: for the newest IFS HRES run (00 / 06 / 12 / 18 UTC, all four go to 144 h at 3-hour steps) it downloads ONLY the
total-precipitation message of each step 0..HORIZON (the step's .index file gives its byte range: ~1.5-2.5 MB of a
~90 MB file, ~65 MB a run), turns the accumulated metres into 3-hour totals in mm, crops the Greece box at the full
0.1 deg and writes wxprecip's "data image":
  out/p<run>.png    NX wide, NY x FRAMES tall; frame k = the rain of the 3 hours ENDING at run + 3k h (frame 0 = 0),
                    each frame north row first, west to east; pixel = mm x 100 in 16 bits (R = high byte, G = low byte)
  out/latest.json   the run, the grid, frames + step (= 3: the page shows 3-hour periods only - the open data have
                    no hourly steps, and splitting them into hours would be an estimate, not model output)

usage:
  python ecm9/build.py --check --published published.json   stdlib only: is a newer complete run out? (GitHub output)
  python ecm9/build.py --published published.json --out out  the build (needs requirements.txt)
  python ecm9/build.py --run 2026102012 --out out             a given run (a local test)
Data: ECMWF open data, CC BY 4.0 - the page credits ECMWF.
"""
import argparse, calendar, json, os, sys, time, urllib.error, urllib.request

# ----------------------------- CONFIG --------------------------------
MIRRORS = [                                    # ECMWF's portal, then the AWS copy - each file falls back to the next one
    "https://data.ecmwf.int/forecasts/",       # (2026-10-21: the AWS copy of the new 0p1 folder answered 503 Slow Down
    "https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com/",   # to every other request; ECMWF limits connections)
]
BOX = dict(lonW=17.5, lonE=31.5, latS=34.6, latN=42.0)   # = wxprecip's grid box (the page's MAX_BOUNDS)
RES = 0.1                                      # deg: the open data's own grid, nothing resampled
STEP = 3                                       # h between the open data's steps (to 144 h)
HORIZON = 96                                   # h of forecast on the map (4 days, as before)
LOOKBACK = 6                                   # runs looked back for the newest complete one (= 36 h)
TRIES = 4                                      # attempts per file and mirror (waits 3, 6, 9 s between them)
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


def tp_message(run, step):
    """the total-precipitation GRIB message of one step (bytes), via the .index byte range"""
    idx = fetch(stem(run, step) + ".index")
    if idx is None: raise RuntimeError("no index for +%dh" % step)
    line = next((json.loads(l) for l in idx.decode().splitlines() if '"param": "tp"' in l), None)
    if not line: raise RuntimeError("no tp in the index of +%dh" % step)
    o, n = line["_offset"], line["_length"]
    msg = fetch(stem(run, step) + ".grib2", (o, o + n - 1))
    if not msg or len(msg) != n: raise RuntimeError("short tp message at +%dh" % step)
    return msg


def build(run, out):
    import numpy as np, eccodes
    from PIL import Image

    lats = lons = None
    prev, frames, got = None, [], 0
    for step in range(0, HORIZON + 1, STEP):
        msg = tp_message(run, step); got += len(msg)
        h = eccodes.codes_new_from_message(msg)
        try:
            ni, nj = eccodes.codes_get(h, "Ni"), eccodes.codes_get(h, "Nj")
            la1 = eccodes.codes_get(h, "latitudeOfFirstGridPointInDegrees")
            lo1 = eccodes.codes_get(h, "longitudeOfFirstGridPointInDegrees")
            d = eccodes.codes_get(h, "iDirectionIncrementInDegrees")
            v = eccodes.codes_get_values(h).reshape(nj, ni)
        finally:
            eccodes.codes_release(h)
        if lats is None:                       # the box's rows (north -> south) + columns (west -> east) in this grid
            la = np.round(np.arange(BOX["latN"], BOX["latS"] - 1e-6, -RES), 4)
            lo = np.round(np.arange(BOX["lonW"], BOX["lonE"] + 1e-6, RES), 4)
            lats, lons = np.rint((la1 - la) / d).astype(int), np.rint(((lo - lo1) % 360) / d).astype(int)
        cur = v[np.ix_(lats, lons)]            # m accumulated since the run's start
        frames.append(np.zeros_like(cur) if prev is None else np.clip(cur - prev, 0, None) * 1000.0)
        prev = cur

    ny, nx = frames[0].shape
    q = np.clip(np.rint(np.stack(frames) * 100), 0, FILL - 1).astype(np.uint16).reshape(len(frames) * ny, nx)
    img = np.zeros((len(frames) * ny, nx, 3), np.uint8)
    img[..., 0], img[..., 1] = q >> 8, q & 255
    os.makedirs(out, exist_ok=True)
    name = "p%d.png" % run
    Image.fromarray(img, "RGB").save(os.path.join(out, name), optimize=True)
    meta = dict(model="ecmwf_ifs_open_0p1", run=run, generated=int(time.time()), file=name, nx=nx, ny=ny,
                frames=len(frames), step=STEP, lonW=BOX["lonW"], latN=BOX["latN"], dx=RES, dy=RES, scale=0.01, fill=FILL,
                source="ECMWF open data (CC BY 4.0)")
    with open(os.path.join(out, "latest.json"), "w") as f: json.dump(meta, f)
    print("run %s: %d frames of %d x %d, tp %.1f MB, png %.1f KB (max 3 h %.1f mm)" % (
        time.strftime("%Y-%m-%d %HZ", time.gmtime(run)), len(frames), nx, ny, got / 1e6,
        os.path.getsize(os.path.join(out, name)) / 1024, float(np.max(frames))))


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
    ap.add_argument("--published", default="")
    ap.add_argument("--run", default="")
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    pub = 0
    if a.published and os.path.exists(a.published):
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
