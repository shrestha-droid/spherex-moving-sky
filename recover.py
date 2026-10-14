#!/usr/bin/env python3
"""Recovery test: are JPL-predicted known asteroids actually visible in SPHEREx?

For each visit of a field, each known asteroid is measured at its predicted pixel,
then the same sky spot is measured in another visit hours later, where the
asteroid has moved away. Detected there and gone later = recovered.

usage: python3 recover.py [RA_DEG] [DEC_DEG] [DETECTOR] [DAYS]  -> web/recovery.json
"""
import json, os, sys, warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from astropy.wcs import WCS
import spherex as S

warnings.filterwarnings("ignore", module="astropy")
R = 10            # stamp half-size, px -> 21x21 stamps
AP, BG = 2.0, (5, 9)  # aperture radius and background annulus, px
SNR_MIN = 5.0
EDGE = 12
# FLAGS bits that make a pixel unusable (bit 19 MP_OUTLIER and 21 MP_SOURCE are kept)
BAD = sum(1 << b for b in (0, 1, 2, 4, 5, 6, 9, 10, 11, 15, 17))


def visits(ra, dec, det, days):
    """First frame of each visit (exposure bursts > ~1 h apart) within `days` of the first."""
    imgs = sorted((i for i in S.images(ra, dec) if i["energy_bandpassname"] == f"SPHEREx-D{det}"),
                  key=lambda i: float(i["t_min"]))
    out = []
    for i in imgs:
        t = float(i["t_min"])
        if out and t - float(out[-1]["t_min"]) < 0.05:
            continue
        if out and t - float(out[0]["t_min"]) > days:
            break
        out.append(i)
    return out


def measure(img, var, dx, dy):
    """Aperture flux and SNR at stamp center + (dx, dy) subpixel offset."""
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    r = np.hypot(xx - dx, yy - dy)
    ann = img[(r >= BG[0]) & (r <= BG[1])]
    ap = r <= AP
    if np.isnan(img[ap]).any() or np.isnan(ann).all():
        return None
    bg = np.nanmedian(ann)
    flux = float(np.sum(img[ap] - bg))
    noise = float(np.sqrt(np.sum(var[ap]) + ap.sum() ** 2 * np.nanvar(ann) / np.isfinite(ann).sum()))
    near = np.where(r <= 2.5, img, -np.inf)
    py, px = np.unravel_index(np.argmax(near), img.shape)
    return {"flux": flux, "snr": flux / noise if noise > 0 else 0.0, "bg": float(bg),
            "peak_off": float(np.hypot(px - R - dx, py - R - dy))}


def to_u8(stamp, lo, hi):
    """asinh stretch to 0..255 for the web viewer."""
    v = np.arcsinh(np.clip((np.nan_to_num(stamp, nan=lo) - lo) / max(hi - lo, 1e-9), 0, None) * 10) / np.arcsinh(10)
    return [int(p) for p in np.clip(v * 255, 0, 255).ravel()]


class Visit:
    def __init__(self, row):
        self.row, self.id = row, row["obs_id"]
        self.mjd = (float(row["t_min"]) + float(row["t_max"])) / 2
        self.fits = S.RemoteFits(row["access_url"])
        self.wcs = WCS(self.fits.header("IMAGE"))
        self.nx, self.ny = self.fits.header("IMAGE")["NAXIS1"], self.fits.header("IMAGE")["NAXIS2"]

    def pixel(self, ra, dec):
        x, y = self.wcs.world_to_pixel_values(ra, dec)
        x, y = float(x), float(y)
        inside = EDGE <= x < self.nx - EDGE and EDGE <= y < self.ny - EDGE
        return (x, y) if inside else None

    def stamp(self, x, y):
        X, Y = int(round(x)), int(round(y))
        img = self.fits.read("IMAGE", X - R, Y - R, 2 * R + 1, 2 * R + 1)
        var = self.fits.read("VARIANCE", X - R, Y - R, 2 * R + 1, 2 * R + 1)
        flags = self.fits.read("FLAGS", X - 1, Y - 1, 3, 3).astype(np.int64)
        img[R - 1:R + 2, R - 1:R + 2][(flags & BAD) != 0] = np.nan  # measure() rejects NaN apertures
        return img, var, x - X, y - Y


def run(ra, dec, det=1, days=2.6):
    rows = visits(ra, dec, det, days)
    print(f"{len(rows)} visits:", [r["obs_id"] for r in rows], flush=True)
    with ThreadPoolExecutor(4) as ex:
        vs = list(ex.map(Visit, rows))
    # every (asteroid, visit) where the asteroid is predicted inside the image
    jobs = []
    for v in vs:
        for a in S.known_asteroids(float(v.row["s_ra"]), float(v.row["s_dec"]), v.mjd, half_width=2.0, vmag=20):
            p = v.pixel(a["ra"], a["dec"])
            if not p:
                continue
            # comparison: the visit furthest in time that still covers this sky spot
            others = [(abs(o.mjd - v.mjd), o) for o in vs if o is not v and abs(o.mjd - v.mjd) > 0.2]
            others = [(dt, o, o.pixel(a["ra"], a["dec"])) for dt, o in others]
            others = [t for t in others if t[2]]
            cmp_v, cmp_p = (max(others, key=lambda t: t[0])[1:] if others else (None, None))
            jobs.append((a, v, p, cmp_v, cmp_p))
    print(f"{len(jobs)} predicted detections to check", flush=True)

    def check(job):
        a, v, p, cv, cp = job
        img, var, dx, dy = v.stamp(*p)
        m = measure(img, var, dx, dy)
        res = {"name": a["name"], "vmag": a["vmag"],
               "visit": v.id, "mjd": v.mjd, "x": round(p[0], 2), "y": round(p[1], 2), "m": m}
        # scale to the asteroid's own peak; the comparison keeps that contrast over its own background
        lo = np.nanmedian(img)
        span = max(np.nanmax(img[R - 1:R + 2, R - 1:R + 2]) - lo, 1e-6) if np.isfinite(img[R - 1:R + 2, R - 1:R + 2]).any() else 1.0
        res["stamp"] = to_u8(img, lo, lo + span)
        if cv and m:
            cimg, cvar, cdx, cdy = cv.stamp(*cp)
            clo = np.nanmedian(cimg)
            res["cmp"] = {"visit": cv.id, "m": measure(cimg, cvar, cdx, cdy), "stamp": to_u8(cimg, clo, clo + span)}
        res["status"] = classify(res)
        return res

    with ThreadPoolExecutor(4) as ex:
        results = list(ex.map(check, jobs))
    return vs, results


def classify(res):
    m = res["m"]
    if not m:
        return "bad pixels"
    if m["snr"] < SNR_MIN or m["peak_off"] > 1.5:
        return "not detected"
    c = res.get("cmp", {}).get("m")
    if c is None:
        return "detected" if "cmp" not in res else "bad pixels"
    return "recovered" if c["flux"] < 0.3 * m["flux"] else "blended"


def summarize(results):
    bins = [(None, 15), (15, 16), (16, 17), (17, 18), (18, 19), (19, 20)]
    out = []
    for lo, hi in bins:
        sel = [r for r in results if (lo is None or r["vmag"] >= lo) and r["vmag"] < hi and r["status"] != "bad pixels"]
        rec = sum(r["status"] in ("recovered", "detected") for r in sel)
        out.append({"label": f"V<{hi}" if lo is None else f"{lo}-{hi}", "n": len(sel), "recovered": rec})
    return out


if __name__ == "__main__":
    assert classify({"m": {"snr": 8, "flux": 10, "peak_off": 0.5}, "cmp": {"m": {"flux": 1}}}) == "recovered"
    assert classify({"m": {"snr": 8, "flux": 10, "peak_off": 0.5}, "cmp": {"m": {"flux": 9}}}) == "blended"
    assert classify({"m": {"snr": 3, "flux": 10, "peak_off": 0.5}}) == "not detected"
    ra = float(sys.argv[1]) if len(sys.argv) > 1 else 180.0
    dec = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    det = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    days = float(sys.argv[4]) if len(sys.argv) > 4 else 2.6
    vs, results = run(ra, dec, det, days)
    summary = summarize(results)
    for s in summary:
        print(f"  {s['label']:>6}: {s['recovered']}/{s['n']}")
    os.makedirs("web", exist_ok=True)
    json.dump({"field": {"ra": ra, "dec": dec, "detector": det, "stamp": 2 * R + 1},
               "visits": [{"id": v.id, "mjd": v.mjd} for v in vs],
               "summary": summary, "results": results}, open("web/recovery.json", "w"))
