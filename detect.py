#!/usr/bin/env python3
"""Blind moving-object search in a SPHEREx sky region, scored against JPL known asteroids.

1. Seeds: local peaks flagged MP_OUTLIER by the SPHEREx pipeline, not on a
   catalogued source (MP_SOURCE), no bad pixels, aperture SNR >= 5.
2. Link seeds from different visits into straight-line tracks (>= 3 visits).
3. Match tracks to JPL sb_ident predictions: known asteroid or new candidate.
Also writes one PNG per visit, resampled onto a common sky grid, for blinking.

usage: python3 detect.py [RA_DEG] [DEC_DEG] [BOX_DEG] [DETECTOR] [DAYS]  -> web/
"""
import itertools, json, os, struct, sys, warnings, zlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from astropy.wcs import WCS
import spherex as S
from recover import BAD, R, SNR_MIN, Visit, measure, to_u8, visits

warnings.filterwarnings("ignore", module="astropy")
PIX = 6.15            # arcsec per SPHEREx pixel
MAX_RATE = 150.0      # arcsec/hour; main belt ~20-60, NEOs faster
# ponytail: rejects slow movers (distant TNOs, Planet X) along with drifting stars; needs
# multi-month baselines (SPHEREx revisits every 6 months) to reach them
MIN_RATE = 2.0
LINK_TOL = 12.0       # arcsec, ~2 px
MATCH_TOL = 12.0      # arcsec, track point vs JPL prediction
OUTLIER, SOURCE = 1 << 19, 1 << 21


def tangent(ra, dec, ra0, dec0):
    """Small-field offsets in arcsec (east, north) from (ra0, dec0)."""
    return (np.asarray(ra) - ra0) * np.cos(np.radians(dec0)) * 3600, (np.asarray(dec) - dec0) * 3600


def grid_wcs(ra0, dec0, n):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [n / 2 + 0.5, n / 2 + 0.5]
    w.wcs.cdelt = [-PIX / 3600, PIX / 3600]
    return w


def region(v, gw, n):
    """Read IMAGE/VARIANCE/FLAGS for the visit-pixel bounding box of the sky grid."""
    c = np.array([[0, 0], [n - 1, 0], [0, n - 1], [n - 1, n - 1]], float)
    ra, dec = gw.pixel_to_world_values(c[:, 0], c[:, 1])
    x, y = v.wcs.world_to_pixel_values(ra, dec)
    x0, y0 = int(np.floor(x.min())) - R, int(np.floor(y.min())) - R
    w, h = int(np.ceil(x.max())) + R - x0 + 1, int(np.ceil(y.max())) + R - y0 + 1
    img, var, flg = (v.fits.read(e, x0, y0, w, h) for e in ("IMAGE", "VARIANCE", "FLAGS"))
    return x0, y0, img, var, flg.astype(np.int64)


def seeds(v, x0, y0, img, var, flg):
    """Outlier-flagged point sources in one visit -> [{ra, dec, x, y, snr, flux, stamp}]."""
    f = np.nan_to_num(img, nan=-np.inf)
    pad = np.pad(f, 1, constant_values=-np.inf)
    peak = np.ones_like(f, bool)
    near_outlier = np.zeros_like(f, bool)
    po = np.pad((flg & OUTLIER) != 0, 1)
    for dy, dx in itertools.product((-1, 0, 1), repeat=2):
        near_outlier |= po[1 + dy:1 + dy + f.shape[0], 1 + dx:1 + dx + f.shape[1]]
        if dy or dx:
            peak &= f > pad[1 + dy:1 + dy + f.shape[0], 1 + dx:1 + dx + f.shape[1]]
    ok = peak & near_outlier & ((flg & (SOURCE | BAD)) == 0) & np.isfinite(img)
    ok[:R, :] = ok[-R:, :] = ok[:, :R] = ok[:, -R:] = False
    out = []
    for yy, xx in zip(*np.nonzero(ok)):
        st, sv = img[yy - R:yy + R + 1, xx - R:xx + R + 1], var[yy - R:yy + R + 1, xx - R:xx + R + 1]
        core = np.clip(np.nan_to_num(st[R - 1:R + 2, R - 1:R + 2] - np.nanmedian(st)), 0, None)
        if core.sum() <= 0:
            continue
        dy, dx = (np.mgrid[-1:2, -1:2] * core).sum(axis=(1, 2)) / core.sum()  # 3x3 centroid
        m = measure(st, sv, dx, dy)
        if not m or m["snr"] < SNR_MIN:
            continue
        x, y = x0 + xx + dx, y0 + yy + dy
        ra, dec = v.wcs.pixel_to_world_values(x, y)
        out.append({"visit": v.id, "mjd": v.mjd, "ra": float(ra), "dec": float(dec), "x": float(x), "y": float(y),
                    "snr": round(m["snr"], 1), "flux": m["flux"],
                    "stamp": to_u8(st, np.nanmedian(st), max(np.nanmax(st), np.nanmedian(st) + 1e-6))})
    return out


def link(dets, ra0, dec0, min_visits=3):
    """Straight-line tracks through detections from >= min_visits different visits."""
    by_visit = {}
    for i, d in enumerate(dets):
        d["e"], d["n"] = (float(v) for v in tangent(d["ra"], d["dec"], ra0, dec0))
        by_visit.setdefault(d["visit"], []).append(i)
    vids = sorted(by_visit, key=lambda k: dets[by_visit[k][0]]["mjd"])
    tracks = set()
    for a, b in itertools.combinations(range(len(vids)), 2):
        for i in by_visit[vids[a]]:
            for j in by_visit[vids[b]]:
                di, dj = dets[i], dets[j]
                dt = (dj["mjd"] - di["mjd"]) * 24
                ve, vn = (dj["e"] - di["e"]) / dt, (dj["n"] - di["n"]) / dt
                if not MIN_RATE <= np.hypot(ve, vn) <= MAX_RATE or np.hypot(ve, vn) * dt < 4 * PIX:  # too fast / not moving
                    continue
                members = {i, j}
                for c in range(len(vids)):
                    if c in (a, b):
                        continue
                    t = (dets[by_visit[vids[c]][0]]["mjd"] - di["mjd"]) * 24
                    pe, pn = di["e"] + ve * t, di["n"] + vn * t
                    best = min(by_visit[vids[c]], key=lambda k: np.hypot(dets[k]["e"] - pe, dets[k]["n"] - pn))
                    if np.hypot(dets[best]["e"] - pe, dets[best]["n"] - pn) < LINK_TOL:
                        members.add(best)
                if len(members) >= min_visits:
                    tracks.add(frozenset(members))
    # keep maximal tracks only: drop any track contained in a bigger one
    tracks = sorted(tracks, key=len, reverse=True)
    kept = [t for k, t in enumerate(tracks) if not any(t < u for u in tracks[:k])]
    out = []
    for t in kept:
        pts = sorted((dets[i] for i in t), key=lambda d: d["mjd"])
        th = np.array([(p["mjd"] - pts[0]["mjd"]) * 24 for p in pts])
        fe, fn = np.polyfit(th, [p["e"] for p in pts], 1), np.polyfit(th, [p["n"] for p in pts], 1)
        res = np.hypot(np.polyval(fe, th) - [p["e"] for p in pts], np.polyval(fn, th) - [p["n"] for p in pts])
        out.append({"points": pts, "rate": round(float(np.hypot(fe[0], fn[0])), 1),
                    "rms": round(float(np.sqrt(np.mean(res ** 2))), 2)})
    return out


def match(tracks, known_by_visit):
    """Label tracks with the JPL known asteroid most of their points sit on."""
    for t in tracks:
        names = []
        for p in t["points"]:
            ks = known_by_visit[p["visit"]]
            d = [np.hypot(*tangent(k["ra"], k["dec"], p["ra"], p["dec"])) for k in ks]
            names.append(ks[int(np.argmin(d))]["name"] if d and min(d) < MATCH_TOL else None)
        best = max(set(names), key=names.count)
        t["known"] = best if best and names.count(best) * 2 > len(names) else None
        t["vmag"] = next((k["vmag"] for ks in known_by_visit.values() for k in ks if k["name"] == t["known"]), None)
        t["status"] = "known" if t["known"] else "candidate"
    return tracks


def png(path, a):
    """Write a uint8 grayscale array as PNG (stdlib only)."""
    h, w = a.shape
    raw = b"".join(b"\x00" + a[y].tobytes() for y in range(h))
    chunk = lambda t, d: struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d))
    open(path, "wb").write(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0))
                           + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def resample(v, x0, y0, img, gw, n):
    """Nearest-neighbour resample of a visit's region onto the common sky grid."""
    gx, gy = np.meshgrid(np.arange(n), np.arange(n))
    ra, dec = gw.pixel_to_world_values(gx.ravel(), gy.ravel())
    x, y = v.wcs.world_to_pixel_values(ra, dec)
    xi, yi = np.round(x).astype(int) - x0, np.round(y).astype(int) - y0
    inside = (xi >= 0) & (yi >= 0) & (xi < img.shape[1]) & (yi < img.shape[0])
    out = np.full(n * n, np.nan)
    out[inside] = img[yi[inside], xi[inside]]
    return out.reshape(n, n)


def run(ra0, dec0, box, det, days):
    n = int(round(box * 3600 / PIX))
    gw = grid_wcs(ra0, dec0, n)
    with ThreadPoolExecutor(4) as ex:
        vs = list(ex.map(Visit, visits(ra0, dec0, det, days)))
    print(f"{len(vs)} visits, {n}x{n} px region", flush=True)
    os.makedirs("web/img", exist_ok=True)
    dets, frames = [], []
    for v in vs:
        x0, y0, img, var, flg = region(v, gw, n)
        s = seeds(v, x0, y0, img, var, flg)
        dets += s
        sky = resample(v, x0, y0, img, gw, n)
        lo, hi = np.nanpercentile(sky, 5), np.nanpercentile(sky, 99.7)
        png(f"web/img/{v.id}.png", np.array(to_u8(sky, lo, hi), np.uint8).reshape(n, n)[::-1])  # north up
        frames.append({"id": v.id, "mjd": v.mjd, "img": f"img/{v.id}.png"})
        print(f"  {v.id}: {len(s)} outlier seeds", flush=True)
    tracks = link(dets, ra0, dec0)
    known = {v.id: S.known_asteroids(float(v.row["s_ra"]), float(v.row["s_dec"]), v.mjd, half_width=2.0, vmag=20)
             for v in vs}
    tracks = match(tracks, known)
    # known asteroids bright enough to see (V<17.5) that sit inside the grid in >= 3 visits: blind completeness
    half = n * PIX / 2
    def in_grid(k):
        e, nn = tangent(k["ra"], k["dec"], ra0, dec0)
        return abs(e) < half - R * PIX and abs(nn) < half - R * PIX
    seen = {}
    for vid, ks in known.items():
        for k in ks:
            if k["vmag"] < 17.5 and in_grid(k):
                seen.setdefault(k["name"], []).append(vid)
    expected = sorted(nm for nm, vids in seen.items() if len(vids) >= 3)
    found = sorted({t["known"] for t in tracks if t["known"]})
    to_px = lambda ra, dec: [round(float(c), 1) for c in gw.world_to_pixel_values(ra, dec)]
    for t in tracks:
        for p in t["points"]:
            p["gx"], p["gy"] = to_px(p["ra"], p["dec"])
            p["gy"] = round(n - 1 - p["gy"], 1)
    overlay = {vid: [dict(name=k["name"], vmag=k["vmag"], gx=to_px(k["ra"], k["dec"])[0],
                          gy=round(n - 1 - to_px(k["ra"], k["dec"])[1], 1)) for k in ks if in_grid(k)]
               for vid, ks in known.items()}
    return {"field": {"ra": ra0, "dec": dec0, "box_deg": box, "detector": det, "n": n, "pix": PIX},
            "frames": frames, "tracks": tracks, "known": overlay,
            "score": {"expected": expected, "found": found,
                      "candidates": sum(t["status"] == "candidate" for t in tracks),
                      "seeds": len(dets)}}


if __name__ == "__main__":
    d = [{"visit": v, "mjd": 1 + h / 24, "ra": 180 + 30 * h / 3600, "dec": 0.0} for v, h in (("a", 0), ("b", 10), ("c", 26))]
    d.append({"visit": "c", "mjd": 1 + 26 / 24, "ra": 180.5, "dec": 0.2})  # noise
    tr = link(d, 180, 0)
    assert len(tr) == 1 and len(tr[0]["points"]) == 3 and abs(tr[0]["rate"] - 30) < 0.1, tr
    ra = float(sys.argv[1]) if len(sys.argv) > 1 else 180.0
    dec = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    box = float(sys.argv[3]) if len(sys.argv) > 3 else 1.3
    det = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    days = float(sys.argv[5]) if len(sys.argv) > 5 else 2.6
    out = run(ra, dec, box, det, days)
    sc = out["score"]
    print(f"{sc['seeds']} seeds -> {len(out['tracks'])} tracks: {len(sc['found'])} known, {sc['candidates']} candidates")
    print(f"blind completeness (V<17.5, in >=3 visits): {len(set(sc['found']) & set(sc['expected']))}/{len(sc['expected'])}")
    for t in out["tracks"]:
        print(f"  {t['status']:9} {t['known'] or '':32} V={t['vmag']}  {len(t['points'])} pts  {t['rate']}\"/h  rms {t['rms']}\"")
    json.dump(out, open("web/detections.json", "w"))
