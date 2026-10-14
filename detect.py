#!/usr/bin/env python3
"""Blind moving-object search in a SPHEREx sky region, scored against JPL known asteroids.

Fast mode (asteroids, 2-150"/h, >= 3 visits over ~2.6 days):
  seeds are peaks the SPHEREx pipeline flagged MP_OUTLIER, not on a catalogued source.
Slow mode (distant worlds, 0.3-3"/h, >= 4 visits over ~a week):
  seeds are any uncatalogued peak that is missing from the median of the other visits
  at the same sky spot (faint stars are in every visit; a slow mover is not).
Both: link seeds into straight-line tracks, reject tracks that sit on the same detector
pixel in different visits (fixed defects that ride the pointing), measure the chance rate
by scrambling, match to JPL, and optionally inject synthetic movers to measure depth.

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
OMEGA = (PIX / 206265) ** 2
FAST = dict(min_rate=2.0, max_rate=150.0, min_visits=3)   # arcsec/hour; main belt ~20-60
SLOW = dict(min_rate=0.3, max_rate=3.0, min_visits=4)     # TNOs / dwarf planets / Planet X
LINK_TOL = 12.0       # arcsec, ~2 px
MATCH_TOL = 12.0      # arcsec, track point vs JPL prediction
DETECTOR_TOL = 2.0    # px: two points this close on the detector = a fixed defect, not a mover
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


class Region:
    """IMAGE/VARIANCE/FLAGS of one visit, for the visit-pixel bounding box of the sky grid."""

    def __init__(self, v, gw, n):
        c = np.array([[0, 0], [n - 1, 0], [0, n - 1], [n - 1, n - 1]], float)
        ra, dec = gw.pixel_to_world_values(c[:, 0], c[:, 1])
        x, y = v.wcs.world_to_pixel_values(ra, dec)
        self.v, self.x0, self.y0 = v, int(np.floor(x.min())) - R, int(np.floor(y.min())) - R
        w, h = int(np.ceil(x.max())) + R - self.x0 + 1, int(np.ceil(y.max())) + R - self.y0 + 1
        self.img, self.var, flg = (v.fits.read(e, self.x0, self.y0, w, h) for e in ("IMAGE", "VARIANCE", "FLAGS"))
        self.flg = flg.astype(np.int64)

    def local(self, ra, dec):
        """Region-array coordinates of sky positions."""
        x, y = self.v.wcs.world_to_pixel_values(ra, dec)
        return np.asarray(x) - self.x0, np.asarray(y) - self.y0

    def stamp(self, ra, dec, img=None):
        img = self.img if img is None else img
        lx, ly = self.local(ra, dec)
        X, Y = int(round(float(lx))), int(round(float(ly)))
        if not (R <= X < img.shape[1] - R and R <= Y < img.shape[0] - R):
            return None
        return img[Y - R:Y + R + 1, X - R:X + R + 1]

    def ap_flux(self, ra, dec, img=None):
        """Quick 3x3 sum minus local median, at sky positions (vectorized). NaN off-region."""
        img = self.img if img is None else img
        lx, ly = self.local(ra, dec)
        X, Y = np.round(lx).astype(int), np.round(ly).astype(int)
        ok = (X >= 6) & (Y >= 6) & (X < img.shape[1] - 6) & (Y < img.shape[0] - 6)
        out = np.full(X.shape, np.nan)
        for k in np.nonzero(ok)[0]:
            box = img[Y[k] - 6:Y[k] + 7, X[k] - 6:X[k] + 7]
            out[k] = np.nansum(box[5:8, 5:8] - np.nanmedian(box))
        return out


def seeds(reg, img, others=(), require_outlier=True):
    """Point sources in one visit -> [{ra, dec, x, y, snr, flux, stamp}].

    Fast mode: must touch an MP_OUTLIER pixel. Slow mode: must be absent (< half as bright)
    in the median of the other visits at the same sky spot."""
    v, flg, var = reg.v, reg.flg, reg.var
    f = np.nan_to_num(img, nan=-np.inf)
    pad = np.pad(f, 1, constant_values=-np.inf)
    peak = np.ones_like(f, bool)
    near_outlier = np.zeros_like(f, bool)
    po = np.pad((flg & OUTLIER) != 0, 1)
    for dy, dx in itertools.product((-1, 0, 1), repeat=2):
        near_outlier |= po[1 + dy:1 + dy + f.shape[0], 1 + dx:1 + dx + f.shape[1]]
        if dy or dx:
            peak &= f > pad[1 + dy:1 + dy + f.shape[0], 1 + dx:1 + dx + f.shape[1]]
    ok = peak & ((flg & (SOURCE | BAD)) == 0) & np.isfinite(img)
    if require_outlier:
        ok &= near_outlier
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
        x, y = reg.x0 + xx + dx, reg.y0 + yy + dy
        ra, dec = (float(c) for c in v.wcs.pixel_to_world_values(x, y))
        out.append({"visit": v.id, "mjd": v.mjd, "ra": ra, "dec": dec, "x": float(x), "y": float(y),
                    "snr": round(m["snr"], 1), "flux": m["flux"]})
    if not require_outlier and others and out:
        here = reg.ap_flux(np.array([d["ra"] for d in out]), np.array([d["dec"] for d in out]), img)
        elsewhere = np.nanmedian(np.array([o.ap_flux(np.array([d["ra"] for d in out]), np.array([d["dec"] for d in out]))
                                           for o in others]), axis=0)
        out = [d for d, h, e in zip(out, here, elsewhere) if np.isfinite(e) and e < 0.5 * h]
    return out


def link(dets, ra0, dec0, min_visits=3, min_rate=2.0, max_rate=150.0, tol=LINK_TOL):
    """Straight-line tracks through detections from >= min_visits different visits.

    Returns (tracks, rejected): rejected tracks reuse one detector pixel across visits."""
    if not dets:
        return [], []
    for d in dets:
        d["e"], d["n"] = (float(c) for c in tangent(d["ra"], d["dec"], ra0, dec0))
    by_visit = {}
    for i, d in enumerate(dets):
        by_visit.setdefault(d["visit"], []).append(i)
    vids = sorted(by_visit, key=lambda k: dets[by_visit[k][0]]["mjd"])
    E = np.array([d["e"] for d in dets]); N = np.array([d["n"] for d in dets])
    T = np.array([d["mjd"] for d in dets]) * 24
    tvis = [T[by_visit[k][0]] for k in vids]
    idx = [np.array(by_visit[k]) for k in vids]
    tracks = set()
    for a, b in itertools.combinations(range(len(vids)), 2):
        ia, ib = idx[a], idx[b]
        dt = tvis[b] - tvis[a]
        de, dn = E[ib][None, :] - E[ia][:, None], N[ib][None, :] - N[ia][:, None]
        disp = np.hypot(de, dn)
        ok = (disp / dt >= min_rate) & (disp / dt <= max_rate) & (disp >= 4 * PIX)
        for pa, pb in zip(*np.nonzero(ok)):
            i, j = ia[pa], ib[pb]
            ve, vn = de[pa, pb] / dt, dn[pa, pb] / dt
            members = {i, j}
            for c in range(len(vids)):
                if c in (a, b):
                    continue
                pe, pn = E[i] + ve * (tvis[c] - tvis[a]), N[i] + vn * (tvis[c] - tvis[a])
                d = np.hypot(E[idx[c]] - pe, N[idx[c]] - pn)
                k = int(np.argmin(d))
                if d[k] < tol:
                    members.add(idx[c][k])
            if len(members) >= min_visits:
                tracks.add(frozenset(members))
    # keep maximal tracks only: drop any track contained in a bigger one
    tracks = sorted(tracks, key=len, reverse=True)
    kept = [t for k, t in enumerate(tracks) if not any(t < u for u in tracks[:k])]
    out, rejected = [], []
    for t in kept:
        pts = sorted((dets[i] for i in t), key=lambda d: d["mjd"])
        th = np.array([(p["mjd"] - pts[0]["mjd"]) * 24 for p in pts])
        fe, fn = np.polyfit(th, [p["e"] for p in pts], 1), np.polyfit(th, [p["n"] for p in pts], 1)
        res = np.hypot(np.polyval(fe, th) - [p["e"] for p in pts], np.polyval(fn, th) - [p["n"] for p in pts])
        tr = {"points": pts, "rate": round(float(np.hypot(fe[0], fn[0])), 2),
              "rms": round(float(np.sqrt(np.mean(res ** 2))), 2)}
        same_pixel = any(np.hypot(p["x"] - q["x"], p["y"] - q["y"]) < DETECTOR_TOL
                         for p, q in itertools.combinations(pts, 2) if "x" in p and "x" in q)
        (rejected if same_pixel else out).append(tr)
    return out, rejected


def false_alarms(dets, ra0, dec0, trials=30, seed=1, **kw):
    """Mean number of tracks linked after rotating each visit's seeds by a random angle.

    Rotation keeps seed density but destroys real motion, so any track left is chance."""
    rng = np.random.default_rng(seed)
    counts = []
    for _ in range(trials):
        fake, ang = [], {}
        for d in dets:
            a = ang.setdefault(d["visit"], rng.uniform(0, 2 * np.pi))
            e, n = tangent(d["ra"], d["dec"], ra0, dec0)
            e, n = e * np.cos(a) - n * np.sin(a), e * np.sin(a) + n * np.cos(a)
            fake.append({k: v for k, v in d.items() if k not in ("x", "y")} |
                        {"ra": ra0 + e / 3600 / np.cos(np.radians(dec0)), "dec": dec0 + n / 3600})
        counts.append(len(link(fake, ra0, dec0, **kw)[0]))
    return float(np.mean(counts))


def match(tracks, known_by_visit):
    """Label tracks with the JPL known object most of their points sit on."""
    for t in tracks:
        names = []
        for p in t["points"]:
            ks = known_by_visit.get(p["visit"], [])
            d = [np.hypot(*tangent(k["ra"], k["dec"], p["ra"], p["dec"])) for k in ks]
            names.append(ks[int(np.argmin(d))]["name"] if d and min(d) < MATCH_TOL else None)
        best = max(set(names), key=names.count)
        t["known"] = best if best and names.count(best) * 2 > len(names) else None
        t["vmag"] = next((k["vmag"] for ks in known_by_visit.values() for k in ks if k["name"] == t["known"]), None)
        t["status"] = "known" if t["known"] else "candidate"
    return tracks


def vet(t, tol=30.0):
    """Deep JPL check for a candidate: no magnitude limit, comets included, around its middle point."""
    p = t["points"][len(t["points"]) // 2]
    def fetch():
        xobs = ",".join(f"{v:.6f}" for v in S.spherex_xyz(p["mjd"]))
        body, _ = S.get("https://ssd-api.jpl.nasa.gov/sb_ident.api", {
            "xobs": xobs, "obs-time": f"{p['mjd'] + 2400000.5:.6f}", "mag-required": "false",
            "fov-ra-center": S.sexa(p["ra"], hours=True), "fov-dec-center": S.sexa(p["dec"]),
            "fov-ra-hwidth": 0.03, "fov-dec-hwidth": 0.03})
        return json.loads(body)
    res = S.cached_json(f"vet_{p['ra']:.5f}_{p['dec']:.5f}_{p['mjd']:.5f}.json", fetch)
    fields = res.get("fields_first", [])
    for row in res.get("data_first_pass", []):
        r = dict(zip(fields, row))
        ra, dec = S.parse_sexa(row[1], hours=True), S.parse_sexa(row[2])
        if np.hypot(*tangent(ra, dec, p["ra"], p["dec"])) < tol:
            t.update(known=row[0], status="known (deep check)", vmag=r.get("Visual magnitude (V)"))
            return t
    t["vetted"] = True
    return t


def vetting_stamps(t, regions):
    """Per point: the visit's stamp, plus the median of every other visit at the same sky spot.

    A real mover is bright in its own stamp and absent from the median; a star is in both."""
    for p in t["points"]:
        own = regions[p["visit"]].stamp(p["ra"], p["dec"])
        rest = [r.stamp(p["ra"], p["dec"]) for vid, r in regions.items() if vid != p["visit"]]
        rest = [s for s in rest if s is not None]
        if own is None:
            continue
        lo = np.nanmedian(own)
        span = max(np.nanmax(own[R - 1:R + 2, R - 1:R + 2]) - lo, 1e-6)
        p["stamp"] = to_u8(own, lo, lo + span)
        if rest:
            med = np.nanmedian(np.array(rest), axis=0)
            p["others"] = to_u8(med, np.nanmedian(med), np.nanmedian(med) + span)


def inject(regions, gw, n, ra0, dec0, rate_range, mags, per_mag, seed=7):
    """Copy each visit's image and add synthetic point sources moving in straight lines.

    The real PSF of the zone each source lands in is used. Returns ({visit: image}, truth)."""
    rng = np.random.default_rng(seed)
    imgs = {vid: r.img.copy() for vid, r in regions.items()}
    t0 = min(r.v.mjd for r in regions.values()) * 24
    truth = []
    half = n * PIX / 2 - 3 * R * PIX
    for m in mags:
        total = 10 ** ((23.9 - m) / 2.5) * 1e-12 / OMEGA   # AB mag -> sum of MJy/sr pixel values
        for _ in range(per_mag):
            e0, n0 = rng.uniform(-half, half, 2)
            rate, ang = rng.uniform(*rate_range), rng.uniform(0, 2 * np.pi)
            pos = {}
            for vid, r in regions.items():
                dt = r.v.mjd * 24 - t0
                e, nn = e0 + rate * dt * np.cos(ang), n0 + rate * dt * np.sin(ang)
                ra, dec = ra0 + e / 3600 / np.cos(np.radians(dec0)), dec0 + nn / 3600
                lx, ly = r.local(ra, dec)
                X, Y = int(np.floor(lx)), int(np.floor(ly))
                if not (R <= X < r.img.shape[1] - R and R <= Y < r.img.shape[0] - R):
                    continue
                k = r.v.fits.psf_kernel(X + r.x0, Y + r.y0, float(lx) - X, float(ly) - Y)
                h = k.shape[0] // 2
                imgs[vid][Y - h:Y + h + 1, X - h:X + h + 1] += total * k
                pos[vid] = (float(ra), float(dec))
            truth.append({"mag": m, "rate": float(rate), "pos": pos})
    return imgs, truth


def recovered(truth, tracks, min_visits):
    """Did some track put >= min_visits points on this synthetic mover's true positions?"""
    for t in tracks:
        hits = sum(1 for p in t["points"] if p["visit"] in truth["pos"] and
                   np.hypot(*tangent(p["ra"], p["dec"], *truth["pos"][p["visit"]])) < LINK_TOL)
        if hits >= min_visits:
            return True
    return False


def png(path, a):
    """Write a uint8 grayscale array as PNG (stdlib only)."""
    h, w = a.shape
    raw = b"".join(b"\x00" + a[y].tobytes() for y in range(h))
    chunk = lambda t, d: struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d))
    open(path, "wb").write(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0))
                           + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def resample(reg, gw, n):
    """Nearest-neighbour resample of a visit's region onto the common sky grid."""
    gx, gy = np.meshgrid(np.arange(n), np.arange(n))
    ra, dec = gw.pixel_to_world_values(gx.ravel(), gy.ravel())
    lx, ly = reg.local(ra, dec)
    xi, yi = np.round(lx).astype(int), np.round(ly).astype(int)
    inside = (xi >= 0) & (yi >= 0) & (xi < reg.img.shape[1]) & (yi < reg.img.shape[0])
    out = np.full(n * n, np.nan)
    out[inside] = reg.img[yi[inside], xi[inside]]
    return out.reshape(n, n)


def run(ra0, dec0, box, det, days, pass_index=0, out_dir="web", rows=None, slow=False,
        injection=None, known_fn=None):
    """Search one sky box. rows: explicit visit rows (else the first frame of each visit).

    injection: dict(mags=[...], per_mag=N) to also measure recovery of synthetic movers.
    known_fn: visit -> [known objects]; defaults to JPL sb_ident asteroids."""
    mode = SLOW if slow else FAST
    n = int(round(box * 3600 / PIX))
    gw = grid_wcs(ra0, dec0, n)
    rows = rows if rows is not None else visits(ra0, dec0, det, days, pass_index)
    until = "WCS-WAVE" if injection else "VARIANCE"
    with ThreadPoolExecutor(4) as ex:
        vs = list(ex.map(lambda r: Visit(r, until=until), rows))
    print(f"{len(vs)} visits, {n}x{n} px region, {'slow' if slow else 'fast'} mode", flush=True)
    os.makedirs(f"{out_dir}/img", exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        regions = dict(zip([v.id for v in vs], ex.map(lambda v: Region(v, gw, n), vs)))
    def search(imgs):
        dets = []
        for vid, r in regions.items():
            others = [o for k, o in regions.items() if k != vid]
            dets += seeds(r, imgs[vid], others, require_outlier=not slow)
        tracks, rejected = link(dets, ra0, dec0, **mode)
        return dets, tracks, rejected
    dets, tracks, rejected = search({vid: r.img for vid, r in regions.items()})
    frames = []
    for vid, r in regions.items():
        sky = resample(r, gw, n)
        lo, hi = np.nanpercentile(sky, 5), np.nanpercentile(sky, 99.7)
        png(f"{out_dir}/img/{vid}.png", (np.array(to_u8(sky, lo, hi), np.uint8).reshape(n, n)[::-1] // 4) * 4)  # north up
        frames.append({"id": vid, "mjd": r.v.mjd, "img": f"img/{vid}.png"})
    print(f"  {len(dets)} seeds -> {len(tracks)} tracks, {len(rejected)} rejected as detector artifacts", flush=True)
    fa = false_alarms(dets, ra0, dec0, **mode)
    print(f"  chance tracks expected (scrambled seeds): {fa:.2f}", flush=True)
    if known_fn is None:
        known_fn = lambda v: S.known_asteroids(float(v.row["s_ra"]), float(v.row["s_dec"]), v.mjd, half_width=2.0, vmag=20)
    known = {v.id: known_fn(v) for v in vs}
    tracks = match(tracks, known)
    if not slow:
        tracks = [vet(t) if t["status"] == "candidate" else t for t in tracks]
    for t in tracks + rejected:
        vetting_stamps(t, regions)
    half = n * PIX / 2
    def in_grid(k):
        e, nn = tangent(k["ra"], k["dec"], ra0, dec0)
        return abs(e) < half - R * PIX and abs(nn) < half - R * PIX
    seen = {}
    for vid, ks in known.items():
        for k in ks:
            if k["vmag"] < 17.5 and in_grid(k):
                seen.setdefault(k["name"], []).append(vid)
    expected = sorted(nm for nm, vids in seen.items() if len(vids) >= mode["min_visits"])
    found = sorted({t["known"] for t in tracks if t["known"]})
    to_px = lambda ra, dec: [round(float(c), 1) for c in gw.world_to_pixel_values(ra, dec)]
    for t in tracks + rejected:
        for p in t["points"]:
            p["gx"], p["gy"] = to_px(p["ra"], p["dec"])
            p["gy"] = round(n - 1 - p["gy"], 1)
            for k in ("e", "n", "flux"):
                p.pop(k, None)
    overlay = {vid: [dict(name=k["name"], vmag=k["vmag"], gx=to_px(k["ra"], k["dec"])[0],
                          gy=round(n - 1 - to_px(k["ra"], k["dec"])[1], 1)) for k in ks if in_grid(k)]
               for vid, ks in known.items()}
    out = {"field": {"ra": ra0, "dec": dec0, "box_deg": box, "detector": det, "n": n, "pix": PIX,
                     "mode": "slow" if slow else "fast", **mode},
           "frames": frames, "tracks": tracks, "known": overlay,
           "rejected": [{"rate": t["rate"], "rms": t["rms"], "points": t["points"]} for t in rejected],
           "score": {"expected": expected, "found": found,
                     "candidates": sum(t["status"] == "candidate" for t in tracks),
                     "seeds": len(dets), "false_alarms": round(fa, 2), "artifacts": len(rejected)}}
    if injection:
        curve = []
        for i, m in enumerate(injection["mags"]):  # one magnitude at a time, so fakes never crowd each other
            imgs, truth = inject(regions, gw, n, ra0, dec0, (mode["min_rate"], mode["max_rate"]), [m],
                                 injection["per_mag"], seed=7 + i)
            _, itracks, _ = search(imgs)
            sel = [tr for tr in truth if len(tr["pos"]) >= mode["min_visits"]]
            curve.append({"mag_ab": m, "n": len(sel), "recovered": sum(recovered(tr, itracks, mode["min_visits"]) for tr in sel)})
        out["injection"] = curve
        print("  injection:", " ".join(f"AB{c['mag_ab']}:{c['recovered']}/{c['n']}" for c in curve), flush=True)
    return out


if __name__ == "__main__":
    d = [{"visit": v, "mjd": 1 + h / 24, "ra": 180 + 30 * h / 3600, "dec": 0.0, "x": 100 + 40 * i, "y": 50.0}
         for i, (v, h) in enumerate((("a", 0), ("b", 10), ("c", 26)))]
    d.append({"visit": "c", "mjd": 1 + 26 / 24, "ra": 180.5, "dec": 0.2, "x": 5.0, "y": 5.0})  # noise
    tr, rej = link(d, 180, 0)
    assert len(tr) == 1 and len(tr[0]["points"]) == 3 and abs(tr[0]["rate"] - 30) < 0.1 and not rej, tr
    for p in d:
        p["x"], p["y"] = 100.0, 50.0  # same detector pixel in every visit: a defect riding the pointing
    tr, rej = link(d, 180, 0)
    assert not tr and len(rej) == 1
    ra = float(sys.argv[1]) if len(sys.argv) > 1 else 180.0
    dec = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    box = float(sys.argv[3]) if len(sys.argv) > 3 else 1.3
    det = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    days = float(sys.argv[5]) if len(sys.argv) > 5 else 2.6
    out = run(ra, dec, box, det, days)
    sc = out["score"]
    print(f"{sc['seeds']} seeds -> {len(out['tracks'])} tracks: {len(sc['found'])} known, {sc['candidates']} candidates")
    for t in out["tracks"]:
        print(f"  {t['status']:9} {t['known'] or '':32} V={t['vmag']}  {len(t['points'])} pts  {t['rate']}\"/h  rms {t['rms']}\"")
    json.dump(out, open("web/detections.json", "w"))
