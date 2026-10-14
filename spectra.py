#!/usr/bin/env python3
"""Asteroid spectra: measure a known asteroid in every SPHEREx exposure it crosses.

SPHEREx's detectors work in pairs (D1+D4, D2+D5, D3+D6), each pair on its own patch
of sky, and each detector row sees a different wavelength. Each visit is a burst of 4
slightly shifted frames, so one visit puts an asteroid on ~8 different wavelengths in
~8 minutes, and more visits add more. JPL's predicted track places the asteroid in
each frame; the WCS-WAVE table gives the exact wavelength of that pixel.
"""
import warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import spherex as S
from recover import Visit, measure

warnings.filterwarnings("ignore", module="astropy")
OMEGA = (6.15 / 206265) ** 2   # sr per pixel
AP_SPEC = 3.0                  # px; wider than detection aperture so the long-wavelength PSF fits
# ponytail: fixed aperture, no PSF-based aperture correction; long-wavelength fluxes read
# somewhat low. Use the PSF extension per image when absolute fluxes matter.
SNR_SPEC = 3.0


def track(points):
    """(mjd, ra, dec) samples -> function t -> (ra, dec), quadratic in time (linear if 2 points)."""
    t = np.array([p[0] for p in points])
    deg = min(2, len(points) - 1)
    fr = np.polyfit(t - t.mean(), [p[1] for p in points], deg)
    fd = np.polyfit(t - t.mean(), [p[2] for p in points], deg)
    return lambda tt: (float(np.polyval(fr, tt - t.mean())), float(np.polyval(fd, tt - t.mean())))


def spectrum(name, points, max_images=64):
    """Measure one asteroid in all exposures (any detector) between its first and last JPL point."""
    points = sorted(points)
    pos = track(points)
    t0, t1 = points[0][0] - 0.02, points[-1][0] + 0.02
    ras, decs = [p[1] for p in points], [p[2] for p in points]
    mid_ra, mid_dec = (min(ras) + max(ras)) / 2, (min(decs) + max(decs)) / 2
    radius = max(np.hypot((max(ras) - min(ras)) * np.cos(np.radians(mid_dec)), max(decs) - min(decs)) / 2, 0.01) + 0.02
    rows = [r for r in S.images(mid_ra, mid_dec, radius) if t0 <= float(r["t_min"]) <= t1]
    rows = sorted(rows, key=lambda r: float(r["t_min"]))[:max_images]

    def one(row):
        v = Visit(row, until="WCS-WAVE")
        ra, dec = pos(v.mjd)
        p = v.pixel(ra, dec)
        if not p:
            return None
        img, var, dx, dy = v.stamp(*p)
        m = measure(img, var, dx, dy, ap_r=AP_SPEC)
        if not m or m["snr"] < SNR_SPEC:
            return None
        wave, bw = v.fits.wavelength(*p)
        if not np.isfinite(wave):
            return None
        f = m["flux"] * OMEGA * 1e12  # MJy -> microJansky
        return {"wave": round(wave, 4), "bw": round(bw, 4), "flux": round(f, 2), "err": round(f / m["snr"], 2),
                "snr": round(m["snr"], 1), "mjd": round(v.mjd, 5), "det": row["energy_bandpassname"][-2:], "obs": v.id}

    with ThreadPoolExecutor(4) as ex:
        pts = [p for p in ex.map(one, rows) if p]
    return {"name": name, "n_images": len(rows), "points": sorted(pts, key=lambda p: p["wave"])}


if __name__ == "__main__":
    f = track([(0.0, 10.0, 1.0), (1.0, 11.0, 1.5), (2.0, 12.0, 2.0)])
    assert np.allclose(f(1.5), (11.5, 1.75))
    # demo: 872 Holda through the May 2025 pass, positions from the recovery test
    import json, sys
    rec = json.load(open("web/recovery.json")) if len(sys.argv) < 2 else json.load(open(sys.argv[1]))
    from recover import visits
    vs = {r["obs_id"]: Visit(r) for r in visits(rec["field"]["ra"], rec["field"]["dec"], 1, 2.6)}
    pts = [(r["mjd"], *map(float, vs[r["visit"]].wcs.pixel_to_world_values(r["x"], r["y"])))
           for r in rec["results"] if "Holda" in r["name"]]
    sp = spectrum("872 Holda", pts)
    print(f"{sp['n_images']} exposures crossed, {len(sp['points'])} measured:")
    for p in sp["points"]:
        print(f"  {p['wave']:.3f} um  {p['flux']:8.1f} +/- {p['err']:6.1f} uJy  SNR {p['snr']:5.1f}  {p['det']}  {p['obs']}")
