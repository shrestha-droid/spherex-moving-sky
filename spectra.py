#!/usr/bin/env python3
"""Asteroid spectra: measure a known asteroid in every SPHEREx exposure it crosses.

SPHEREx's detectors work in pairs (D1+D4, D2+D5, D3+D6), each pair on its own patch
of sky, and each detector row sees a different wavelength. Each visit is a burst of 4
slightly shifted frames, so one visit puts an asteroid on ~8 different wavelengths in
~8 minutes, and more visits add more. JPL's predicted track places the asteroid in
each frame; the WCS-WAVE table gives the exact wavelength of that pixel.
"""
import json, warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import spherex as S
from recover import Visit, measure, measure_psf

warnings.filterwarnings("ignore", module="astropy")
OMEGA = (6.15 / 206265) ** 2   # sr per pixel
AP_SPEC = 3.0                  # px; the real PSF puts > 99.6% of the light inside this at every wavelength
SNR_SPEC = 3.0
ROT_MAG = 0.15                 # mag; typical asteroid spin brightness swing between visits hours apart
KNEE = 2.5                     # um; reflectance slope fitted below this, held flat beyond


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
        kernel = v.fits.psf_kernel(*p, dx, dy)     # this image's own PSF, for this detector zone
        m = measure_psf(img, var, dx, dy, kernel)  # optimal (PSF-weighted) extraction
        if not m or m["snr"] < SNR_SPEC:
            return None
        wave, bw = v.fits.wavelength(*p)
        if not np.isfinite(wave):
            return None
        omega = v.fits.pixel_area_sr(v.wcs, *p)    # local pixel solid angle from the TAN-SIP Jacobian
        f = m["flux"] * omega * 1e12               # sum of MJy/sr x sr -> MJy -> microJansky
        return {"wave": round(wave, 4), "bw": round(bw, 4), "flux": round(f, 2), "err": round(f / m["snr"], 2),
                "snr": round(m["snr"], 1), "mjd": round(v.mjd, 5), "det": row["energy_bandpassname"][-2:], "obs": v.id,
                "omega_rel": round(omega / OMEGA, 4), "visit": v.id.rsplit("_", 1)[0]}

    with ThreadPoolExecutor(4) as ex:
        pts = [p for p in ex.map(one, rows) if p]
    r_au, phase = geometry(name, float(np.mean([p[0] for p in points])))
    sp = {"name": name, "n_images": len(rows), "points": sorted(pts, key=lambda p: p["wave"]),
          "r_au": r_au, "phase_deg": phase, "photometry": "psf"}
    sp["model"] = fit_model(sp["points"], sp["r_au"])
    sp["jpl"] = sbdb(name)
    return sp


def slope_albedo_test(specs, trials=1000, seed=3):
    """Spearman rho between fitted slope and JPL albedo, with Monte Carlo over photometric noise
    plus a per-visit spin offset (ROT_MAG), refitting every slope each trial.

    Selection (stated, objective): spectra with a model fit (>= 6 points spanning >= 2.5 um) and a
    published SBDB albedo."""
    rng = np.random.default_rng(seed)
    use = [sp for sp in specs if sp.get("model") and sp.get("jpl", {}).get("albedo")]
    if len(use) < 5:
        return None
    alb = np.array([float(sp["jpl"]["albedo"]) for sp in use])
    rank = lambda a: np.argsort(np.argsort(a))
    rho = lambda a, b: float(np.corrcoef(rank(a), rank(b))[0, 1])
    nominal = np.array([sp["model"]["slope_per_um"] for sp in use])
    r0 = rho(nominal, alb)
    perm_p = float(np.mean([abs(rho(nominal, rng.permutation(alb))) >= abs(r0) for _ in range(20000)]))
    rhos = []
    for _ in range(trials):
        slopes = []
        for sp in use:
            pts = sp["points"]
            offs = {v: 10 ** (-0.4 * rng.normal(0, ROT_MAG)) for v in {p.get("visit") for p in pts}}
            pert = [{**p, "flux": p["flux"] * offs[p.get("visit")] + rng.normal(0, p["err"])} for p in pts]
            m = fit_model(pert, sp["r_au"])
            slopes.append(m["slope_per_um"] if m else np.nan)
        slopes = np.array(slopes)
        ok = np.isfinite(slopes)
        if ok.sum() >= 5:
            rhos.append(rho(slopes[ok], alb[ok]))
    lo, med, hi = np.percentile(rhos, [16, 50, 84])
    return {"n": len(use), "n_spectra": sum(1 for sp in specs if sp.get("points")), "rho": round(r0, 2),
            "perm_p": round(perm_p, 3), "mc_median": round(float(med), 2), "mc_16": round(float(lo), 2),
            "mc_84": round(float(hi), 2), "mc_frac_positive": round(float(np.mean(np.array(rhos) > 0)), 3),
            "trials": len(rhos), "selection": "model fit (>= 6 points over >= 2.5 um) and published SBDB albedo",
            "phase_deg": [sp.get("phase_deg") for sp in use]}


def sbdb(name):
    """Published physical data from JPL SBDB: taxonomic class, albedo, diameter (for comparison only)."""
    num = name.split()[0]
    def fetch():
        body, _ = S.get("https://ssd-api.jpl.nasa.gov/sbdb.api", {"sstr": num, "phys-par": "1"})
        return json.loads(body)
    try:
        phys = {p["name"]: p["value"] for p in S.cached_json(f"sbdb_{num}.json", fetch).get("phys_par", [])}
    except Exception:
        return {}
    return {k: phys[v] for k, v in (("tholen", "spec_T"), ("smass", "spec_B"), ("albedo", "albedo"), ("diameter_km", "diameter")) if phys.get(v)}


def sun_fnu(um):
    """Solar spectrum shape in F_nu (5778 K blackbody), arbitrary units."""
    nu = 2.998e14 / np.asarray(um)
    return nu ** 3 / np.expm1(6.626e-34 * nu / (1.381e-23 * 5778))


def bb_fnu(um, T):
    nu = 2.998e14 / np.asarray(um)
    return nu ** 3 / np.expm1(6.626e-34 * nu / (1.381e-23 * T))


def geometry(name, mjd):
    """(Sun distance au, solar phase angle deg) of a numbered body at mjd, seen from SPHEREx (JPL Horizons)."""
    num = name.split()[0]
    def fetch():
        body, _ = S.get("https://ssd.jpl.nasa.gov/api/horizons.api", {
            "format": "json", "COMMAND": f"'{num};'", "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'OBSERVER'",
            "CENTER": "'500@-163182'", "TLIST": f"'{mjd + 2400000.5}'", "TLIST_TYPE": "'JD'",
            "QUANTITIES": "'19,24'", "CSV_FORMAT": "'YES'", "OBJ_DATA": "'NO'"})
        return json.loads(body)
    try:
        res = S.cached_json(f"geom_{num}_{mjd:.3f}.json", fetch)["result"]
        head = res.split("$$SOE")[0].strip().splitlines()[-2]   # column names line
        cols = [c.strip() for c in head.split(",")]
        row = [c.strip() for c in res.split("$$SOE")[1].split("$$EOE")[0].strip().split(",")]
        return round(float(row[cols.index("r")]), 4), round(float(row[cols.index("S-T-O")]), 2)
    except Exception:
        return None, None


def helio_r(name, mjd):
    return geometry(name, mjd)[0]


def fit_model(pts, r_au):
    """Reflected sunlight with a linear spectral slope + one blackbody for the asteroid's heat.

    F(l) = A * sun(l) * (1 + s * (min(l, 2.5 um) - 1 um)) + B * BB(l, T). A, A*s, B are linear; T is
    grid-searched. Errors include a ROT_MAG spin term, since points come from different visits.
    The free amplitude B absorbs the unknown diameter (D^2/Delta^2) and the hemisphere-integration
    factor of a T_ss*cos(theta)^(1/4) surface, so the fixed-T_ss curve adds no degenerate parameter.
    ponytail: a single-temperature blackbody, not full NEATM; enough to separate reflected light
    from heat, not to measure sizes or albedos."""
    if len(pts) < 6 or max(p["wave"] for p in pts) - min(p["wave"] for p in pts) < 2.5:
        return None
    lam = np.array([p["wave"] for p in pts]); f = np.array([p["flux"] for p in pts])
    sig = np.hypot([p["err"] for p in pts], f * (10 ** (0.4 * ROT_MAG) - 1))
    best = None
    # T from the Sun distance (subsolar temperature, beaming eta=1.2); a free T is degenerate with
    # the heat amplitude when only the last detector constrains it. Grid-search only if r is unknown.
    temps = [394 / np.sqrt(r_au) * 1.2 ** -0.25] if r_au else np.arange(120, 451, 5)
    lc = np.minimum(lam, KNEE) - 1
    for T in temps:
        X = np.stack([sun_fnu(lam), sun_fnu(lam) * lc, bb_fnu(lam, T)], 1) / sig[:, None]
        coef, *_ = np.linalg.lstsq(X, f / sig, rcond=None)
        if coef[0] <= 0 or coef[2] < 0:
            continue
        chi2 = float(np.sum((X @ coef - f / sig) ** 2))
        if best is None or chi2 < best[0]:
            best = (chi2, T, coef)
    if best is None:
        return None
    chi2, T, (a, a_s, b) = best
    grid = np.geomspace(0.75, 5.0, 60)
    refl, heat = a * sun_fnu(grid) * (1 + a_s / a * (np.minimum(grid, KNEE) - 1)), b * bb_fnu(grid, T)
    i48 = int(np.argmin(abs(grid - 4.8)))
    return {"T": round(float(T)), "T_fixed": bool(r_au), "slope_per_um": round(float(a_s / a), 3), "chi2_dof": round(chi2 / max(len(pts) - 4, 1), 2),
            "thermal_frac_4p8": round(float(heat[i48] / (heat[i48] + refl[i48])), 2),
            "T_expected": round(394 / np.sqrt(r_au), 0) if r_au else None,
            "grid": [round(float(g), 3) for g in grid], "refl": [round(float(v), 1) for v in refl],
            "heat": [round(float(v), 1) for v in heat]}


if __name__ == "__main__":
    # the reflectance continuum is continuous at the 2.5 um knee (flat part pinned to the slope's value)
    refl = lambda lam, s: 1 + s * (np.minimum(lam, KNEE) - 1)
    assert abs(refl(KNEE - 1e-9, 0.2) - refl(KNEE + 1e-9, 0.2)) < 1e-9 and refl(4.0, 0.2) == refl(KNEE, 0.2)
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
