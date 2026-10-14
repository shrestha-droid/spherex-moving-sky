#!/usr/bin/env python3
"""Planet X mode: can SPHEREx find a slow, distant world blind?

For a known dwarf planet (Makemake, Haumea: ~50 au, V~17, ~1"/h), find the SPHEREx
visits that caught it, run the slow-mover search with no hints, then check the result
against JPL. Synthetic movers planted in the same real images (with SPHEREx's own PSF)
measure how faint a Planet-X-like object could be and still be found.

usage: python3 distant.py            # Makemake and Haumea -> web/planetx/, web/planetx.json
"""
import json, os, re, sys
import numpy as np
import spherex as S
import detect
from recover import Visit, measure, BAD, R

TARGETS = [("136472", "Makemake"), ("136108", "Haumea")]
WINDOW_DAYS = 10
BOX = 0.35
INJECT = dict(mags=[14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5], per_mag=8)


def ephem(cmd, mjds):
    """(ra, dec, V, r_au) of a body seen from SPHEREx at each MJD, from JPL Horizons."""
    out = []
    for i in range(0, len(mjds), 40):
        chunk = mjds[i:i + 40]
        def fetch(chunk=chunk):
            body, _ = S.get("https://ssd.jpl.nasa.gov/api/horizons.api", {
                "format": "json", "COMMAND": f"'{cmd};'", "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'OBSERVER'",
                "CENTER": "'500@-163182'", "TLIST": " ".join(f"'{m + 2400000.5:.6f}'" for m in chunk),
                "TLIST_TYPE": "'JD'", "QUANTITIES": "'1,9,19'", "ANG_FORMAT": "'DEG'", "CSV_FORMAT": "'YES'",
                "OBJ_DATA": "'NO'"})
            return json.loads(body)
        res = S.cached_json(f"eph_{cmd}_{chunk[0]:.5f}_{len(chunk)}.json", fetch)["result"]
        for line in res.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
            p = [x.strip() for x in line.split(",")]
            out.append((float(p[3]), float(p[4]), float(p[5]), float(p[7])))
    return out


def inside(region, ra, dec):
    """Point in an IRSA s_region 'POLYGON ICRS ...' (small-field tangent approximation)."""
    v = [float(x) for x in re.findall(r"-?\d+\.?\d*", region)]
    xs, ys = (np.array(v[0::2]) - ra) * np.cos(np.radians(dec)), np.array(v[1::2]) - dec
    xs = (xs + 180) % 360 - 180
    hit = False
    for i in range(len(xs)):
        j = i - 1
        if (ys[i] > 0) != (ys[j] > 0) and 0 < (xs[j] - xs[i]) * (0 - ys[i]) / (ys[j] - ys[i]) + xs[i]:
            hit = not hit
    return hit


def coverage(cmd, det=1):
    """Visits (first frame of each exposure burst) whose image contains the body, with its positions."""
    mean = ephem(cmd, [61000.0])[0]
    rows = [r for r in S.images(mean[0], mean[1], 2.0) if r["energy_bandpassname"] == f"SPHEREx-D{det}"]
    rows.sort(key=lambda r: float(r["t_min"]))
    firsts = []
    for r in rows:
        if not firsts or float(r["t_min"]) - float(firsts[-1]["t_min"]) > 0.05:
            firsts.append(r)
    eph = ephem(cmd, [(float(r["t_min"]) + float(r["t_max"])) / 2 for r in firsts])
    return [(r, e) for r, e in zip(firsts, eph) if inside(r["s_region"], e[0], e[1])]


def best_window(cov):
    """The run of visits within WINDOW_DAYS that has the most visits."""
    best = []
    for i in range(len(cov)):
        w = [c for c in cov[i:] if float(c[0]["t_min"]) - float(cov[i][0]["t_min"]) <= WINDOW_DAYS]
        if len(w) > len(best):
            best = w
    return best


def forced(cov_window):
    """Measure the body at its predicted pixel in each visit: is it there at all?"""
    out = []
    for r, (ra, dec, V, _) in cov_window:
        v = Visit(r)
        p = v.pixel(ra, dec)
        if not p:
            continue
        img, var, dx, dy = v.stamp(*p)
        m = measure(img, var, dx, dy)
        out.append({"visit": v.id, "mjd": v.mjd, "V": V, "snr": round(m["snr"], 1) if m else None})
    return out


def run_target(cmd, name):
    print(f"\n=== {name}", flush=True)
    cov = coverage(cmd)
    win = best_window(cov)
    print(f"  {len(cov)} visits contain {name}; best {WINDOW_DAYS}-day window has {len(win)}", flush=True)
    if len(win) < detect.SLOW["min_visits"]:
        return {"name": name, "visits": len(win), "skipped": True}
    ras, decs = [w[1][0] for w in win], [w[1][1] for w in win]
    ra0, dec0 = round(float(np.mean(ras)), 4), round(float(np.mean(decs)), 4)
    eph = {w[0]["obs_id"]: w[1] for w in win}
    known_fn = lambda v: [{"name": name, "ra": eph[v.id][0], "dec": eph[v.id][1], "vmag": eph[v.id][2]}] if v.id in eph else []
    out_dir = f"web/planetx/{name.lower()}"
    det = detect.run(ra0, dec0, BOX, 1, WINDOW_DAYS, out_dir=out_dir, rows=[w[0] for w in win], slow=True,
                     injection=INJECT, known_fn=known_fn)
    det["forced"] = forced(win)
    det["target"] = {"name": name, "id": cmd, "r_au": round(float(np.mean([w[1][3] for w in win])), 2),
                     "V": round(float(np.mean([w[1][2] for w in win])), 2),
                     "motion_arcsec_h": round(float(np.hypot(*detect.tangent(ras[-1], decs[-1], ras[0], decs[0])) /
                                              ((float(win[-1][0]["t_min"]) - float(win[0][0]["t_min"])) * 24)), 2)}
    json.dump(det, open(f"{out_dir}/detections.json", "w"))
    found = name in det["score"]["found"]
    print(f"  blind search found {name}: {found}; forced SNR per visit: {[f['snr'] for f in det['forced']]}", flush=True)
    return {"name": name, "id": cmd, "dir": f"planetx/{name.lower()}/", "found": found, "visits": len(win),
            "target": det["target"], "tracks": len(det["tracks"]), "candidates": det["score"]["candidates"],
            "false_alarms": det["score"]["false_alarms"], "artifacts": det["score"]["artifacts"],
            "injection": det.get("injection"), "forced": det["forced"]}


if __name__ == "__main__":
    assert inside("POLYGON ICRS 10 10 12 10 12 12 10 12", 11, 11) and not inside("POLYGON ICRS 10 10 12 10 12 12 10 12", 13, 11)
    results = [run_target(c, n) for c, n in (TARGETS if len(sys.argv) < 2 else [t for t in TARGETS if t[1] in sys.argv[1:]])]
    os.makedirs("web", exist_ok=True)
    json.dump(results, open("web/planetx.json", "w"))
    print(json.dumps([{k: r.get(k) for k in ("name", "found", "visits", "candidates", "false_alarms", "injection")} for r in results], indent=1))
