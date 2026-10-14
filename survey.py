#!/usr/bin/env python3
"""Survey: run recovery, blind search and spectra over many ecliptic fields and sky passes.

Each field-pass goes to web/fields/<id>/. A finished field writes done.json, so a re-run
resumes where it stopped. web/survey.json indexes everything with combined statistics.

usage: python3 survey.py            # all FIELDS
       python3 survey.py 3          # only the first 3
"""
import json, os, sys, time
import numpy as np
import detect, recover, spectra

# (ra, dec, pass) along the ecliptic around RA 180, where SPHEREx revisits every ~6 months
FIELDS = [(180, 0.0, 0), (180, 0.0, 1), (180, 0.0, 2), (170, 4.3, 0), (190, -4.3, 0),
          (170, 4.3, 1), (190, -4.3, 1), (160, 8.4, 0), (200, -8.4, 0), (210, -12.2, 0)]
BOX = 0.6          # degrees; blind-search region
BRIGHT = 17.5      # V limit for "should be detectable" (single-exposure limit is ~17.5)
SPEC_V = 16.5      # V limit for spectra
MAX_SPEC = 5


def fid(ra, dec, p):
    return f"ra{ra:03.0f}_dec{dec:+05.1f}_p{p}"


def paths(vs, results):
    """name -> [(mjd, ra, dec), ...] for asteroids from the recovery test, plus their V."""
    by_id = {v.id: v for v in vs}
    out, mag = {}, {}
    for r in results:
        ra, dec = by_id[r["visit"]].wcs.pixel_to_world_values(r["x"], r["y"])
        out.setdefault(r["name"], []).append((r["mjd"], float(ra), float(dec)))
        mag[r["name"]] = r["vmag"]
    return out, mag


def choose_box(path, mag, ra0, dec0, box=BOX):
    """Center the search box where the most bright multi-visit asteroids fit entirely."""
    half = box / 2 - 0.03
    bright = [p for n, p in path.items() if mag[n] < BRIGHT and len(p) >= 3]
    best, center = 0, (ra0, dec0)
    for p in bright:
        cra, cdec = np.mean([q[1] for q in p]), np.mean([q[2] for q in p])
        cos = np.cos(np.radians(cdec))
        k = sum(all(abs((q[1] - cra) * cos) < half and abs(q[2] - cdec) < half for q in other) for other in bright)
        if k > best:
            best, center = k, (round(float(cra), 4), round(float(cdec), 4))
    return center, best


def run_field(ra, dec, p):
    out = f"web/fields/{fid(ra, dec, p)}"
    if os.path.exists(f"{out}/done.json"):
        return json.load(open(f"{out}/done.json"))
    os.makedirs(out, exist_ok=True)
    t = time.time()
    print(f"\n=== {fid(ra, dec, p)}", flush=True)
    vs, results = recover.run(ra, dec, 1, 2.6, p)
    if len(vs) < 3:
        print("  fewer than 3 visits, skipped", flush=True)
        info = {"id": fid(ra, dec, p), "skipped": True}
        json.dump(info, open(f"{out}/done.json", "w"))
        return info
    summary = recover.summarize(results)
    path, mag = paths(vs, results)
    (cra, cdec), nbright = choose_box(path, mag, ra, dec)
    print(f"  search box at {cra}, {cdec} covering {nbright} bright asteroid paths", flush=True)
    det = detect.run(cra, cdec, BOX, 1, 2.6, pass_index=p, out_dir=out)
    targets = sorted((n for n in path if mag[n] < SPEC_V and len(path[n]) >= 2), key=lambda n: mag[n])[:MAX_SPEC]
    specs = []
    for n in targets:
        sp = spectra.spectrum(n, path[n])
        sp["vmag"] = mag[n]
        specs.append(sp)
        print(f"  spectrum {n}: {len(sp['points'])} points from {sp['n_images']} exposures", flush=True)
    recover.slim(results)
    json.dump({"field": {"ra": ra, "dec": dec, "pass": p, "detector": 1, "stamp": 2 * recover.R + 1},
               "visits": [{"id": v.id, "mjd": v.mjd} for v in vs], "summary": summary, "results": results},
              open(f"{out}/recovery.json", "w"))
    json.dump(det, open(f"{out}/detections.json", "w"))
    json.dump(specs, open(f"{out}/spectra.json", "w"))
    sc = det["score"]
    info = {"id": fid(ra, dec, p), "ra": ra, "dec": dec, "pass": p, "box": [cra, cdec, BOX],
            "start_mjd": min(v.mjd for v in vs), "visits": len(vs), "summary": summary,
            "predicted": len(results), "tracks": len(det["tracks"]),
            "found": sc["found"], "expected": sc["expected"], "candidates": sc["candidates"],
            "false_alarms": sc["false_alarms"], "spectra": [s["name"] for s in specs if s["points"]],
            "minutes": round((time.time() - t) / 60, 1)}
    json.dump(info, open(f"{out}/done.json", "w"))
    return info


def combine(infos):
    done = [i for i in infos if not i.get("skipped")]
    bins = {}
    for i in done:
        for b in i["summary"]:
            c = bins.setdefault(b["label"], {"label": b["label"], "n": 0, "recovered": 0})
            c["n"] += b["n"]
            c["recovered"] += b["recovered"]
    exp = sum(len(i["expected"]) for i in done)
    hit = sum(len(set(i["found"]) & set(i["expected"])) for i in done)
    return {"fields": len(done), "area_deg2": round(len(done) * BOX ** 2, 2),
            "predicted": sum(i["predicted"] for i in done), "completeness": list(bins.values()),
            "blind_expected": exp, "blind_found": hit,
            "known_tracks": sum(len(i["found"]) for i in done),
            "candidates": sum(i["candidates"] for i in done),
            "false_alarms": round(sum(i["false_alarms"] for i in done), 2),
            "spectra": sum(len(i["spectra"]) for i in done)}


if __name__ == "__main__":
    c, k = choose_box({"a": [(0, 10, 0), (1, 10.1, 0), (2, 10.2, 0)], "b": [(0, 50, 0)] * 3}, {"a": 15, "b": 15}, 0, 0)
    assert k == 1 and abs(c[0] - 10.1) < 1e-6
    fields = FIELDS[:int(sys.argv[1])] if len(sys.argv) > 1 else FIELDS
    infos = []
    for ra, dec, p in fields:
        try:
            infos.append(run_field(ra, dec, p))
        except Exception as e:  # one bad field must not sink a multi-hour run
            print(f"  FAILED {fid(ra, dec, p)}: {e!r}", flush=True)
        json.dump({"combined": combine(infos), "fields": infos}, open("web/survey.json", "w"))
    print("\n", json.dumps(combine(infos), indent=1))
