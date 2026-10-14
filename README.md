# SPHEREx Moving Sky

NASA Space Apps Challenge 2026: [Planet X and SPHEREx](https://www.spaceappschallenge.org/2026/challenges/planet-x-and-spherex/).

This public web tool blinks the SPHEREx sky across visits and finds moving objects automatically. It then proves the detector works by recovering known asteroids from JPL's catalogue.

- How it works (read this before judging): [HOW_IT_WORKS.md](HOW_IT_WORKS.md)
- Pitch, video script and judge Q&A: [PITCH.md](PITCH.md)

## Run
```bash
python3 -m venv .venv && .venv/bin/pip install numpy astropy
.venv/bin/python recover.py 180 0 1 2.6          # known-asteroid recovery -> web/recovery.json
.venv/bin/python detect.py 179.66 0.80 0.7 1 2.6 # blind search -> web/detections.json + web/img/
python3 -m http.server 8000 -d web              # open http://localhost:8000
```
The arguments are RA and Dec in degrees, the box size in degrees (for `detect.py` only), the detector (1–6), and the number of days of visits.

Everything we download is cached in `data/`, which git ignores. A re-run makes no network calls.

## Files
| File | What it does |
|---|---|
| `spherex.py` | Data access: IRSA image search, SPHEREx position from JPL Horizons, known asteroids from JPL sb_ident, and `RemoteFits`, which reads pixel rectangles from remote FITS files with HTTP range requests |
| `recover.py` | Measures each JPL-predicted asteroid at its predicted pixel and confirms it is gone in a later visit, giving a completeness curve |
| `detect.py` | Blind search: outlier-flagged seeds, straight-line linking across 3 or more visits, then matching to JPL (known or candidate). Also writes the blink frames |
| `web/index.html` | The viewer: blink player, track overlays, recovery chart and stamp gallery. No dependencies |

Each script runs its own small self-checks (asserts) at the start.

## Data credits
SPHEREx Quick Release 2, NASA/IPAC Infrared Science Archive, [DOI 10.26131/IRSA652](https://doi.org/10.26131/IRSA652). JPL Horizons. JPL Small-Body Identification API.
