# Planet X and SPHEREx

NASA Space Apps Challenge 2026: [Planet X and SPHEREx](https://www.spaceappschallenge.org/2026/challenges/planet-x-and-spherex/).

The goal is a public web tool that shows how the SPHEREx sky changes over time. It finds moving objects and proves the detector works by recovering known asteroids.

## Data sources
| What | Where |
|---|---|
| SPHEREx QR2 spectral images (2040² px, 6.15″/px, 3.5° field, detectors D1–D6 cover 0.75–5 µm) | [IRSA SIA](https://irsa.ipac.caltech.edu/SIA) collection `spherex_qr2`. Cutouts come from `access_url?center=RA,DECdeg&size=Ndeg` |
| SPHEREx position (Horizons ID `-163182`) | [JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html) |
| Known asteroids in a field, as seen from SPHEREx | [JPL sb_ident API](https://ssd-api.jpl.nasa.gov/doc/sb_ident.html) with `xobs` |

Please cite SPHEREx QR2 as DOI [10.26131/IRSA652](https://doi.org/10.26131/IRSA652).

## Spike (metadata only, no downloads)
```bash
python3 spike.py 180 0 1   # RA, Dec, detector
```
For the ecliptic field at RA 180°, Dec 0°, the spike found 55 D1 images spread over 3 sky passes. JPL predicts 21 known asteroids brighter than V=19 in the first image. Exposures come in bursts of 4 frames about 2 minutes apart. The same field is seen again hours later and again 6 months later.

## Plan
1. **Recovery test.** Take cutouts at JPL-predicted positions in two images about 10 hours apart, then check whether the asteroid is present in one and absent in the other. Report recovered / predicted.
2. **Blind detector.** Difference images from the same epoch, with the wavelength mismatch handled. SPHEREx's linear variable filter means each pixel sees a different wavelength. Then link detections into tracks.
3. **Score.** Cross-match tracks against sb_ident. Unmatched tracks become *candidates*, never claimed discoveries.
4. **Web tool.** A blink and slider view of the same sky across epochs, with known asteroids and candidates overlaid and a recovery-rate scoreboard.
