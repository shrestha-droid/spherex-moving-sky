# How it works

Every team member should be able to explain this page on a whiteboard. Judges will ask.

## 1. The data
- **SPHEREx** is a NASA infrared telescope in a polar orbit, launched in March 2025. It maps the whole sky every **6 months** in 102 colours (0.75–5 µm).
- Each exposure is a 2040×2040-pixel image covering **3.5°×3.5°**, with **6.15″ per pixel**. There are 6 detectors, D1–D6, one per wavelength range. We use **D1** (0.75–1.1 µm). It sees reflected sunlight, which is what asteroids give off.
- The telescope has **no filter wheel**. A *linear variable filter* means each row of the detector sees a slightly different wavelength. The same star therefore looks a bit brighter or fainter depending on where it lands. That is why a naive "subtract image A from image B" fails, and why we don't do it.
- Images come from NASA/IPAC **IRSA** (collection `spherex_qr2`). We never download whole 69 MB files. We send **HTTP range requests** and read only the pixel rows we need, a few KB per asteroid.

## 2. Each image file has 4 useful layers
| Layer | What it is |
|---|---|
| IMAGE | brightness (MJy/sr) |
| VARIANCE | the noise per pixel, so we can compute signal-to-noise (SNR) |
| FLAGS | 32 bit-flags per pixel from the SPHEREx pipeline |
| WCS | the mapping from pixel to sky coordinates (RA/Dec) |

Two flags matter most:
- **Bit 19, `MP_OUTLIER`**: the pipeline compared overlapping exposures and this pixel didn't agree. **Moving objects trigger it.** We confirmed this on real known asteroids.
- **Bit 21, `MP_SOURCE`**: the pixel belongs to a catalogued star or galaxy. A real asteroid sitting on one of these is *blended*.

Other bits mark cosmic rays, hot pixels and similar defects. We treat those pixels as bad.

## 3. Step one: prove it works on known asteroids (`recover.py`)
1. Find every SPHEREx visit to a field. A visit is a burst of 4 frames; visits are hours apart.
2. Ask **JPL Horizons** where SPHEREx was at the moment of each exposure. Its position matters: SPHEREx is about 650 km up, and that parallax shifts nearby asteroids.
3. Ask **JPL's Small-Body Identification API** (`sb_ident`) which known asteroids should be in that image, *as seen from SPHEREx*.
4. At each predicted pixel, measure the light inside a 2-pixel circle (*aperture photometry*) minus the local background. Compute SNR from the VARIANCE layer.
5. Look at the **same patch of sky in another visit**, hours later. The asteroid should have moved on.
   - Bright then, empty later → **recovered**
   - Still bright later → a star is in the way → **blended**
   - Nothing there → too faint → **not detected**

Grouped by brightness (V magnitude), this gives the **completeness curve**: how faint SPHEREx can see in one exposure.

## 4. Step two: find movers with no hints (`detect.py`)
1. **Seeds.** Within a sky region, find bright local peaks with SNR ≥ 5 near a `MP_OUTLIER` pixel, not on a catalogued source, and with no bad pixels.
2. **Linking.** For every pair of seeds from two different visits, compute the speed. Keep it if it is between 2 and 150″/hour, which is the asteroid range. Then predict where that object would be in the other visits. If seeds are waiting there (within 12″), it's a **track**. We require **3 or more visits**. Noise and cosmic rays don't line up in a straight line three times.
3. **Fit.** Fit straight-line motion to the track. The residual (RMS, in arcseconds) shows how clean the track is.
4. **Identify.** Compare each track with JPL's predictions. A match makes it a **known** asteroid. If nothing matches, it is a **candidate**.
5. **Deep check.** Each unmatched track gets a second JPL query with no brightness limit and comets included.
6. **False-alarm test.** Rotate each visit's seeds by a random angle around the field center, 30 times, and count how many tracks still link. Rotation keeps the density of noise peaks but destroys real motion, so whatever links is pure chance. If we find N candidates and chance alone predicts about N, the candidates are not significant. We say so openly.
7. **Score.** Of the bright known asteroids that crossed the region, how many did we find blind?

## 5. What a candidate is, and is not
A candidate is a moving source with no match in JPL's catalogue at our search limits. It could be:
- an asteroid fainter than our search limit,
- a comet,
- an object with an uncertain orbit,
- a coincidence of noise, or
- something genuinely unknown.

We call it a candidate and show the evidence (stamps, rate, fit residual). **We never claim a discovery.** Real discoveries go through the Minor Planet Center with follow-up observations.

## 6. Known limits
- **Slow movers.** Objects slower than 2″/h are rejected, because slowly drifting stars can fake them. Very distant objects such as TNOs or a "Planet X" move slower than that. Catching them needs SPHEREx's 6-month revisits, not a few days. That is the natural next step and a good answer when judges ask "why Planet X?".
- **Straight lines.** Over 2–3 days asteroid paths are almost straight; longer gaps need real orbit fitting.
- **Regions.** We search small regions (about 1°) to keep downloads small. The pipeline itself has no size limit.

## Glossary
- **V magnitude**: brightness on a backwards scale. Bigger means fainter; each +5 is 100× fainter.
- **″ (arcsecond)**: 1/3600 of a degree. The full Moon is about 1800″ across.
- **SNR**: signal-to-noise ratio. 5 means the signal is 5× the noise, a standard detection threshold.
- **MJD**: Modified Julian Date, a day count used in astronomy.
