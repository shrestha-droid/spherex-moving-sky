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

## 4b. Planet X mode: slow movers (`detect.py` slow mode, `distant.py`)
- A body at 50–500 AU moves only about **0.3–3″/h**. Over a few hours that's less than a pixel, so SPHEREx's outlier flags don't fire. Over a week, though, it moves **8–55 pixels**, which is enough to link.
- **Seeds:** any uncatalogued peak (SNR ≥ 5) that is **missing from the median of every other visit** at the same sky spot. A faint star is present in every visit, so it drops out; a slow mover is only there briefly.
- **Linking:** straight lines at 0.3–3″/h through **4 or more visits** within about 10 days, followed by the same scramble test and same-pixel artifact filter as fast mode.
- **Why not across 6 months?** In half a year, Earth's own orbit around the Sun shifts a distant object's apparent position by degrees. Linking across passes needs a full orbit search. Within one pass, the motion is a straight line.
- **The real test:** the dwarf planets **Makemake** (52.7 AU, V 17.1, 1.4″/h) and **Haumea** (49.9 AU, V 17.3, 1.5″/h). Both were **found blind**: 0 unidentified tracks, against 0.8 and 0.03 expected by chance.
- **Depth:** we plant synthetic slow movers into the same real images, one brightness at a time. Each is drawn with **SPHEREx's own PSF** for the detector zone it lands in. We then search with the identical pipeline. The 50% recovery limit is about **AB 17.7**, roughly V 18. Synthetic sources never enter the real results.

## 4c. Defences against fake tracks
- **Detector-fixed defects.** A warm pixel or persistence spot stays at the same detector (x, y). Because SPHEREx points somewhere new each visit, that fixed pixel shows up at a different place on the sky each time and could fake a track. The scramble test rotates detections on the *sky*, so it can't catch this. So we **reject any track that uses the same detector pixel (within 2 px) in two visits**. Rejected tracks are shown in the ledger, not hidden.
- **Robust noise.** The sky noise comes from the median absolute deviation of the background ring. A plain variance let a single bright star in the ring inflate the noise 50×. That bug once hid Makemake in 4 of its 8 visits; with the robust estimate it measures SNR 5–12 in 7 of 8.

## 5. Step three: a spectrum for every bright asteroid (`spectra.py`)
- SPHEREx's six detectors work in **pairs**: D1+D4, D2+D5 and D3+D6. Each pair looks at its own patch of sky, so at any instant an asteroid is in exactly two detectors: one short-wavelength, one long-wavelength. It is *not* imaged in all six bands at once.
- Each visit is a burst of **4 slightly shifted frames**. Because every detector row sees a different wavelength, those shifts put the asteroid on about **8 different wavelengths in about 8 minutes**. Later visits, which land on other detectors, add more.
- JPL's predicted track places the asteroid in every exposure it crosses, in any detector. We measure its brightness there, using a slightly wider aperture so the longer-wavelength blur fits. The file's per-pixel `WCS-WAVE` table gives the exact wavelength.
- Converting to µJy: pixel values are MJy/sr. Multiply by the pixel's solid angle, (6.15″)² in steradians, then by 10¹².
- **Reflectance** is the asteroid's flux divided by the Sun's spectrum, normalised near 1 µm. Flat or reddish means reflected sunlight. A steep rise past about 3.5 µm is **heat** from the asteroid's surface.
- Example: 872 Holda gave 22 points from 0.87 to 4.85 µm, covering all six detectors, with a clear thermal rise at 4.8 µm.
- **PSF correction:** every file stores SPHEREx's real PSF for 121 detector zones. A 3-pixel aperture captures more than **99.6% of the light at every wavelength from 0.9 to 4.7 µm**, because the 6.15″ pixels dominate the blur. We still apply each point's own correction.
- **Spin:** points from different visits carry an extra ±0.15 mag uncertainty in the fit, and the chart colours points by visit so any offset between visits is visible.
- **Model:** reflected sunlight with a linear spectral slope, plus one blackbody for the asteroid's heat. The temperature is *fixed* from its Sun distance (JPL Horizons; subsolar temperature with the standard beaming factor η = 1.2). A free temperature is degenerate with the heat amplitude when only the last detector constrains it.
- **Validation against JPL SBDB, not used in the fit:** across 9 asteroids our slope rises with published albedo (Spearman ρ = 0.68, permutation p ≈ 0.05). Dark, primitive asteroids come out flat; bright S-types come out red.

## 6. Step four: a survey, not one field (`survey.py`)
We run ten field-passes along the ecliptic, the band where most asteroids are. Some are the same patch of sky six months apart, which means a completely different set of asteroids. Each search box is centred automatically where bright known asteroids cross. The combined numbers on the site (completeness, blind finds, candidates against chance) come from all fields together.

## 6b. Live mode (`web/api/live.py`)
The website is not only a recording. The **Query NASA, live** panel calls a Vercel serverless function, standard-library Python only, which:
1. asks IRSA's SIA service which SPHEREx exposures cover a position;
2. reads the FITS headers and **only the requested pixel rows** of the 69 MB file with HTTP range requests (for example, about 40–50 KB for a 49×49 cutout, in about 1 s);
3. asks JPL which known asteroids are in that cutout as seen from SPHEREx at that moment, or where any named body is.

Every request's byte range, size and time is shown in the on-page network log. The function only contacts IRSA and JPL and only accepts SPHEREx file URLs (other URLs return HTTP 400), so it can't be used as an open proxy. Sky-to-pixel uses the header's TAN-SIP inverse distortion, written from scratch, and agrees with astropy to 0.003 px.

## 7. What a candidate is, and is not
A candidate is a moving source with no match in JPL's catalogue at our search limits. It could be:
- an asteroid fainter than our search limit,
- a comet,
- an object with an uncertain orbit,
- a coincidence of noise, or
- something genuinely unknown.

We call it a candidate and show the evidence (stamps, rate, fit residual). **We never claim a discovery.** Real discoveries go through the Minor Planet Center with follow-up observations.

## 8. Known limits
- **Slow movers below 0.3″/h** (beyond about 500 AU) need more than a week; across 6-month passes they need an orbit search. Planet X mode covers 0.3–3″/h.
- **Straight lines.** Over 2–3 days asteroid paths are almost straight; longer gaps need real orbit fitting.
- **Regions.** We search small regions (about 1°) to keep downloads small. The pipeline itself has no size limit.

## Glossary
- **V magnitude**: brightness on a backwards scale. Bigger means fainter; each +5 is 100× fainter.
- **″ (arcsecond)**: 1/3600 of a degree. The full Moon is about 1800″ across.
- **SNR**: signal-to-noise ratio. 5 means the signal is 5× the noise, a standard detection threshold.
- **MJD**: Modified Julian Date, a day count used in astronomy.
