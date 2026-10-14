# Pitch & demo script

## One-liner
**"We built a tool that watches NASA's newest all-sky telescope for things that move, and proves it works by finding known asteroids blind."**

## 30-second video (shot list)
| Time | Visual | Voice-over |
|---|---|---|
| 0–5s | Blink viewer playing: stars still, one dot hopping | "This is the same patch of sky, photographed by NASA's SPHEREx telescope hours apart. One of these dots is moving." |
| 5–12s | Track overlay appears, then the JPL label | "Our pipeline found it with no hints. Then it checked JPL's catalogue: it's asteroid 872 Holda." |
| 12–20s | Completeness chart, then the "Here, then gone" stamps | "We tested it on every known asteroid in the field. Every one brighter than magnitude 17 was recovered." |
| 20–27s | Scoreboard tile: candidates vs. expected by chance | "And it can't fool itself: we measure how many fake tracks noise alone would make, so a real discovery has to beat the odds." |
| 24–27s | Spectrum card: 872 Holda, its curve rising at 4.8 µm | "And because SPHEREx sees every colour, each asteroid leaves a spectrum. This one is glowing with its own heat." |
| 27–30s | Title card | "SPHEREx Moving Sky. Anyone can watch the sky change." |

## Judge Q&A (practise these)
- **"Why not just subtract images?"** SPHEREx has no filter wheel: each detector row sees a different wavelength, so stars change brightness between exposures. We use the pipeline's outlier flags plus motion linking across 3 or more visits instead.
- **"How do you know it works?"** Known-asteroid recovery. We predict positions with JPL's sb_ident *from SPHEREx's own orbit position* (from JPL Horizons), then measure. Bright asteroids are recovered at 100%.
- **"Did you discover anything?"** Not yet, and we can prove it. Our false-alarm test scrambles the detections to measure how many tracks noise alone produces. In our first field, chance predicted 1.0 and we found 1, so that candidate is noise. Any real candidate has to beat the chance rate, then get follow-up and a Minor Planet Center report.
- **"What about Planet X?"** Distant objects move under 2″/h, too slowly for a few days of data. SPHEREx revisits every 6 months, and linking across those revisits is our next step.
- **"Is the spectrum real?"** Yes. Every point is a real exposure, and its wavelength comes from SPHEREx's own per-pixel calibration table. Points come from the detector pairs as the asteroid moves across them. The shape is indicative, because we don't correct for how the blur grows with wavelength and the asteroid spins between visits.
- **"How much data do you download?"** A few KB per asteroid. We read only the pixels we need from the full images using HTTP range requests.
- **"Did you use AI?"** Be honest and follow the Space Apps disclosure rules. Say what you built, what the AI built, and what you verified.

## Space Apps checklist
- [ ] Check the 2026 rules on AI-tool disclosure and add the required statement to the project page
- [ ] Project page: summary, link to the repo, link to the demo, 30-second video
- [ ] Credit the data: SPHEREx QR2, DOI 10.26131/IRSA652; JPL Horizons; JPL SBDB/sb_ident
- [ ] Deploy `web/` as a static site (any static host) so judges can click it
