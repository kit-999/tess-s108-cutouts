# TESS Sectors 108 and 109: first-look light curves for 30 TIC IDs

Fetched on 2026-10-26, on request, for a community vetting game. **A first look, not a vet.** Sector 108 had no SPOC or QLP
light curves in MAST when these were fetched, so each light curve here is my own simple aperture photometry on a TESSCut
cutout of the full-frame images. 30 of the 30 have one.

## What each star has (in `stars/`)

- `TIC<id>_s108_lc.csv`: time (BTJD = BJD − 2457000), normalized flux, and flux in e⁻/s (the aperture sum minus a background:
  the median of the cutout's border pixels in each frame, times the aperture's pixel count). Frames with a nonzero quality
  flag are dropped.
- `TIC<id>_s108_lc.png`: the light curve with one-hour bins; a Lomb-Scargle amplitude spectrum; and the 20-minute bins folded
  on the box search's best period.
- `TIC<id>_s108_meta.json`: the TIC catalog row used (RA, Dec, Tmag, Teff, contamination ratio), the exact TESSCut request,
  the fetch times in UTC, both hashes (below), and the numbers in the table.

## Hashes

TESSCut stamps each cutout with the time it was made, in the comments of its CHECKSUM and DATASUM cards, so a re-fetched
file's sha256 differs from mine even though its data are identical. `fits_sha256` names the copy I fetched;
`data_sha256_time_flux_quality` survives a re-fetch: sha256 over the TIME, FLUX and QUALITY columns' bytes as stored in the
file (big-endian), in that order.

## Method (`method/batch.py`)

- Position: the MAST TIC catalog. Cutout: TESSCut, 15 × 15 pixels, Sector 108, followed only to MAST's own hosts.
- Aperture: pixels connected to the brightest pixel near the centre whose median flux is above the background by more than
  the largest of 1% of that pixel's, 5 σ of the border pixels, or 50 e⁻/s.
- Periodogram: Lomb-Scargle from 0.1 cycles/day to the Nyquist frequency, as semi-amplitudes in parts per thousand.
- Box search: astropy's BoxLeastSquares on 20-minute bins, periods from 0.5 days to half the time span, durations of 1, 2, 3
  and 4.8 hours. SDE is the best peak's height above the mean of the search's power, in standard deviations.

## The 30, in the manifest's order

"What the curve shows" is my own reading of each plot, with its numbers filled in from the files. Every star used 10,999 of
its 11,083 frames (the other 84 carry a quality flag).

| rank | TIC (plot) | Tmag | what the curve shows | scatter, 1-h bins (ppm) | box search, best | fetched (UTC) | data hash (start) |
|---|---|---|---|---|---|---|---|
| 1 | [589826](stars/TIC589826_s108_lc.png) | 6.02 | pulsates fast: 41.5 and 44.9 cycles a day at 1.4 and 1.3 ppt; no dips stand out | 291 | 5.167 d, 0.87 ppt, SDE 3.8 | 2026-10-26T17:32:52Z | `81d5ce202cc0` |
| 2 | [589906](stars/TIC589906_s108_lc.png) | 8.51 | wanders by about 4 ppt over a day or so; no dips stand out | 999 | 10.485 d, 2.21 ppt, SDE 2.5 | 2026-10-26T17:33:26Z | `e48ee8ea520c` |
| 3 | [590274](stars/TIC590274_s108_lc.png) | 7.47 | a small wave near one day (0.99 cycles a day, 0.36 ppt); no dips stand out | 442 | 1.010 d, 0.45 ppt, SDE 4.7 | 2026-10-26T17:34:05Z | `a40cdff14dca` |
| 4 | [592638](stars/TIC592638_s108_lc.png) | 8.21 | small wiggles, 1.5 ppt top to bottom, with sharp periodogram peaks at 0.65, 7.07, 11.11 cycles a day (all under 0.3 ppt) | 350 | 1.555 d, 0.32 ppt, SDE 4.4 | 2026-10-26T17:34:11Z | `bb85d0cf9b1f` |
| 5 | [593228](stars/TIC593228_s108_lc.png) | 9.22 | the same light as 593230 (Tmag 8.3), 8.8″ away, inside one TESS pixel; the TIC puts 2.3 times this star's own light into its aperture: a swing every 6.5 days (4.3 ppt) and flares | 3,857 | 6.485 d, 4.99 ppt, SDE 2.5 | 2026-10-26T17:34:18Z | `05cf919b03bf` |
| 6 | [593230](stars/TIC593230_s108_lc.png) | 8.27 | the same blended light as 593228 (rank 5) | 3,888 | 6.485 d, 4.99 ppt, SDE 2.5 | 2026-10-26T17:34:25Z | `96ab3b0948ad` |
| 7 | [595279](stars/TIC595279_s108_lc.png) | 9.08 | slow drifts and a step near the second orbit's start; no dips stand out | 431 | 10.177 d, 0.48 ppt, SDE 3.0 | 2026-10-26T17:34:31Z | `c52ae4e24c5c` |
| 8 | [607907](stars/TIC607907_s108_lc.png) | 8.41 | quiet: 1.0 ppt top to bottom in hour bins | 220 | 6.944 d, 0.33 ppt, SDE 3.7 | 2026-10-26T17:34:37Z | `1b4929ee89c8` |
| 9 | [608932](stars/TIC608932_s108_lc.png) | 8.30 | small wiggles, 2.7 ppt top to bottom; no dips stand out | 544 | 13.468 d, 0.98 ppt, SDE 3.1 | 2026-10-26T17:34:43Z | `3434683ecead` |
| 10 | [610232](stars/TIC610232_s108_lc.png) | 8.81 | a smooth swing, strongest every 5.2 days, 17 ppt top to bottom | 4,560 | 10.267 d, 6.49 ppt, SDE 2.4 | 2026-10-26T17:34:50Z | `8ec2f5c043e6` |
| 11 | [612908](stars/TIC612908_s108_lc.png) | 7.08 | bumps at the start of each orbit, where TESS's scattered light usually sits; otherwise quiet | 366 | 7.029 d, 0.71 ppt, SDE 3.3 | 2026-10-26T17:34:57Z | `73f4fa90349b` |
| 12 | [614129](stars/TIC614129_s108_lc.png) | 8.75 | a slow wave, 3.6 ppt top to bottom; no dips stand out | 1,052 | 6.270 d, 1.14 ppt, SDE 2.7 | 2026-10-26T17:35:03Z | `124798970dca` |
| 13 | [614556](stars/TIC614556_s108_lc.png) | 8.04 | bumps at the start of each orbit, where scattered light usually sits; otherwise quiet | 253 | 13.196 d, 0.43 ppt, SDE 2.7 | 2026-10-26T17:35:11Z | `af061de650a7` |
| 14 | [615133](stars/TIC615133_s108_lc.png) | 8.10 | a strong wave every 0.89 days, 35 ppt top to bottom, its height changing through the sector | 10,655 | 0.889 d, 13.38 ppt, SDE 6.0 | 2026-10-26T17:35:18Z | `0758f53ed1ea` |
| 15 | [627296](stars/TIC627296_s108_lc.png) | 8.52 | wiggles, 3.4 ppt top to bottom; no dips stand out | 740 | 6.883 d, 0.99 ppt, SDE 2.6 | 2026-10-26T17:35:24Z | `747e11283bb9` |
| 16 | [664269](stars/TIC664269_s108_lc.png) | 7.98 | a ramp at the sector's start, then quiet | 268 | 13.670 d, 0.61 ppt, SDE 4.0 | 2026-10-26T17:35:30Z | `0971238d91fa` |
| 17 | [672060](stars/TIC672060_s108_lc.png) | 8.91 | quiet: 1.2 ppt top to bottom in hour bins | 239 | 3.418 d, 0.28 ppt, SDE 2.9 | 2026-10-26T17:35:37Z | `c326faaaf02e` |
| 18 | [672624](stars/TIC672624_s108_lc.png) | 8.12 | a wave every 0.55 days, 11 ppt top to bottom; its high box-search score is this smooth wave, not a box-shaped dip | 2,785 | 0.548 d, 4.24 ppt, SDE 8.4 | 2026-10-26T17:35:44Z | `8d8be985ce90` |
| 19 | [675445](stars/TIC675445_s108_lc.png) | 8.55 | a big smooth swing, 155 ppt top to bottom, about two cycles in the sector | 71,645 | 10.970 d, 83.41 ppt, SDE 2.6 | 2026-10-26T17:35:52Z | `55b5b4609b3b` |
| 20 | [682256](stars/TIC682256_s108_lc.png) | 9.79 | faint and quiet: 2.2 ppt top to bottom in hour bins | 490 | 7.033 d, 0.73 ppt, SDE 3.4 | 2026-10-26T17:36:26Z | `61f596ae5966` |
| 21 | [682491](stars/TIC682491_s108_lc.png) | 7.11 | one slow hump across the sector, 5.5 ppt | 1,719 | 10.375 d, 1.70 ppt, SDE 2.6 | 2026-10-26T17:36:34Z | `eb96e9fd8ce6` |
| 22 | [683071](stars/TIC683071_s108_lc.png) | 8.38 | wiggles over hours to days, 5.1 ppt top to bottom | 1,025 | 7.652 d, 1.33 ppt, SDE 2.5 | 2026-10-26T17:36:41Z | `68ca949f9481` |
| 23 | [709603](stars/TIC709603_s108_lc.png) | 8.25 | quiet: 1.0 ppt top to bottom in hour bins | 241 | 12.937 d, 0.33 ppt, SDE 2.9 | 2026-10-26T17:36:49Z | `a2fadac414e2` |
| 24 | [710187](stars/TIC710187_s108_lc.png) | 8.67 | quiet, a slow drift of 2.1 ppt | 466 | 6.615 d, 0.74 ppt, SDE 3.1 | 2026-10-26T17:36:55Z | `9f9f066617fa` |
| 25 | [711970](stars/TIC711970_s108_lc.png) | 8.95 | quiet, a slow drift of 1.9 ppt | 472 | 6.841 d, 0.70 ppt, SDE 3.9 | 2026-10-26T17:37:02Z | `bd60bbb9fc1d` |
| 26 | [737327](stars/TIC737327_s108_lc.png) | 9.15 | a wave every 1.64 days that grows from the first orbit to the second, 45 ppt top to bottom | 14,516 | 1.625 d, 17.90 ppt, SDE 4.7 | 2026-10-26T17:37:10Z | `11c58d4b897d` |
| 27 | [737546](stars/TIC737546_s108_lc.png) | 10.54 | eclipses: a primary 32% deep every 3.069 days (8 seen, about 10 hours each) and a secondary 6.5% deep at half phase | 21,386 | 3.070 d, 266.65 ppt, SDE 5.8 | 2026-10-26T17:37:17Z | `d0428db133c1` |
| 28 | [740069](stars/TIC740069_s108_lc.png) | 8.19 | wiggles over days, 10 ppt top to bottom | 1,778 | 5.236 d, 2.51 ppt, SDE 2.3 | 2026-10-26T17:37:23Z | `e2429a357cf8` |
| 29 | [767526](stars/TIC767526_s108_lc.png) | 5.70 | bright and red (Teff 3884 K): a slow swing of 14 ppt across the sector | 6,458 | 13.678 d, 4.96 ppt, SDE 2.9 | 2026-10-26T17:37:29Z | `a9167d31ba04` |
| 30 | [767943](stars/TIC767943_s108_lc.png) | 8.07 | wiggles over a day or two, 4.7 ppt top to bottom | 1,076 | 5.164 d, 1.49 ppt, SDE 2.2 | 2026-10-26T17:37:36Z | `16cdee00e25c` |

## Sector 109, first half

The same method on Sector 109's first half (2026-10-05 to 2026-10-18, camera 2, CCD 1):
every star used 5,570 of its 5,583 frames. Files: `stars/TIC<id>_s109_*`.

| rank | TIC (plot) | what the curve shows | scatter, 1-h bins (ppm) | fetched (UTC) | data hash (start) |
|---|---|---|---|---|---|
| 1 | [589826](stars/TIC589826_s109_lc.png) | the same fast pulsation, 41.5 and 44.9 cycles a day; no dips stand out | 283 | 2026-10-26T17:43:41Z | `31ffc183bd5a` |
| 2 | [589906](stars/TIC589906_s109_lc.png) | wanders by about 5 ppt, as in 108; no dips stand out | 764 | 2026-10-26T17:43:46Z | `8c35904ed568` |
| 3 | [590274](stars/TIC590274_s109_lc.png) | the same small wave near one day; no dips stand out | 393 | 2026-10-26T17:43:52Z | `6542ce90076b` |
| 4 | [592638](stars/TIC592638_s109_lc.png) | the same small wiggles and sharp periodogram peaks | 305 | 2026-10-26T17:43:57Z | `7074909e3fe4` |
| 5 | [593228](stars/TIC593228_s109_lc.png) | the blended pair again: the slow swing (every 6.4 days), and a flare-shaped jump of 28% at BTJD 4325.44, just before the mid-sector gap, rising in one frame and fading over hours | 4,562 | 2026-10-26T17:44:19Z | `d509333aee8e` |
| 6 | [593230](stars/TIC593230_s109_lc.png) | the same blended light (rank 5) | 4,562 | 2026-10-26T17:44:41Z | `d509333aee8e` |
| 7 | [595279](stars/TIC595279_s109_lc.png) | a slow drift; no dips stand out | 324 | 2026-10-26T17:44:46Z | `6097e1d4d818` |
| 8 | [607907](stars/TIC607907_s109_lc.png) | quiet: 0.8 ppt top to bottom | 178 | 2026-10-26T17:44:51Z | `4b5aaefffb92` |
| 9 | [608932](stars/TIC608932_s109_lc.png) | small wiggles, 1.9 ppt; no dips stand out | 402 | 2026-10-26T17:44:57Z | `a53679703410` |
| 10 | [610232](stars/TIC610232_s109_lc.png) | the same smooth swing, 13 ppt top to bottom | 5,341 | 2026-10-26T17:45:16Z | `337df105e072` |
| 11 | [612908](stars/TIC612908_s109_lc.png) | quiet; no bumps at this half's start | 307 | 2026-10-26T17:45:21Z | `90c003d697d4` |
| 12 | [614129](stars/TIC614129_s109_lc.png) | a dip at the start, then a slow rise, 5.6 ppt | 1,915 | 2026-10-26T17:45:26Z | `3a79f108133b` |
| 13 | [614556](stars/TIC614556_s109_lc.png) | quiet; its box search's best period, 6.56 days, is shared by 4 different stars, so it is the instrument | 177 | 2026-10-26T17:45:45Z | `b758b7530392` |
| 14 | [615133](stars/TIC615133_s109_lc.png) | the same wave, every 0.881 days here (108: 0.888), 29 ppt top to bottom | 11,628 | 2026-10-26T17:45:49Z | `87f739421fcd` |
| 15 | [627296](stars/TIC627296_s109_lc.png) | a hump, then a slow decline, 7.8 ppt | 1,460 | 2026-10-26T17:45:55Z | `c0782349c776` |
| 16 | [664269](stars/TIC664269_s109_lc.png) | a ramp at the start, then quiet, as in 108 | 900 | 2026-10-26T17:46:00Z | `5232e127aa9c` |
| 17 | [672060](stars/TIC672060_s109_lc.png) | quiet; its box search's best period, 6.56 days, is shared by 4 different stars, so it is the instrument | 283 | 2026-10-26T17:46:06Z | `3eeabe23fe62` |
| 18 | [672624](stars/TIC672624_s109_lc.png) | the same wave, every 0.552 days, 14 ppt top to bottom | 3,005 | 2026-10-26T17:46:11Z | `e0505091508b` |
| 19 | [675445](stars/TIC675445_s109_lc.png) | the big slow swing, and a flare-shaped jump of 39% above the level before it at BTJD 4328.63, rising in under an hour and fading over a day or more; no other star shows anything then | 57,267 | 2026-10-26T17:46:17Z | `70f6eea76ccb` |
| 20 | [682256](stars/TIC682256_s109_lc.png) | faint and quiet: 2.6 ppt top to bottom | 731 | 2026-10-26T17:46:23Z | `5b32dca0961d` |
| 21 | [682491](stars/TIC682491_s109_lc.png) | the slow hump again, 2.4 ppt | 754 | 2026-10-26T17:46:30Z | `df24552d984f` |
| 22 | [683071](stars/TIC683071_s109_lc.png) | wiggles, 3.6 ppt top to bottom | 983 | 2026-10-26T17:46:50Z | `89b4d19eebc0` |
| 23 | [709603](stars/TIC709603_s109_lc.png) | quiet; its box search's best period, 6.56 days, is shared by 4 different stars, so it is the instrument | 201 | 2026-10-26T17:46:54Z | `6c714aa5a0c5` |
| 24 | [710187](stars/TIC710187_s109_lc.png) | quiet, a slow drift of 2.1 ppt | 536 | 2026-10-26T17:47:01Z | `cfe498922854` |
| 25 | [711970](stars/TIC711970_s109_lc.png) | a ramp at the start, then a slow rise, 3.2 ppt | 1,321 | 2026-10-26T17:47:06Z | `235ec2f83fa7` |
| 26 | [737327](stars/TIC737327_s109_lc.png) | the same wave, every 1.65 days, shrinking through this half, 32 ppt top to bottom | 9,945 | 2026-10-26T17:47:12Z | `e3fbdcdf6dbc` |
| 27 | [737546](stars/TIC737546_s109_lc.png) | the same eclipses: 5 primaries, each within 33 minutes of the time 108's period predicts, at the same 32% depth, with the secondaries between | 24,139 | 2026-10-26T17:47:21Z | `7425f753d0dd` |
| 28 | [740069](stars/TIC740069_s109_lc.png) | wiggles, 7 ppt, with a deeper dip of 5.4 ppt near BTJD 4324.3 (about a day wide) | 1,413 | 2026-10-26T17:47:26Z | `840a85ea79ec` |
| 29 | [767526](stars/TIC767526_s109_lc.png) | the slow swing, 14 ppt; its box search's best period, 6.56 days, is shared by 4 different stars, so it is the instrument | 4,636 | 2026-10-26T17:47:30Z | `05a9fce0e487` |
| 30 | [767943](stars/TIC767943_s109_lc.png) | wiggles, 3.6 ppt top to bottom | 810 | 2026-10-26T17:47:36Z | `d373e338ab4f` |

- A box-search period shared by different stars belongs to the instrument: four stars here land on the same 6.56 days.
- The two flare-shaped jumps (ranks 5 and 6 at BTJD 4325.44, rank 19 at 4328.63) happen at different times, and no other star
  of the 30 shows anything at either moment.

## Read with care

- A box search always returns a best period, and its SDE is not a verdict either way: a smooth wave can score high (rank 18),
  and clear eclipses can score low when their aliases fill the search (rank 27, eclipses 32% deep, SDE 5.8). Look at the plot.
- A strongest sinusoid at 0.1 to 0.2 cycles a day, the bottom of the search, is usually a slow drift across the sector, not a
  period.
- Bumps at the start of each orbit (about BTJD 4291 and 4305) are where TESS's scattered light usually sits.
- Ranks 5 and 6 (TIC 593228 and 593230) are 8.8″ apart, inside one TESS pixel: their two light curves are one blended light.
- Bright stars saturate and bleed along columns; `aperture_touches_edge` in a star's meta file says when the 15 × 15 cutout
  may have cut some flux off.
- Neighbours inside the aperture dilute a dip and can make one: `contratio` is the TIC catalog's contamination ratio.

TESS data are public. The files derived here are dedicated to the public domain (CC0).
