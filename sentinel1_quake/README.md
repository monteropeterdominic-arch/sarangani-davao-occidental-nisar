# Sentinel-1 radar and Sentinel-2 coast check over the 8 June 2026 M7.8 earthquake

Before-and-after radar interferometry (InSAR) and a tide-matched optical check of the coast, for Sarangani and Davao Occidental. Everything is on the team grid: EPSG:32651, 30 m, 5147 × 4810 pixels, origin 648150 E / 735870 N.

**In short:**
- This data **can't show the reported coastal uplift**. The coast nearest the epicentre has no reliable measurements in any pair spanning the quake.
- Inland, between General Santos and Maasim, the pairs from before the quake change by about as much as the pairs spanning it. What's visible there is mostly water vapour in the atmosphere.
- The one quake-linked pattern is a **loss of radar coherence along the Glan coast** that shows up in both viewing directions.
- The tide-matched Sentinel-2 images show **no visible uplift** of reef flats on Sarangani Island or at the Glan peninsula tip.

## Files

| Folder / file | What it holds |
|---|---|
| `geotiff/` | 18 Cloud-Optimized GeoTIFFs (Int16, nodata −32768, DEFLATE, 4 internal overviews). Multiply by the scale factor stored in each file; QGIS and GDAL apply it automatically. |
| `figures/1_fringes_all_pairs.png` | Interferogram colour cycles (fringes) for all 7 pairs, 400 m cells |
| `figures/2_los_change.png` | Line-of-sight change in mm for the 4 exported pairs, reliable areas only |
| `figures/3_coherence.png` | Coherence for the 4 exported pairs |
| `figures/4_profile_gensan_maasim.png` | Change along the 34 km General Santos to Maasim line, all 7 pairs, with elevation |
| `figures/5_coherence_loss.png` | Stable pixels lost between the pre-quake pair and the pair spanning the quake |
| `optical/` | Sentinel-2 true-colour panels of three coasts, before and after the quake, matched by tide (captions below) |
| `pairs.csv` | All 7 pairs with dates, satellites, baselines, the share of land that unwrapped reliably, and the General Santos → Maasim change |

| GeoTIFF pattern | Value | Scale |
|---|---|---|
| `*_los_mm.tif` | Line-of-sight change, mm, positive = toward the satellite | 0.1 |
| `*_coherence.tif` | Coherence, 0 to 1 | 0.0001 |
| `*_reliable.tif` | 1 = unwrapping reliable, 0 = not | 1 |
| `*_wrapped_phase.tif` | Wrapped phase, radians; one cycle = 27.7 mm line of sight | 0.0001 |
| `t163_coherence_loss.tif`, `t069_coherence_loss.tif` | Figure 5 values: percentage points of stable pixels lost after the quake (positive = lost), 390 m blocks | 0.01 |

File names are `t<track>_<reference date>_<secondary date>_<quake|control>_*.tif`. The four exported pairs are the shortest pair spanning the quake on each track, plus its pre-quake control.

## Acquisition dates

| Track | Direction | Before | After | Days | Type |
|---|---|---|---|---|---|
| 163 | descending | 2026-05-28 (S1C) | 2026-06-03 (S1A) | 6 | control, exported |
| 163 | descending | 2026-06-03 (S1A) | 2026-06-10 (S1D) | 7 | spans quake, exported |
| 163 | descending | 2026-05-28 (S1C) | 2026-06-10 (S1D) | 13 | spans quake |
| 163 | descending | 2026-06-03 (S1A) | 2026-06-15 (S1A) | 12 | spans quake |
| 69 | ascending | 2026-05-22 (S1C) | 2026-06-03 (S1C) | 12 | control, exported |
| 69 | ascending | 2026-06-04 (S1D) | 2026-06-16 (S1D) | 12 | spans quake, exported |
| 69 | ascending | 2026-06-03 (S1C) | 2026-06-16 (S1D) | 13 | spans quake |

Dates are UTC.
- **Track 163** passes at 21:24 UTC (05:24 PHT the next day) and looks west.
- **Track 69** passes at 09:55 UTC (17:55 PHT the same day) and looks east.

NISAR has nothing public over this area before 2026-06-19, so it can't give a before-and-after pair.

## Processing steps

1. Picked Sentinel-1 dates on both tracks just before and just after 8 June, plus one pair per track entirely before the quake as a control.
2. Checked with an ASF search that every burst existed on every date.
3. Submitted 7 ASF HyP3 `INSAR_ISCE_MULTI_BURST` jobs: 10 × 2 looks (40 m pixels), water mask on, phase filter 0.5, SNAPHU unwrapping.
4. Converted unwrapped phase to line-of-sight change: mm = −phase × 55.465 / (4π). Positive means the ground moved toward the satellite.
5. Kept only the largest reliably unwrapped area of each pair (SNAPHU connected components). Referenced it to its own median north of 713000 N.
6. Averaged the wrapped phase into 400 m cells to show fringes where single pixels are too noisy.
7. Resampled everything to the team grid with nearest-neighbour resampling and wrote Int16 COGs.

## Observations

1. **Coherence is low.**
   - Only 10–18% of land unwrapped reliably in the 12- and 13-day pairs.
   - The 6-day control on track 163 reached 51%, and the 7-day pair spanning the quake reached 31%.
   - Most of the loss comes from vegetation changing between dates, which is normal for C-band radar in the tropics.
2. **The coast nearest the epicentre** (Glan, Malapatan, Jose Abad Santos, Balut Island) has no reliably unwrapped pixels in any pair spanning the quake.
3. **General Santos is stable** in radar terms.
   - Coherence there is 0.71–0.84 in the pairs spanning the quake.
   - It moved less than 8 mm relative to the area north of it in every pair spanning the quake.
   - The ascending control pair shows 27 mm there with no quake, another sign of atmosphere.
4. **Between General Santos and Maasim**, the change is atmosphere, not ground movement.
   - The ascending pairs spanning the quake show +4.3 and +7.8 cm toward the satellite, but the ascending control shows −3.9 cm with no quake.
   - Adding the control to the 06-03 → 06-16 pair gives 05-22 → 06-16, which also spans the quake: +0.4 cm.
   - The descending pairs show +1.0 to +3.4 cm, and their control +0.7 cm.
   - The phase follows elevation along this line (correlation −0.44 to −0.70), as atmospheric water vapour does.
   - No quake movement above about 4 cm is visible on this line.
5. **Coherence near Glan dropped after the quake in both directions.**
   - In a 4 km box, the share of pixels with coherence above 0.5 fell from 22% to 2% (descending) and from 23% to 4% (ascending).
   - Maasim, about the same distance from the epicentre (29 km), lost less: 26% to 22% (descending) and 26% to 15% (ascending).
   - Possible causes: ground disturbance, or movement too steep for the radar to follow. Rain or harvests could also do this.
   - On track 163 coherence fell over most of the area, so that track is weak evidence on its own.

## Cross-check with the landslide spots and NISAR

We tested whether the optical landslide spots in `../earthquake/` lost Sentinel-1 coherence more than other steep land. Steep land here means slopes of 20° or more that were checked and showed no slide. Coherence drop = control-pair coherence minus quake-pair coherence.

| | Spots | Other steep land |
|---|---|---|
| Track 163: median coherence drop | 0.024 | 0.025 |
| Track 69: median coherence drop | 0.000 | −0.003 |
| Track 163: coherence above 0.5 **before** the quake | 9% | 19% |
| Track 69: coherence above 0.5 **before** the quake | 2% | 4% |
| **NISAR** (L-band): flagged as changed land | **27.9%** | **11.9%** |

**Result: Sentinel-1 can't confirm or rule out the landslides.** Even before the quake, its coherence on these forested slopes was at the noise level (median about 0.3), so there was nothing to lose. Its 5.5 cm waves scatter off the leaves.

NISAR's 24 cm L-band waves reach further into the canopy. On the same spots, NISAR flagged changed land 2.3 times as often as on other steep land. This is why NISAR is the main radar in this project. The NISAR figure comes from the land-cover change classifier in `../nisar_ml/landcover/`, which compares 2026 with the 2020 map, so part of it may predate the quake.

## Optical check for coastal uplift (`optical/`)

If the coast rose by a metre or more, reef flats that used to stay under water should show up as a new pale strip along the shore. Sentinel-2 true-colour images (10 m) from before and after the quake were compared, matched by tide height, along three coasts:
- the east coast of Sarangani Island
- the south tip of the Glan peninsula
- the Davao Occidental coast around Jose Abad Santos

**Result:**
- **Sarangani Island and the Glan tip:** no visible uplift. A pale strip appears on 2026-06-14, but that morning had a very low tide, and pre-quake images at the same low tide show the same strip. At high tide the coast looks the same before and after the quake.
- **Jose Abad Santos:** the coast is too steep for this check to say much. A 1–2 m rise moves the waterline by less than one 10 m pixel.
- **Muddy water near the Jose Abad Santos river mouths** on 2026-06-14 and 2026-07-14 doesn't appear in the pre-quake images. It could come from rain or from landslide sediment upstream; this hasn't been tested. About 63 ha of the optical landslide spots are in Jose Abad Santos.

| File | What it shows (tide from the GOT4.10c model) |
|---|---|
| `lowtide_sarangani.png` | Sarangani Island east coast at low tide: 2024-05-07 (−90 cm), 2025-05-25 (−93 cm), 2026-04-30 (−76 cm) before; 2026-06-14 (−92 cm) after |
| `hightide_sarangani.png` | Same coast at high tide: 2026-04-22 (+71 cm), 2026-05-20 (+42 cm) before; 2026-06-19 (+63 cm) after |
| `lowtide_tip.png` | Glan peninsula south tip at low tide: 2025-05-25 (−93 cm), 2025-05-27 (−81 cm), 2026-04-30 (−76 cm) before; 2026-06-14 (−92 cm) after |
| `hightide_tip.png` | Same coast at high tide: 2026-04-22 (+71 cm), 2026-05-20 (+42 cm) before; 2026-06-19 (+63 cm) after |
| `jas_low.png` | Jose Abad Santos town coast at low tide: 2024-05-07 (−90 cm), 2024-08-18 (−82 cm), 2026-04-30 (−76 cm) before; 2026-06-14 (−95 cm), 2026-07-14 (−75 cm) after |
| `jas_high.png` | Same coast at high tide: 2026-04-09 (+53 cm), 2026-04-22 (+70 cm) before; 2026-06-08 (+43 cm, 10:06 PHT, after the quake), 2026-08-05 (+50 cm), 2026-08-18 (+42 cm) after |
| `jasN_low.png` | Coast north of Jose Abad Santos town at low tide, same dates as `jas_low.png` |

All images were taken between 10:05 and 10:15 PHT.

**Method:**
1. Listed every Sentinel-2 L2A scene over the coast and measured how much of a 2 km shoreline strip was cloud-free, using the SCL band.
2. Computed the tide at each image time with the NASA GSFC GOT4.10c model (pyTMD).
3. Matched post-quake images with pre-quake images at a similar tide.

## Limitations

- **Line of sight only.** Each number is movement toward or away from the satellite, a mix of vertical and mostly east–west horizontal motion. Separating them needs reliable data from both tracks at the same place, which we don't have near the coast.
- **No reliable radar data where the uplift was reported.** Sentinel-1's 5.5 cm wavelength loses coherence over vegetation, and large or steep movement near the epicentre makes it worse.
- **Atmosphere.** Single pairs carry centimetre-level water-vapour signals, as the control pairs show.
- **Relative values only.** `*_los_mm.tif` covers only each pair's largest reliable area (2–16% of the product). The values are relative to that area's northern part, not absolute.
- **Aftershocks.** The pairs spanning the quake also include up to 8.4 days of aftershocks and after-slip.
- **Optical check.**
  - Only three coasts were checked, and the Glan town coast was not compared by tide.
  - The half-degree tide model may be off by 10–20 cm near this coast.
  - 10 m pixels catch a metre-scale rise of a reef flat, not a few tens of centimetres.
  - Haze and sun glint differ between dates.

## What could show the uplift instead

NISAR L-band pairs after 19 June could show post-quake movement, since L-band keeps coherence through vegetation much better; this hasn't been tested.

## Sources and credit

- **Radar:** Copernicus Sentinel-1 IW SLC, VV, 2026, from the Alaska Satellite Facility (ASF DAAC).
  - Track 163 bursts 349004–349009 and track 69 bursts 146086–146091, sub-swaths IW2 and IW3.
  - Processed with ASF HyP3 `INSAR_ISCE_MULTI_BURST` (hyp3_isce2 4.1.4, ISCE 2.6.5).
- **Elevation:** Copernicus GLO-30 DEM.
- **Optical:** Copernicus Sentinel-2 L2A, read from the Earth Search STAC catalogue (Element 84).
- **Tides:** Ray, R. D., NASA GSFC GOT4.10c ocean tide model, computed with pyTMD.
- **Earthquake:** USGS event us7000srb1, M7.8, 2026-06-07 23:37 UTC (06-08 07:37 PHT), 5.599 N 125.056 E, 57 km deep.

Credit lines for anything published:
- "Contains modified Copernicus Sentinel data 2026, processed by ASF DAAC HyP3 2026."
- "Contains modified Copernicus Sentinel data 2024–2026."
