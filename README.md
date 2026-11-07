# Sarangani and Davao Occidental: NISAR, Landsat and MODIS layers (2025–2027)

GIS layers for the provinces of **Sarangani** and **Davao Occidental**, Mindanao, Philippines. They cover:
- the **M7.8 earthquake of 8 June 2026**
- the **2026 El Niño**
- a **12-month outlook** (Nov 2026 – Oct 2027)

NISAR L-band radar is the main input. Landsat HLS, MODIS, GPM rainfall and SRTM terrain support it. NASA Space Apps Challenge 2026.

The interactive web map that uses these layers lives in a separate repository: [nisar-mindanao](https://github.com/monteropeterdominic-arch/nisar-mindanao) ([live map](https://monteropeterdominic-arch.github.io/nisar-mindanao/)).

## Format

| | |
|---|---|
| CRS | **EPSG:32651** (WGS 84 / UTM 51N) for every raster and the boundary (the boundary also comes in EPSG:4326) |
| Grid | 30 m, 5147 × 4810 px, origin 648150 E, 735870 N. MODIS layers are 250 m or 1 km. |
| Extent | Both provinces. Pixels outside the boundary are nodata. |
| Storage | Compressed GeoTIFF. Most layers are **Int16 with a scale factor** in the file metadata. QGIS and ArcGIS apply the scale automatically. In Python, multiply by the scale: `rasterio` → `src.scales`, GDAL → `band.GetScale()`. Nodata is −32768 (Int16), 255 (Byte) or NaN (Float32). |

| Kind of value | Scale | Example |
|---|---|---|
| Radar backscatter, change (dB) | 0.01 | stored −1191 = −11.91 dB |
| EVI, NDVI, NDMI, NDWI, NBR | 0.0001 | stored 5440 = 0.544 |
| Land temperature (°C), slope and aspect (°) | 0.01 | |
| Elevation (m) | 0.1 | |
| Classes, counts | none (Byte) | |

## Folders

### `nisar/`: NISAR L2 GCOV radar (NASA–ISRO, provisional)
HH and HV backscatter in dB. The original 10 m data are averaged to 30 m in linear power. HV (cross-pol) responds to leafy and woody vegetation. HH (co-pol) responds to soil moisture, roughness and flooding.
- `backscatter_by_month/`: **September 2026** (4 passes) and **October 2026** (3 dates, 7 frames). Each covers 100% of both provinces. Every pass in the month is averaged in power. June 2026 is left out because only one June pass has been downloaded so far (28 June, track 39, about a third of the area). ASF lists four June passes over the area: 19 Jun (track 75), 21 Jun (111), 26 Jun (3) and 28 Jun (39).
- `backscatter_by_date/`: the **8 dates downloaded so far**, out of the 38 that ASF lists over the area up to 19 October: 28 Jun, 20/23/25/30 Sep, 14/17/19 Oct 2026. Each file merges every track and frame flown that day (frames on the same date are averaged), so these are not yet one file per track. HH coverage of the two provinces ranges from 19% (30 Sep) to 89% (23 Sep and 17 Oct). Per track, the downloads are: track 39 – 3 dates (28 Jun, 20 Sep, 14 Oct), track 75 – 2, track 111 – 2, track 3 – 1.
- `change_by_month/`: **September → October 2026** change in dB, for HH and HV. Each value is the per-pixel average of every same-track Sep→Oct pair: 5 pairs for HH (covers 99.9%, median +0.46 dB) and 3 pairs for HV (covers 87%, median +0.34 dB).
- `change_same_track/`: 11 pairs of passes on the **same orbit track and frame**, so the viewing geometry is identical. The value is the later pass minus the earlier one, in dB. The 3 km nearest each frame edge is removed.
  - 28 Jun → 20 Sep: HH median −0.55 dB (drying).
  - 20–25 Sep → 14–19 Oct: HH up to +1.05 dB (rain returned).

### `nisar_ml/`: machine learning with NISAR at the centre
- `vegetation_by_month/`, `vegetation_by_date/`: **EVI and NDMI estimated from NISAR radar plus terrain.**
  - Gradient-boosted trees trained on 168,000 pixels where NISAR and Landsat saw the ground within 8 days.
  - 5 km spatial block cross-validation: EVI R² 0.37, NDMI R² 0.32.
  - This fills Landsat's cloud gaps.
- `landcover/`: **2026 land cover from NISAR.**
  - 8 classes, learned from the 2020 land-cover map. Block CV: 93% forest vs non-forest.
  - Includes classifier confidence (%) and possible change versus 2020: 1 agree, 2 forest → non-forest, 3 vegetation → built/bare, 4 land → water, 5 other.
  - Also includes the 2020 map in the same 8 classes.

Classes: 1 closed forest, 2 open forest and mangrove, 3 brush/shrubs, 4 grassland, 5 annual crop, 6 perennial crop, 7 built-up and barren, 8 water and fishpond.

### `earthquake/`: M7.8, 8 June 2026
The event: 7:37 am PHT, epicentre 5.60°N 125.06°E, about 32 km south of Maasim. Depth 33 km (PHIVOLCS) or 57 km (USGS). Thrust faulting at the Cotabato Trench.
- Landsat HLS EVI and NBR, before (Mar – 7 Jun) and after (9 Jun – Sep), for 2026 and the same season of 2025. Also their differences, and the **2026 change minus the normal 2025 seasonal change**.
- `Landslide_Candidates_after_8Jun2026_30m.tif` plus polygons: slope ≥ 20°, EVI dropped from ≥ 0.35 to ≤ 0.20, and the drop exceeds the 2025 seasonal change by at least 0.25. Result: **305 patches, 165 ha**, two thirds around Glan and Jose Abad Santos. PhilSA estimated about 138 ha.
- Raster classes: 1 checked, no slide · 2 possible landslide · 4 not checked (cloudy in every view) · 5 checked but no elevation data.
- `Possible_landslides_after_8Jun2026_ranked.csv`: 152 clusters ranked by area, with nearest place names.
- Not field-checked. Heavy rain or land clearing can look like a slide.

### `landsat_hls/`: Landsat HLS L30 v2.0, quarterly medians 2025 Q1 – 2026 Q3
EVI, NDMI (plant moisture) and NDWI (surface water), plus 2026-minus-2025 anomalies. `ValidObs` gives the number of clear views per pixel.

### `modis/`
Monthly maximum NDVI at 250 m (MOD13Q1) and mean daytime land temperature at 1 km (MOD11A2), Dec 2024 – Oct 2026. Includes 2026-minus-2025 anomalies.

### `outlook/`: 12-month El Niño scenario, Nov 2026 – Oct 2027
The scenario assumes a strong El Niño (ONI +2.5 °C, range +2.0 to +3.0), as NOAA CPC forecast on 8 Oct 2026.
- **Rainfall:** GPM IMERG V07, 1998 – Oct 2026. Nov–Jan El Niño strength explains Nov–Apr rain: r = −0.81, and leave-one-year-out r = 0.83. The expected total is about **41% of normal** for Nov 2026 – Apr 2027 (341 vs 825 mm).
- **Vegetation and heat:** for each MODIS pixel and each month ahead, the anomaly is modelled as a + b·ONI + c·(Sep–Oct anomaly). The pixel's own seasonal cycle and 2000–2026 trend are removed first.
- 12-band rasters, one band per month, with the band description naming the month: expected value, anomaly versus normal, chance of below normal, and sensitivity per +1 °C of ONI.
- Expected results:
  - EVI is lowest in **April 2027**: −0.073 on average, −0.084 in Sarangani.
  - About **71% of the land** has a 2-in-3 chance of below-normal greenness that month.
  - Daytime land temperature is about **+2.2 °C** above normal in April 2027.
- `Outlook_Nov2026_Oct2027_ElNino_scenario.csv`: month-by-month figures with rainfall ranges.
- `MODIS_AOI_monthly_2000_2026.csv`: the long-term area series.
- **Hindcasts** (each El Niño hidden in turn):
  - Timing and size were roughly right for 2002–03, 2009–10 and 2023–24.
  - **2015–16, the strongest event, was underestimated**: EVI −0.063 expected vs −0.119 actual. A very strong event may hit harder than shown.
  - This is a scenario, not a weather forecast.

### `terrain/`, `boundary/`, `tables/`, `scripts/`
- `terrain/`: SRTM 30 m elevation, slope and aspect.
- `boundary/`: the two provinces, dissolved from the 2020 regional land-cover map. Davao Occidental 2,003 km², Sarangani 3,261 km².
- `tables/`: monthly area averages for 2025–2026.
- `scripts/`: the PyQGIS/GDAL/scikit-learn scripts that produced everything here. Paths point to the original working folder.

## Limits
- No public NISAR data over the area covers the time before the quake yet: the earliest pass is 19 June, and ASF is still adding older NISAR data through the end of 2026. So NISAR backscatter shows the months after the quake, not the quake itself. Ground uplift from the quake needs Sentinel-1 interferometry spanning 8 June (for example path 69, 4 → 16 June, Sentinel-1D on both dates).
- About 30% of the area could not be checked for landslides, because every Landsat view in one of the windows was cloudy.
- The radar → vegetation estimates explain about a third of the variance. Use them for broad patterns, not field-level values.

## Sources
NASA/ISRO NISAR L2 GCOV (provisional) · NASA HLS L30 v2.0 · MODIS MOD13Q1 and MOD11A2 v6.1 (via AppEEARS) · GPM IMERG V07 (via GIOVANNI) · NASA SRTMGL1 · NOAA CPC Oceanic Niño Index and ENSO discussion (8 Oct 2026) · IRI ENSO forecast · PHIVOLCS, USGS, NDRRMC and PhilSA reports on the 8 June 2026 earthquake · 2020 regional land cover (NAMRIA).
