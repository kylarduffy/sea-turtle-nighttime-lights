# Nighttime radiance near Queensland sea turtle nesting habitat (2025)

An exploratory spatial analysis of satellite-observed nighttime radiance near Queensland sea turtle nesting coastlines, grouped by Queensland Government nesting significance classification.

**Research question:** How does artificial nighttime light exposure vary among sea turtle nesting areas of different conservation significance along the Queensland coast?

## Overview

This project compares two descriptive aspects of 2025 VIIRS nighttime radiance: **occurrence** (whether the mean radiance in a coastal sampling buffer is greater than zero) and **intensity** (the median radiance among positive observations). Satellite-observed upward radiance is used as a *proxy for nearby nighttime light*, not a direct measure of illumination experienced by turtles on the beach.

## Data sources

- **Queensland Government, QSpatial — Sea turtle nesting areas:** mapped coastline features and the `significan` field (`Present`, `Significant`, `Very significant`). Record the exact dataset URL, download date, version and reuse conditions before publication.
- **Earth Observation Group — VIIRS annual nighttime lights, 2025:** `VNL_npp_2025_global_vcmslcfg_v2_c202604011200.average_masked.dat.tif` (annual average masked radiance). Record the dataset's official citation, licence and download URL before publication.
- **Natural Earth:** medium-resolution Australian land boundary, accessed through the `rnaturalearth` R package for map backgrounds.

## Methods

1. Import nesting coastline data; drop the Z coordinate; split multipart lines into single lines; remove exact duplicate geometries (retaining the first associated record).
2. Divide the resulting **691 cleaned coastline lines** into **3,221 approximately 1 km or shorter segments**.
3. Crop the 2025 VIIRS raster to Queensland and clamp negative radiance estimates to zero.
4. Create **1 km buffers** around each segment in a projected coordinate system; transform buffers to the raster's geographic CRS.
5. Extract area-weighted mean VIIRS radiance for each buffer using `terra::extract(..., fun = mean, exact = TRUE)`.
6. Summarise positive-radiance occurrence and intensity by nesting significance; check sensitivity to segment-length weighting, original-line aggregation and **500 m buffers**.
7. Export summary tables and four presentation figures.

## Key descriptive results

| Nesting significance | Segments | Positive segments (1 km) | Median radiance among positive segments |
|---|---:|---:|---:|
| Present | 1,361 | 24.8% | 0.421 |
| Significant | 1,573 | 11.2% | 1.12 |
| Very significant | 287 | 12.9% | 0.826 |

Present habitat had the highest *proportion* of positive-radiance segments, whereas Significant habitat had the highest *median intensity among positive segments*. The ranking of categories remained unchanged when 500 m buffers were used, although individual values changed.

These are **descriptive**, not causal or inferential, comparisons.

## Reproduce the analysis

Use R and the packages `sf`, `terra`, `dplyr`, `lwgeom`, `ggplot2` and `rnaturalearth` (including the medium-scale Natural Earth data dependency, if prompted).

Open the RStudio project at the repository root and arrange the downloaded input data as follows:

```text
project-root/
├── data/
│   ├── nesting_areas/
│   │   └── Sea_turtle_nesting_data/
│   │       ├── Sea_turtle_nesting_areas.shp
│   │       ├── Sea_turtle_nesting_areas.shx
│   │       ├── Sea_turtle_nesting_areas.dbf
│   │       └── ... other shapefile components
│   └── nighttime_lights/
│       └── VNL_npp_2025_global_vcmslcfg_v2_c202604011200.average_masked.dat.tif
├── code/
│   └── clean_turtle_nighttime_lights.R
├── outputs/
├── README.md
└── .gitignore
```

Run `source("code/clean_turtle_nighttime_lights.R")` from the project root, or open the script in RStudio and run it from top to bottom. The script creates `outputs/` if necessary. The full-resolution raster may be large and extraction may take time.

**Important:** The input datasets are not bundled in this example repository structure. Obtain them from their original providers, check their licence terms, and place them at the paths shown. If your files have different names, update `nesting_path` and `viirs_path` near the top of the script.

## Outputs

- `figure1_nesting_significance.png` — mapped nesting significance categories.
- `figure2_positive_radiance.png` — percentage of segments with positive radiance.
- `figure3_positive_radiance_intensity.png` — median radiance among positive segments.
- `figure4_spatial_radiance.png` — spatial distribution of extracted radiance.
- Six `.csv` files with primary results and sensitivity summaries.

## Limitations

- VIIRS measures satellite-observed upward radiance, **not light levels at turtle nesting beaches**, spectral composition, or ecological impacts.
- Zero or masked radiance does **not** demonstrate the absence of artificial light.
- Buffers overlap, so adjacent segments are spatially dependent and should not be treated as independent samples.
- Buffer areas can include land and water; nesting coastline classifications are not measures of turtle abundance or nesting success.
- Exact duplicate geometries were removed, but partially overlapping geometries and differences in original coastline lengths may remain.
- Results are sensitive to spatial resolution and analytical definitions; the 500 m comparison is a limited sensitivity check.

## Repository notes

Keep the original data source credits, check redistribution permissions, and avoid committing large raw rasters or generated intermediate TIFFs. The supplied `.gitignore` excludes the `data/` directory and typical temporary files. Presentation figures and summary CSVs can be committed to `outputs/`.
