# Queensland sea turtle nesting habitat and nighttime radiance (2025)
# Research question: How does artificial nighttime light exposure vary among
# sea turtle nesting areas of different conservation significance in Queensland?
# Run from the RStudio project root. See README.md for data requirements.
 
# 1. Setup ---------------------------------------------------------------
library(sf)            # Spatial vector data and geometry operations
library(terra)         # VIIRS raster processing and zonal extraction
library(dplyr)         # Data manipulation and summaries
library(lwgeom)        # Subdivide coastlines into equal-length segments
library(ggplot2)       # Figures
library(rnaturalearth) # Coastline basemap

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

nesting_path <- "data/nesting_areas/Sea_turtle_nesting_data/Sea_turtle_nesting_areas.shp"
viirs_path <- "data/nighttime_lights/VNL_npp_2025_global_vcmslcfg_v2_c202604011200.average_masked.dat.tif"
stopifnot(file.exists(nesting_path), file.exists(viirs_path))

# 2. Read and clean nesting coastlines -----------------------------------
nesting_areas <- st_read(nesting_path, quiet = TRUE)
# Remove the Z dimension, then separate multipart coastlines into lines.
nesting_lines <- nesting_areas |>
  st_zm() |>
  st_cast("LINESTRING", warn = FALSE)

# Remove geometrically identical lines; keep the first associated record.
# This does not eliminate partial overlaps or reconcile conflicting metadata.
nesting_lines_clean <- nesting_lines[
  !duplicated(st_as_binary(st_geometry(nesting_lines), EWKB = TRUE)),
]

# The source CRS is projected in metres; check before buffering or measuring.
stopifnot(!st_is_longlat(nesting_lines_clean))

# 3. Prepare the 2025 VIIRS raster ---------------------------------------
night_lights <- rast(viirs_path)
qld_extent <- ext(137.5, 154.1, -28.5, -8.6) # Longitude, latitude (WGS84)
# The supplied VIIRS raster uses geographic coordinates; crop to the study area.
stopifnot(terra::is.lonlat(night_lights))
night_lights_qld <- crop(night_lights, qld_extent)
# Treat negative radiance estimates as zero, as in the exploratory analysis.
night_lights_analysis <- terra::clamp(
  night_lights_qld, lower = 0, values = TRUE
)

# 4. Create approximately 1 km coastline segments -----------------------
# Each original line is split into equal fractions, with length <= 1 km.
nesting_lines_clean <- nesting_lines_clean |>
  mutate(
    length_m = as.numeric(st_length(geometry)),
    n_segments = pmax(1L, ceiling(length_m / 1000))
  )

all_segments <- lapply(seq_len(nrow(nesting_lines_clean)), function(i) {
  line <- nesting_lines_clean[i, ]
  n <- line$n_segments
  pieces <- lapply(seq_len(n), function(j) {
    lwgeom::st_linesubstring(
      st_geometry(line), from = (j - 1) / n, to = j / n
    )
  })
  st_sf(
    beach_id = rep(line$beach_id, n),
    locality = rep(line$locality, n),
    significan = rep(line$significan, n),
    source_line = rep(i, n),
    geometry = st_sfc(lapply(pieces, function(x) x[[1]]), crs = st_crs(line))
  )
})
nesting_segments <- do.call(rbind, all_segments) |>
  mutate(segment_length_m = as.numeric(st_length(geometry)))

# 5. Extract mean radiance within 1 km buffers ---------------------------
# Buffer in the projected CRS (metres), then transform to raster CRS.
nesting_buffers <- st_buffer(nesting_segments, dist = 1000)
nesting_buffers_wgs84 <- st_transform(nesting_buffers, crs = 4326)
# exact=TRUE weights partially covered pixels by their covered fraction.
all_radiance <- terra::extract(
  night_lights_analysis, terra::vect(nesting_buffers_wgs84),
  fun = mean, na.rm = TRUE, exact = TRUE
)
stopifnot(nrow(all_radiance) == nrow(nesting_segments),
          identical(all_radiance$ID, seq_len(nrow(nesting_segments))))
nesting_segments$mean_radiance <- all_radiance[[2]]
stopifnot(!anyNA(nesting_segments$mean_radiance))

# 6. Descriptive summaries -----------------------------------------------
# Positive = extracted mean radiance > 0, not necessarily light on the beach.
radiance_by_significance <- nesting_segments |>
  st_drop_geometry() |>
  group_by(significan) |>
  summarise(
    n_segments = n(),
    positive_segments = sum(mean_radiance > 0),
    percent_positive = mean(mean_radiance > 0) * 100,
    average_radiance = mean(mean_radiance),
    median_radiance = median(mean_radiance),
    max_radiance = max(mean_radiance),
    .groups = "drop"
  )

positive_radiance_summary <- nesting_segments |>
  st_drop_geometry() |>
  filter(mean_radiance > 0) |>
  group_by(significan) |>
  summarise(
    n_positive = n(),
    mean_positive = mean(mean_radiance),
    median_positive = median(mean_radiance),
    max_positive = max(mean_radiance),
    .groups = "drop"
  )

# Segment-length weighting checks whether short segments distort occurrence.
radiance_by_length <- nesting_segments |>
  st_drop_geometry() |>
  group_by(significan) |>
  summarise(
    total_length_km = sum(segment_length_m) / 1000,
    positive_length_km = sum(segment_length_m[mean_radiance > 0]) / 1000,
    percent_positive_length = 100 * positive_length_km / total_length_km,
    .groups = "drop"
  )

# Original-line summaries are supplementary: lines differ in length.
radiance_by_line <- nesting_segments |>
  st_drop_geometry() |>
  group_by(source_line, significan) |>
  summarise(
    n_segments = n(),
    line_length_km = sum(segment_length_m) / 1000,
    mean_line_radiance = weighted.mean(mean_radiance, w = segment_length_m),
    any_positive = any(mean_radiance > 0),
    .groups = "drop"
  )
line_level_summary <- radiance_by_line |>
  group_by(significan) |>
  summarise(
    n_lines = n(),
    positive_lines = sum(any_positive),
    percent_positive_lines = mean(any_positive) * 100,
    average_line_radiance = mean(mean_line_radiance),
    median_line_radiance = median(mean_line_radiance),
    .groups = "drop"
  )

# 7. Buffer-distance sensitivity: 500 m versus 1 km ----------------------
# Repeat the same extraction with narrower buffers around the same segments.
nesting_buffers_500m <- st_transform(
  st_buffer(nesting_segments, dist = 500), crs = 4326
)
radiance_500m <- terra::extract(
  night_lights_analysis, terra::vect(nesting_buffers_500m),
  fun = mean, na.rm = TRUE, exact = TRUE
)
stopifnot(nrow(radiance_500m) == nrow(nesting_segments),
          identical(radiance_500m$ID, seq_len(nrow(nesting_segments))))
nesting_segments$mean_radiance_500m <- radiance_500m[[2]]
stopifnot(!anyNA(nesting_segments$mean_radiance_500m))

buffer_comparison <- nesting_segments |>
  st_drop_geometry() |>
  group_by(significan) |>
  summarise(
    n_segments = n(),
    positive_1km = mean(mean_radiance > 0) * 100,
    positive_500m = mean(mean_radiance_500m > 0) * 100,
    difference = positive_500m - positive_1km,
    .groups = "drop"
  )
intensity_buffer_comparison <- nesting_segments |>
  st_drop_geometry() |>
  group_by(significan) |>
  summarise(
    average_1km = mean(mean_radiance),
    average_500m = mean(mean_radiance_500m),
    median_positive_1km = median(mean_radiance[mean_radiance > 0]),
    median_positive_500m = median(mean_radiance_500m[mean_radiance_500m > 0]),
    .groups = "drop"
  )

# 8. Save summary tables -------------------------------------------------
write.csv(radiance_by_significance, "outputs/radiance_by_significance.csv", row.names = FALSE)
write.csv(positive_radiance_summary, "outputs/positive_radiance_summary.csv", row.names = FALSE)
write.csv(radiance_by_length, "outputs/radiance_by_length.csv", row.names = FALSE)
write.csv(line_level_summary, "outputs/line_level_summary.csv", row.names = FALSE)
write.csv(buffer_comparison, "outputs/buffer_comparison.csv", row.names = FALSE)
write.csv(intensity_buffer_comparison, "outputs/intensity_buffer_comparison.csv", row.names = FALSE)

# 9. Presentation figures ------------------------------------------------
# Same three category colours throughout the study-area map and bar charts.
australia_detailed <- ne_countries(
  country = "Australia", scale = "medium", returnclass = "sf"
)

# Figure 1: mapped nesting significance
nesting_significance_map <- ggplot() +
  geom_sf(
    data = australia_detailed,
    fill = "grey95",
    colour = "grey65",
    linewidth = 0.2
  ) +
  geom_sf(
    data = nesting_lines_clean,
    aes(colour = significan),
    linewidth = 0.8,
    alpha = 0.9
  ) +
  scale_colour_manual(
    name = "Nesting significance",
    values = c(
      "Present" = "#DDAA33",
      "Significant" = "#3388AA",
      "Very significant" = "#BB4477"
    ),
    breaks = c(
      "Present",
      "Significant",
      "Very significant"
    )
  ) +
  coord_sf(
    xlim = c(137.5, 154.1),
    ylim = c(-28.5, -8.6),
    expand = FALSE
  ) +
  labs(
    title = "Sea turtle nesting habitat along the Queensland coast",
    subtitle = "Queensland Government nesting significance classifications",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.2
    ),
    axis.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.background = element_rect(
      fill = "#EAF3F8",
      colour = NA
    )
  )

ggsave(
  filename = "outputs/figure1_nesting_significance.png",
  plot = nesting_significance_map,
  width = 7,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

# Figure 2: positive-radiance occurrence
figure2_data <- nesting_segments |>
  sf::st_drop_geometry() |>
  dplyr::group_by(significan) |>
  dplyr::summarise(
    n_segments = dplyr::n(),
    percent_positive = mean(mean_radiance > 0) * 100,
    .groups = "drop"
  )

figure2_plot <- ggplot(
  figure2_data,
  aes(
    x = percent_positive,
    y = factor(
      significan,
      levels = c("Very significant", "Significant", "Present")
    ),
    fill = significan
  )
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent_positive)),
    hjust = -0.2,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Present" = "#DDAA33",
      "Significant" = "#3388AA",
      "Very significant" = "#BB4477"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 30),
    breaks = seq(0, 30, 5),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Occurrence of positive nighttime radiance",
    subtitle = "Queensland sea turtle nesting habitat, 2025",
    x = "Coastline segments with positive radiance (%)",
    y = NULL,
    caption = "Mean VIIRS radiance within 1 km coastline buffers"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.caption = element_text(
      size = 9,
      hjust = 0
    ),
    plot.margin = margin(15, 25, 15, 15)
  )

figure2_plot <- figure2_plot +
  labs(
    caption = "Positive mean VIIRS radiance (>0) within 1 km coastline buffers; 2025 annual composite"
  )

ggsave(
  filename = "outputs/figure2_positive_radiance.png",
  plot = figure2_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

# Figure 3: positive-radiance intensity
figure3_data <- nesting_segments |>
  sf::st_drop_geometry() |>
  dplyr::filter(mean_radiance > 0) |>
  dplyr::group_by(significan) |>
  dplyr::summarise(
    n_positive = dplyr::n(),
    median_positive = median(mean_radiance),
    .groups = "drop"
  )

figure3_plot <- ggplot(
  figure3_data,
  aes(
    x = median_positive,
    y = factor(
      significan,
      levels = c("Very significant", "Significant", "Present")
    ),
    fill = significan
  )
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f", median_positive)),
    hjust = -0.2,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Present" = "#DDAA33",
      "Significant" = "#3388AA",
      "Very significant" = "#BB4477"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1.5),
    breaks = seq(0, 1.5, 0.25),
    expand = c(0, 0)
  ) +
  labs(
    title = "Intensity of positive nighttime radiance",
    subtitle = "Queensland sea turtle nesting habitat, 2025",
    x = "Median VIIRS radiance (nW/cm²/sr)",
    y = NULL,
    caption = "Positive radiance segments only (>0); mean radiance within 1 km coastline buffers"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.caption = element_text(
      size = 9,
      hjust = 0
    ),
    plot.margin = margin(15, 25, 15, 15)
  )

ggsave(
  filename = "outputs/figure3_positive_radiance_intensity.png",
  plot = figure3_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

# Figure 4: geographic pattern in radiance
radiance_map_presentation <- ggplot() +
  geom_sf(
    data = australia_detailed,
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.2
  ) +
  geom_sf(
    data = nesting_segments |>
      dplyr::filter(mean_radiance == 0),
    colour = "grey65",
    linewidth = 0.7
  ) +
  geom_sf(
    data = nesting_segments |>
      dplyr::filter(mean_radiance > 0),
    aes(colour = log1p(mean_radiance)),
    linewidth = 1.2
  ) +
  scale_colour_viridis_c(
    option = "magma",
    name = "Mean VIIRS radiance",
    breaks = log1p(c(0.1, 0.5, 1, 5, 10, 30)),
    labels = c("0.1", "0.5", "1", "5", "10", "30"),
    guide = guide_colourbar(
      barheight = grid::unit(6, "cm")
    )
  ) +
  coord_sf(
    xlim = c(137.5, 154.1),
    ylim = c(-28.5, -8.6),
    expand = FALSE
  ) +
  labs(
    title = "Nighttime radiance near sea turtle nesting habitat",
    subtitle = "Queensland coastline, 2025",
    x = NULL,
    y = NULL,
    caption = "Mean VIIRS radiance within 1 km coastline buffers; grey = zero radiance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.2
    ),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.background = element_rect(
      fill = "#EAF3F8",
      colour = NA
    ),
    plot.caption = element_text(size = 8)
  )

radiance_map_presentation <- radiance_map_presentation +
  labs(
    colour = "Mean VIIRS radiance\n(nW/cm²/sr)"
  )

ggsave(
  filename = "outputs/figure4_spatial_radiance.png",
  plot = radiance_map_presentation,
  width = 8,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)
