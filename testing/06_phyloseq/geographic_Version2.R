# geographic.R
# maps, isolation by distance, dbMEMs, spatial patterns

source("setup.R")

outdir <- file.path(output_root, "geographic")

# ---- control block ----
use_ps      <- ps_all_combined
use_lakes   <- NULL

# ---- build geo metadata ----
# use station-renamed morph object so IDs match the combined dataset
meta_geo <- data.frame(sample_data(ps_morph_by_station)) %>%
  rownames_to_column("Sample_ID")
meta_geo <- filter_geo_metadata(meta_geo, c("Latitude", "Longitude"))
meta_geo$Lake <- factor(meta_geo$Lake, levels = lake_order, ordered = TRUE)

if (!is.null(use_lakes))
  meta_geo <- meta_geo %>% filter(Lake %in% use_lakes)

# attach richness from combined data
shared_geo <- intersect(meta_geo$Sample_ID, sample_names(use_ps))
if (length(shared_geo) > 0) {
  alpha_geo <- compute_alpha(prune_samples(shared_geo, use_ps),
                              "Combined", alpha_metrics, lake_order)
  meta_geo <- left_join(meta_geo,
                         alpha_geo %>% select(Sample_ID, Observed, InvSimpson),
                         by = "Sample_ID")
} else {
  cat("no shared samples between geo metadata and combined data\n")
  meta_geo$Observed <- NA_real_
  meta_geo$InvSimpson <- NA_real_
}

basemap <- ne_states(country = c("United States of America", "Canada"),
                     returnclass = "sf")

# ============================================================
# sampling map
# ============================================================

p_map <- ggplot() +
  geom_sf(data = basemap, fill = "grey92", color = "grey60") +
  geom_point(data = meta_geo,
             aes(Longitude, Latitude, color = Lake, size = Observed),
             alpha = 0.85) +
  scale_color_manual(values = lake_colors) +
  coord_sf(xlim = c(-93, -75), ylim = c(41, 49)) +
  theme_minimal(base_size = 10) +
  labs(title = "Sampling sites (size = observed richness)", size = "Richness")
save_plot(p_map, file.path(outdir, "figures", "map_richness.pdf"))

# ============================================================
# isolation by distance
# ============================================================

ps_geo <- prune_samples(meta_geo$Sample_ID, use_ps) %>% to_pa()

ibd_res <- map(lake_order, ~ tryCatch(
  run_ibd(ps_geo, meta_geo, .x),
  error = function(e) { cat("  IBD failed for", .x, ":", e$message, "\n"); NULL }))
names(ibd_res) <- lake_order

mantel_tbl <- imap_dfr(ibd_res, function(r, lk) {
  if (is.null(r)) return(tibble(Lake = lk, rho = NA, p = NA))
  tibble(Lake = lk, rho = r$r_stat, p = r$p_val)
}) %>% mutate(sig = sig_stars(p))

save_stats(mantel_tbl,
           file.path(outdir, "stats", "mantel_ibd"),
           caption = "Mantel tests - isolation by distance (within-lake)")

# distance-decay plots
dd_all <- bind_rows(compact(map(ibd_res, "dd")))
if (nrow(dd_all) > 0) {
  p_dd <- ggplot(dd_all, aes(geo_km, dissim, color = Lake)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~ Lake, scales = "free") +
    scale_color_manual(values = lake_colors) +
    theme_minimal(base_size = 10) +
    labs(title = "Distance-decay of community similarity",
         x = "Geographic distance (km)", y = "Jaccard dissimilarity")
  save_plot(p_dd, file.path(outdir, "figures", "ibd_distance_decay.pdf"),
            width = 14, height = 10)
}

# ============================================================
# distance from lake centroid
# ============================================================

meta_geo <- meta_geo %>%
  group_by(Lake) %>%
  mutate(
    cent_lon = mean(Longitude, na.rm = TRUE),
    cent_lat = mean(Latitude, na.rm = TRUE),
    dist_center_km = distHaversine(
      cbind(Longitude, Latitude),
      cbind(cent_lon, cent_lat)) / 1000
  ) %>% ungroup()

p_cent <- ggplot(meta_geo, aes(dist_center_km, Observed, color = Lake)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ Lake, scales = "free") +
  scale_color_manual(values = lake_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Richness vs distance from lake centroid",
       x = "Distance (km)", y = "Observed richness")
save_plot(p_cent, file.path(outdir, "figures", "richness_vs_centroid_distance.pdf"),
          width = 14, height = 10)

# ============================================================
# dbMEMs
# ============================================================

ps_geo_pa <- prune_samples(meta_geo$Sample_ID, use_ps) %>% to_pa()

coords_df <- meta_geo %>%
  filter(Sample_ID %in% sample_names(ps_geo_pa)) %>%
  arrange(match(Sample_ID, sample_names(ps_geo_pa))) %>%
  select(Longitude, Latitude)

geo_dm <- as.dist(distm(coords_df, fun = distHaversine) / 1000)
dbmem_out <- dbmem(geo_dm, silent = TRUE)

dbmem_txt <- c()
if (!is.null(dbmem_out) && ncol(as.data.frame(dbmem_out)) > 0) {
  n_mem <- ncol(as.data.frame(dbmem_out))
  cat("positive dbMEM eigenvectors:", n_mem, "\n")

  otu_geo <- as(otu_table(ps_geo_pa), "matrix")
  if (taxa_are_rows(ps_geo_pa)) otu_geo <- t(otu_geo)

  rda_global <- rda(otu_geo ~ ., data = as.data.frame(dbmem_out))
  anova_global <- anova(rda_global, permutations = 999)
  rda_p <- anova_global$`Pr(>F)`[1]
  cat("global RDA p-value:", rda_p, "\n")

  dbmem_txt <- c(sprintf("dbMEM eigenvectors: %d", n_mem),
                  sprintf("Global RDA p-value: %.4f %s", rda_p, sig_stars(rda_p)))

  if (rda_p < 0.05) {
    fwd <- forward.sel(otu_geo, as.data.frame(dbmem_out),
                       alpha = 0.05, nperm = 999)
    save_stats(as.data.frame(fwd),
               file.path(outdir, "stats", "dbmem_forward_selection"),
               caption = "Forward-selected dbMEM eigenvectors")
    dbmem_txt <- c(dbmem_txt, sprintf("Selected dbMEMs: %d", nrow(fwd)))
  }
} else {
  dbmem_txt <- "No positive dbMEM eigenvectors found"
  cat(dbmem_txt, "\n")
}

# ============================================================
# map NMDS scores onto geography
# ============================================================

ord_geo <- ordinate(ps_geo_pa, "NMDS", "jaccard")
scores_df <- as.data.frame(scores(ord_geo, display = "sites"))
scores_df$Sample_ID <- rownames(scores_df)
map_df <- left_join(meta_geo, scores_df, by = "Sample_ID")

p_nmds1 <- ggplot() +
  geom_sf(data = basemap, fill = "grey92", color = "grey60") +
  geom_point(data = map_df, aes(Longitude, Latitude, color = NMDS1), size = 4) +
  scale_color_viridis_c() +
  coord_sf(xlim = c(-93, -75), ylim = c(41, 49)) +
  theme_minimal() + labs(title = "NMDS Axis 1")

p_nmds2 <- ggplot() +
  geom_sf(data = basemap, fill = "grey92", color = "grey60") +
  geom_point(data = map_df, aes(Longitude, Latitude, color = NMDS2), size = 4) +
  scale_color_viridis_c(option = "magma") +
  coord_sf(xlim = c(-93, -75), ylim = c(41, 49)) +
  theme_minimal() + labs(title = "NMDS Axis 2")

p_nmds_map <- p_nmds1 + p_nmds2
save_plot(p_nmds_map, file.path(outdir, "figures", "map_nmds_axes.pdf"), width = 14, height = 8)

# ============================================================
# text summary
# ============================================================

txt <- c(
  "=== GEOGRAPHIC ANALYSIS RESULTS ===",
  paste("Generated:", Sys.Date()), "",
  sprintf("Sites with geo data: %d", nrow(meta_geo)), "",
  "--- Mantel tests (IBD within-lake) ---")
for (i in seq_len(nrow(mantel_tbl))) {
  r <- mantel_tbl[i, ]
  if (is.na(r$rho)) {
    txt <- c(txt, sprintf("  %s: insufficient data", r$Lake))
  } else {
    txt <- c(txt, sprintf("  %s: rho = %.3f, p = %.4f %s",
                           r$Lake, r$rho, r$p, r$sig))
  }
}
txt <- c(txt, "", "--- dbMEMs ---", dbmem_txt)

save_summary(txt, file.path(outdir, "stats", "geographic_results_summary.txt"))