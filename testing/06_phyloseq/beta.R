# beta.R
# ordination, PERMANOVA, betadisper

source("setup.R")

outdir <- file.path(output_root, "beta")

# ---- control block ----
use_ps_list  <- ps_all_methods
use_markers  <- NULL
use_lakes    <- NULL
use_distance <- "jaccard"
use_binary   <- TRUE
use_method   <- "NMDS"
color_var    <- "Lake"
formula_str  <- "~ Lake"

# ---- filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)

# ============================================================
# within-marker ordinations (Jaccard P/A)
# ============================================================

ord_jac <- imap(ps_filt, ~ run_ordination(
  .x, "NMDS", "jaccard", TRUE, "Lake",
  title = paste(.y, "- NMDS Jaccard"),
  color_palette = lake_colors))

p_jac <- wrap_plots(map(ord_jac, "plot"), ncol = 2) +
  plot_annotation(title = "NMDS Jaccard (P/A)")
save_plot(p_jac, file.path(outdir, "figures", "nmds_jaccard_lake_all_markers.pdf"),
          width = 16, height = 14)

# PERMANOVA - Jaccard
perm_jac <- imap_dfr(ps_filt, function(ps, m) {
  res <- run_permanova(ps, "~ Lake", "jaccard", TRUE)
  as.data.frame(res) %>% rownames_to_column("Term") %>% mutate(Marker = m)
})
perm_jac_clean <- perm_jac %>% filter(!is.na(`Pr(>F)`))
save_stats(perm_jac_clean,
           file.path(outdir, "stats", "permanova_jaccard_by_marker"),
           caption = "PERMANOVA - lake effect (Jaccard P/A)")

# ============================================================
# within-marker ordinations (Bray-Curtis abundance)
# ============================================================

ord_bray <- imap(ps_filt, ~ run_ordination(
  .x, "NMDS", "bray", FALSE, "Lake",
  title = paste(.y, "- NMDS Bray-Curtis"),
  color_palette = lake_colors))

p_bray <- wrap_plots(map(ord_bray, "plot"), ncol = 2) +
  plot_annotation(title = "NMDS Bray-Curtis (abundance)")
save_plot(p_bray, file.path(outdir, "figures", "nmds_bray_lake_all_markers.pdf"),
          width = 16, height = 14)

# PERMANOVA - Bray-Curtis
perm_bray <- imap_dfr(ps_filt, function(ps, m) {
  res <- run_permanova(ps, "~ Lake", "bray", FALSE)
  as.data.frame(res) %>% rownames_to_column("Term") %>% mutate(Marker = m)
})
perm_bray_clean <- perm_bray %>% filter(!is.na(`Pr(>F)`))
save_stats(perm_bray_clean,
           file.path(outdir, "stats", "permanova_bray_by_marker"),
           caption = "PERMANOVA - lake effect (Bray-Curtis)")

# ============================================================
# betadisper
# ============================================================

bd_results <- imap_dfr(ps_filt, function(ps, m) {
  bd <- run_betadisper(ps, "Lake", "jaccard", TRUE)
  pt <- bd$permtest
  tibble(Marker = m, F_stat = pt$tab$F[1],
         p_value = pt$tab$`Pr(>F)`[1])
}) %>% mutate(sig = sig_stars(p_value))

save_stats(bd_results,
           file.path(outdir, "stats", "betadisper_jaccard"),
           caption = "Betadisper - homogeneity of dispersion (Jaccard)")

# ============================================================
# between-marker comparison
# ============================================================

shared_samps <- Reduce(intersect, map(ps_markers, sample_names))
ps_mk_cmp <- build_marker_ps(ps_markers, "Species", shared_samps)

ord_mk <- run_ordination(ps_mk_cmp, "PCoA", "jaccard", TRUE, "Marker",
                          title = "PCoA - between markers (Jaccard P/A)",
                          color_palette = marker_colors)
save_plot(ord_mk$plot, file.path(outdir, "figures", "pcoa_jaccard_between_markers.pdf"))

perm_mk <- run_permanova(ps_mk_cmp, "~ Marker", "jaccard", TRUE)
perm_mk_df <- as.data.frame(perm_mk) %>% rownames_to_column("Term")
save_stats(perm_mk_df,
           file.path(outdir, "stats", "permanova_marker_effect"),
           caption = "PERMANOVA - marker effect (Jaccard P/A)")

# ============================================================
# combined data ordination
# ============================================================

p_comb1 <- run_ordination(ps_markers_combined, "NMDS", "jaccard", TRUE, "Lake",
                           title = "Combined markers",
                           color_palette = lake_colors)$plot
p_comb2 <- run_ordination(ps_all_combined, "NMDS", "jaccard", TRUE, "Lake",
                           title = "Combined + morphology",
                           color_palette = lake_colors)$plot

p_combined <- p_comb1 + p_comb2 +
  plot_annotation(title = "Beta diversity - combined data (P/A)")
save_plot(p_combined, file.path(outdir, "figures", "nmds_jaccard_combined.pdf"),
          width = 14, height = 8)

# ============================================================
# text summary
# ============================================================

txt <- c(
  "=== BETA DIVERSITY RESULTS ===",
  paste("Generated:", Sys.Date()), "",
  "--- PERMANOVA: Lake effect (Jaccard P/A) ---")
for (i in seq_len(nrow(perm_jac_clean))) {
  row <- perm_jac_clean[i, ]
  txt <- c(txt, sprintf("  %s: R2 = %.3f, F = %.2f, p = %s %s",
                         row$Marker, row$R2, row$F, row$`Pr(>F)`,
                         sig_stars(row$`Pr(>F)`)))
}
txt <- c(txt, "", "--- PERMANOVA: Lake effect (Bray-Curtis) ---")
for (i in seq_len(nrow(perm_bray_clean))) {
  row <- perm_bray_clean[i, ]
  txt <- c(txt, sprintf("  %s: R2 = %.3f, F = %.2f, p = %s %s",
                         row$Marker, row$R2, row$F, row$`Pr(>F)`,
                         sig_stars(row$`Pr(>F)`)))
}
txt <- c(txt, "", "--- PERMANOVA: Marker effect (Jaccard P/A) ---")
mk_row <- perm_mk_df %>% filter(Term == "Model")
if (nrow(mk_row) > 0) {
  txt <- c(txt, sprintf("  R2 = %.3f, F = %.1f, p = %s %s",
                         mk_row$R2, mk_row$F, mk_row$`Pr(>F)`,
                         sig_stars(mk_row$`Pr(>F)`)))
}
txt <- c(txt, "", "--- Betadisper (Jaccard P/A) ---")
for (i in seq_len(nrow(bd_results))) {
  txt <- c(txt, sprintf("  %s: F = %.2f, p = %.4f %s",
                         bd_results$Marker[i], bd_results$F_stat[i],
                         bd_results$p_value[i], bd_results$sig[i]))
}
save_summary(txt, file.path(outdir, "stats", "beta_results_summary.txt"))