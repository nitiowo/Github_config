# alpha.R
# alpha diversity: compute, compare between markers and lakes, concordance

source("setup.R")

outdir <- file.path(output_root, "alpha")

# ---- control block ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL
use_metrics <- c("Observed", "InvSimpson")
use_rank    <- "ASV"

# ---- filter and aggregate ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)
if (use_rank != "ASV") ps_filt <- lapply(ps_filt, agg_rank, use_rank)

# ---- compute ----
alpha_all <- compute_alpha_all(ps_filt, use_metrics, lake_order)

alpha_long <- alpha_all %>%
  pivot_longer(cols = all_of(use_metrics), names_to = "Metric", values_to = "Value")

# ---- summary table ----
alpha_summary <- alpha_all %>%
  group_by(Marker, Lake) %>%
  summarise(across(all_of(use_metrics),
                   list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
            n = n(), .groups = "drop")

save_stats(alpha_summary,
           file.path(outdir, "stats", "alpha_summary_by_marker_lake"),
           caption = "Alpha diversity by marker and lake (mean +/- SD)")

# ---- between markers ----
kw_mk <- run_kruskal(alpha_long, "Marker")
pw_mk <- run_pairwise_wilcox(alpha_long, "Marker")

save_stats(kw_mk, file.path(outdir, "stats", "alpha_kruskal_markers"),
           caption = "Kruskal-Wallis: alpha diversity between markers")
save_stats(pw_mk, file.path(outdir, "stats", "alpha_pairwise_markers"),
           caption = "Pairwise Wilcoxon: alpha diversity between markers")

# ---- between lakes per marker ----
kw_lk <- run_kruskal(alpha_long, "Lake", group_by_vars = "Marker")
pw_lk <- run_pairwise_wilcox(alpha_long, "Lake", group_by_vars = "Marker")

save_stats(kw_lk, file.path(outdir, "stats", "alpha_kruskal_lakes_per_marker"),
           caption = "Kruskal-Wallis: alpha diversity between lakes (per marker)")
save_stats(pw_lk %>% filter(p.adj < 0.05),
           file.path(outdir, "stats", "alpha_pairwise_lakes_significant"),
           caption = "Significant pairwise lake comparisons")

# ---- concordance ----
alpha_wide <- alpha_all %>%
  select(Sample_ID, Marker, Station_ID, Lake, all_of(use_metrics)) %>%
  pivot_longer(cols = all_of(use_metrics), names_to = "Metric", values_to = "Value") %>%
  pivot_wider(names_from = Marker, values_from = Value)

marker_pairs <- combn(names(ps_filt), 2, simplify = FALSE)

cor_table <- alpha_wide %>%
  group_by(Metric) %>%
  reframe({
    map_dfr(marker_pairs, function(pr) {
      vals <- cur_data() %>% select(all_of(pr)) %>% drop_na()
      if (nrow(vals) < 5) return(NULL)
      ct <- cor.test(vals[[1]], vals[[2]], method = "spearman")
      tibble(Marker_A = pr[1], Marker_B = pr[2],
             rho = ct$estimate, p = ct$p.value, n = nrow(vals))
    })
  }) %>%
  mutate(sig = sig_stars(p))

save_stats(cor_table, file.path(outdir, "stats", "alpha_concordance_spearman"),
           caption = "Spearman correlations of site-level alpha diversity")

# ---- text summary ----
txt <- c(
  "=== ALPHA DIVERSITY RESULTS ===",
  paste("Generated:", Sys.Date()), "",
  "--- Between markers ---")
for (i in seq_len(nrow(kw_mk))) {
  txt <- c(txt, sprintf("  %s: H = %.1f, df = %.0f, p = %.2e %s",
                         kw_mk$Metric[i], kw_mk$H[i], kw_mk$df[i],
                         kw_mk$p_value[i], kw_mk$sig[i]))
}
txt <- c(txt, "", "--- Between lakes (per marker) ---")
for (i in seq_len(nrow(kw_lk))) {
  txt <- c(txt, sprintf("  %s %s: H = %.1f, df = %.0f, p = %.2e %s",
                         kw_lk$Marker[i], kw_lk$Metric[i], kw_lk$H[i],
                         kw_lk$df[i], kw_lk$p_value[i], kw_lk$sig[i]))
}
txt <- c(txt, "", "--- Concordance (Spearman) ---")
for (i in seq_len(nrow(cor_table))) {
  txt <- c(txt, sprintf("  %s vs %s (%s): rho = %.3f, p = %.2e %s",
                         cor_table$Marker_A[i], cor_table$Marker_B[i],
                         cor_table$Metric[i], cor_table$rho[i],
                         cor_table$p[i], cor_table$sig[i]))
}
save_summary(txt, file.path(outdir, "stats", "alpha_results_summary.txt"))

# ---- plots ----

# between markers
p_mk <- ggplot(alpha_long, aes(x = Marker, y = Value, fill = Marker)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = marker_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Alpha diversity - between markers", y = "Diversity value") +
  theme(legend.position = "bottom")
save_plot(p_mk, file.path(outdir, "figures", "alpha_boxplot_between_markers.pdf"))

# between lakes, faceted by metric and marker
p_lk <- ggplot(alpha_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_grid(Metric ~ Marker, scales = "free_y") +
  scale_fill_manual(values = lake_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Alpha diversity - between lakes (per marker)", y = "Diversity value") +
  theme(legend.position = "bottom")
save_plot(p_lk, file.path(outdir, "figures", "alpha_boxplot_lake_by_marker.pdf"),
          width = 16, height = 10)

# concordance scatterplots
scatter_plots <- purrr::compact(map(marker_pairs, function(pr) {
  df <- alpha_wide %>% filter(Metric == "Observed") %>%
    select(Sample_ID, Lake, all_of(pr)) %>% drop_na()
  if (nrow(df) < 3) return(NULL)
  ggplot(df, aes(x = .data[[pr[1]]], y = .data[[pr[2]]], color = Lake)) +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
    scale_color_manual(values = lake_colors) +
    theme_minimal(base_size = 10) +
    labs(title = paste(pr[1], "vs", pr[2]))
}))

if (length(scatter_plots) > 0) {
  p_conc <- wrap_plots(scatter_plots, ncol = 3) +
    plot_annotation(title = "Site-level richness concordance (Observed)")
  save_plot(p_conc, file.path(outdir, "figures", "alpha_concordance_scatter.pdf"),
            width = 16, height = 10)
}

# combined markers alpha
alpha_comb <- compute_alpha(ps_markers_combined, "Combined", use_metrics, lake_order)
alpha_comb_long <- alpha_comb %>%
  pivot_longer(cols = all_of(use_metrics), names_to = "Metric", values_to = "Value")

p_comb <- ggplot(alpha_comb_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = lake_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Alpha diversity - combined markers (P/A)", y = "Diversity value")
save_plot(p_comb, file.path(outdir, "figures", "alpha_boxplot_combined_markers.pdf"),
          width = 10, height = 6)