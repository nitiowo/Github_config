# test_all.R
# Run each analysis script with multiple parameter configurations
# Descriptive output filenames so you can tell what each file is at a glance

source("setup.R")

test_output <- file.path(output_root, "test_runs")
dir.create(test_output, showWarnings = FALSE, recursive = TRUE)

cat("=== TEST RUNNER ===\n")
cat("Output root:", test_output, "\n\n")

# Helper: run a code block, catch errors, log result
run_test <- function(test_name, expr) {
  cat("\n--- TEST:", test_name, "---\n")
  t0 <- Sys.time()
  result <- tryCatch(
    { force(expr); "PASS" },
    error = function(e) paste("FAIL:", e$message)
  )
  elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
  cat("  ", result, "(", elapsed, "s)\n")
  data.frame(test = test_name, result = result, seconds = as.numeric(elapsed),
             stringsAsFactors = FALSE)
}

results <- list()

# ============================================================
# 1. COMPOSITION: multiple parameter combos
# ============================================================

cat("\n========== COMPOSITION TESTS ==========\n")

comp_dir <- file.path(test_output, "composition", "figures")
dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

# -- Test 1a: Genus, facet by Marker, all methods
results[[length(results) + 1]] <- run_test("comp: genus facet-marker all", {
  ps_filt <- filter_ps_list(ps_all_methods, NULL, NULL)
  rank <- "Genus"
  df_all <- imap_dfr(ps_filt, function(ps, m) {
    ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(NULL)
    ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))
    df <- psmelt(ps_rel)
    df$Marker <- m
    df
  })
  df_all$Lake <- factor(df_all$Lake, levels = lake_order, ordered = TRUE)
  df_all$Marker <- factor(df_all$Marker, levels = names(ps_filt))
  top_taxa <- df_all %>% group_by(.data[[rank]]) %>%
    summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>% slice_head(n = 15) %>% pull(.data[[rank]])
  df_all[[rank]] <- if_else(df_all[[rank]] %in% top_taxa, as.character(df_all[[rank]]), "Other")
  df_all[[rank]] <- factor(df_all[[rank]], levels = c(as.character(top_taxa), "Other"))
  fill_pal <- c(colorRampPalette(brewer.pal(12, "Set3"))(15), "grey70")
  names(fill_pal) <- levels(df_all[[rank]])
  p <- ggplot(df_all, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    facet_wrap(~ Marker, scales = "free_x") +
    scale_fill_manual(values = fill_pal) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
          legend.position = "right") +
    labs(title = "All markers - Genus (faceted by Marker)")
  save_plot(p, file.path(comp_dir, "barplot_genus_all_facet-marker.pdf"), width = 16, height = 8)
})

# -- Test 1b: Phylum, facet by Lake, COI only
results[[length(results) + 1]] <- run_test("comp: phylum facet-lake coi-only", {
  ps_filt <- filter_ps_list(ps_coi, NULL, NULL)
  rank <- "Phylum"
  for (m in names(ps_filt)) {
    ps <- ps_filt[[m]]
    ps_agg <- agg_rank(ps, rank)
    ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))
    df <- psmelt(ps_rel)
    df$Marker <- m
    df$Lake <- factor(df$Lake, levels = lake_order, ordered = TRUE)
    top_taxa <- df %>% group_by(.data[[rank]]) %>%
      summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>% slice_head(n = 10) %>% pull(.data[[rank]])
    df[[rank]] <- if_else(df[[rank]] %in% top_taxa, as.character(df[[rank]]), "Other")
    df[[rank]] <- factor(df[[rank]], levels = c(as.character(top_taxa), "Other"))
    fill_pal <- c(colorRampPalette(brewer.pal(12, "Set3"))(10), "grey70")
    names(fill_pal) <- levels(df[[rank]])
    p <- ggplot(df, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
      geom_bar(stat = "identity", position = "stack", width = 1) +
      facet_wrap(~ Lake, scales = "free_x") +
      scale_fill_manual(values = fill_pal) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4)) +
      labs(title = paste(m, "- Phylum (faceted by Lake)"))
    save_plot(p, file.path(comp_dir,
      paste0("barplot_phylum_", tolower(m), "_facet-lake_coi-only.pdf")), width = 14, height = 8)
  }
})

# -- Test 1c: Order, no rotifers, Erie+Ontario only
results[[length(results) + 1]] <- run_test("comp: order no-rotifera erie-ontario", {
  tsub <- list(rank = "Phylum", exclude = "Rotifera")
  ps_filt <- filter_ps_list(ps_markers, NULL, c("Erie", "Ontario"), tsub = tsub)
  rank <- "Order"
  for (m in names(ps_filt)) {
    ps <- ps_filt[[m]]
    if (nsamples(ps) == 0 || ntaxa(ps) == 0) next
    ps_agg <- agg_rank(ps, rank)
    ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))
    df <- psmelt(ps_rel)
    df$Marker <- m
    df$Lake <- factor(df$Lake, levels = lake_order, ordered = TRUE)
    top_taxa <- df %>% group_by(.data[[rank]]) %>%
      summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>% slice_head(n = 15) %>% pull(.data[[rank]])
    df[[rank]] <- if_else(df[[rank]] %in% top_taxa, as.character(df[[rank]]), "Other")
    df[[rank]] <- factor(df[[rank]], levels = c(as.character(top_taxa), "Other"))
    fill_pal <- c(colorRampPalette(brewer.pal(12, "Set3"))(15), "grey70")
    names(fill_pal) <- levels(df[[rank]])
    p <- ggplot(df, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
      geom_bar(stat = "identity", position = "stack", width = 1) +
      scale_fill_manual(values = fill_pal) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
      labs(title = paste(m, "- Order (no Rotifera, Erie+Ontario)"))
    save_plot(p, file.path(comp_dir,
      paste0("barplot_order_", tolower(m), "_no-rotifera_erie-ontario.pdf")))
  }
})

# -- Test 1d: Genus, only Arthropoda, facet by Marker
results[[length(results) + 1]] <- run_test("comp: genus only-arthropoda facet-marker", {
  tsub <- list(rank = "Phylum", include = "Arthropoda")
  ps_filt <- filter_ps_list(ps_markers, NULL, NULL, tsub = tsub)
  rank <- "Genus"
  df_all <- imap_dfr(ps_filt, function(ps, m) {
    if (ntaxa(ps) == 0) return(data.frame())
    ps_agg <- agg_rank(ps, rank)
    ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))
    df <- psmelt(ps_rel)
    df$Marker <- m
    df
  })
  if (nrow(df_all) == 0) stop("no Arthropoda data")
  df_all$Lake <- factor(df_all$Lake, levels = lake_order, ordered = TRUE)
  df_all$Marker <- factor(df_all$Marker, levels = names(ps_filt))
  top_taxa <- df_all %>% group_by(.data[[rank]]) %>%
    summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>% slice_head(n = 15) %>% pull(.data[[rank]])
  df_all[[rank]] <- if_else(df_all[[rank]] %in% top_taxa, as.character(df_all[[rank]]), "Other")
  df_all[[rank]] <- factor(df_all[[rank]], levels = c(as.character(top_taxa), "Other"))
  fill_pal <- c(colorRampPalette(brewer.pal(12, "Set3"))(15), "grey70")
  names(fill_pal) <- levels(df_all[[rank]])
  p <- ggplot(df_all, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    facet_wrap(~ Marker, scales = "free_x") +
    scale_fill_manual(values = fill_pal) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4)) +
    labs(title = "Arthropoda only - Genus (faceted by Marker)")
  save_plot(p, file.path(comp_dir,
    "barplot_genus_all_facet-marker_only-arthropoda.pdf"), width = 16, height = 8)
})

# ============================================================
# 2. ALPHA DIVERSITY: multiple parameter combos
# ============================================================

cat("\n========== ALPHA DIVERSITY TESTS ==========\n")

alpha_dir <- file.path(test_output, "alpha")
dir.create(file.path(alpha_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(alpha_dir, "stats"), showWarnings = FALSE, recursive = TRUE)

# -- Test 2a: all methods, boxplot by marker
results[[length(results) + 1]] <- run_test("alpha: all methods by marker", {
  ps_filt <- filter_ps_list(ps_all_methods, NULL, NULL)
  alpha_all <- compute_alpha_all(ps_filt, alpha_metrics, lake_order)
  alpha_long <- alpha_all %>%
    pivot_longer(cols = all_of(alpha_metrics), names_to = "Metric", values_to = "Value")
  p <- ggplot(alpha_long, aes(x = Marker, y = Value, fill = Marker)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    facet_wrap(~ Metric, scales = "free_y") +
    scale_fill_manual(values = marker_colors) +
    theme_minimal(base_size = 10) +
    labs(title = "Alpha diversity - all methods by marker")
  save_plot(p, file.path(alpha_dir, "figures",
    "alpha_boxplot_by-marker_all-methods.pdf"))
})

# -- Test 2b: COI only, boxplot by lake
results[[length(results) + 1]] <- run_test("alpha: coi by lake", {
  ps_filt <- filter_ps_list(ps_coi, NULL, NULL)
  alpha_all <- compute_alpha_all(ps_filt, alpha_metrics, lake_order)
  alpha_long <- alpha_all %>%
    pivot_longer(cols = all_of(alpha_metrics), names_to = "Metric", values_to = "Value")
  p <- ggplot(alpha_long, aes(x = Lake, y = Value, fill = Lake)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    facet_grid(Metric ~ Marker, scales = "free_y") +
    scale_fill_manual(values = lake_colors) +
    theme_minimal(base_size = 10) +
    labs(title = "Alpha diversity - COI only, by lake")
  save_plot(p, file.path(alpha_dir, "figures",
    "alpha_boxplot_by-lake_coi-only.pdf"), width = 12, height = 8)
})

# -- Test 2c: only Rotifera, all markers, by lake faceted by marker
results[[length(results) + 1]] <- run_test("alpha: rotifera-only by lake facet marker", {
  tsub <- list(rank = "Phylum", include = "Rotifera")
  ps_filt <- filter_ps_list(ps_markers, NULL, NULL, tsub = tsub)
  # Drop any markers with 0 taxa
  ps_filt <- ps_filt[sapply(ps_filt, ntaxa) > 0]
  if (length(ps_filt) == 0) stop("no Rotifera in any marker")
  alpha_all <- compute_alpha_all(ps_filt, alpha_metrics, lake_order)
  alpha_long <- alpha_all %>%
    pivot_longer(cols = all_of(alpha_metrics), names_to = "Metric", values_to = "Value")
  p <- ggplot(alpha_long, aes(x = Lake, y = Value, fill = Lake)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    facet_grid(Metric ~ Marker, scales = "free_y") +
    scale_fill_manual(values = lake_colors) +
    theme_minimal(base_size = 10) +
    labs(title = "Alpha diversity - Rotifera only, by lake")
  save_plot(p, file.path(alpha_dir, "figures",
    "alpha_boxplot_by-lake_facet-marker_only-rotifera.pdf"), width = 14, height = 8)
})

# -- Test 2d: no rotifers, Superior+Michigan, KW stats
results[[length(results) + 1]] <- run_test("alpha: no-rotifera sup-mich stats", {
  tsub <- list(rank = "Phylum", exclude = "Rotifera")
  ps_filt <- filter_ps_list(ps_markers, NULL, c("Superior", "Michigan"), tsub = tsub)
  ps_filt <- ps_filt[sapply(ps_filt, ntaxa) > 0]
  alpha_all <- compute_alpha_all(ps_filt, alpha_metrics, lake_order)
  alpha_long <- alpha_all %>%
    pivot_longer(cols = all_of(alpha_metrics), names_to = "Metric", values_to = "Value")
  kw <- run_kruskal(alpha_long, "Lake", group_by_vars = "Marker")
  save_stats(kw, file.path(alpha_dir, "stats",
    "alpha_kw_lake_no-rotifera_superior-michigan"),
    caption = "KW: Lake effect, no Rotifera, Superior+Michigan")
  p <- ggplot(alpha_long, aes(x = Lake, y = Value, fill = Lake)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    facet_grid(Metric ~ Marker, scales = "free_y") +
    scale_fill_manual(values = lake_colors) +
    theme_minimal(base_size = 10) +
    labs(title = "Alpha - no Rotifera, Superior+Michigan")
  save_plot(p, file.path(alpha_dir, "figures",
    "alpha_boxplot_by-lake_no-rotifera_superior-michigan.pdf"), width = 14, height = 8)
})

# ============================================================
# 3. BETA DIVERSITY: multiple parameter combos
# ============================================================

cat("\n========== BETA DIVERSITY TESTS ==========\n")

beta_dir <- file.path(test_output, "beta")
dir.create(file.path(beta_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(beta_dir, "stats"), showWarnings = FALSE, recursive = TRUE)

# -- Test 3a: Jaccard NMDS, all markers, color by lake
results[[length(results) + 1]] <- run_test("beta: jaccard nmds all by lake", {
  ps_filt <- filter_ps_list(ps_all_methods, NULL, NULL)
  ord_plots <- imap(ps_filt, ~ run_ordination(
    .x, "NMDS", "jaccard", TRUE, "Lake",
    title = paste(.y, "- NMDS Jaccard"),
    color_palette = lake_colors))
  p <- wrap_plots(map(ord_plots, "plot"), ncol = 2) +
    plot_annotation(title = "NMDS Jaccard - all markers by lake")
  save_plot(p, file.path(beta_dir, "figures",
    "nmds_jaccard_by-lake_all-methods.pdf"), width = 16, height = 14)
})

# -- Test 3b: Bray NMDS, COI only, color by lake
results[[length(results) + 1]] <- run_test("beta: bray nmds coi by lake", {
  ps_filt <- filter_ps_list(ps_coi, NULL, NULL)
  ord_plots <- imap(ps_filt, ~ run_ordination(
    .x, "NMDS", "bray", FALSE, "Lake",
    title = paste(.y, "- NMDS Bray"),
    color_palette = lake_colors))
  p <- wrap_plots(map(ord_plots, "plot"), ncol = 2) +
    plot_annotation(title = "NMDS Bray-Curtis - COI only by lake")
  save_plot(p, file.path(beta_dir, "figures",
    "nmds_bray_by-lake_coi-only.pdf"), width = 14, height = 8)
})

# -- Test 3c: Jaccard NMDS, no rotifers, all markers
results[[length(results) + 1]] <- run_test("beta: jaccard no-rotifera all markers", {
  tsub <- list(rank = "Phylum", exclude = "Rotifera")
  ps_filt <- filter_ps_list(ps_markers, NULL, NULL, tsub = tsub)
  ps_filt <- ps_filt[sapply(ps_filt, ntaxa) > 0]
  ord_plots <- imap(ps_filt, ~ run_ordination(
    .x, "NMDS", "jaccard", TRUE, "Lake",
    title = paste(.y, "- NMDS Jaccard (no Rotifera)"),
    color_palette = lake_colors))
  p <- wrap_plots(map(ord_plots, "plot"), ncol = 2) +
    plot_annotation(title = "NMDS Jaccard - no Rotifera")
  save_plot(p, file.path(beta_dir, "figures",
    "nmds_jaccard_by-lake_markers_no-rotifera.pdf"), width = 16, height = 10)
})

# -- Test 3d: PERMANOVA for Erie+Ontario, Jaccard
results[[length(results) + 1]] <- run_test("beta: permanova erie-ontario jaccard", {
  ps_filt <- filter_ps_list(ps_markers, NULL, c("Erie", "Ontario"))
  perm_results <- imap_dfr(ps_filt, function(ps, m) {
    if (nsamples(ps) < 4) return(tibble())
    res <- run_permanova(ps, "~ Lake", "jaccard", TRUE)
    as.data.frame(res) %>% rownames_to_column("Term") %>% mutate(Marker = m)
  })
  perm_clean <- perm_results %>% filter(!is.na(`Pr(>F)`))
  save_stats(perm_clean, file.path(beta_dir, "stats",
    "permanova_jaccard_erie-ontario"),
    caption = "PERMANOVA Jaccard - Erie vs Ontario")
})

# ============================================================
# 4. HEATMAPS: multiple combos
# ============================================================

cat("\n========== HEATMAP TESTS ==========\n")

heat_dir <- file.path(test_output, "heatmaps", "figures")
dir.create(heat_dir, showWarnings = FALSE, recursive = TRUE)

# -- Test 4a: Genus, all methods
results[[length(results) + 1]] <- run_test("heatmap: genus all methods", {
  ps_filt <- filter_ps_list(ps_all_methods, NULL, NULL)
  for (m in names(ps_filt)) {
    fname <- paste0("heatmap_genus_", tolower(m), ".pdf")
    pdf(file.path(heat_dir, fname), width = 14, height = 10)
    plot_top_heatmap(ps_filt[[m]], rank = "Genus", top_n = 30,
                     annotation_var = "Lake",
                     lake_colors = lake_colors, lake_order = lake_order)
    dev.off()
    cat("  saved:", fname, "\n")
  }
})

# -- Test 4b: Order, no rotifers, COI only
results[[length(results) + 1]] <- run_test("heatmap: order no-rotifera coi", {
  tsub <- list(rank = "Phylum", exclude = "Rotifera")
  ps_filt <- filter_ps_list(ps_coi, NULL, NULL, tsub = tsub)
  ps_filt <- ps_filt[sapply(ps_filt, ntaxa) > 0]
  for (m in names(ps_filt)) {
    fname <- paste0("heatmap_order_", tolower(m), "_no-rotifera.pdf")
    pdf(file.path(heat_dir, fname), width = 14, height = 10)
    plot_top_heatmap(ps_filt[[m]], rank = "Order", top_n = 20,
                     annotation_var = "Lake",
                     lake_colors = lake_colors, lake_order = lake_order)
    dev.off()
    cat("  saved:", fname, "\n")
  }
})

# -- Test 4c: Genus, Erie only
results[[length(results) + 1]] <- run_test("heatmap: genus erie-only folmer", {
  ps_filt <- filter_ps_list(list(Folmer = ps_folmer), NULL, c("Erie"))
  for (m in names(ps_filt)) {
    fname <- paste0("heatmap_genus_", tolower(m), "_erie-only.pdf")
    pdf(file.path(heat_dir, fname), width = 14, height = 10)
    plot_top_heatmap(ps_filt[[m]], rank = "Genus", top_n = 25,
                     annotation_var = "Lake",
                     lake_colors = lake_colors, lake_order = lake_order)
    dev.off()
    cat("  saved:", fname, "\n")
  }
})

# ============================================================
# 5. OVERLAP: different subsets
# ============================================================

cat("\n========== OVERLAP TESTS ==========\n")

olap_dir <- file.path(test_output, "overlap")
dir.create(file.path(olap_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(olap_dir, "stats"), showWarnings = FALSE, recursive = TRUE)

# -- Test 5a: Species, all methods
results[[length(results) + 1]] <- run_test("overlap: species all methods", {
  ps_filt <- filter_ps_list(ps_all_methods, NULL, NULL)
  rank <- "Species"
  taxa_sets <- imap(ps_filt, ~ get_taxa_set(.x, rank))
  taxa_sets <- taxa_sets[lengths(taxa_sets) > 0]
  if (length(taxa_sets) < 2) stop("too few methods")
  # Upset plot
  upset_mat <- make_upset_matrix(taxa_sets)
  pdf(file.path(olap_dir, "figures", "upset_species_all-methods.pdf"), width = 10, height = 6)
  print(upset(upset_mat, sets = names(taxa_sets), order.by = "freq", text.scale = 1.3))
  dev.off()
  cat("  saved upset\n")
})

# -- Test 5b: Genus, COI only
results[[length(results) + 1]] <- run_test("overlap: genus coi-only", {
  ps_filt <- filter_ps_list(ps_coi, NULL, NULL)
  rank <- "Genus"
  taxa_sets <- imap(ps_filt, ~ get_taxa_set(.x, rank))
  taxa_sets <- taxa_sets[lengths(taxa_sets) > 0]
  if (length(taxa_sets) < 2) stop("too few methods")
  vc <- unname(marker_colors[names(taxa_sets)])
  venn_obj <- venn.diagram(
    x = taxa_sets, category.names = names(taxa_sets),
    filename = NULL, output = TRUE, imagetype = "none",
    col = vc, fill = adjustcolor(vc, alpha.f = 0.3),
    cat.col = vc, cat.cex = 1.2, margin = 0.1,
    main = paste(rank, "overlap - COI only"))
  pdf(file.path(olap_dir, "figures", "venn_genus_coi-only.pdf"), width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(venn_obj)
  dev.off()
  cat("  saved venn\n")
})

# -- Test 5c: Species, no rotifers
results[[length(results) + 1]] <- run_test("overlap: species no-rotifera", {
  tsub <- list(rank = "Phylum", exclude = "Rotifera")
  ps_filt <- filter_ps_list(ps_all_methods, NULL, NULL, tsub = tsub)
  ps_filt <- ps_filt[sapply(ps_filt, ntaxa) > 0]
  rank <- "Species"
  taxa_sets <- imap(ps_filt, ~ get_taxa_set(.x, rank))
  taxa_sets <- taxa_sets[lengths(taxa_sets) > 0]
  if (length(taxa_sets) < 2) stop("too few methods with species")
  upset_mat <- make_upset_matrix(taxa_sets)
  pdf(file.path(olap_dir, "figures", "upset_species_no-rotifera.pdf"), width = 10, height = 6)
  print(upset(upset_mat, sets = names(taxa_sets), order.by = "freq", text.scale = 1.3))
  dev.off()
  cat("  saved upset\n")
})

# ============================================================
# 6. EXPLORATORY: basic run
# ============================================================

cat("\n========== EXPLORATORY TEST ==========\n")

expl_dir <- file.path(test_output, "exploratory", "stats")
dir.create(expl_dir, showWarnings = FALSE, recursive = TRUE)

results[[length(results) + 1]] <- run_test("exploratory: all datasets", {
  datasets <- c(ps_all_methods,
                list(`Combined markers` = ps_markers_combined,
                     `Combined all` = ps_all_combined))
  summary_df <- imap_dfr(datasets, ~ summarise_ps(.x, .y, tax_ranks))
  summary_df <- summary_df %>%
    select(dataset, samples, total_asvs, rank, unique_taxa,
           unassigned_asvs, pct_reads_unassigned)
  save_stats(summary_df, file.path(expl_dir, "dataset_summary_all"),
             caption = "Dataset summary - all methods")
})

# ============================================================
# 7. RAREFACTION: skipped (slow)
# ============================================================

cat("\n========== RAREFACTION TESTS (SKIPPED) ==========\n")

# ============================================================
# 8. FOCAL TAXON: multiple groups
# ============================================================

cat("\n========== FOCAL TAXON TESTS ==========\n")

focal_dir <- file.path(test_output, "focal_taxon")
dir.create(file.path(focal_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(focal_dir, "stats"), showWarnings = FALSE, recursive = TRUE)

# -- Test 8a: Calanoida default
results[[length(results) + 1]] <- run_test("focal: calanoida all methods", {
  focal_summary <- focal_detection_summary(ps_all_methods, "Order", "Calanoida", "Species")
  save_stats(focal_summary, file.path(focal_dir, "stats", "calanoida_detection_all-methods"),
             caption = "Calanoida detection by method")
  # Richness by lake
  focal_lake <- imap_dfr(ps_all_methods, function(ps, m) {
    ps_sub <- subset_focal_taxon(ps, "Order", "Calanoida")
    if (is.null(ps_sub)) return(tibble())
    ps_agg <- agg_rank(ps_sub, "Species") %>% to_pa()
    sd <- data.frame(sample_data(ps_agg))
    sd$richness <- sample_sums(ps_agg)
    sd$Marker <- m
    sd
  })
  if (nrow(focal_lake) > 0) {
    focal_lake$Lake <- factor(focal_lake$Lake, levels = lake_order, ordered = TRUE)
    p <- ggplot(focal_lake, aes(Lake, richness, fill = Lake)) +
      geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.4) +
      facet_wrap(~ Marker) +
      scale_fill_manual(values = lake_colors) +
      theme_minimal(base_size = 10) +
      labs(title = "Calanoida species richness by lake")
    save_plot(p, file.path(focal_dir, "figures",
      "focal_calanoida_richness_by-lake_all-methods.pdf"))
  }
})

# -- Test 8b: Rotifera focal
results[[length(results) + 1]] <- run_test("focal: rotifera all methods", {
  focal_summary <- focal_detection_summary(ps_all_methods, "Phylum", "Rotifera", "Species")
  save_stats(focal_summary, file.path(focal_dir, "stats", "rotifera_detection_all-methods"),
             caption = "Rotifera detection by method")
  focal_lake <- imap_dfr(ps_all_methods, function(ps, m) {
    ps_sub <- subset_focal_taxon(ps, "Phylum", "Rotifera")
    if (is.null(ps_sub)) return(tibble())
    ps_agg <- agg_rank(ps_sub, "Species") %>% to_pa()
    sd <- data.frame(sample_data(ps_agg))
    sd$richness <- sample_sums(ps_agg)
    sd$Marker <- m
    sd
  })
  if (nrow(focal_lake) > 0) {
    focal_lake$Lake <- factor(focal_lake$Lake, levels = lake_order, ordered = TRUE)
    p <- ggplot(focal_lake, aes(Lake, richness, fill = Lake)) +
      geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.4) +
      facet_wrap(~ Marker) +
      scale_fill_manual(values = lake_colors) +
      theme_minimal(base_size = 10) +
      labs(title = "Rotifera species richness by lake")
    save_plot(p, file.path(focal_dir, "figures",
      "focal_rotifera_richness_by-lake_all-methods.pdf"))
  }
})

# ============================================================
# 9. DIFFERENTIAL: COI subset
# ============================================================

cat("\n========== DIFFERENTIAL TESTS ==========\n")

diff_dir <- file.path(test_output, "differential")
dir.create(file.path(diff_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(diff_dir, "stats"), showWarnings = FALSE, recursive = TRUE)

# -- Test 9a: SIMPER, COI only, genus, lake
results[[length(results) + 1]] <- run_test("diff: simper coi genus lake", {
  ps_filt <- filter_ps_list(ps_coi, NULL, NULL)
  simper_all <- list()
  for (m in names(ps_filt)) {
    res <- tryCatch(
      run_simper_analysis(ps_filt[[m]], "Lake", "Genus", 15, NULL),
      error = function(e) { cat("  simper error:", e$message, "\n"); NULL })
    if (!is.null(res)) {
      top <- res$top %>%
        group_by(taxon) %>%
        summarise(mean_contrib = mean(average, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_contrib)) %>% head(15)
      top$Marker <- m
      simper_all[[m]] <- top
    }
  }
  if (length(simper_all) > 0)
    save_stats(bind_rows(simper_all), file.path(diff_dir, "stats",
      "simper_genus_lake_coi-only"),
      caption = "SIMPER top genera - COI only")
})

# -- Test 9b: IndVal, all markers, genus, lake
results[[length(results) + 1]] <- run_test("diff: indval markers genus lake", {
  ps_filt <- filter_ps_list(ps_markers, NULL, NULL)
  indval_all <- list()
  for (m in names(ps_filt)) {
    res <- tryCatch(
      run_indicator(ps_filt[[m]], "Lake", "Genus", NULL),
      error = function(e) { cat("  indval error:", e$message, "\n"); NULL })
    if (!is.null(res)) {
      sig <- res$sign
      sig <- sig[sig$p.value < 0.05, ]
      if (nrow(sig) > 0) {
        sig$taxon <- rownames(sig)
        sig$Marker <- m
        indval_all[[m]] <- sig %>% arrange(p.value) %>% head(20)
      }
    }
  }
  if (length(indval_all) > 0)
    save_stats(bind_rows(indval_all), file.path(diff_dir, "stats",
      "indval_genus_lake_markers"),
      caption = "IndVal significant genera - all markers")
})

# ============================================================
# 10. GEOGRAPHIC: basic run
# ============================================================

cat("\n========== GEOGRAPHIC TEST ==========\n")

geo_dir <- file.path(test_output, "geographic")
dir.create(file.path(geo_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(geo_dir, "stats"), showWarnings = FALSE, recursive = TRUE)

results[[length(results) + 1]] <- run_test("geo: map + ibd", {
  meta_geo <- data.frame(sample_data(ps_morph_by_station))
  meta_geo$Sample_ID <- rownames(meta_geo)
  meta_geo <- filter_geo_metadata(meta_geo, c("Latitude", "Longitude"))
  meta_geo$Lake <- factor(meta_geo$Lake, levels = lake_order, ordered = TRUE)

  shared_geo <- intersect(meta_geo$Sample_ID, sample_names(ps_all_combined))
  alpha_geo <- compute_alpha(prune_samples(shared_geo, ps_all_combined),
                              "Combined", alpha_metrics, lake_order)
  meta_geo <- left_join(meta_geo,
                         alpha_geo %>% select(Sample_ID, Observed, InvSimpson),
                         by = "Sample_ID")

  basemap <- ne_states(country = c("United States of America", "Canada"),
                       returnclass = "sf")

  p_map <- ggplot() +
    geom_sf(data = basemap, fill = "grey92", color = "grey60") +
    geom_point(data = meta_geo, aes(Longitude, Latitude, color = Lake, size = Observed),
               alpha = 0.85) +
    scale_color_manual(values = lake_colors) +
    coord_sf(xlim = c(-93, -75), ylim = c(41, 49)) +
    theme_minimal(base_size = 10) +
    labs(title = "Sampling sites - combined richness")
  save_plot(p_map, file.path(geo_dir, "figures", "map_richness_combined.pdf"))

  # IBD
  ps_geo <- prune_samples(meta_geo$Sample_ID, ps_all_combined) %>% to_pa()
  ibd_res <- map(lake_order, ~ tryCatch(
    run_ibd(ps_geo, meta_geo, .x),
    error = function(e) NULL))
  names(ibd_res) <- lake_order
  mantel_tbl <- imap_dfr(ibd_res, function(r, lk) {
    if (is.null(r)) return(tibble(Lake = lk, rho = NA, p = NA))
    tibble(Lake = lk, rho = r$r_stat, p = r$p_val)
  }) %>% mutate(sig = sig_stars(p))
  save_stats(mantel_tbl, file.path(geo_dir, "stats", "mantel_ibd_combined"),
             caption = "Mantel IBD - combined data")
})

# ============================================================
# 11. VARPART: basic run
# ============================================================

cat("\n========== VARPART TEST ==========\n")

vp_dir <- file.path(test_output, "varpart")
dir.create(file.path(vp_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(vp_dir, "stats"), showWarnings = FALSE, recursive = TRUE)

results[[length(results) + 1]] <- run_test("varpart: markers lake+mesh vs lat/lon", {
  ps_filt <- filter_ps_list(ps_markers, NULL, NULL)
  vp_results <- list()
  for (m in names(ps_filt)) {
    vp <- tryCatch(
      run_varpart(ps_filt[[m]], c("Lake", "Mesh"), c("Latitude", "Longitude"), TRUE),
      error = function(e) { cat("  vp error:", e$message, "\n"); NULL })
    if (!is.null(vp)) {
      vp_results[[m]] <- vp
      pdf(file.path(vp_dir, "figures", paste0("varpart_", tolower(m), ".pdf")), width = 8, height = 6)
      plot(vp, bg = c("steelblue", "tomato"), Xnames = c("Lake+Mesh", "Lat/Lon"))
      title(main = paste(m, "- variance partitioning"))
      dev.off()
    }
  }
  if (length(vp_results) > 0) {
    vp_table <- imap_dfr(vp_results, function(vp, m) {
      fracs <- vp$part$indfract
      tibble(Dataset = m, Env_only = fracs$Adj.R.squared[1],
             Shared = fracs$Adj.R.squared[2],
             Spatial_only = fracs$Adj.R.squared[3],
             Residual = fracs$Adj.R.squared[4])
    })
    save_stats(vp_table, file.path(vp_dir, "stats", "varpart_fractions_markers"),
               caption = "Varpart fractions - markers")
  }
})

# ============================================================
# 12. GROUND TRUTH: basic run
# ============================================================

cat("\n========== GROUND TRUTH TEST ==========\n")

gt_dir <- file.path(test_output, "ground_truth", "stats")
dir.create(gt_dir, showWarnings = FALSE, recursive = TRUE)

results[[length(results) + 1]] <- run_test("ground_truth: basic", {
  trebitz_file <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/data/trebitz_lists/Trebitz_Zoops_2026_overall_taxfixed.csv"
  if (!file.exists(trebitz_file)) {
    # Try alternate paths
    alt <- "/Volumes/Samsung_1TB/Zooplankton/Metagenomics/data/trebitz_lists/Trebitz_Zoops_2026_overall.txt"
    if (file.exists(alt)) trebitz_file <- alt
    else stop("trebitz file not found")
  }
  raw <- read.csv(trebitz_file, stringsAsFactors = FALSE, check.names = FALSE)
  raw <- raw[, colnames(raw) != "" & !is.na(colnames(raw)), drop = FALSE]
  trebitz_overall <- raw %>% select(Species) %>% distinct()

  unexpected_all <- list()
  for (m in names(ps_all_methods)) {
    unexpected <- compare_to_ground_truth(ps_all_methods[[m]], trebitz_overall, "Species")
    if (nrow(unexpected) > 0) {
      unexpected$Marker <- m
      unexpected_all[[m]] <- unexpected
    }
  }
  if (length(unexpected_all) > 0) {
    df <- bind_rows(unexpected_all)
    save_stats(df, file.path(gt_dir, "unexpected_overall_all-methods"),
               caption = "Unexpected detections vs Trebitz overall list")
  }
})

# ============================================================
# RESULTS SUMMARY
# ============================================================

cat("\n\n============================\n")
cat("  TEST RESULTS SUMMARY\n")
cat("============================\n\n")

results_df <- bind_rows(results)
results_df$status <- ifelse(grepl("^PASS", results_df$result), "PASS", "FAIL")

for (i in seq_len(nrow(results_df))) {
  cat(sprintf("  [%s] %s (%.1fs)\n",
              results_df$status[i], results_df$test[i], results_df$seconds[i]))
  if (results_df$status[i] == "FAIL")
    cat("        ", results_df$result[i], "\n")
}

cat(sprintf("\n  %d/%d tests passed\n",
            sum(results_df$status == "PASS"), nrow(results_df)))
cat(sprintf("  Total time: %.0f seconds\n", sum(results_df$seconds)))

# Save results table
save_stats(results_df, file.path(test_output, "test_results"),
           caption = "Test run results")

cat("\n  Output in:", test_output, "\n")
