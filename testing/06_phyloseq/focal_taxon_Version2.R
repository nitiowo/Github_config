# focal_taxon.R
# deep dive on a single taxon group (default: Calanoida)

source("setup.R")

outdir <- file.path(output_root, "focal_taxon")

# ---- control block ----
focal_rank <- "Order"
focal_name <- "Calanoida"
use_ps_list <- ps_all_methods

# ---- detection summary ----
focal_summary <- focal_detection_summary(use_ps_list, focal_rank, focal_name, "Species")
cat("\n", focal_name, "detection:\n")
print(focal_summary)

save_stats(focal_summary,
           file.path(outdir, "stats", paste0(tolower(focal_name), "_detection_summary")),
           caption = paste(focal_name, "detection by method"))

# ---- distribution across lakes ----
focal_lake_data <- imap_dfr(use_ps_list, function(ps, m) {
  ps_sub <- subset_focal_taxon(ps, focal_rank, focal_name)
  if (is.null(ps_sub)) return(tibble())
  ps_agg <- agg_rank(ps_sub, "Species") %>% to_pa()
  sd <- data.frame(sample_data(ps_agg))
  sd$richness <- sample_sums(ps_agg)
  sd$Marker <- m
  sd$Sample_ID <- rownames(sd)
  sd
})

if (nrow(focal_lake_data) > 0) {
  focal_lake_data$Lake <- factor(focal_lake_data$Lake,
                                  levels = lake_order, ordered = TRUE)

  p_focal_lake <- ggplot(focal_lake_data, aes(Lake, richness, fill = Lake)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4) +
    facet_wrap(~ Marker) +
    scale_fill_manual(values = lake_colors) +
    theme_minimal(base_size = 10) +
    labs(title = paste(focal_name, "richness by lake and marker"),
         y = paste(focal_name, "species per sample"))
  save_plot(p_focal_lake,
            file.path(outdir, "figures",
                      paste0(tolower(focal_name), "_richness_by_lake.pdf")))
}

# ---- venn: which markers detect which species ----
focal_sets <- imap(use_ps_list, function(ps, m) {
  ps_sub <- subset_focal_taxon(ps, focal_rank, focal_name)
  if (is.null(ps_sub)) return(character(0))
  get_taxa_set(ps_sub, "Species")
})
focal_sets <- focal_sets[lengths(focal_sets) > 0]

if (length(focal_sets) >= 2) {
  fc <- unname(marker_colors[names(focal_sets)])
  venn_focal <- venn.diagram(
    x = focal_sets, category.names = names(focal_sets),
    filename = NULL, output = TRUE, imagetype = "none",
    col = fc, fill = adjustcolor(fc, alpha.f = 0.3),
    cat.col = fc, cat.cex = 1.2, margin = 0.1,
    main = paste(focal_name, "- species overlap"))

  venn_path <- file.path(outdir, "figures",
                          paste0(tolower(focal_name), "_venn.pdf"))
  pdf(venn_path, width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(venn_focal)
  dev.off()
  cat("saved:", venn_path, "\n")
}

# ---- composition within focal group ----
focal_ps_list <- list()
for (m in names(use_ps_list)) {
  ps_sub <- subset_focal_taxon(use_ps_list[[m]], focal_rank, focal_name)
  if (!is.null(ps_sub) && ntaxa(ps_sub) > 0 &&
      !is.null(access(ps_sub, "tax_table", errorIfNULL = FALSE))) {
    focal_ps_list[[m]] <- ps_sub
  }
}

if (length(focal_ps_list) > 0) {
  for (m in names(focal_ps_list)) {
    ps <- focal_ps_list[[m]]
    ps_agg <- agg_rank(ps, "Species")
    ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))
    df <- psmelt(ps_rel)
    df$Lake <- factor(df$Lake, levels = lake_order, ordered = TRUE)

    p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Species)) +
      geom_bar(stat = "identity", position = "stack", width = 1) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
            legend.position = "right") +
      labs(title = paste(focal_name, "-", m, "species composition"),
           y = "Relative abundance", x = "")

    save_plot(p, file.path(outdir, "figures",
                            paste0(tolower(focal_name), "_composition_",
                                   tolower(m), ".pdf")))
  }
}

# ---- text summary ----
txt <- c(
  paste("=== FOCAL TAXON:", toupper(focal_name), "==="),
  paste("Generated:", Sys.Date()), "")

txt <- c(txt, "--- Detection summary ---")
for (i in seq_len(nrow(focal_summary))) {
  txt <- c(txt, sprintf("  %s: %d species in %d samples",
                         focal_summary$Marker[i],
                         focal_summary$n_taxa[i],
                         focal_summary$n_samples_with[i]))
}

if (length(focal_sets) > 0) {
  all_focal <- unique(unlist(focal_sets))
  inter_focal <- Reduce(intersect, focal_sets)
  txt <- c(txt, "",
           sprintf("Total %s species (union): %d", focal_name, length(all_focal)),
           sprintf("Shared across all detecting methods: %d", length(inter_focal)))
}

save_summary(txt, file.path(outdir, "stats",
                             paste0(tolower(focal_name), "_results_summary.txt")))