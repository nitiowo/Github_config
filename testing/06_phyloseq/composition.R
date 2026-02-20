# composition.R
# stacked barplots of relative abundance

source("setup.R")

outdir <- file.path(output_root, "composition")

# ---- control block ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL
use_ranks   <- c("Phylum", "Order", "Genus")
use_tsub    <- NULL
facet_var   <- NULL
top_n       <- 15

# ---- filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes, tsub = use_tsub)

# ---- generate one plot per rank x marker ----
for (rank in use_ranks) {
  for (m in names(ps_filt)) {
    ps <- ps_filt[[m]]
    ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(use_tsub)
    ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))

    df <- psmelt(ps_rel)
    df$Lake <- factor(df$Lake, levels = lake_order, ordered = TRUE)

    # top N taxa
    top_taxa <- df %>%
      group_by(.data[[rank]]) %>%
      summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = top_n) %>%
      pull(.data[[rank]])

    df[[rank]] <- if_else(df[[rank]] %in% top_taxa,
                          as.character(df[[rank]]), "Other")
    df[[rank]] <- factor(df[[rank]], levels = c(as.character(top_taxa), "Other"))

    fill_pal <- c(colorRampPalette(brewer.pal(12, "Set3"))(top_n), "grey70")
    names(fill_pal) <- levels(df[[rank]])

    p <- ggplot(df, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
      geom_bar(stat = "identity", position = "stack", width = 1) +
      scale_fill_manual(values = fill_pal) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
            legend.position = "right") +
      labs(title = paste(m, "-", rank), y = "Relative abundance", x = "")

    if (!is.null(facet_var))
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free_x")

    fname <- paste0("barplot_", tolower(rank), "_", tolower(m), ".pdf")
    save_plot(p, file.path(outdir, "figures", fname))
  }
}