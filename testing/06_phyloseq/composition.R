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

# ---- build filename suffix from current settings ----
comp_suffix <- ""
if (!is.null(use_tsub)) {
  if (!is.null(use_tsub$include))
    comp_suffix <- paste0(comp_suffix, "_only-", tolower(paste(use_tsub$include, collapse = "-")))
  if (!is.null(use_tsub$exclude))
    comp_suffix <- paste0(comp_suffix, "_no-", tolower(paste(use_tsub$exclude, collapse = "-")))
}
if (!is.null(use_lakes))
  comp_suffix <- paste0(comp_suffix, "_", tolower(paste(use_lakes, collapse = "-")))

# ---- generate barplots ----
# When facet_var = "Marker", combine all markers into one figure per rank
# When facet_var = "Lake" or NULL, produce one figure per rank x marker

for (rank in use_ranks) {

  if (!is.null(facet_var) && facet_var == "Marker") {
    # Combined plot: melt all markers into one data frame
    df_all <- imap_dfr(ps_filt, function(ps, m) {
      ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(use_tsub)
      ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))
      df <- psmelt(ps_rel)
      df$Marker <- m
      df
    })
    df_all$Lake <- factor(df_all$Lake, levels = lake_order, ordered = TRUE)
    df_all$Marker <- factor(df_all$Marker, levels = names(ps_filt))

    # Top N taxa across all markers combined
    top_taxa <- df_all %>%
      group_by(.data[[rank]]) %>%
      summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = top_n) %>%
      pull(.data[[rank]])

    df_all[[rank]] <- if_else(df_all[[rank]] %in% top_taxa,
                              as.character(df_all[[rank]]), "Other")
    df_all[[rank]] <- factor(df_all[[rank]],
                             levels = c(as.character(top_taxa), "Other"))

    fill_pal <- c(colorRampPalette(brewer.pal(12, "Set3"))(top_n), "grey70")
    names(fill_pal) <- levels(df_all[[rank]])

    p <- ggplot(df_all, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
      geom_bar(stat = "identity", position = "stack", width = 1) +
      facet_wrap(~ Marker, scales = "free_x") +
      scale_fill_manual(values = fill_pal) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
            legend.position = "right") +
      labs(title = paste("All markers -", rank, "(faceted by Marker)"),
           y = "Relative abundance", x = "")

    fname <- paste0("barplot_", tolower(rank), "_all_facet-marker",
                    comp_suffix, ".pdf")
    save_plot(p, file.path(outdir, "figures", fname), width = 16, height = 8)

  } else {
    # Per-marker plots (optionally faceted by Lake or other variable)
    for (m in names(ps_filt)) {
      ps <- ps_filt[[m]]
      ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(use_tsub)
      ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))

      df <- psmelt(ps_rel)
      df$Marker <- m
      df$Lake <- factor(df$Lake, levels = lake_order, ordered = TRUE)

      # Top N taxa
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

      facet_tag <- if (!is.null(facet_var)) paste0("_facet-", tolower(facet_var)) else ""
      fname <- paste0("barplot_", tolower(rank), "_", tolower(m),
                      facet_tag, comp_suffix, ".pdf")
      save_plot(p, file.path(outdir, "figures", fname))
    }
  }
}