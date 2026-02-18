# ground_truth.R
# compare detected taxa against trebitz ground-truth lists

source("setup.R")

outdir <- file.path(output_root, "ground_truth")

# ---- control block ----
use_ps_list <- ps_all_methods
use_rank    <- "Species"

# placeholder paths -- update when CSVs are ready
trebitz_files <- list(
  Superior = "data/trebitz_taxa_superior.csv",
  Michigan = "data/trebitz_taxa_michigan.csv",
  Huron    = "data/trebitz_taxa_huron.csv",
  Erie     = "data/trebitz_taxa_erie.csv",
  Ontario  = "data/trebitz_taxa_ontario.csv",
  Overall  = "data/trebitz_taxa_overall.csv"
)

# ---- load (skip missing files) ----
trebitz <- list()
for (nm in names(trebitz_files)) {
  fp <- trebitz_files[[nm]]
  if (file.exists(fp)) {
    trebitz[[nm]] <- load_trebitz(fp)
    cat("loaded", nm, ":", nrow(trebitz[[nm]]), "taxa\n")
  } else {
    cat("not found:", fp, "- skipping\n")
  }
}

if (length(trebitz) == 0) {
  cat("\nno trebitz files found. create the CSVs and update paths above.\n")
} else {

  # ============================================================
  # per-lake: unexpected detections
  # ============================================================

  unexpected_all <- list()
  for (lk in intersect(lake_order, names(trebitz))) {
    for (m in names(use_ps_list)) {
      unexpected <- compare_to_ground_truth(
        use_ps_list[[m]], trebitz[[lk]], rank = use_rank, lake = lk)
      if (nrow(unexpected) > 0) {
        unexpected$Lake <- lk
        unexpected$Marker <- m
        unexpected_all[[paste(lk, m)]] <- unexpected
      }
    }
  }

  if (length(unexpected_all) > 0) {
    unexpected_df <- bind_rows(unexpected_all)
    save_stats(unexpected_df,
               file.path(outdir, "stats", "unexpected_detections_per_lake"),
               caption = "Taxa detected but absent from Trebitz ground-truth (per lake)")

    # summary count
    unexpected_summary <- unexpected_df %>%
      group_by(Lake, Marker) %>%
      summarise(n_unexpected = n(), .groups = "drop")
    save_stats(unexpected_summary,
               file.path(outdir, "stats", "unexpected_detections_summary"),
               caption = "Count of unexpected detections per lake and marker")
  }

  # ============================================================
  # overall: unexpected detections
  # ============================================================

  if ("Overall" %in% names(trebitz)) {
    overall_unexpected <- list()
    for (m in names(use_ps_list)) {
      unexpected <- compare_to_ground_truth(
        use_ps_list[[m]], trebitz[["Overall"]], rank = use_rank)
      if (nrow(unexpected) > 0) {
        unexpected$Marker <- m
        overall_unexpected[[m]] <- unexpected
      }
    }

    if (length(overall_unexpected) > 0) {
      overall_df <- bind_rows(overall_unexpected)
      save_stats(overall_df,
                 file.path(outdir, "stats", "unexpected_detections_overall"),
                 caption = "Taxa detected but absent from overall Trebitz list")
    }
  }

  # ============================================================
  # lake richness: our data vs trebitz
  # ============================================================

  if (all(lake_order %in% names(trebitz))) {
    trebitz_richness <- imap_dfr(trebitz[lake_order], function(df, lk) {
      tibble(Lake = lk,
             trebitz_species = length(unique(na.omit(df$Species))))
    })
    trebitz_richness$Lake <- factor(trebitz_richness$Lake,
                                     levels = lake_order, ordered = TRUE)

    alpha_comb_lake <- compute_alpha(ps_markers_combined, "Combined",
                                      "Observed", lake_order) %>%
      group_by(Lake) %>%
      summarise(mean_observed = mean(Observed, na.rm = TRUE), .groups = "drop")

    compare_df <- left_join(trebitz_richness, alpha_comb_lake, by = "Lake")
    save_stats(compare_df,
               file.path(outdir, "stats", "trebitz_vs_markers_lake_richness"),
               caption = "Lake richness: Trebitz vs combined markers")

    # correlation
    if (nrow(compare_df) >= 4) {
      ct <- cor.test(compare_df$trebitz_species, compare_df$mean_observed,
                     method = "spearman")
      cat(sprintf("richness correlation: rho = %.3f, p = %.4f\n",
                  ct$estimate, ct$p.value))
    }

    p_treb <- ggplot(compare_df, aes(trebitz_species, mean_observed, color = Lake)) +
      geom_point(size = 4) +
      geom_text_repel(aes(label = Lake)) +
      scale_color_manual(values = lake_colors) +
      theme_minimal(base_size = 10) +
      labs(title = "Lake richness: Trebitz vs our combined markers",
           x = "Trebitz species count", y = "Mean observed richness (markers)")
    save_plot(p_treb,
              file.path(outdir, "figures", "trebitz_richness_correlation.pdf"))
  }

  # ============================================================
  # text summary
  # ============================================================

  txt <- c("=== GROUND TRUTH COMPARISON ===",
           paste("Generated:", Sys.Date()), "")

  if (length(unexpected_all) > 0) {
    txt <- c(txt, "--- Unexpected detections per lake ---")
    for (lk in lake_order) {
      lk_rows <- bind_rows(unexpected_all[grepl(paste0("^", lk), names(unexpected_all))])
      if (nrow(lk_rows) > 0) {
        per_marker <- lk_rows %>% group_by(Marker) %>% summarise(n = n(), .groups = "drop")
        for (j in seq_len(nrow(per_marker))) {
          txt <- c(txt, sprintf("  %s - %s: %d unexpected species",
                                 lk, per_marker$Marker[j], per_marker$n[j]))
        }
      }
    }
  }

  if ("Overall" %in% names(trebitz) && length(overall_unexpected) > 0) {
    txt <- c(txt, "", "--- Unexpected detections (overall) ---")
    for (m in names(overall_unexpected)) {
      txt <- c(txt, sprintf("  %s: %d species not in overall list",
                             m, nrow(overall_unexpected[[m]])))
    }
  }

  save_summary(txt, file.path(outdir, "stats", "ground_truth_results_summary.txt"))
}