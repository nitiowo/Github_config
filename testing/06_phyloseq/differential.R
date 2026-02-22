# differential.R
# ANCOM-BC2, SIMPER, indicator species analysis

source("setup.R")

outdir <- file.path(output_root, "differential")

# ---- control block ----
use_ps_list  <- ps_markers
use_markers  <- NULL
use_lakes    <- NULL
use_rank     <- "Genus"
compare_var  <- "Lake"
use_tsub     <- NULL
ancom_prev   <- 0.10
simper_top_n <- 15

# ---- filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes, tsub = use_tsub)

# ============================================================
# ANCOM-BC2
# ============================================================

ancom_all_sig <- list()
for (m in names(ps_filt)) {
  cat("ANCOM-BC2:", m, "\n")
  res <- tryCatch(
    run_ancombc(ps_filt[[m]], compare_var, use_rank, ancom_prev),
    error = function(e) { cat("  error:", e$message, "\n"); NULL })

  if (!is.null(res)) {
    df <- res$res
    diff_cols <- grep("^diff_", colnames(df), value = TRUE)
    sig <- df %>% filter(if_any(all_of(diff_cols), ~ . == TRUE))
    sig$Marker <- m
    ancom_all_sig[[m]] <- sig
    cat("  ", nrow(sig), "significant taxa\n")

    # volcano for first contrast
    lfc_cols <- grep("^lfc_", colnames(df), value = TRUE)
    q_cols   <- grep("^q_",   colnames(df), value = TRUE)

    if (length(lfc_cols) > 0 && length(q_cols) > 0) {
      vdf <- tibble(taxon = df$taxon,
                     lfc = df[[lfc_cols[1]]],
                     q   = df[[q_cols[1]]]) %>%
        mutate(neg_log_q = -log10(q), is_sig = q < 0.05)

      p <- ggplot(vdf, aes(lfc, neg_log_q, color = is_sig)) +
        geom_point(alpha = 0.7) +
        geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
        geom_text_repel(data = filter(vdf, is_sig), aes(label = taxon),
                        size = 2.5, max.overlaps = 15) +
        scale_color_manual(values = c("grey60", "red")) +
        theme_minimal(base_size = 10) +
        labs(title = paste(m, "-", lfc_cols[1]),
             x = "Log fold change", y = "-log10(q)")
      save_plot(p, file.path(outdir, "figures",
                              paste0("volcano_", tolower(m), "_", tolower(compare_var), ".pdf")))
    }
  }
}

if (length(ancom_all_sig) > 0) {
  ancom_combined <- bind_rows(ancom_all_sig)
  save_stats(ancom_combined,
             file.path(outdir, "stats", "ancombc_significant_taxa"),
             caption = "ANCOM-BC2 significant taxa")
}

# ============================================================
# SIMPER
# ============================================================

simper_all <- list()
for (m in names(ps_filt)) {
  cat("SIMPER:", m, "\n")
  res <- tryCatch(
    run_simper_analysis(ps_filt[[m]], compare_var, use_rank, simper_top_n, use_tsub),
    error = function(e) { cat("  error:", e$message, "\n"); NULL })

  if (!is.null(res)) {
    top <- res$top %>%
      group_by(taxon) %>%
      summarise(mean_contrib = mean(average, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_contrib)) %>%
      head(simper_top_n)
    top$Marker <- m
    simper_all[[m]] <- top
  }
}

if (length(simper_all) > 0) {
  simper_combined <- bind_rows(simper_all)
  save_stats(simper_combined,
             file.path(outdir, "stats", "simper_top_genera"),
             caption = "SIMPER top contributing genera")
}

# ============================================================
# indicator species
# ============================================================

indval_all <- list()
for (m in names(ps_filt)) {
  cat("IndVal:", m, "\n")
  res <- tryCatch(
    run_indicator(ps_filt[[m]], compare_var, use_rank, use_tsub),
    error = function(e) { cat("  error:", e$message, "\n"); NULL })

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

if (length(indval_all) > 0) {
  indval_combined <- bind_rows(indval_all)
  save_stats(indval_combined,
             file.path(outdir, "stats", "indval_significant"),
             caption = "Significant indicator species (IndVal)")
}

# ============================================================
# text summary
# ============================================================

txt <- c(
  "=== DIFFERENTIAL ABUNDANCE RESULTS ===",
  paste("Generated:", Sys.Date()), "")

txt <- c(txt, "--- ANCOM-BC2 ---")
for (m in names(ancom_all_sig)) {
  txt <- c(txt, sprintf("  %s: %d significant genera", m, nrow(ancom_all_sig[[m]])))
}

txt <- c(txt, "", "--- SIMPER (top contributors to between-lake dissimilarity) ---")
for (m in names(simper_all)) {
  top3 <- head(simper_all[[m]]$taxon, 3)
  txt <- c(txt, sprintf("  %s: %s", m, paste(top3, collapse = ", ")))
}

txt <- c(txt, "", "--- Indicator Species ---")
for (m in names(indval_all)) {
  txt <- c(txt, sprintf("  %s: %d significant indicators", m, nrow(indval_all[[m]])))
}

save_summary(txt, file.path(outdir, "stats", "differential_results_summary.txt"))