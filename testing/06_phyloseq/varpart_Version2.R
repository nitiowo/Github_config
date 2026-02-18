# varpart.R
# variance partitioning: Lake + Mesh vs spatial (lat/lon)

source("setup.R")

outdir <- file.path(output_root, "varpart")

# ---- control block ----
use_ps_list   <- ps_markers
use_markers   <- NULL
env_vars      <- c("Lake", "Mesh")
spatial_vars  <- c("Latitude", "Longitude")
use_binary    <- TRUE

# ---- filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers)

# ============================================================
# per-marker variance partitioning
# ============================================================

vp_results <- list()
for (m in names(ps_filt)) {
  cat("varpart:", m, "\n")
  vp <- tryCatch(
    run_varpart(ps_filt[[m]], env_vars, spatial_vars, use_binary),
    error = function(e) { cat("  error:", e$message, "\n"); NULL })

  if (!is.null(vp)) {
    vp_results[[m]] <- vp

    pdf(file.path(outdir, "figures", paste0("varpart_", tolower(m), ".pdf")),
        width = 8, height = 6)
    plot(vp, bg = c("steelblue", "tomato"),
         Xnames = c("Lake + Mesh", "Lat/Lon"))
    title(main = paste(m, "- variance partitioning"))
    dev.off()
    cat("  saved figure\n")
  }
}

# ============================================================
# combined data
# ============================================================

cat("varpart: combined data\n")
vp_comb <- tryCatch(
  run_varpart(ps_all_combined, env_vars, spatial_vars, use_binary),
  error = function(e) { cat("  error:", e$message, "\n"); NULL })

if (!is.null(vp_comb)) {
  vp_results[["Combined"]] <- vp_comb

  pdf(file.path(outdir, "figures", "varpart_combined.pdf"), width = 8, height = 6)
  plot(vp_comb, bg = c("steelblue", "tomato"),
       Xnames = c("Lake + Mesh", "Lat/Lon"))
  title(main = "Combined data - variance partitioning")
  dev.off()
}

# ============================================================
# extract fractions into a table
# ============================================================

vp_table <- imap_dfr(vp_results, function(vp, m) {
  # varpart returns: [a] env only, [b] shared, [c] spatial only, [d] residual
  # these are in vp$part$indfract
  fracs <- vp$part$indfract
  tibble(
    Dataset = m,
    Env_only = fracs$Adj.R.squared[1],
    Shared = fracs$Adj.R.squared[2],
    Spatial_only = fracs$Adj.R.squared[3],
    Residual = fracs$Adj.R.squared[4]
  )
})

save_stats(vp_table,
           file.path(outdir, "stats", "varpart_fractions"),
           caption = "Variance partitioning: adjusted R-squared fractions")

# ============================================================
# text summary
# ============================================================

txt <- c(
  "=== VARIANCE PARTITIONING RESULTS ===",
  paste("Generated:", Sys.Date()),
  sprintf("Env variables: %s", paste(env_vars, collapse = ", ")),
  sprintf("Spatial variables: %s", paste(spatial_vars, collapse = ", ")),
  "")

for (i in seq_len(nrow(vp_table))) {
  r <- vp_table[i, ]
  txt <- c(txt, sprintf("--- %s ---", r$Dataset),
           sprintf("  Env only (adj R2): %.3f", r$Env_only),
           sprintf("  Shared: %.3f", r$Shared),
           sprintf("  Spatial only: %.3f", r$Spatial_only),
           sprintf("  Residual: %.3f", r$Residual), "")
}

save_summary(txt, file.path(outdir, "stats", "varpart_results_summary.txt"))