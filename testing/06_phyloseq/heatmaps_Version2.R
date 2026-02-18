# heatmaps.R
# top-taxa heatmaps per marker

source("setup.R")

outdir <- file.path(output_root, "heatmaps")

# ---- control block ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL
use_rank    <- "Genus"
use_tsub    <- NULL
top_n       <- 30

# ---- filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes, tsub = use_tsub)

# ---- one heatmap per marker ----
for (m in names(ps_filt)) {
  fname <- paste0("heatmap_", tolower(use_rank), "_", tolower(m), ".pdf")
  pdf(file.path(outdir, "figures", fname), width = 14, height = 10)
  plot_top_heatmap(ps_filt[[m]], rank = use_rank, top_n = top_n,
                   annotation_var = "Lake",
                   lake_colors = lake_colors, lake_order = lake_order,
                   tsub = use_tsub)
  dev.off()
  cat("saved:", fname, "\n")
}