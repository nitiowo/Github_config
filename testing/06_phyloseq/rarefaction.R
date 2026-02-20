# rarefaction.R
# iNEXT rarefaction curves

source("setup.R")

outdir <- file.path(output_root, "rarefaction")

# ---- control block ----
use_ps_list <- ps_all_methods
use_markers <- NULL

# ---- filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers)

# ---- run iNEXT and plot ----
run_inext <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  iNEXT(as.list(data.frame(otu)), q = 0, datatype = "abundance", nboot = 50)
}

rare_plots <- list()
for (m in names(ps_filt)) {
  cat("  iNEXT:", m, "\n")
  res <- run_inext(ps_filt[[m]])
  rare_plots[[m]] <- ggiNEXT(res, type = 1) +
    theme_minimal(base_size = 10) +
    labs(title = paste(m, "- rarefaction")) +
    theme(legend.position = "none")
}

p_all <- wrap_plots(rare_plots, ncol = 2) +
  plot_annotation(title = "Rarefaction curves")
save_plot(p_all, file.path(outdir, "figures", "rarefaction_all.pdf"),
          width = 16, height = 10)