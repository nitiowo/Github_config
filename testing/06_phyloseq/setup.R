# setup.R
# load all data, build phyloseq objects, define palettes
# source this at the top of every analysis script

source("zoop_functions.R")
set.seed(42)

# -- output root (change this to redirect all output) --
output_root <- "output"

# -- palettes and constants --
lake_order <- c("Superior", "Michigan", "Huron", "Erie", "Ontario")

lake_colors <- c(
  Superior = "#1b9e77", Michigan = "#d95f02", Huron = "#7570b3",
  Erie     = "#e7298a", Ontario  = "#66a61e"
)

marker_colors <- c(
  Folmer = "#e41a1c", Leray = "#377eb8",
  `18S`  = "#4daf4a", Morphology = "#984ea3"
)

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order",
               "Family", "Genus", "Species")

alpha_metrics <- c("Observed", "InvSimpson")

# -- load phyloseq objects --
ps_folmer <- readRDS("setup_output/folmer_ps.RDS")
ps_leray  <- readRDS("setup_output/leray_ps.RDS")
ps_18S    <- readRDS("setup_output/ssu_ps.RDS")
ps_morph  <- readRDS("setup_output/morph_ps.RDS")

# set lake factor ordering on everything
ps_folmer <- set_lake_order(ps_folmer, lake_order)
ps_leray  <- set_lake_order(ps_leray, lake_order)
ps_18S    <- set_lake_order(ps_18S, lake_order)
ps_morph  <- set_lake_order(ps_morph, lake_order)

# -- named lists --
ps_markers     <- list(Folmer = ps_folmer, Leray = ps_leray, `18S` = ps_18S)
ps_all_methods <- c(ps_markers, list(Morphology = ps_morph))
ps_coi         <- list(Folmer = ps_folmer, Leray = ps_leray)

# -- P/A versions --
ps_folmer_pa <- to_pa(ps_folmer)
ps_leray_pa  <- to_pa(ps_leray)
ps_18S_pa    <- to_pa(ps_18S)
ps_morph_pa  <- to_pa(ps_morph)

ps_markers_pa     <- list(Folmer = ps_folmer_pa, Leray = ps_leray_pa, `18S` = ps_18S_pa)
ps_all_methods_pa <- c(ps_markers_pa, list(Morphology = ps_morph_pa))
ps_coi_pa         <- list(Folmer = ps_folmer_pa, Leray = ps_leray_pa)

# -- station-aggregated for morph comparisons --
ps_markers_by_station <- lapply(ps_markers, aggregate_to_station)
ps_morph_by_station   <- rename_samples_to_station(ps_morph, "Station_ID")

# -- combined datasets --
ps_markers_combined <- combine_ps_pa(ps_markers, rank = "Species")
ps_coi_combined     <- combine_ps_pa(ps_coi, rank = "Species")

ps_station_all  <- c(ps_markers_by_station, list(Morphology = ps_morph_by_station))
ps_all_combined <- combine_ps_pa(ps_station_all, rank = "Species")

# -- morph info --
morph_sample_names <- sample_names(ps_morph)
morph_station_ids  <- data.frame(sample_data(ps_morph))$Station_ID
shared_stations    <- intersect(
  Reduce(intersect, lapply(ps_markers_by_station, sample_names)),
  sample_names(ps_morph_by_station)
)

# -- create output directories --
analysis_dirs <- c("exploratory", "composition", "rarefaction", "alpha",
                   "beta", "differential", "heatmaps", "overlap",
                   "geographic", "focal_taxon", "varpart", "ground_truth")
for (d in analysis_dirs) {
  dir.create(file.path(output_root, d, "stats"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_root, d, "figures"), recursive = TRUE, showWarnings = FALSE)
}

cat("Setup complete.\n")
cat("  Markers:", paste(names(ps_markers), collapse = ", "), "\n")
cat("  Samples per marker:", paste(sapply(ps_markers, nsamples), collapse = ", "), "\n")
cat("  Morph samples:", nsamples(ps_morph), "\n")
cat("  Shared stations (all methods):", length(shared_stations), "\n")
cat("  Combined markers:", ntaxa(ps_markers_combined), "taxa\n")
cat("  Combined all:", ntaxa(ps_all_combined), "taxa\n")
cat("  Output root:", output_root, "\n")