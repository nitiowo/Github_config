# run_all.R
# runs all analysis scripts in order and saves all outputs
# change output_root here to redirect everything to a different directory

output_root <- "output"

cat("============================\n")
cat("  running full analysis\n")
cat("  output ->", output_root, "\n")
cat("============================\n\n")

source("setup.R")

scripts <- c(
  "exploratory.R",
  "composition.R",
  "rarefaction.R",
  "alpha.R",
  "beta.R",
  "differential.R",
  "heatmaps.R",
  "overlap.R",
  "geographic.R",
  "focal_taxon.R",
  "varpart.R",
  "ground_truth.R"
)

for (s in scripts) {
  if (file.exists(s)) {
    cat("\n--- running:", s, "---\n")
    tryCatch(source(s), error = function(e) cat("  ERROR:", e$message, "\n"))
  } else {
    cat("  skipping (not found):", s, "\n")
  }
}

cat("\n============================\n")
cat("  done. output in:", output_root, "\n")
cat("============================\n")