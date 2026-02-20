# exploratory.R
# dataset summaries: counts, unique taxa, % unassigned at each rank

source("setup.R")

outdir <- file.path(output_root, "exploratory")

# ---- what to summarize ----
datasets <- c(ps_all_methods,
              list(`Combined markers` = ps_markers_combined,
                   `Combined all` = ps_all_combined))

# ---- build summary table ----
summary_df <- imap_dfr(datasets, ~ summarise_ps(.x, .y, tax_ranks))

# reorder columns
summary_df <- summary_df %>%
  select(dataset, samples, total_asvs, rank, unique_taxa,
         unassigned_asvs, pct_reads_unassigned)

# ---- print to console ----
for (nm in names(datasets)) {
  summarise_ps_print(datasets[[nm]], nm, tax_ranks)
}

# ---- save ----
save_stats(summary_df,
           file.path(outdir, "stats", "dataset_summary"),
           caption = "Dataset summary: taxa counts and % reads unassigned")

# text summary
txt <- capture.output({
  for (nm in names(datasets))
    summarise_ps_print(datasets[[nm]], nm, tax_ranks)
})
save_summary(txt, file.path(outdir, "stats", "dataset_summary.txt"))