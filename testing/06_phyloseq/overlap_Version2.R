# overlap.R
# venn diagrams and upset plots of taxa overlap between methods

source("setup.R")

outdir <- file.path(output_root, "overlap")

# ---- control block ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_ranks   <- c("Species", "Genus", "Order")

# ---- filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers)

for (rank in use_ranks) {
  cat("overlap at", rank, "level\n")

  # get taxa sets
  taxa_sets <- imap(ps_filt, ~ get_taxa_set(.x, rank))

  # drop any methods with zero taxa at this rank
  taxa_sets <- taxa_sets[lengths(taxa_sets) > 0]
  if (length(taxa_sets) < 2) {
    cat("  fewer than 2 methods with data, skipping\n")
    next
  }

  # stats
  all_taxa <- unique(unlist(taxa_sets))
  intersect_all <- Reduce(intersect, taxa_sets)

  set_summary <- imap_dfr(taxa_sets, ~ tibble(
    Method = .y, n_taxa = length(.x),
    n_unique = length(setdiff(.x, unlist(taxa_sets[names(taxa_sets) != .y])))
  ))
  set_summary <- bind_rows(
    set_summary,
    tibble(Method = "Union (all)", n_taxa = length(all_taxa), n_unique = NA_integer_),
    tibble(Method = "Intersection (all)", n_taxa = length(intersect_all), n_unique = NA_integer_)
  )

  save_stats(set_summary,
             file.path(outdir, "stats", paste0(tolower(rank), "_set_sizes")),
             caption = paste(rank, "level set sizes"))

  # pairwise overlap matrix
  methods <- names(taxa_sets)
  pw_mat <- expand.grid(A = methods, B = methods, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(shared = length(intersect(taxa_sets[[A]], taxa_sets[[B]]))) %>%
    ungroup()
  save_stats(pw_mat,
             file.path(outdir, "stats", paste0(tolower(rank), "_pairwise_overlap")),
             caption = paste(rank, "level pairwise overlap"))

  # venn diagram
  vc <- unname(marker_colors[names(taxa_sets)])
  venn_obj <- venn.diagram(
    x = taxa_sets,
    category.names = names(taxa_sets),
    filename = NULL, output = TRUE, imagetype = "none",
    col = vc, fill = adjustcolor(vc, alpha.f = 0.3),
    cat.col = vc, cat.cex = 1.2, margin = 0.1,
    main = paste(rank, "overlap"))

  venn_path <- file.path(outdir, "figures", paste0("venn_", tolower(rank), ".pdf"))
  pdf(venn_path, width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(venn_obj)
  dev.off()
  cat("saved:", venn_path, "\n")

  # upset plot
  upset_mat <- make_upset_matrix(taxa_sets)

  upset_path <- file.path(outdir, "figures", paste0("upset_", tolower(rank), ".pdf"))
  pdf(upset_path, width = 10, height = 6)
  print(upset(upset_mat,
              sets = names(taxa_sets),
              order.by = "freq",
              mainbar.y.label = paste("Shared", tolower(rank)),
              sets.x.label = paste(rank, "per method"),
              text.scale = 1.3,
              mb.ratio = c(0.6, 0.4)))
  dev.off()
  cat("saved:", upset_path, "\n")
}

# ---- text summary ----
txt <- c("=== TAXA OVERLAP RESULTS ===",
         paste("Generated:", Sys.Date()), "")

for (rank in use_ranks) {
  taxa_sets <- imap(ps_filt, ~ get_taxa_set(.x, rank))
  taxa_sets <- taxa_sets[lengths(taxa_sets) > 0]
  if (length(taxa_sets) < 2) next

  all_taxa <- unique(unlist(taxa_sets))
  inter <- Reduce(intersect, taxa_sets)

  txt <- c(txt, paste0("--- ", rank, " ---"))
  for (nm in names(taxa_sets))
    txt <- c(txt, sprintf("  %s: %d", nm, length(taxa_sets[[nm]])))
  txt <- c(txt,
           sprintf("  Union: %d", length(all_taxa)),
           sprintf("  Intersection (all with data): %d", length(inter)),
           "")
}

save_summary(txt, file.path(outdir, "stats", "overlap_results_summary.txt"))