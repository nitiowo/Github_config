# zoop_functions.R
# helper functions for zooplankton phyloseq analysis
# source this via setup.R (not directly)

library(phyloseq)
library(tidyverse)
library(vegan)
library(ape)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(ggpubr)
library(VennDiagram)
library(UpSetR)
library(eulerr)
library(indicspecies)
library(adespatial)
library(geosphere)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(reshape2)
library(knitr)
library(kableExtra)
library(DT)
library(iNEXT)
library(ggtree)
library(DECIPHER)
library(Biostrings)
library(pheatmap)
library(flextable)


# ============================================================
# data prep
# ============================================================

set_lake_order <- function(ps, lake_order) {
  sd <- data.frame(sample_data(ps))
  sd$Lake <- factor(sd$Lake, levels = lake_order, ordered = TRUE)
  if ("Mesh" %in% colnames(sd)) sd$Mesh <- factor(sd$Mesh)
  sample_data(ps) <- sample_data(sd)
  ps
}

agg_rank <- function(ps, rank = "Species") {
  if (rank == "ASV") return(ps)
  tt <- access(ps, "tax_table", errorIfNULL = FALSE)
  if (is.null(tt)) {
    warning("no tax_table, skipping aggregation")
    return(ps)
  }
  if (!(rank %in% rank_names(ps))) {
    warning("rank '", rank, "' not in tax_table, skipping")
    return(ps)
  }
  tax_glom(ps, taxrank = rank, NArm = FALSE)
}

to_pa <- function(ps) {
  ot <- as(otu_table(ps), "matrix")
  ot[ot > 0] <- 1
  otu_table(ps) <- otu_table(ot, taxa_are_rows = taxa_are_rows(ps))
  ps
}

subset_taxa_custom <- function(ps, tsub = NULL) {
  if (is.null(tsub)) return(ps)
  tt <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
  rc <- tsub$rank
  if (!is.null(tsub$include)) {
    keep <- which(tt[[rc]] %in% tsub$include)
    ps <- prune_taxa(taxa_names(ps)[keep], ps)
  }
  if (!is.null(tsub$exclude)) {
    keep <- which(!(tt[[rc]] %in% tsub$exclude))
    ps <- prune_taxa(taxa_names(ps)[keep], ps)
  }
  prune_samples(sample_sums(ps) > 0, ps)
}

prepare_ps <- function(ps, rank = "Species", tsub = NULL) {
  ps %>% agg_rank(rank) %>% subset_taxa_custom(tsub)
}

rename_samples_to_station <- function(ps, station_col = "Station_ID") {
  sd <- data.frame(sample_data(ps))
  if (!(station_col %in% colnames(sd)))
    stop("column '", station_col, "' not in sample_data")
  new_names <- sd[[station_col]]
  if (any(duplicated(new_names)))
    warning("duplicate station IDs found")
  sample_names(ps) <- new_names
  ps
}

aggregate_to_station <- function(ps, station_col = "Station_ID") {
  sd <- data.frame(sample_data(ps))
  if (!(station_col %in% colnames(sd)))
    stop("column '", station_col, "' not in sample_data")

  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)

  stations <- sd[[station_col]]
  unique_st <- unique(stations)

  agg_mat <- matrix(0L, nrow = length(unique_st), ncol = ncol(otu),
                    dimnames = list(unique_st, colnames(otu)))
  for (st in unique_st) {
    idx <- which(stations == st)
    if (length(idx) == 1) {
      agg_mat[st, ] <- as.integer(otu[idx, ] > 0)
    } else {
      agg_mat[st, ] <- as.integer(colSums(otu[idx, , drop = FALSE]) > 0)
    }
  }

  new_sd <- sd %>%
    mutate(.station = stations) %>%
    group_by(.station) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    as.data.frame()
  rownames(new_sd) <- new_sd$.station
  new_sd$.station <- NULL

  new_ps <- phyloseq(otu_table(agg_mat, taxa_are_rows = FALSE),
                     sample_data(new_sd), tax_table(ps))
  prune_taxa(taxa_sums(new_ps) > 0, new_ps)
}


# ============================================================
# filtering
# ============================================================

filter_ps_list <- function(ps_list, markers = NULL, lakes = NULL,
                           mesh = NULL, tsub = NULL) {
  if (!is.null(markers))
    ps_list <- ps_list[intersect(names(ps_list), markers)]

  lapply(ps_list, function(ps) {
    if (!is.null(lakes)) {
      sd <- data.frame(sample_data(ps))
      keep <- rownames(sd)[sd$Lake %in% lakes]
      if (length(keep) > 0) ps <- prune_samples(keep, ps)
    }
    if (!is.null(mesh)) {
      sd <- data.frame(sample_data(ps))
      keep <- rownames(sd)[sd$Mesh %in% mesh]
      if (length(keep) > 0) ps <- prune_samples(keep, ps)
    }
    if (!is.null(tsub)) ps <- subset_taxa_custom(ps, tsub)
    ps
  })
}


# ============================================================
# combining markers
# ============================================================

combine_ps_pa <- function(ps_list, rank = "Species", morph_only = FALSE,
                          morph_samples = NULL) {
  ps_prepped <- lapply(ps_list, function(ps) {
    ps <- agg_rank(ps, rank)
    ps <- to_pa(ps)
    if (morph_only && !is.null(morph_samples))
      ps <- prune_samples(intersect(sample_names(ps), morph_samples), ps)
    ps
  })

  get_ids <- function(ps) {
    if (rank == "ASV") return(taxa_names(ps))
    tt <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
    ids <- tt[[rank]]
    ids[is.na(ids)] <- paste0("Unassigned_", seq_along(ids[is.na(ids)]))
    ids
  }

  all_taxa    <- unique(unlist(lapply(ps_prepped, get_ids)))
  all_samples <- Reduce(intersect, lapply(ps_prepped, sample_names))

  mat <- matrix(0L, nrow = length(all_taxa), ncol = length(all_samples),
                dimnames = list(all_taxa, all_samples))

  for (ps in ps_prepped) {
    otu <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) otu <- t(otu)
    ids <- get_ids(ps)
    for (i in seq_along(ids)) {
      tid <- ids[i]
      shared_s <- intersect(colnames(otu), all_samples)
      mat[tid, shared_s] <- pmax(mat[tid, shared_s],
                                  as.integer(otu[i, shared_s] > 0))
    }
  }

  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  sd <- sample_data(ps_prepped[[1]])[all_samples, ]
  tax_mat <- matrix(rownames(mat), ncol = 1,
                    dimnames = list(rownames(mat), rank))

  phyloseq(otu_table(mat, taxa_are_rows = TRUE),
           sd, tax_table(tax_mat))
}


# ============================================================
# long data frame builders (for ggplot)
# ============================================================

build_long_df <- function(ps_list, rank = "Genus", relative = TRUE,
                          tsub = NULL, lakes = NULL) {
  # melt one or more ps objects into a single long df for plotting
  imap_dfr(ps_list, function(ps, m) {
    ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(tsub)
    if (!is.null(lakes)) {
      sd <- data.frame(sample_data(ps_agg))
      keep <- rownames(sd)[sd$Lake %in% lakes]
      ps_agg <- prune_samples(keep, ps_agg)
    }
    if (relative)
      ps_agg <- transform_sample_counts(ps_agg, function(x) x / sum(x))
    df <- psmelt(ps_agg)
    df$Marker <- m
    df
  })
}


# ============================================================
# output helpers
# ============================================================

save_plot <- function(p, filepath, width = 12, height = 8) {
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  ggsave(filepath, p, width = width, height = height)
  cat("saved:", filepath, "\n")
}

save_stats <- function(df, filepath_base, caption = NULL) {
  # writes .csv and .docx versions of a data frame
  dir.create(dirname(filepath_base), recursive = TRUE, showWarnings = FALSE)

  csv_path <- paste0(filepath_base, ".csv")
  write.csv(df, csv_path, row.names = FALSE)
  cat("saved:", csv_path, "\n")

  # word table via flextable
  docx_path <- paste0(filepath_base, ".docx")
  ft <- flextable(as.data.frame(df))
  ft <- autofit(ft)
  ft <- theme_booktabs(ft)
  if (!is.null(caption)) ft <- set_caption(ft, caption)
  flextable::save_as_docx(ft, path = docx_path)
  cat("saved:", docx_path, "\n")
}

save_summary <- function(text, filepath) {
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, filepath)
  cat("saved:", filepath, "\n")
}

sig_stars <- function(p) {
  case_when(p < 0.001 ~ "***",
            p < 0.01  ~ "**",
            p < 0.05  ~ "*",
            TRUE      ~ "ns")
}


# ============================================================
# exploratory summary
# ============================================================

summarise_ps <- function(ps, name = "dataset",
                         tax_ranks = c("Phylum","Class","Order",
                                       "Family","Genus","Species")) {
  tt <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  taxon_totals <- rowSums(otu)
  grand_total <- sum(taxon_totals)

  rows <- list()
  for (r in tax_ranks) {
    if (r %in% colnames(tt)) {
      n_unique <- length(unique(na.omit(tt[[r]])))
      n_na <- sum(is.na(tt[[r]]))
      pct_un <- if (grand_total > 0) round(100 * sum(taxon_totals[is.na(tt[[r]])]) / grand_total, 1) else 0
      rows[[r]] <- tibble(rank = r, unique_taxa = n_unique,
                          unassigned_asvs = n_na, pct_reads_unassigned = pct_un)
    }
  }
  summary_df <- bind_rows(rows)
  summary_df$dataset <- name
  summary_df$samples <- nsamples(ps)
  summary_df$total_asvs <- ntaxa(ps)
  summary_df
}

summarise_ps_print <- function(ps, name, tax_ranks) {
  # console-friendly version
  cat("===", name, "===\n")
  cat("  Samples:", nsamples(ps), "| ASVs:", ntaxa(ps), "\n")
  df <- summarise_ps(ps, name, tax_ranks)
  for (i in seq_len(nrow(df))) {
    cat("  ", df$rank[i], ": ", df$unique_taxa[i], " unique, ",
        df$pct_reads_unassigned[i], "% reads unassigned\n", sep = "")
  }
  cat("\n")
}


# ============================================================
# alpha diversity
# ============================================================

compute_alpha <- function(ps, marker_name, metrics = c("Observed","InvSimpson"),
                          lake_order = NULL) {
  ad <- estimate_richness(ps, measures = metrics)
  ad$Sample_ID <- sample_names(ps)
  ad$Marker <- marker_name
  meta <- data.frame(sample_data(ps))
  meta$Sample_ID <- rownames(meta)
  ad <- left_join(ad, meta, by = "Sample_ID")
  if (!is.null(lake_order))
    ad$Lake <- factor(ad$Lake, levels = lake_order, ordered = TRUE)
  ad
}

compute_alpha_all <- function(ps_list, metrics = c("Observed","InvSimpson"),
                              lake_order = NULL) {
  bind_rows(imap(ps_list, ~ compute_alpha(.x, .y, metrics, lake_order)))
}

run_kruskal <- function(alpha_long, group_var, group_by_vars = NULL) {
  grp <- c(group_by_vars, "Metric")
  alpha_long %>%
    group_by(across(all_of(grp))) %>%
    summarise(
      H = tryCatch(kruskal.test(Value ~ .data[[group_var]])$statistic, error = function(e) NA_real_),
      df = tryCatch(kruskal.test(Value ~ .data[[group_var]])$parameter, error = function(e) NA_real_),
      p_value = tryCatch(kruskal.test(Value ~ .data[[group_var]])$p.value, error = function(e) NA_real_),
      .groups = "drop") %>%
    mutate(sig = sig_stars(p_value))
}

run_pairwise_wilcox <- function(alpha_long, group_var, group_by_vars = NULL) {
  grp <- c(group_by_vars, "Metric")
  alpha_long %>%
    group_by(across(all_of(grp))) %>%
    reframe({
      pw <- tryCatch(
        pairwise.wilcox.test(Value, .data[[group_var]], p.adjust.method = "BH"),
        error = function(e) NULL)
      if (is.null(pw)) return(tibble())
      as.data.frame(pw$p.value) %>%
        rownames_to_column("Group1") %>%
        pivot_longer(-Group1, names_to = "Group2", values_to = "p.adj") %>%
        filter(!is.na(p.adj))
    }) %>%
    mutate(sig = sig_stars(p.adj))
}


# ============================================================
# beta diversity
# ============================================================

run_ordination <- function(ps, method = "NMDS", distance = "jaccard",
                           binary = TRUE, color_var = "Lake",
                           shape_var = NULL, title = "",
                           color_palette = NULL,
                           point_size = 2.5, text_size = 10) {
  ps_use <- if (binary) to_pa(ps) else ps
  ord <- ordinate(ps_use, method = method, distance = distance)

  p <- plot_ordination(ps_use, ord, color = color_var, shape = shape_var) +
    geom_point(size = point_size, alpha = 0.8) +
    theme_minimal(base_size = text_size) +
    labs(title = title)

  if (!is.null(color_palette))
    p <- p + scale_color_manual(values = color_palette)

  # 95% t-ellipses if >= 2 groups with >= 3 points
  grps <- data.frame(sample_data(ps_use))[[color_var]]
  grp_n <- table(grps)
  if (sum(grp_n >= 3) >= 2)
    p <- p + stat_ellipse(aes(group = .data[[color_var]]),
                          type = "t", level = 0.95, linetype = 2)

  list(plot = p, ordination = ord, ps = ps_use)
}

run_permanova <- function(ps, formula_str, distance = "jaccard",
                          binary = TRUE, nperm = 999) {
  ps_use <- if (binary) to_pa(ps) else ps
  otu <- as(otu_table(ps_use), "matrix")
  if (taxa_are_rows(ps_use)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_use))
  dm <- vegdist(otu, method = distance, binary = binary)
  adonis2(as.formula(paste("dm", formula_str)), data = meta,
          permutations = nperm)
}

run_betadisper <- function(ps, group_var = "Lake",
                           distance = "jaccard", binary = TRUE) {
  ps_use <- if (binary) to_pa(ps) else ps
  otu <- as(otu_table(ps_use), "matrix")
  if (taxa_are_rows(ps_use)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_use))
  dm <- vegdist(otu, method = distance, binary = binary)
  bd <- betadisper(dm, meta[[group_var]])
  list(betadisper = bd, permtest = permutest(bd, permutations = 999))
}

build_marker_ps <- function(ps_list, rank = "Species", shared_samps = NULL) {
  if (is.null(shared_samps))
    shared_samps <- Reduce(intersect, map(ps_list, sample_names))

  combined <- imap(ps_list, function(ps, m) {
    ps <- prune_samples(shared_samps, ps) %>% agg_rank(rank) %>% to_pa()
    sample_names(ps) <- paste0(m, "__", sample_names(ps))
    sd <- data.frame(sample_data(ps))
    sd$Marker <- m
    sd$Original_Sample <- gsub(paste0(m, "__"), "", rownames(sd))
    sample_data(ps) <- sample_data(sd)
    ps
  })
  Reduce(merge_phyloseq, combined)
}


# ============================================================
# differential abundance
# ============================================================

run_ancombc <- function(ps, group_var = "Lake", rank = "Genus",
                        prev_cut = 0.10) {
  ps_agg <- agg_rank(ps, rank)
  prev <- apply(as(otu_table(ps_agg), "matrix"),
                if (taxa_are_rows(ps_agg)) 1 else 2,
                function(x) mean(x > 0))
  ps_filt <- prune_taxa(names(prev[prev >= prev_cut]), ps_agg)
  if (ntaxa(ps_filt) < 3) { cat("too few taxa after filtering\n"); return(NULL) }

  ANCOMBC::ancombc2(data = ps_filt, fix_formula = group_var,
                    p_adj_method = "BH", prv_cut = 0,
                    group = group_var, struc_zero = TRUE,
                    neg_lb = TRUE, verbose = FALSE)
}

run_simper_analysis <- function(ps, group_var = "Lake",
                                rank = "Genus", top_n = 15, tsub = NULL) {
  ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(tsub)
  otu <- as(otu_table(ps_agg), "matrix")
  if (taxa_are_rows(ps_agg)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_agg))
  sim <- simper(otu, meta[[group_var]], permutations = 99)
  summ <- lapply(names(summary(sim)), function(comp) {
    s <- summary(sim)[[comp]]
    s$taxon <- rownames(s)
    s$comparison <- comp
    head(s, top_n)
  })
  list(simper = sim, top = bind_rows(summ))
}

run_indicator <- function(ps, group_var = "Lake", rank = "Genus", tsub = NULL) {
  ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(tsub)
  otu <- as(otu_table(ps_agg), "matrix")
  if (taxa_are_rows(ps_agg)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_agg))
  multipatt(otu, meta[[group_var]], func = "IndVal.g",
            control = how(nperm = 999))
}


# ============================================================
# heatmaps
# ============================================================

plot_top_heatmap <- function(ps, rank = "Genus", top_n = 30,
                             annotation_var = "Lake",
                             lake_colors = NULL, lake_order = NULL,
                             tsub = NULL) {
  ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(tsub)
  ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))

  otu <- as(otu_table(ps_rel), "matrix")
  if (!taxa_are_rows(ps_rel)) otu <- t(otu)

  tt <- data.frame(tax_table(ps_rel), stringsAsFactors = FALSE)
  rn <- tt[[rank]]
  rn[is.na(rn)] <- paste0("Unassigned_", seq_along(rn[is.na(rn)]))
  rownames(otu) <- make.unique(rn)

  idx <- order(rowMeans(otu), decreasing = TRUE)[1:min(top_n, nrow(otu))]
  otu_top <- otu[idx, , drop = FALSE]

  ann <- data.frame(sample_data(ps_rel))[, annotation_var, drop = FALSE]
  ann_colors <- list()
  if ("Lake" %in% colnames(ann) && !is.null(lake_colors) && !is.null(lake_order))
    ann_colors[["Lake"]] <- lake_colors[levels(factor(ann$Lake, lake_order))]

  pheatmap(otu_top, annotation_col = ann, annotation_colors = ann_colors,
           cluster_cols = TRUE, cluster_rows = TRUE,
           fontsize_row = 7, fontsize_col = 5,
           color = viridis(100), border_color = NA)
}


# ============================================================
# geographic
# ============================================================

filter_geo_metadata <- function(meta, required_cols = c("Latitude","Longitude")) {
  complete <- complete.cases(meta[, required_cols, drop = FALSE])
  if (sum(!complete) > 0)
    cat("  dropped", sum(!complete), "samples missing", paste(required_cols, collapse="/"), "\n")
  meta[complete, , drop = FALSE]
}

run_ibd <- function(ps, meta, lake) {
  samps <- meta %>% filter(Lake == lake) %>% pull(Sample_ID)
  ps_lk <- prune_samples(samps, ps) %>% to_pa()
  if (nsamples(ps_lk) < 4) return(NULL)

  otu <- as(otu_table(ps_lk), "matrix")
  if (taxa_are_rows(ps_lk)) otu <- t(otu)
  comm_d <- vegdist(otu, "jaccard", binary = TRUE)

  coords <- meta %>%
    filter(Sample_ID %in% sample_names(ps_lk)) %>%
    arrange(match(Sample_ID, sample_names(ps_lk)))
  geo_d <- as.dist(distm(coords[, c("Longitude","Latitude")],
                          fun = distHaversine) / 1000)

  mt <- mantel(comm_d, geo_d, method = "spearman", permutations = 999)

  dd <- tibble(geo_km = as.vector(geo_d),
               dissim = as.vector(comm_d),
               Lake   = lake)
  list(r_stat = mt$statistic, p_val = mt$signif, dd = dd)
}


# ============================================================
# trees
# ============================================================

build_nj_tree <- function(ps, rank = "Species", max_seqs = 500) {
  if (!is.null(ps@refseq)) {
    seqs <- refseq(ps)
  } else {
    seqs <- Biostrings::DNAStringSet(taxa_names(ps))
    names(seqs) <- taxa_names(ps)
  }

  if (rank != "ASV") {
    ps_agg <- agg_rank(ps, rank)
    if (!is.null(ps@refseq)) {
      tt_orig <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
      tt_orig$asv <- taxa_names(ps)
      reps <- tt_orig %>% group_by(.data[[rank]]) %>%
        slice_head(n = 1) %>% pull(asv)
      seqs <- seqs[reps]
    } else {
      seqs <- Biostrings::DNAStringSet(taxa_names(ps_agg))
      names(seqs) <- taxa_names(ps_agg)
    }
  }

  if (length(seqs) > max_seqs) {
    cat("  subsampling to", max_seqs, "sequences\n")
    seqs <- seqs[sample(length(seqs), max_seqs)]
  }

  aligned <- DECIPHER::AlignSeqs(seqs, verbose = FALSE)
  dm <- DECIPHER::DistanceMatrix(aligned, verbose = FALSE)
  tree <- ape::nj(dm)
  ape::ladderize(tree)
}


# ============================================================
# venn / upset
# ============================================================

get_taxa_set <- function(ps, rank = "Species") {
  ps_agg <- agg_rank(ps, rank) %>% to_pa()
  tt <- data.frame(tax_table(ps_agg), stringsAsFactors = FALSE)
  if (rank == "ASV") return(taxa_names(ps_agg))
  ids <- tt[[rank]]
  unique(ids[!is.na(ids)])
}

make_upset_matrix <- function(taxa_sets) {
  all_taxa <- unique(unlist(taxa_sets))
  mat <- data.frame(row.names = all_taxa)
  for (nm in names(taxa_sets))
    mat[[nm]] <- as.integer(all_taxa %in% taxa_sets[[nm]])
  mat
}


# ============================================================
# ground truth
# ============================================================

load_trebitz <- function(filepath) {
  read.csv(filepath, stringsAsFactors = FALSE)
}

compare_to_ground_truth <- function(ps, trebitz_df, rank = "Species", lake = NULL) {
  ps_agg <- agg_rank(ps, rank) %>% to_pa()
  if (!is.null(lake)) {
    sd <- data.frame(sample_data(ps_agg))
    keep <- rownames(sd)[sd$Lake == lake]
    if (length(keep) == 0) return(tibble(taxon = character(), rank = character()))
    ps_agg <- prune_samples(keep, ps_agg)
    ps_agg <- prune_taxa(taxa_sums(ps_agg) > 0, ps_agg)
  }
  tt <- data.frame(tax_table(ps_agg), stringsAsFactors = FALSE)
  detected <- unique(na.omit(tt[[rank]]))
  known <- unique(na.omit(trebitz_df[[rank]]))
  unexpected <- setdiff(detected, known)
  if (length(unexpected) == 0) return(tibble(taxon = character(), rank = character()))
  tibble(taxon = unexpected, rank = rank)
}


# ============================================================
# variance partitioning
# ============================================================

run_varpart <- function(ps, env_vars, spatial_vars = NULL, binary = TRUE) {
  ps_use <- if (binary) to_pa(ps) else ps
  otu <- as(otu_table(ps_use), "matrix")
  if (taxa_are_rows(ps_use)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_use))

  all_vars <- env_vars
  if (!is.null(spatial_vars) && is.character(spatial_vars))
    all_vars <- c(all_vars, spatial_vars)

  # only use columns that actually exist
  available <- intersect(all_vars, colnames(meta))
  if (length(available) < length(all_vars)) {
    missing_cols <- setdiff(all_vars, colnames(meta))
    cat("  columns not found:", paste(missing_cols, collapse = ", "), "\n")
    env_vars <- intersect(env_vars, colnames(meta))
    if (is.character(spatial_vars))
      spatial_vars <- intersect(spatial_vars, colnames(meta))
    all_vars <- c(env_vars, if (is.character(spatial_vars)) spatial_vars)
  }
  if (length(all_vars) == 0) { warning("no valid variables"); return(NULL) }

  complete <- complete.cases(meta[, all_vars, drop = FALSE])
  if (sum(complete) < 4) {
    cat("  only", sum(complete), "complete samples, need >= 4\n")
    return(NULL)
  }
  if (sum(!complete) > 0) cat("  dropped", sum(!complete), "incomplete samples\n")

  meta <- meta[complete, , drop = FALSE]
  otu <- otu[complete, , drop = FALSE]
  otu <- otu[, colSums(otu) > 0, drop = FALSE]
  row_ok <- rowSums(otu) > 0
  otu <- otu[row_ok, , drop = FALSE]
  meta <- meta[row_ok, , drop = FALSE]

  if (nrow(otu) < 4 || ncol(otu) < 2) {
    cat("  too few samples/taxa after cleaning\n")
    return(NULL)
  }
  otu[!is.finite(otu)] <- 0

  env_df <- data.frame(lapply(meta[, env_vars, drop = FALSE], function(x) {
    if (is.factor(x) || is.character(x)) as.numeric(as.factor(x))
    else as.numeric(x)
  }))
  env_df[is.na(env_df)] <- 0

  if (!is.null(spatial_vars) && length(spatial_vars) > 0) {
    if (is.character(spatial_vars)) {
      spat_df <- data.frame(lapply(meta[, spatial_vars, drop = FALSE], as.numeric))
    } else {
      spat_df <- spatial_vars[complete, , drop = FALSE]
      spat_df <- spat_df[row_ok, , drop = FALSE]
      spat_df <- data.frame(lapply(spat_df, as.numeric))
    }
    spat_df[!is.finite(as.matrix(spat_df))] <- 0
    if (ncol(spat_df) > 0) varpart(otu, env_df, spat_df)
    else varpart(otu, env_df)
  } else {
    varpart(otu, env_df)
  }
}


# ============================================================
# focal taxon
# ============================================================

subset_focal_taxon <- function(ps, focal_rank, focal_name) {
  tt <- access(ps, "tax_table", errorIfNULL = FALSE)
  if (is.null(tt)) return(NULL)
  if (!(focal_rank %in% colnames(tt))) return(NULL)
  ps_sub <- subset_taxa_custom(ps, list(rank = focal_rank, include = focal_name))
  if (ntaxa(ps_sub) == 0) return(NULL)
  ps_sub
}

focal_detection_summary <- function(ps_list, focal_rank, focal_name,
                                    rank = "Species") {
  imap_dfr(ps_list, function(ps, m) {
    ps_sub <- subset_focal_taxon(ps, focal_rank, focal_name)
    if (is.null(ps_sub)) return(tibble(Marker = m, n_taxa = 0, n_samples_with = 0))
    ps_agg <- agg_rank(ps_sub, rank)
    n_taxa <- ntaxa(ps_agg)
    sums <- sample_sums(to_pa(ps_agg))
    tibble(Marker = m, n_taxa = n_taxa, n_samples_with = sum(sums > 0))
  })
}