#!/usr/bin/env Rscript
#
# Step 3: Resolve taxonomy using rgbif::name_backbone()
#
# Can be run in two ways:
#   1. Pipeline:     Rscript 03_resolve_taxonomy.R
#   2. Interactive:  source("scripts/03_resolve_taxonomy.R") in RStudio
#
# RESUME: If interrupted, re-run and it will skip already-processed taxa IDs
#
# Resolution strategy:
#   1. Try name_backbone() - accept if matchType EXACT or HIGHERRANK AND status ACCEPTED
#   2. If not, try name_backbone(verbose=TRUE) - pick best ACCEPTED candidate
#   3. If still no match, go up one rank and retry
#   4. Taxonomy pulled directly from $kingdom, $phylum, etc. columns
#
# MASTER LIST: Checks a cross-project CSV before any API call so each name
#   is only resolved once across all projects
#
# Outputs:
#   03_taxonomy_fixed/all_taxa_with_ids.csv       (output file 1)
#   03_taxonomy_fixed/unique_taxa.csv             (output file 2)
#   03_taxonomy_fixed/resolved_taxonomy.csv       (output file 3)
#   03_taxonomy_fixed/resolution_log.csv          (output file 4)
#
# Requires: rgbif, dplyr, stringr, Biostrings, yaml

# ---- Setup ----

suppressPackageStartupMessages({
  library(rgbif)
  library(dplyr)
  library(stringr)
  library(Biostrings)
  library(yaml)
})

cfg <- yaml::read_yaml("config/params.yaml")

input_fasta   <- "02_outgroups/combined_with_outgroups.fasta"
out_dir       <- "03_taxonomy_fixed"
file_all_taxa <- file.path(out_dir, "all_taxa_with_ids.csv")   # Output 1
file_unique   <- file.path(out_dir, "unique_taxa.csv")          # Output 2
file_resolved <- file.path(out_dir, "resolved_taxonomy.csv")    # Output 3
file_log      <- file.path(out_dir, "resolution_log.csv")       # Output 4
dir.create(out_dir, showWarnings = FALSE)

# Master taxonomy list - shared across projects
master_path  <- cfg$master_taxonomy_path
TARGET_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ---- Parse FASTA ----

message("=== Step 1: Parsing FASTA headers into table ===")

fa      <- readDNAStringSet(input_fasta)
headers <- names(fa)

# Detect if headers have accession numbers
# An accession token is short, alphanumeric, contains digits, no semicolons
has_accession <- function(h) {
  tok <- str_split_fixed(h, "\\s+", 2)[1]
  grepl("^[A-Za-z0-9_.]+$", tok) &&
    !grepl(";", tok) &&
    nchar(tok) <= 20 &&
    grepl("[0-9]", tok)
}

sample_check  <- head(headers, min(20, length(headers)))
use_accession <- mean(sapply(sample_check, has_accession)) >= 0.5
message(sprintf("Header format: %s",
                if (use_accession) "accession + taxonomy" else "taxonomy only"))

parsed_list <- lapply(headers, function(h) {
  if (use_accession) {
    parts     <- str_split_fixed(h, "\\s+", 2)
    accession <- parts[1]
    taxonomy  <- parts[2]
  } else {
    accession <- NA_character_
    taxonomy  <- h
  }
  ranks <- str_trim(str_split(str_remove(taxonomy, ";\\s*$"), ";")[[1]])
  ranks <- ranks[nzchar(ranks)]
  c(Accession = accession, setNames(ranks, paste0("V", seq_along(ranks))))
})

max_ranks   <- max(sapply(parsed_list, length)) - 1
padded_list <- lapply(parsed_list, function(x) {
  n <- length(x) - 1
  if (n < max_ranks) {
    x <- c(x, setNames(rep(NA_character_, max_ranks - n),
                       paste0("V", (n + 1):max_ranks)))
  }
  x
})

all_taxa_df  <- as.data.frame(do.call(rbind, padded_list), stringsAsFactors = FALSE)
rownames(all_taxa_df) <- NULL

write.csv(all_taxa_df, file = file_all_taxa, row.names = FALSE)
message(sprintf("Output 1 written: %s (%d sequences)", file_all_taxa, nrow(all_taxa_df)))

# ---- Unique taxonomy table ----

message("\n=== Step 2: Building unique taxonomy table ===")

rank_cols <- setdiff(colnames(all_taxa_df), "Accession")

all_taxa_df$tax_signature <- apply(all_taxa_df[, rank_cols], 1, function(row) {
  paste(ifelse(is.na(row), "NA", row), collapse = ";")
})

unique_groups <- all_taxa_df %>%
  group_by(tax_signature) %>%
  summarise(
    accessions = paste(Accession, collapse = "|"),
    n_seqs     = n(),
    .groups    = "drop"
  )

unique_groups$taxa_id <- paste0("TAX_", sprintf("%05d", seq_len(nrow(unique_groups))))

rank_matrix <- do.call(rbind, str_split(unique_groups$tax_signature, ";"))
rank_matrix[rank_matrix == "NA"] <- NA_character_
colnames(rank_matrix) <- rank_cols[1:ncol(rank_matrix)]

if (ncol(rank_matrix) < length(rank_cols)) {
  extra <- matrix(NA_character_, nrow = nrow(rank_matrix),
                  ncol = length(rank_cols) - ncol(rank_matrix))
  colnames(extra) <- rank_cols[(ncol(rank_matrix) + 1):length(rank_cols)]
  rank_matrix <- cbind(rank_matrix, extra)
}

unique_taxa_df <- cbind(
  data.frame(taxa_id = unique_groups$taxa_id, stringsAsFactors = FALSE),
  as.data.frame(rank_matrix, stringsAsFactors = FALSE),
  data.frame(accessions = unique_groups$accessions,
             n_seqs     = unique_groups$n_seqs,
             stringsAsFactors = FALSE)
)

write.csv(unique_taxa_df, file = file_unique, row.names = FALSE)
message(sprintf("Output 2 written: %s (%d unique taxonomies from %d sequences)",
                file_unique, nrow(unique_taxa_df), nrow(all_taxa_df)))

# ---- Helper: extract query name from a rank field ----

#' Extract a query string from a raw rank field.
#' Takes the first capitalized word (genus), stripping SILVA _X suffixes.
#' If the field looks like a binomial, passes both words as the query.
#' If the field starts with a lowercase letter (species epithet only, e.g. "spiniger"),
#'   returns needs_genus=TRUE so the loop can prepend the genus from the rank above.
#' Returns go_up=TRUE for problematic terms that should skip the rank.
#' @param s  Raw rank string from FASTA header
#' @return list(query, epithet, exception, go_up, needs_genus)
extract_query <- function(s) {
  base <- list(query = NA_character_, epithet = NA_character_,
               exception = "none", go_up = FALSE, needs_genus = FALSE)

  if (is.na(s) || !nzchar(str_squish(s))) {
    return(modifyList(base, list(exception = "NA_or_empty", go_up = TRUE)))
  }

  words <- str_split(str_squish(s), "\\s+")[[1]]
  word  <- str_remove(words[1], "_X{1,3}$")

  # Species epithet only (starts with lowercase) - need genus from rank above
  if (str_detect(word, "^[a-z]")) {
    return(modifyList(base, list(epithet = word, exception = "epithet_only",
                                 needs_genus = TRUE)))
  }

  # Skip "Candidatus" prefix
  if (word == "Candidatus" && length(words) > 1) {
    word <- str_remove(words[2], "_X{1,3}$")
  }

  if (!str_detect(word, "^[A-Z]")) {
    return(modifyList(base, list(exception = "not_capitalized", go_up = TRUE)))
  }

  if (str_detect(word, "(?i)^(uncultured|unidentified|unknown|metagenome|metagenomic|environmental)")) {
    return(modifyList(base, list(query = word, exception = "uncultured_environmental",
                                 go_up = TRUE)))
  }

  if (str_detect(word, "^[A-Z]{0,2}[0-9]{3,}$") || str_detect(word, "^[0-9]+$")) {
    return(modifyList(base, list(exception = "accession_numeric", go_up = TRUE)))
  }

  # Full binomial: first word capitalised, second word lowercase epithet
  if (length(words) >= 2 && str_detect(words[2], "^[a-z]")) {
    return(modifyList(base, list(query = paste(word, words[2]))))
  }

  modifyList(base, list(query = word))
}

# ---- Core rgbif resolution function ----

RANK_COLS_GBIF <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

#' Resolve a single name via rgbif::name_backbone().
#' Tries non-verbose first, then verbose with filtering.
#' @param query_name  Name to resolve
#' @param taxa_id     For console logging
#' @return list(ranks, matchType, confidence, usageKey) or NULL on failure
resolve_gbif <- function(query_name, taxa_id) {

  # Non-verbose attempt
  result <- tryCatch(
    name_backbone(name = query_name),
    error = function(e) {
      message(sprintf("    [%s] rgbif error: %s", taxa_id, e$message))
      NULL
    }
  )

  if (!is.null(result) && is.data.frame(result) && nrow(result) > 0) {
    mt <- result$matchType[1]
    st <- if ("status" %in% colnames(result)) result$status[1] else NA_character_

    if (!is.na(mt) && mt %in% c("EXACT", "HIGHERRANK") &&
        !is.na(st) && st == "ACCEPTED") {
      message(sprintf("    [%s] GBIF: %s / %s (conf %s)", taxa_id, mt, st, result$confidence[1]))
      return(extract_gbif_ranks(result[1, ], mt, result$confidence[1], result$usageKey[1]))
    }
  }

  # Verbose attempt - filter to best ACCEPTED candidate
  verbose_result <- tryCatch(
    name_backbone(name = query_name, verbose = TRUE),
    error = function(e) {
      message(sprintf("    [%s] rgbif verbose error: %s", taxa_id, e$message))
      NULL
    }
  )

  if (is.null(verbose_result) || !is.data.frame(verbose_result) ||
      nrow(verbose_result) == 0) {
    return(NULL)
  }

  # Ensure status column exists before filtering - some results omit it
  if (!"status" %in% colnames(verbose_result)) {
    message(sprintf("    [%s] GBIF verbose: no status column in result", taxa_id))
    return(NULL)
  }

  accepted <- verbose_result %>%
    filter(!is.na(.data$matchType),
           .data$matchType %in% c("EXACT", "HIGHERRANK"),
           !is.na(.data$status),
           .data$status == "ACCEPTED") %>%
    arrange(desc(.data$matchType == "EXACT"), desc(.data$confidence))

  if (nrow(accepted) == 0) {
    message(sprintf("    [%s] GBIF verbose: no ACCEPTED + EXACT/HIGHERRANK candidates", taxa_id))
    return(NULL)
  }

  best <- accepted[1, ]
  message(sprintf("    [%s] GBIF verbose: %s / %s (conf %s)",
                  taxa_id, best$matchType, best$status, best$confidence))
  extract_gbif_ranks(best, best$matchType, best$confidence, best$usageKey)
}

#' Pull 7 rank fields from a name_backbone() result row.
#' @param row        One-row data.frame
#' @param match_type matchType string
#' @param conf       Confidence value
#' @param usage_key  usageKey value
#' @return list(ranks, matchType, confidence, usageKey)
extract_gbif_ranks <- function(row, match_type, conf, usage_key) {
  ranks <- setNames(rep(NA_character_, 7), TARGET_RANKS)
  for (i in seq_along(TARGET_RANKS)) {
    col <- RANK_COLS_GBIF[i]
    if (col %in% colnames(row) && !is.na(row[[col]]) && nzchar(row[[col]])) {
      ranks[TARGET_RANKS[i]] <- row[[col]]
    }
  }
  list(ranks      = ranks,
       matchType  = match_type,
       confidence = as.numeric(conf),
       usageKey   = as.character(usage_key))
}

# ---- Master taxonomy list helpers ----

MASTER_COLS <- c("query_name", TARGET_RANKS, "matchType", "confidence",
                 "usageKey", "resolved_at")

#' Load the master taxonomy list.
#' Returns empty data.frame with correct column types if file doesn't exist.
load_master <- function(path) {
  if (!is.null(path) && nzchar(path) && file.exists(path)) {
    df <- read.csv(path, stringsAsFactors = FALSE,
                   colClasses = c(usageKey = "character", confidence = "numeric"))
    missing_cols <- setdiff(MASTER_COLS, colnames(df))
    if (length(missing_cols) > 0) df[missing_cols] <- NA_character_
    return(df)
  }
  # Build empty data.frame with explicit types to avoid bind_rows() coercion errors
  char_cols <- setdiff(MASTER_COLS, "confidence")
  df <- as.data.frame(
    lapply(setNames(char_cols, char_cols), function(x) character(0)),
    stringsAsFactors = FALSE
  )
  df$confidence <- numeric(0)
  df[, MASTER_COLS]  # Return columns in canonical order
}

#' Append one resolved entry to the master list and save it.
#' Returns the updated master data.frame.
append_master <- function(path, master_df, query_name, gbif_result) {
  if (is.null(path) || !nzchar(path)) return(master_df)

  new_row <- as.data.frame(
    c(list(query_name  = query_name),
      as.list(gbif_result$ranks),
      list(matchType   = gbif_result$matchType,
           confidence  = gbif_result$confidence,
           usageKey    = gbif_result$usageKey,
           resolved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
    stringsAsFactors = FALSE
  )

  master_df <- bind_rows(master_df, new_row)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(master_df, file = path, row.names = FALSE)
  master_df
}

# ---- Log and output helpers ----

#' Append one row to the resolution log (Output 4).
append_log <- function(taxa_id, attempt_rank, raw_name, query_name,
                       exception, match_type, confidence, final_outcome) {
  row <- data.frame(
    taxa_id       = taxa_id,
    attempt_rank  = attempt_rank,
    raw_name      = ifelse(is.na(raw_name), "", raw_name),
    query_name    = ifelse(is.na(query_name), "", query_name),
    exception     = exception,
    matchType     = ifelse(is.na(match_type), "", match_type),
    confidence    = ifelse(is.na(confidence), NA_real_, confidence),
    final_outcome = final_outcome,
    stringsAsFactors = FALSE
  )
  if (!file.exists(file_log)) {
    write.csv(row, file = file_log, row.names = FALSE)
  } else {
    write.table(row, file = file_log, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE, quote = TRUE)
  }
}

#' Append one resolved taxonomy to Output 3.
append_resolved <- function(taxa_id, ranks, match_type, confidence,
                             usage_key, accessions) {
  row <- data.frame(
    taxa_id    = taxa_id,
    Kingdom    = ranks["Kingdom"],
    Phylum     = ranks["Phylum"],
    Class      = ranks["Class"],
    Order      = ranks["Order"],
    Family     = ranks["Family"],
    Genus      = ranks["Genus"],
    Species    = ranks["Species"],
    matchType  = ifelse(is.na(match_type), "UNRESOLVED", match_type),
    confidence = ifelse(is.na(confidence), NA_real_, confidence),
    usageKey   = ifelse(is.na(usage_key), "", usage_key),
    accessions = accessions,
    stringsAsFactors = FALSE
  )
  if (!file.exists(file_resolved)) {
    write.csv(row, file = file_resolved, row.names = FALSE)
  } else {
    write.table(row, file = file_resolved, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE, quote = TRUE)
  }
}

# ---- Main resolution loop ----

message("\n=== Step 3: Resolving taxonomy ===")

master_df <- load_master(master_path)
if (nrow(master_df) > 0) {
  message(sprintf("Master taxonomy list loaded: %d entries", nrow(master_df)))
} else {
  message("No master taxonomy list found - all names will be queried via API")
}

# Resume: skip already-completed taxa IDs
completed_ids <- character()
if (file.exists(file_resolved)) {
  prev          <- read.csv(file_resolved, stringsAsFactors = FALSE)
  completed_ids <- unique(prev$taxa_id)
  message(sprintf("Resuming: %d taxa already resolved", length(completed_ids)))
}

todo_df <- unique_taxa_df %>% filter(!taxa_id %in% completed_ids)
message(sprintf("Taxa to resolve: %d / %d", nrow(todo_df), nrow(unique_taxa_df)))

if (nrow(todo_df) == 0) {
  message("All taxa already resolved!")
} else {

  for (i in seq_len(nrow(todo_df))) {
    row     <- todo_df[i, ]
    tid     <- row$taxa_id
    acc_str <- row$accessions

    message(sprintf("\n== [%d/%d] %s ==", i, nrow(todo_df), tid))

    rank_values <- as.character(row[, rank_cols])
    non_na_idx  <- which(!is.na(rank_values) & nzchar(str_squish(rank_values)))

    # Nothing to work with
    if (length(non_na_idx) == 0) {
      message(sprintf("  [%s] All ranks NA - unresolved", tid))
      append_log(tid, "NONE", NA, NA, "all_NA", NA, NA, "unresolved")
      append_resolved(tid, setNames(rep(NA_character_, 7), TARGET_RANKS),
                      NA, NA, NA, acc_str)
      next
    }

    last_idx    <- max(non_na_idx)
    resolved    <- FALSE
    final_ranks <- setNames(rep(NA_character_, 7), TARGET_RANKS)
    final_mt    <- NA_character_
    final_conf  <- NA_real_
    final_key   <- NA_character_

    # Try ranks from lowest upward
    for (col_idx in last_idx:1) {
      raw_name <- rank_values[col_idx]
      col_name <- rank_cols[col_idx]

      message(sprintf("  [%s] Trying rank %s: '%s'", tid, col_name, raw_name))

      q <- extract_query(raw_name)

      if (q$go_up) {
        message(sprintf("  [%s] Skipping (%s) - going up", tid, q$exception))
        append_log(tid, col_name, raw_name, q$query, q$exception, NA, NA, "go_up")
        next
      }

      # Epithet-only field (e.g. "spiniger") - combine with genus from rank above
      if (q$needs_genus) {
        genus_from_above <- NA_character_
        for (above_idx in (col_idx - 1):1) {
          above_name <- rank_values[above_idx]
          if (!is.na(above_name) && nzchar(str_squish(above_name))) {
            qa <- extract_query(above_name)
            if (!qa$go_up && !qa$needs_genus && !is.na(qa$query)) {
              # Use only the first word (genus) from the rank above
              genus_from_above <- str_split(qa$query, "\\s+")[[1]][1]
              break
            }
          }
        }

        if (is.na(genus_from_above)) {
          message(sprintf("  [%s] Epithet '%s' but no usable genus found above - going up",
                          tid, q$epithet))
          append_log(tid, col_name, raw_name, NA, "epithet_no_genus_above", NA, NA, "go_up")
          next
        }

        raw_name  <- paste(genus_from_above, q$epithet)
        q         <- extract_query(raw_name)  # Re-parse as binomial
        message(sprintf("  [%s] Epithet-only: combined with genus above -> '%s'",
                        tid, raw_name))
      }

      query_name <- q$query
      message(sprintf("  [%s] Query: '%s'", tid, query_name))

      # Check master list before hitting the API
      master_hit <- NULL
      if (nrow(master_df) > 0) {
        idx <- match(query_name, master_df$query_name)
        if (!is.na(idx)) {
          message(sprintf("  [%s] Found '%s' in master list", tid, query_name))
          mr         <- master_df[idx, ]
          master_hit <- list(
            ranks      = setNames(as.character(unlist(mr[, TARGET_RANKS])), TARGET_RANKS),
            matchType  = mr$matchType,
            confidence = as.numeric(mr$confidence),
            usageKey   = mr$usageKey
          )
        }
      }

      gbif_result <- if (!is.null(master_hit)) master_hit else resolve_gbif(query_name, tid)

      if (is.null(gbif_result)) {
        message(sprintf("  [%s] No result - going up", tid))
        append_log(tid, col_name, raw_name, query_name, q$exception, NA, NA, "go_up")
        next
      }

      # Persist new results to master list
      if (is.null(master_hit) && !is.null(master_path) && nzchar(master_path)) {
        master_df <- append_master(master_path, master_df, query_name, gbif_result)
      }

      final_ranks <- gbif_result$ranks
      final_mt    <- gbif_result$matchType
      final_conf  <- gbif_result$confidence
      final_key   <- gbif_result$usageKey
      resolved    <- TRUE

      message(sprintf("  [%s] RESOLVED: %s / conf %s", tid, final_mt, final_conf))
      append_log(tid, col_name, raw_name, query_name, q$exception,
                 final_mt, final_conf, "resolved")
      break
    }

    if (!resolved) {
      message(sprintf("  [%s] UNRESOLVED after all ranks", tid))
      append_log(tid, "EXHAUSTED", NA, NA, "all_ranks_failed", NA, NA, "unresolved")

      # Retain original taxonomy if it fits within the 7-rank schema
      # Fill from Kingdom downward - first original field -> Kingdom, etc.
      orig_values <- rank_values[non_na_idx]
      if (length(orig_values) <= 7) {
        final_ranks[seq_along(orig_values)] <- orig_values
        message(sprintf("  [%s] Retaining %d original field(s) in resolved output",
                        tid, length(orig_values)))
      } else {
        message(sprintf("  [%s] >7 original fields - leaving resolved output as NA", tid))
      }
    }

    # Save result immediately (incremental - supports resume)
    append_resolved(tid, final_ranks, final_mt, final_conf, final_key, acc_str)

    # Gentle throttle every 10 queries
    if (i %% 10 == 0) Sys.sleep(1)
  }
}

# ---- Summary ----

message("\n=== Resolution Summary ===")

resolved_all <- read.csv(file_resolved, stringsAsFactors = FALSE)
log_all      <- read.csv(file_log, stringsAsFactors = FALSE)

n_total      <- nrow(resolved_all)
n_resolved   <- sum(resolved_all$matchType != "UNRESOLVED", na.rm = TRUE)
n_unresolved <- n_total - n_resolved

mt_counts <- resolved_all %>%
  filter(matchType != "UNRESOLVED") %>%
  count(matchType) %>%
  arrange(desc(n))

msg <- paste0(
  sprintf("Total unique taxa:  %d\n", n_total),
  sprintf("Resolved:           %d\n", n_resolved),
  sprintf("Unresolved:         %d\n", n_unresolved),
  "\nMatch type breakdown (resolved):\n",
  paste(sprintf("  %-15s %d", mt_counts$matchType, mt_counts$n), collapse = "\n"),
  "\n"
)
message(msg)

message("Output files:")
message(sprintf("  1. All taxa:          %s", file_all_taxa))
message(sprintf("  2. Unique taxa:       %s", file_unique))
message(sprintf("  3. Resolved taxonomy: %s", file_resolved))
message(sprintf("  4. Resolution log:    %s", file_log))
if (!is.null(master_path) && nzchar(master_path)) {
  message(sprintf("  Master list:          %s (%d entries)", master_path, nrow(master_df)))
}
