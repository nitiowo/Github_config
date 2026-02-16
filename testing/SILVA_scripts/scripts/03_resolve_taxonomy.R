#!/usr/bin/env Rscript
#
# Step 3: Resolve taxonomy using taxize.
#
# Can be run in two ways:
#   1. Pipeline:     Rscript 03_resolve_taxonomy.R
#   2. Interactive:  source("scripts/03_resolve_taxonomy.R") in RStudio
#
# RESUME: If interrupted, re-run and it will skip already-processed taxa IDs.
#
# Resolution strategy:
#   1. Extract genus from species field (first word)
#   2. Use genus to resolve taxonomy via WORMS → NCBI → ITIS
#   3. If genus fails, try GBIF canonical name
#   4. If still fails, go up one rank
#   5. Once resolved, insert CLEANED species string back into Species field
#      (with exceptions: remove parentheticals, check for valid format, etc.)
#
# NCBI note: If NCBI returns "Metazoa" at Kingdom level, replace with "Animalia"
#
# Outputs:
#   03_taxonomy_fixed/all_taxa_with_ids.csv       (output file 1)
#   03_taxonomy_fixed/unique_taxa.csv             (output file 2)
#   03_taxonomy_fixed/resolved_taxonomy.csv       (output file 3)
#   03_taxonomy_fixed/resolution_log.csv          (output file 4)
#
# Requires: taxize, dplyr, stringr, Biostrings, yaml

# ═══════════════════════════════════════════════════════════════
# 0.  SETUP
# ═══════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(taxize)
  library(dplyr)
  library(stringr)
  library(Biostrings)
  library(yaml)
})

cfg <- yaml::read_yaml("config/params.yaml")

input_fasta    <- "02_outgroups/combined_with_outgroups.fasta"
out_dir        <- "03_taxonomy_fixed"
file_all_taxa  <- file.path(out_dir, "all_taxa_with_ids.csv")       # Output 1
file_unique    <- file.path(out_dir, "unique_taxa.csv")              # Output 2
file_resolved  <- file.path(out_dir, "resolved_taxonomy.csv")        # Output 3
file_log       <- file.path(out_dir, "resolution_log.csv")           # Output 4
dir.create(out_dir, showWarnings = FALSE)

SLEEP_SEC     <- cfg$sleep_between_batches
TARGET_RANKS  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ── Set NCBI API key if provided ─────────────────────────────
if (!is.null(cfg$ncbi_api_key) && nzchar(cfg$ncbi_api_key)) {
  Sys.setenv(ENTREZ_KEY = cfg$ncbi_api_key)
  message("NCBI Entrez API key set from config.")
}

# ═══════════════════════════════════════════════════════════════
# 1.  PARSE FASTA INTO TABLE (Output File 1)
# ═══════════════════════════════════════════════════════════════

message("═══ Step 1: Parsing FASTA headers into table ═══")

fa      <- readDNAStringSet(input_fasta)
headers <- names(fa)

# Parse each header into accession + variable-length ranks
parsed_list <- lapply(headers, function(h) {
  parts     <- str_split_fixed(h, "\\s+", 2)
  accession <- parts[1]
  taxonomy  <- parts[2]
  ranks     <- str_trim(str_split(str_remove(taxonomy, ";\\s*$"), ";")[[1]])
  ranks     <- ranks[ranks != ""]
  c(Accession = accession, setNames(ranks, paste0("V", seq_along(ranks))))
})

# Pad all rows to the same width
max_ranks   <- max(sapply(parsed_list, length)) - 1
padded_list <- lapply(parsed_list, function(x) {
  n <- length(x) - 1
  if (n < max_ranks) {
    x <- c(x, setNames(rep(NA_character_, max_ranks - n), paste0("V", (n + 1):max_ranks)))
  }
  x
})

all_taxa_df <- as.data.frame(do.call(rbind, padded_list), stringsAsFactors = FALSE)
rownames(all_taxa_df) <- NULL

# Write Output File 1
write.csv(all_taxa_df, file = file_all_taxa, row.names = FALSE)
message(sprintf("Output 1 written: %s (%d sequences, %d columns)", 
                file_all_taxa, nrow(all_taxa_df), ncol(all_taxa_df)))

# ═══════════════════════════════════════════════════════════════
# 2.  CREATE UNIQUE TAXONOMY TABLE (Output File 2)
# ═══════════════════════════════════════════════════════════════

message("\n═══ Step 2: Building unique taxonomy table ═══")

# Columns that hold taxonomy ranks (everything except Accession)
rank_cols <- setdiff(colnames(all_taxa_df), "Accession")

# Create a taxonomy signature string for grouping
all_taxa_df$tax_signature <- apply(all_taxa_df[, rank_cols], 1, function(row) {
  paste(ifelse(is.na(row), "NA", row), collapse = ";")
})

# Group by taxonomy signature, collect accessions
unique_groups <- all_taxa_df %>%
  group_by(tax_signature) %>%
  summarise(
    accessions = paste(Accession, collapse = "|"),
    n_seqs     = n(),
    .groups    = "drop"
  )

# Assign unique taxa IDs
unique_groups$taxa_id <- paste0("TAX_", sprintf("%05d", seq_len(nrow(unique_groups))))

# Reconstruct rank columns from signature
rank_matrix <- do.call(rbind, str_split(unique_groups$tax_signature, ";"))
rank_matrix[rank_matrix == "NA"] <- NA_character_

# Ensure column count matches
colnames(rank_matrix) <- rank_cols[1:ncol(rank_matrix)]
if (ncol(rank_matrix) < length(rank_cols)) {
  extra <- matrix(NA, nrow = nrow(rank_matrix), 
                  ncol = length(rank_cols) - ncol(rank_matrix))
  colnames(extra) <- rank_cols[(ncol(rank_matrix) + 1):length(rank_cols)]
  rank_matrix <- cbind(rank_matrix, extra)
}

unique_taxa_df <- cbind(
  data.frame(taxa_id = unique_groups$taxa_id, stringsAsFactors = FALSE),
  as.data.frame(rank_matrix, stringsAsFactors = FALSE),
  data.frame(
    accessions = unique_groups$accessions,
    n_seqs     = unique_groups$n_seqs,
    stringsAsFactors = FALSE
  )
)

# Write Output File 2
write.csv(unique_taxa_df, file = file_unique, row.names = FALSE)
message(sprintf("Output 2 written: %s (%d unique taxonomies from %d sequences)",
                file_unique, nrow(unique_taxa_df), nrow(all_taxa_df)))

# ═══════════════════════════════════════════════════════════════
# 3.  NAME EXTRACTION AND CLEANING
# ═══════════════════════════════════════════════════════════════

#' Extract genus name from a species field.
#' Takes the first word for genus extraction (database query).
#' Now allows problematic terms - they just trigger going up one rank
#' @param species_str  The species string
#' @return list(genus = genus name or NA, exception = reason, go_up = should we skip this rank)
extract_genus <- function(species_str) {
  if (is.na(species_str) || !nzchar(str_squish(species_str))) {
    return(list(genus = NA_character_, exception = "NA_or_empty", go_up = TRUE))
  }
  
  # Take first word
  words <- str_split(str_squish(species_str), "\\s+")[[1]]
  genus <- words[1]
  
  # Remove SILVA _X suffixes
  genus <- str_remove(genus, "_X{1,3}$")
  
  # Remove "Candidatus" if it's the first word
  if (genus == "Candidatus" && length(words) > 1) {
    genus <- words[2]
    genus <- str_remove(genus, "_X{1,3}$")
  }
  
  # Check if first letter is capitalized (RULE 1)
  if (!str_detect(genus, "^[A-Z]")) {
    return(list(genus = NA_character_, exception = "not_capitalized", go_up = TRUE))
  }
  
  # Check for problematic patterns (RULE 2) - but DON'T set genus to NA
  # Just flag to go up one rank
  if (str_detect(genus, "(?i)^(uncultured|unidentified|unknown|metagenome|metagenomic|environmental)")) {
    return(list(genus = genus, exception = "uncultured_environmental", go_up = TRUE))
  }
  
  # Check for accession-like or pure numeric
  if (str_detect(genus, "^[A-Z]{0,2}[0-9]{3,}$") || str_detect(genus, "^[0-9]+$")) {
    return(list(genus = NA_character_, exception = "accession_numeric", go_up = TRUE))
  }
  
  list(genus = genus, exception = "none", go_up = FALSE)
}

#' Clean species string for final output.
#' Applies rules for what to return as Species field.
#' @param species_str  The original species string
#' @param genus        The extracted genus (for validation)
#' @return Cleaned species string or NA
clean_species_for_output <- function(species_str, genus) {
  if (is.na(species_str) || !nzchar(str_squish(species_str))) {
    return(NA_character_)
  }
  
  cleaned <- str_squish(species_str)
  
  # RULE 5: Remove SILVA _X suffixes
  cleaned <- str_remove_all(cleaned, "_X{1,3}")
  
  # RULE 3: Strip parenthetical common names
  cleaned <- str_remove_all(cleaned, "\\s*\\([^)]*\\)")
  cleaned <- str_squish(cleaned)
  
  # RULE 1: Check if first letter is capitalized
  first_char <- substr(cleaned, 1, 1)
  if (!grepl("^[A-Z]", first_char)) {
    return(NA_character_)
  }
  
  # RULE 2: Do not return if contains problematic terms
  if (str_detect(cleaned, "(?i)(uncultured|unidentified|unknown|metagenome|metagenomic|environmental)")) {
    return(NA_character_)
  }
  
  # RULE 4: Handle sp./cf./aff. qualifiers
  # Check if there's something after the qualifier
  if (str_detect(cleaned, "\\b(sp\\.|cf\\.|aff\\.)")) {
    # Extract everything after the qualifier
    after_qualifier <- str_extract(cleaned, "(?<=sp\\.|cf\\.|aff\\.).*$")
    after_qualifier <- str_squish(after_qualifier)
    
    # If nothing after qualifier (or just whitespace), return NA
    if (is.na(after_qualifier) || !nzchar(after_qualifier)) {
      return(NA_character_)
    }
    
    # Otherwise return the full cleaned string
    return(cleaned)
  }
  
  # Return cleaned species string
  return(cleaned)
}

# ═══════════════════════════════════════════════════════════════
# 4.  GBIF LOOKUP: get canonical name (used as fallback)
# ═══════════════════════════════════════════════════════════════

#' Look up a name in GBIF and return the best canonical name.
#' Prefers HIGHERRANK matches over EXACT matches among ACCEPTED results.
#' @param query_name  The genus name to look up
#' @return list(canonical_name, success, details)
gbif_lookup <- function(query_name) {
  db_info <- tryCatch({
    result <- get_gbifid_(query_name)
    result[[1]]
  }, error = function(e) {
    return(data.frame())
  })
  
  if (!is.data.frame(db_info) || nrow(db_info) == 0) {
    return(list(canonical_name = NA_character_, success = FALSE, details = "no_results"))
  }
  
  # Filter to ACCEPTED status with EXACT or HIGHERRANK match
  db_filtered <- db_info %>%
    filter(
      status == "ACCEPTED",
      matchtype %in% c("EXACT", "HIGHERRANK")
    )
  
  if (nrow(db_filtered) == 0) {
    # Fallback: just ACCEPTED, any matchtype
    db_filtered <- db_info %>% filter(status == "ACCEPTED")
    if (nrow(db_filtered) == 0) {
      return(list(canonical_name = NA_character_, success = FALSE, details = "no_accepted"))
    }
  }
  
  # Prefer HIGHERRANK over EXACT
  higherrank <- db_filtered %>% filter(matchtype == "HIGHERRANK")
  if (nrow(higherrank) > 0) {
    canonical <- higherrank$canonicalname[1]
  } else {
    canonical <- db_filtered$canonicalname[1]
  }
  
  if (is.na(canonical) || !nzchar(canonical)) {
    return(list(canonical_name = NA_character_, success = FALSE, details = "empty_canonical"))
  }
  
  return(list(canonical_name = canonical, success = TRUE, details = "ok"))
}

# ═══════════════════════════════════════════════════════════════
# 5.  CLASSIFICATION LOOKUP: single DB attempt
# ═══════════════════════════════════════════════════════════════

#' Extract the 7 target ranks from a classification result.
#' Uses the unlist + matrix approach for reliable extraction.
#' @param query_name  The name to classify (genus)
#' @param db          Database to query ("worms", "ncbi", or "itis")
#' @param taxa_id     For console logging
#' @return list(ranks = named vector, success = logical) or NULL on error
classify_single_db <- function(query_name, db, taxa_id) {
  message(sprintf("    [%s] %s: classifying '%s'...", taxa_id, toupper(db), query_name))
  
  classy <- tryCatch({
    suppressWarnings({
      raw <- classification(query_name, db = db)
      
      # Convert to dataframe using the reliable unlist approach
      if (is.null(raw) || length(raw) == 0) return(NULL)
      
      # Check if the result is NA (classification not found)
      first_result <- raw[[1]]
      if (length(first_result) == 1 && is.na(first_result)) {
        message(sprintf("    [%s] %s: no classification found for '%s'", taxa_id, toupper(db), query_name))
        return(NULL)
      }
      
      # Convert to dataframe: unlist into matrix with 3 columns (name, rank, id)
      cls_df <- as.data.frame(matrix(unlist(first_result), ncol = 3), stringsAsFactors = FALSE)
      colnames(cls_df) <- c("name", "rank", "id")
      
      if (nrow(cls_df) == 0) return(NULL)
      
      cls_df
    })
  }, error = function(e) {
    message(sprintf("    [%s] %s ERROR: %s", taxa_id, toupper(db), e$message))
    return(NULL)
  })
  
  if (is.null(classy)) return(NULL)
  
  # Extract the 7 target ranks (but leave Species as NA - we'll fill it later)
  rank_map <- setNames(rep(NA_character_, 7), TARGET_RANKS)
  for (r in TARGET_RANKS[1:6]) {  # Kingdom through Genus only
    matching_rows <- classy[tolower(classy$rank) == tolower(r), ]
    if (nrow(matching_rows) > 0) {
      rank_map[r] <- matching_rows$name[1]
    }
  }
  # Species will be filled from cleaned original data later
  
  # NCBI-specific fix: Replace "Metazoa" with "Animalia" at Kingdom level
  if (db == "ncbi" && !is.na(rank_map["Kingdom"]) && rank_map["Kingdom"] == "Metazoa") {
    message(sprintf("    [%s] NCBI: replacing 'Metazoa' with 'Animalia' at Kingdom level", taxa_id))
    rank_map["Kingdom"] <- "Animalia"
  }
  
  # Must have at least Kingdom + one other rank to count as success
  n_filled <- sum(!is.na(rank_map[1:6]))  # Don't count Species
  if (n_filled >= 2) {
    message(sprintf("    [%s] %s: SUCCESS (%d/6 ranks filled, Species will be added from cleaned original)", 
                    taxa_id, toupper(db), n_filled))
    return(list(ranks = rank_map, success = TRUE))
  } else {
    message(sprintf("    [%s] %s: only %d/6 ranks found, treating as failure", taxa_id, toupper(db), n_filled))
    return(NULL)
  }
}

# ═══════════════════════════════════════════════════════════════
# 6.  CLASSIFICATION: WORMS → NCBI → ITIS cascade
# ═══════════════════════════════════════════════════════════════

#' Look up full classification: try WORMS → NCBI → ITIS in order.
#' Stop as soon as one succeeds.
#' @param query_name  The genus name to classify
#' @param taxa_id     For console logging
#' @return list(ranks, source_db, worms_tried, worms_success, ncbi_tried, ncbi_success, itis_tried, itis_success)
classify_name <- function(query_name, taxa_id) {
  
  results <- list(
    ranks         = setNames(rep(NA_character_, 7), TARGET_RANKS),
    source_db     = NA_character_,
    worms_tried   = FALSE,
    worms_success = FALSE,
    ncbi_tried    = FALSE,
    ncbi_success  = FALSE,
    itis_tried    = FALSE,
    itis_success  = FALSE
  )
  
  # ── Try WORMS ───────────────────────────────────────────────
  results$worms_tried <- TRUE
  worms_result <- classify_single_db(query_name, "worms", taxa_id)
  
  if (!is.null(worms_result) && worms_result$success) {
    results$ranks         <- worms_result$ranks
    results$source_db     <- "worms"
    results$worms_success <- TRUE
    return(results)
  }
  
  # ── WORMS failed → try NCBI ─────────────────────────────────
  message(sprintf("    [%s] WORMS failed for '%s', trying NCBI...", taxa_id, query_name))
  results$ncbi_tried <- TRUE
  ncbi_result <- classify_single_db(query_name, "ncbi", taxa_id)
  
  if (!is.null(ncbi_result) && ncbi_result$success) {
    results$ranks        <- ncbi_result$ranks
    results$source_db    <- "ncbi"
    results$ncbi_success <- TRUE
    return(results)
  }
  
  # ── NCBI failed → try ITIS ──────────────────────────────────
  message(sprintf("    [%s] NCBI failed for '%s', trying ITIS...", taxa_id, query_name))
  results$itis_tried <- TRUE
  itis_result <- classify_single_db(query_name, "itis", taxa_id)
  
  if (!is.null(itis_result) && itis_result$success) {
    results$ranks        <- itis_result$ranks
    results$source_db    <- "itis"
    results$itis_success <- TRUE
    return(results)
  }
  
  # ── All three failed ────────────────────────────────────────
  message(sprintf("    [%s] All databases (WORMS, NCBI, ITIS) failed for '%s'", taxa_id, query_name))
  return(results)
}

# ═══════════════════════════════════════════════════════════════
# 7.  GAP DETECTION IN TAXONOMY
# ═══════════════════════════════════════════════════════════════

#' Check if there are missing intermediate ranks in a taxonomy.
#' @param ranks  Named vector of 7 ranks
#' @return list(has_gap = logical, gap_description = string)
check_taxonomy_gaps <- function(ranks) {
  # Only check Kingdom through Genus (not Species)
  ranks_to_check <- ranks[1:6]
  is_filled <- !is.na(ranks_to_check) & ranks_to_check != "NA"
  
  if (sum(is_filled) == 0) {
    return(list(has_gap = FALSE, gap_description = ""))
  }
  
  filled_positions <- which(is_filled)
  first_filled <- min(filled_positions)
  last_filled <- max(filled_positions)
  
  positions_in_range <- first_filled:last_filled
  na_in_range <- positions_in_range[!is_filled[positions_in_range]]
  
  if (length(na_in_range) > 0) {
    gap_ranks <- names(ranks_to_check)[na_in_range]
    gap_desc <- paste0("Missing intermediate rank(s): ", paste(gap_ranks, collapse = ", "))
    return(list(has_gap = TRUE, gap_description = gap_desc))
  }
  
  return(list(has_gap = FALSE, gap_description = ""))
}

# ═══════════════════════════════════════════════════════════════
# 8.  LOG TRACKING
# ═══════════════════════════════════════════════════════════════

#' Append one log row to the resolution log (Output File 4)
append_log <- function(taxa_id, attempt_rank, raw_name, extracted_genus,
                       exception_hit, gbif_tried, gbif_success, canonical_name,
                       worms_tried, worms_success, ncbi_tried, ncbi_success,
                       itis_tried, itis_success, final_outcome, 
                       cleaned_species = "", has_gap = FALSE, gap_description = "") {
  log_row <- data.frame(
    taxa_id        = taxa_id,
    attempt_rank   = attempt_rank,
    raw_name       = ifelse(is.na(raw_name), "", raw_name),
    extracted_genus = ifelse(is.na(extracted_genus), "", extracted_genus),
    cleaned_species = cleaned_species,
    exception      = exception_hit,
    gbif_tried     = gbif_tried,
    gbif_success   = gbif_success,
    canonical_name = ifelse(is.na(canonical_name), "", canonical_name),
    worms_tried    = worms_tried,
    worms_success  = worms_success,
    ncbi_tried     = ncbi_tried,
    ncbi_success   = ncbi_success,
    itis_tried     = itis_tried,
    itis_success   = itis_success,
    final_outcome  = final_outcome,
    has_gap        = has_gap,
    gap_description = gap_description,
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(file_log)) {
    write.csv(log_row, file = file_log, row.names = FALSE)
  } else {
    write.table(log_row, file = file_log, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE, quote = TRUE)
  }
}

# ═══════════════════════════════════════════════════════════════
# 9.  SAVE SINGLE RESOLVED RESULT (Output File 3)
# ═══════════════════════════════════════════════════════════════

#' Append one resolved taxonomy to the resolved file
append_resolved <- function(taxa_id, ranks, source_db, accessions) {
  row <- data.frame(
    taxa_id    = taxa_id,
    Kingdom    = ranks["Kingdom"],
    Phylum     = ranks["Phylum"],
    Class      = ranks["Class"],
    Order      = ranks["Order"],
    Family     = ranks["Family"],
    Genus      = ranks["Genus"],
    Species    = ranks["Species"],
    source_db  = ifelse(is.na(source_db), "UNRESOLVED", source_db),
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

# ═══════════════════════════════════════════════════════════════
# 10.  MAIN RESOLUTION LOOP
# ═══════════════════════════════════════════════════════════════

message("\n═══ Step 3: Resolving taxonomy ═══")

# Check for already-completed taxa IDs (resume support)
completed_ids <- character()
if (file.exists(file_resolved)) {
  prev <- read.csv(file_resolved, stringsAsFactors = FALSE)
  completed_ids <- unique(prev$taxa_id)
  message(sprintf("Resuming: %d taxa IDs already resolved", length(completed_ids)))
}

# Filter to unprocessed taxa
todo_df <- unique_taxa_df %>% filter(!taxa_id %in% completed_ids)
message(sprintf("Taxa to resolve: %d / %d", nrow(todo_df), nrow(unique_taxa_df)))

if (nrow(todo_df) == 0) {
  message("All taxa already resolved!")
} else {
  
  for (i in seq_len(nrow(todo_df))) {
    row <- todo_df[i, ]
    tid <- row$taxa_id
    accessions_str <- row$accessions
    
    # Progress header for each taxon
    message(sprintf("\n══ [%d/%d] %s ══════════════════════════════", 
                    i, nrow(todo_df), tid))
    
    # Get the rank columns for this row (V1, V2, ... Vn)
    rank_values <- as.character(row[, rank_cols])
    non_na_idx <- which(!is.na(rank_values))
    
    if (length(non_na_idx) == 0) {
      message(sprintf("  [%s] All ranks are NA - marking as unresolved", tid))
      append_log(tid, "NONE", NA_character_, NA_character_, "all_NA",
                 FALSE, FALSE, NA_character_,
                 FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, "unresolved")
      append_resolved(tid, setNames(rep(NA_character_, 7), TARGET_RANKS), 
                      NA_character_, accessions_str)
      next
    }
    
    last_non_na <- max(non_na_idx)
    
    resolved <- FALSE
    final_ranks <- setNames(rep(NA_character_, 7), TARGET_RANKS)
    final_source <- NA_character_
    original_species <- NA_character_
    cleaned_species_output <- NA_character_
    
    # Iterate from lowest rank upward
    for (col_idx in last_non_na:1) {
      raw_name  <- rank_values[col_idx]
      col_name  <- rank_cols[col_idx]
      
      message(sprintf("  [%s] Trying rank %s: '%s'", tid, col_name, raw_name))
      
      # Store original species if this is the lowest rank
      if (col_idx == last_non_na) {
        original_species <- raw_name
      }
      
      # ── Extract genus from this rank ───────────────────────
      genus_extract <- extract_genus(raw_name)
      genus         <- genus_extract$genus
      exception     <- genus_extract$exception
      should_go_up  <- genus_extract$go_up
      
      # If extraction says to go up (problematic patterns), skip this rank
      if (should_go_up) {
        message(sprintf("  [%s] Skipping this rank (exception: %s) - going up", tid, exception))
        append_log(tid, col_name, raw_name, genus, exception,
                   FALSE, FALSE, NA_character_,
                   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, "go_up")
        next
      }
      
      # If genus is NA (capitalization or accession issue), also go up
      if (is.na(genus)) {
        message(sprintf("  [%s] Could not extract valid genus (exception: %s) - going up", tid, exception))
        append_log(tid, col_name, raw_name, NA_character_, exception,
                   FALSE, FALSE, NA_character_,
                   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, "go_up")
        next
      }
      
      message(sprintf("  [%s] Extracted genus: '%s' (exception: %s)", tid, genus, exception))
      
      # ── STEP 1: Try classification with extracted genus ───────
      message(sprintf("  [%s] Attempting classification with genus (WORMS → NCBI → ITIS)", tid))
      class_result <- classify_name(genus, tid)
      
      gbif_tried <- FALSE
      gbif_success <- FALSE
      canonical_name <- NA_character_
      
      if (class_result$worms_success || class_result$ncbi_success || class_result$itis_success) {
        # SUCCESS with genus
        final_ranks  <- class_result$ranks
        final_source <- class_result$source_db
        
        # Clean the species for output
        cleaned_species_output <- clean_species_for_output(original_species, genus)
        final_ranks["Species"] <- cleaned_species_output
        
        resolved <- TRUE
        
        # Check for gaps (only Kingdom through Genus)
        gap_check <- check_taxonomy_gaps(final_ranks)
        
        if (gap_check$has_gap) {
          message(sprintf("  [%s] ⚠ GAP DETECTED: %s", tid, gap_check$gap_description))
        }
        
        message(sprintf("  [%s] ✓ RESOLVED via %s using genus '%s'", tid, toupper(final_source), genus))
        if (!is.na(cleaned_species_output)) {
          message(sprintf("  [%s]   Cleaned species: '%s' (from '%s')", tid, cleaned_species_output, original_species))
        } else {
          message(sprintf("  [%s]   Species set to NA (invalid format or problematic content)", tid))
        }
        
        append_log(tid, col_name, raw_name, genus, exception,
                   gbif_tried, gbif_success, canonical_name,
                   class_result$worms_tried, class_result$worms_success,
                   class_result$ncbi_tried, class_result$ncbi_success,
                   class_result$itis_tried, class_result$itis_success, "resolved",
                   ifelse(is.na(cleaned_species_output), "", cleaned_species_output),
                   gap_check$has_gap, gap_check$gap_description)
        break
      }
      
      # ── STEP 2: All DBs failed → try GBIF canonical name ──────
      message(sprintf("  [%s] All DBs failed with genus, trying GBIF canonical...", tid))
      gbif_tried <- TRUE
      gbif_result <- gbif_lookup(genus)
      
      if (gbif_result$success) {
        gbif_success <- TRUE
        canonical_name <- gbif_result$canonical_name
        message(sprintf("  [%s] GBIF: canonical name = '%s'", tid, canonical_name))
        
        # Retry classification with canonical name
        message(sprintf("  [%s] Retrying classification with canonical name (WORMS → NCBI → ITIS)", tid))
        class_result2 <- classify_name(canonical_name, tid)
        
        if (class_result2$worms_success || class_result2$ncbi_success || class_result2$itis_success) {
          # SUCCESS with canonical name
          final_ranks  <- class_result2$ranks
          final_source <- class_result2$source_db
          
          # Clean the species for output
          cleaned_species_output <- clean_species_for_output(original_species, genus)
          final_ranks["Species"] <- cleaned_species_output
          
          resolved <- TRUE
          
          gap_check <- check_taxonomy_gaps(final_ranks)
          
          if (gap_check$has_gap) {
            message(sprintf("  [%s] ⚠ GAP DETECTED: %s", tid, gap_check$gap_description))
          }
          
          message(sprintf("  [%s] ✓ RESOLVED via %s using GBIF canonical name", tid, toupper(final_source)))
          if (!is.na(cleaned_species_output)) {
            message(sprintf("  [%s]   Cleaned species: '%s' (from '%s')", tid, cleaned_species_output, original_species))
          } else {
            message(sprintf("  [%s]   Species set to NA (invalid format or problematic content)", tid))
          }
          
          append_log(tid, col_name, raw_name, genus, exception,
                     gbif_tried, gbif_success, canonical_name,
                     class_result2$worms_tried, class_result2$worms_success,
                     class_result2$ncbi_tried, class_result2$ncbi_success,
                     class_result2$itis_tried, class_result2$itis_success, "resolved",
                     ifelse(is.na(cleaned_species_output), "", cleaned_species_output),
                     gap_check$has_gap, gap_check$gap_description)
          break
        }
        
        # Canonical name also failed
        message(sprintf("  [%s] ✗ All DBs failed even with GBIF canonical name - going up", tid))
        append_log(tid, col_name, raw_name, genus, exception,
                   gbif_tried, gbif_success, canonical_name,
                   class_result2$worms_tried, class_result2$worms_success,
                   class_result2$ncbi_tried, class_result2$ncbi_success,
                   class_result2$itis_tried, class_result2$itis_success, "go_up")
        
      } else {
        # GBIF also failed
        message(sprintf("  [%s] GBIF failed (%s) - going up", tid, gbif_result$details))
        append_log(tid, col_name, raw_name, genus, exception,
                   gbif_tried, FALSE, NA_character_,
                   class_result$worms_tried, class_result$worms_success,
                   class_result$ncbi_tried, class_result$ncbi_success,
                   class_result$itis_tried, class_result$itis_success, "go_up")
      }
    }
    
    # If nothing resolved after trying all ranks
    if (!resolved) {
      message(sprintf("  [%s] ✗✗ UNRESOLVED after exhausting all ranks", tid))
      append_log(tid, "EXHAUSTED", NA_character_, NA_character_, "all_ranks_failed",
                 FALSE, FALSE, NA_character_,
                 FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, "unresolved")
    }
    
    # Save result (resolved or unresolved)
    append_resolved(tid, final_ranks, final_source, accessions_str)
    
    # Brief sleep to avoid API throttling
    if (i %% 10 == 0) Sys.sleep(1)
  }
}

# ═══════════════════════════════════════════════════════════════
# 11.  SUMMARY
# ═══════════════════════════════════════════════════════════════

message("\n═══ Resolution Summary ═══")

resolved_all <- read.csv(file_resolved, stringsAsFactors = FALSE)
log_all <- read.csv(file_log, stringsAsFactors = FALSE)

n_total      <- nrow(resolved_all)
n_resolved   <- sum(resolved_all$source_db != "UNRESOLVED")
n_unresolved <- sum(resolved_all$source_db == "UNRESOLVED")

# Count rows with gaps
n_gaps <- sum(log_all$final_outcome == "resolved" & log_all$has_gap == TRUE, na.rm = TRUE)

# Count species set to NA after cleaning
n_species_na <- sum(resolved_all$source_db != "UNRESOLVED" & 
                      (is.na(resolved_all$Species) | resolved_all$Species == "NA"), na.rm = TRUE)

db_counts <- resolved_all %>%
  count(source_db) %>%
  arrange(desc(n))

msg <- paste0(
  sprintf("Total unique taxa:      %d\n", n_total),
  sprintf("Resolved:               %d\n", n_resolved),
  sprintf("  With taxonomy gaps:   %d (⚠ check log)\n", n_gaps),
  sprintf("  Species set to NA:    %d (invalid format/qualifier only)\n", n_species_na),
  sprintf("Unresolved:             %d\n", n_unresolved),
  "\nSource DB breakdown:\n",
  paste(sprintf("  %-15s %d", db_counts$source_db, db_counts$n), collapse = "\n"),
  "\n"
)
message(msg)

if (n_gaps > 0) {
  message(sprintf("\n⚠ WARNING: %d resolved taxa have intermediate rank gaps.", n_gaps))
  message("  Check resolution_log.csv where has_gap = TRUE for manual inspection.")
}

message("\nOutput files:")
message(sprintf("  1. All taxa:          %s", file_all_taxa))
message(sprintf("  2. Unique taxa:       %s", file_unique))
message(sprintf("  3. Resolved taxonomy: %s", file_resolved))
message(sprintf("  4. Resolution log:    %s", file_log))