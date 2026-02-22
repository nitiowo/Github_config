# Zooplankton Metabarcoding Analysis - Conversation Context

## Project Overview

This is a Great Lakes zooplankton metabarcoding project analyzing biodiversity across 96 sampling sites in the five Great Lakes (Superior, Michigan, Huron, Erie, Ontario). The project uses three metabarcoding markers plus morphological identification.

### Markers
- **COI-Folmer** (LCO1490 primers) → `ps_folmer` — 94 samples, 5,979 ASVs, 31 species
- **COI-Leray** (mICOintF primers) → `ps_leray` — 96 samples, 3,317 ASVs, 32 species
- **18S rRNA** (SSUF04/SSUR22 primers) → `ps_18S` — 95 samples, 1,195 ASVs. **Currently has near-zero taxonomic resolution below Class level** (100% unassigned at genus/species). User stated this is a mistake they will fix — treat 18S as having similar taxonomy to COI markers going forward, but allow excluding it via control block.
- **Morphology** → `ps_morph` — 29 samples (M01–M29), 87 taxa, 65 species. Biomass-adjusted abundance data.

### Key Data Details
- Marker samples are named Z01–Z96 (with 1–2 dropouts per marker)
- Morph samples are named M01–M29, one per station
- Marker samples have 2 per station (linked by `Station_ID` in metadata)
- All morph station IDs exist in the marker data, but not vice versa
- **Only sites with morphological data currently have lat/long/depth info** — this may be filled in later, so code must gracefully handle missing geo data

### Metadata Columns
`Sample_ID`, `Station_ID`, `Mesh` (53 or 100, categorical), `Year_Sampled`, `Month_Sampled`, `Lake`, `Morph_avail` (YES/NO), `Latitude`, `Longitude`, `Site_Depth`

### Comparison Rules
- **Between markers**: use presence/absence data, can use all samples
- **Between markers and morphology**: use presence/absence, only stations where morph data is available
- **Combined datasets**: merge by presence/absence only
- **Within a single marker**: abundance data is fine
- Lakes always ordered West to East: Superior, Michigan, Huron, Erie, Ontario
- Colors must be consistent throughout all plots

---

## Repository Structure

The code lives in `nitiowo/Metagenomics` repo. The phyloseq setup script is at `06_phyloseq/ps_setup.R`. Phyloseq objects are saved as RDS files in `06_phyloseq/setup_output/`.

Analysis scripts were developed in a separate repo `nitiowo/Github_config` under `testing/06_phyloseq/`.

---

## Architecture: Separate R Scripts (Not Rmd)

The user moved away from a monolithic Rmd approach to **separate R scripts** for each analysis type. The reasoning:
- Faster iteration (run individual lines in RStudio without knitting)
- Easier to tweak plots (ggplot portions are exposed, not wrapped in functions)
- No caching headaches
- Can combine into a final Rmd later once figures are finalized

### File Structure

```
06_phyloseq/
  zoop_functions.R      # all helper functions
  setup.R               # load data, build all ps objects, define palettes
  cheatsheet.md         # quick reference for objects and options
  run_all.R             # sources everything, saves all outputs
  exploratory.R         # summary tables, taxa counts per rank
  composition.R         # stacked barplots
  rarefaction.R         # iNEXT curves
  alpha.R               # alpha diversity + stats
  beta.R                # ordination, PERMANOVA, betadisper
  differential.R        # ANCOM-BC2, SIMPER, indicator species
  heatmaps.R            # top-taxa heatmaps
  overlap.R             # venn diagrams, upset plots
  geographic.R          # maps, IBD, dbMEMs, spatial
  focal_taxon.R         # calanoida deep-dive
  varpart.R             # variance partitioning
  ground_truth.R        # trebitz comparison
  output/
    <analysis_name>/
      stats/            # CSV + .docx tables
      figures/          # PDF plots
```

### Design Principles

1. **Each script has a "control block" at the top** where the user sets what to analyze (which markers, lakes, ranks, taxa subsets, etc.) then runs the rest of the script.

2. **ggplot code is exposed** — not wrapped in plotting functions. Functions handle data preparation only; the actual `ggplot()` call is visible so the user can freely edit titles, colors, sizes, facets, etc.

3. **Pre-built ps objects** in `setup.R` so the user picks the right object rather than subsetting on the fly:
   - Raw abundance: `ps_folmer`, `ps_leray`, `ps_18S`, `ps_morph`
   - Named lists: `ps_markers`, `ps_all_methods`, `ps_coi`
   - P/A versions: `ps_folmer_pa`, etc., `ps_markers_pa`, `ps_all_methods_pa`, `ps_coi_pa`
   - Station-aggregated: `ps_markers_by_station`, `ps_morph_by_station`
   - Combined: `ps_markers_combined`, `ps_coi_combined`, `ps_all_combined`

4. **Station-level aggregation** handles the M01–M29 vs Z01–Z96 mismatch. `aggregate_to_station()` collapses marker replicates to P/A per station. `rename_samples_to_station()` renames morph samples from M## to their Station_ID. Both sides then share Station_ID as sample names.

5. **Output strategy**: each script saves CSV + .docx (via flextable) for tables, PDF for figures, and a human-readable .txt summary. The .docx files produce editable Word tables for the manuscript.

6. **`run_all.R`** sources everything in order. Has an `output_root` variable at the top that can be changed to redirect all output to a different directory.

---

## Control Block Options (from cheatsheet.md)

```r
# which methods
use_ps_list <- ps_all_methods       # all 4
use_ps_list <- ps_markers           # 3 molecular only
use_ps_list <- ps_coi               # COI only

# subset markers within a list
use_markers <- NULL                 # all
use_markers <- c("Folmer", "Leray")

# filter by lake
use_lakes <- NULL                   # all 5
use_lakes <- c("Erie", "Ontario")

# taxonomic rank
use_rank <- "ASV" / "Species" / "Genus" / "Family" / "Order" / "Class" / "Phylum"

# taxa subsetting
use_tsub <- NULL
use_tsub <- list(rank = "Phylum", include = "Arthropoda")
use_tsub <- list(rank = "Phylum", exclude = "Rotifera")
use_tsub <- list(rank = "Order", include = "Calanoida")

# P/A vs abundance
use_pa <- TRUE / FALSE

# distance metric
use_distance <- "jaccard" / "bray"

# ordination
use_method <- "NMDS" / "PCoA"

# comparison variable
compare_var <- "Lake" / "Marker" / "Mesh"

# mesh filter
use_mesh <- NULL / c("53") / c("100")
```

---

## Key Functions in zoop_functions.R

### Data Prep
- `set_lake_order(ps, lake_order)` — factor Lake W→E, Mesh to factor
- `agg_rank(ps, rank)` — `tax_glom` wrapper, skips if rank="ASV" or no tax_table
- `to_pa(ps)` — convert to presence/absence
- `subset_taxa_custom(ps, tsub)` — filter by rank include/exclude
- `prepare_ps(ps, rank, tsub)` — aggregate + subset in one step
- `rename_samples_to_station(ps, station_col)` — rename sample names to Station_ID
- `aggregate_to_station(ps, station_col)` — collapse replicates to P/A per station
- `filter_ps_list(ps_list, markers, lakes, mesh, tsub)` — subset a named list
- `combine_ps_pa(ps_list, rank)` — merge multiple ps objects into one P/A matrix
- `build_long_df(ps_list, rank, relative, tsub, lakes)` — melt for ggplot

### Output
- `save_plot(p, filepath, width, height)` — ggsave wrapper
- `save_stats(df, filepath_base, caption)` — writes .csv + .docx (flextable)
- `save_summary(text, filepath)` — writeLines for .txt summaries
- `sig_stars(p)` — returns "***", "**", "*", or "ns"

### Analysis
- `summarise_ps(ps, name, tax_ranks)` — returns tibble with unique taxa, % unassigned per rank
- `compute_alpha(ps, marker_name, metrics, lake_order)` — estimate_richness wrapper
- `run_kruskal(alpha_long, group_var, group_by_vars)` — grouped KW test
- `run_pairwise_wilcox(alpha_long, group_var, group_by_vars)` — grouped pairwise Wilcoxon
- `run_ordination(ps, method, distance, binary, color_var, ...)` — ordinate + plot
- `run_permanova(ps, formula_str, distance, binary, nperm)` — adonis2 wrapper
- `run_betadisper(ps, group_var, distance, binary)` — betadisper + permutest
- `build_marker_ps(ps_list, rank, shared_samps)` — combine for between-marker comparison
- `run_ancombc(ps, group_var, rank, prev_cut)` — ANCOM-BC2
- `run_simper_analysis(ps, group_var, rank, top_n, tsub)` — SIMPER
- `run_indicator(ps, group_var, rank, tsub)` — IndVal
- `plot_top_heatmap(ps, rank, top_n, ...)` — pheatmap of top taxa
- `run_ibd(ps, meta, lake)` — Mantel test for isolation by distance. Returns `list(r_stat, p_val, dd)`
- `filter_geo_metadata(meta, required_cols)` — drops rows with missing lat/lon
- `build_nj_tree(ps, rank, max_seqs)` — NJ tree from refseq
- `get_taxa_set(ps, rank)` — unique non-NA taxa names at a rank
- `make_upset_matrix(taxa_sets)` — binary membership matrix
- `compare_to_ground_truth(ps, trebitz_df, rank, lake)` — flag unexpected detections
- `run_varpart(ps, env_vars, spatial_vars, binary)` — variance partitioning (handles missing data)
- `subset_focal_taxon(ps, focal_rank, focal_name)` — extract one taxon group
- `focal_detection_summary(ps_list, focal_rank, focal_name, rank)` — count detections per method

---

## Known Issues and Fixes Applied

### Sample name mismatch (morph M## vs marker Z##)
- **Problem**: morph samples are M01–M29, marker samples are Z01–Z96, linked by Station_ID
- **Fix**: `rename_samples_to_station()` renames morph to Station_ID; `aggregate_to_station()` collapses marker data to station-level P/A. Both now share Station_ID as sample names.
- This fix must be applied **anywhere** morph data is combined with marker data (combining, geographic analysis, etc.)

### `combine_ps_pa` zero-dimension error
- **Cause**: sample names didn't match between station-aggregated markers and morph
- **Fix**: build `ps_morph_by_station` via `rename_samples_to_station()` before combining

### Geographic metadata mismatch
- **Cause**: `meta_geo` was built from `ps_morph` (M## names) but `ps_all_combined` uses Station_IDs
- **Fix**: build `meta_geo` from `ps_morph_by_station` instead

### `tax_glom` error on objects without tax_table
- **Cause**: combined P/A objects or subsetted objects can lose their tax_table
- **Fix**: `agg_rank()` checks for tax_table existence and valid rank before calling `tax_glom`

### Mantel `$` operator error
- **Cause**: accessing mantel results incorrectly
- **Fix**: `run_ibd()` returns `list(r_stat = mt$statistic, p_val = mt$signif, dd = ...)`. Callers use `r$r_stat` and `r$p_val`.

### Varpart NA/NaN/Inf error
- **Cause**: missing values in metadata or zero-sum taxa/samples after subsetting
- **Fix**: `run_varpart()` does extensive cleaning — drops missing rows, removes zero-sum columns/rows, replaces non-finite values, checks for minimum sample count, handles missing columns gracefully

---

## Actual Results Data (from results_dump.txt)

### Dataset Sizes
| Marker | Samples | ASVs | Species (unique) | Species (% unassigned reads) |
|--------|---------|------|-----------------|------------------------------|
| Folmer | 94 | 5,979 | 31 | 84.4% |
| Leray | 96 | 3,317 | 32 | 77.8% |
| 18S | 95 | 1,195 | 0 | 100% |
| Morphology | 29 | 87 | 65 | 26.8% |
| Combined markers | 93 | — | 120 | 0% |
| Combined all | 29 | — | 176 | 0% |

### Alpha Diversity Summary (mean by marker × lake)
| Marker | Lake | Observed (mean) | InvSimpson (mean) | n |
|--------|------|----------------|-------------------|---|
| 18S | Superior | 38.3 | 3.76 | 20 |
| 18S | Michigan | 56.2 | 6.97 | 20 |
| 18S | Huron | 27.3 | 4.12 | 16 |
| 18S | Erie | 52.7 | 5.08 | 19 |
| 18S | Ontario | 42.2 | 4.96 | 20 |
| Folmer | Superior | 124 | 10.5 | 20 |
| Folmer | Michigan | 143 | 10.5 | 20 |
| Folmer | Huron | 59 | 5.09 | 15 |
| Folmer | Erie | 113 | 10.0 | 19 |
| Folmer | Ontario | 155 | 9.73 | 20 |
| Leray | Superior | 94.9 | 13.7 | 20 |
| Leray | Michigan | 87.1 | 12.1 | 20 |
| Leray | Huron | 99.5 | 10.1 | 16 |
| Leray | Erie | 135 | 17.5 | 20 |
| Leray | Ontario | 78.4 | 10.7 | 20 |
| Morphology | Superior | 17.3 | 4.61 | 6 |
| Morphology | Michigan | 28 | 5.44 | 4 |
| Morphology | Huron | 26.4 | 5.81 | 5 |
| Morphology | Erie | 38.1 | 5.59 | 9 |
| Morphology | Ontario | 31 | 6.93 | 5 |

### Key Statistical Results

**KW between markers:**
- Observed: H = 139, p = 4.86e-30
- InvSimpson: H = 86.2, p = 1.42e-18

**Pairwise (Observed):** Folmer vs Leray p = 0.244 (ns); all other pairs significant

**KW between lakes (per marker):**
- Folmer Observed: H = 27.7, p = 1.45e-05
- Leray Observed: H = 11.4, p = 0.022
- 18S Observed: H = 16.1, p = 0.003
- Morphology Observed: H = 22.7, p = 1.48e-04
- Leray InvSimpson: H = 5.41, p = 0.247 (ns)
- Morphology InvSimpson: H = 2.39, p = 0.665 (ns)

**PERMANOVA Jaccard (lake effect):**
| Marker | R² | F | p |
|--------|-----|------|------|
| Folmer | 0.126 | 3.22 | 0.001 |
| Leray | 0.153 | 4.12 | 0.001 |
| 18S | 0.160 | 4.29 | 0.001 |
| Morphology | 0.546 | 7.22 | 0.001 |

**PERMANOVA Bray-Curtis (lake effect):**
| Marker | R² | F | p |
|--------|-----|------|------|
| Folmer | 0.274 | 8.42 | 0.001 |
| Leray | 0.248 | 7.51 | 0.001 |
| 18S | 0.239 | 7.06 | 0.001 |
| Morphology | 0.483 | 5.60 | 0.001 |

**PERMANOVA marker effect (Jaccard):** R² = 0.519, F = 149.1, p = 0.001

**Venn (species level):** Folmer 31, Leray 32, 18S 0, Morphology 65. Union = 104. Intersection of all 4 = 0.

**Focal taxon (Calanoida):**
| Marker | Species | Samples detected |
|--------|---------|-----------------|
| Folmer | 14 | 92/94 |
| Leray | 16 | 96/96 |
| 18S | 0 | 0 |
| Morphology | 15 | 29/29 |

---

## Paper Context

### Scope
- Main story: "here's the biodiversity of the Great Lakes" (ecological focus)
- Marker comparison is secondary but included (not totally novel, but extensive 96-site sampling is the strength)
- Exploratory/descriptive, no specific hypothesis
- Includes management/monitoring angle

### Marker Names in Paper
- COI-Folmer, COI-Leray, 18S

### Ground Truth
- Trebitz 2019 (living online dataset) — comparison framed as "contextualizing results"
- Identify misidentifications/range expansions (taxa in our data not in Trebitz list)
- Detection gaps not needed
- Compare between-lake richness patterns with Trebitz

### Morphology Role
- Traditional benchmark being compared against
- But all data sources included when assessing biodiversity

### Results Section
- A filled-in results skeleton exists (see `results_skeleton_filled.md` in conversation)
- Many gaps remain for: composition descriptions, SIMPER/IndVal specifics, geographic results, Trebitz comparison, variance partitioning fractions, concordance Spearman values, betadisper values
- Occupancy modeling excluded from paper
- Variance partitioning with depth excluded (too much missing data)
- Focal taxon (Calanoida) gets its own brief section

---

## User Preferences

### Code Style
- Write code and comments like an **intermediate programmer**
- Avoid AI "tells" in comments and code structure
- Comment from "this is what this code does and why" perspective, not "instructions for someone else"
- No over-engineered fallback chains or excessive error handling patterns

### Workflow
- Prefers R scripts over Rmd for exploration phase
- Wants exposed ggplot code for easy tweaking
- Control blocks at top of each script
- Separate output directories per analysis type (not one giant pile)
- No numbered script prefixes (just descriptive names)
- Verbose file naming: `alpha_boxplot_lake_by_marker_genus.pdf`

### Tables
- CSV always
- .docx via flextable for manuscript (editable Word tables)
- Human-readable .txt summaries for quick review

### Figures
- PDF output
- Default dimensions are fine

### Month/Year
- Not important for this dataset, don't include in analyses

### Missing Data
- Gracefully exclude missing lat/long/depth rows so code doesn't need changing when data is filled in later