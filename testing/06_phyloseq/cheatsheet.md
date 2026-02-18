# Phyloseq Analysis Cheatsheet

Quick reference for all available objects, options, and common patterns.

## Phyloseq Objects (available after `source("setup.R")`)

### Raw abundance data (one object per method)

| Object | Samples | Description |
|--------|---------|-------------|
| `ps_folmer` | Z01-Z96 (minus dropouts) | COI-Folmer ASVs, abundance counts |
| `ps_leray` | Z01-Z96 (minus dropouts) | COI-Leray ASVs, abundance counts |
| `ps_18S` | Z01-Z96 (minus dropouts) | 18S ASVs, abundance counts |
| `ps_morph` | M01-M29 | Morphological IDs, biomass-adjusted abundance |

### Named lists (for looping / imap)

| Object | Contents |
|--------|----------|
| `ps_markers` | `list(Folmer=ps_folmer, Leray=ps_leray, "18S"=ps_18S)` |
| `ps_all_methods` | markers + `Morphology=ps_morph` |
| `ps_coi` | `list(Folmer=ps_folmer, Leray=ps_leray)` — COI only |

### Presence/absence versions

| Object | Description |
|--------|-------------|
| `ps_folmer_pa`, `ps_leray_pa`, `ps_18S_pa`, `ps_morph_pa` | P/A of each method |
| `ps_markers_pa` | P/A named list (3 markers) |
| `ps_all_methods_pa` | P/A named list (3 markers + morph) |
| `ps_coi_pa` | P/A named list (COI only) |

### Station-aggregated (for morph comparisons)

Marker samples (Z01-Z96, 2 per station) collapsed to 1 row per station by
P/A. Morph samples renamed from M01-M29 to their Station_ID values so
sample names match.

| Object | Description |
|--------|-------------|
| `ps_markers_by_station` | List of station-level P/A per marker |
| `ps_morph_by_station` | Morph with Station_ID as sample names |

### Combined datasets

| Object | Description |
|--------|-------------|
| `ps_markers_combined` | 3 markers merged, sample-level P/A, shared samples |
| `ps_coi_combined` | COI-Folmer + COI-Leray merged, sample-level P/A |
| `ps_all_combined` | 3 markers + morph, station-level P/A |

---

## Control Block Options

### `use_ps_list` — which methods to include

```r
use_ps_list <- ps_all_methods       # all 4 methods
use_ps_list <- ps_markers           # 3 molecular markers only
use_ps_list <- ps_coi               # COI-Folmer + COI-Leray only
use_ps_list <- ps_all_methods_pa    # all 4, already P/A
use_ps_list <- ps_markers_pa        # 3 markers, already P/A
use_ps_list <- list(Folmer = ps_folmer)  # single marker
```

### `use_markers` — subset within a list

```r
use_markers <- NULL                        # use everything in use_ps_list
use_markers <- c("Folmer", "Leray")        # only these two
use_markers <- c("Folmer")                 # just one
```

### `use_lakes` — filter samples by lake

```r
use_lakes <- NULL                                          # all 5 lakes
use_lakes <- c("Erie", "Ontario")                          # just these
use_lakes <- c("Superior", "Michigan", "Huron")            # upper lakes
```

### `use_rank` — taxonomic aggregation level

```r
use_rank <- "ASV"       # no aggregation
use_rank <- "Species"
use_rank <- "Genus"
use_rank <- "Family"
use_rank <- "Order"
use_rank <- "Class"
use_rank <- "Phylum"
```

### `use_tsub` — taxa subsetting

```r
use_tsub <- NULL                                           # all taxa
use_tsub <- list(rank = "Phylum", include = "Arthropoda")  # only arthropods
use_tsub <- list(rank = "Phylum", exclude = "Rotifera")    # exclude rotifers
use_tsub <- list(rank = "Order",  include = "Calanoida")   # only calanoids
use_tsub <- list(rank = "Class",  include = c("Copepoda", "Branchiopoda"))
```

### `use_pa` / `use_binary` — presence/absence vs abundance

```r
use_pa <- TRUE     # convert to P/A before analysis
use_pa <- FALSE    # use raw abundance counts
```

### `use_distance` — dissimilarity metric

```r
use_distance <- "jaccard"    # P/A (set use_binary = TRUE)
use_distance <- "bray"       # abundance (set use_binary = FALSE)
```

### `use_method` — ordination method

```r
use_method <- "NMDS"
use_method <- "PCoA"
use_method <- "DCA"
```

### `compare_var` — what to compare between

```r
compare_var <- "Lake"          # between lakes
compare_var <- "Marker"        # between markers (needs combined ps)
compare_var <- "Mesh"          # between mesh sizes
compare_var <- "Station_ID"    # between stations
```

### `facet_var` — what to facet/panel by

```r
facet_var <- NULL              # no faceting
facet_var <- "Lake"
facet_var <- "Marker"
facet_var <- "Mesh"
```

### `use_mesh` — filter by mesh size

```r
use_mesh <- NULL               # both mesh sizes
use_mesh <- c("53")            # 53 um only
use_mesh <- c("100")           # 100 um only
```

---

## Palettes

| Name | Description |
|------|-------------|
| `lake_order` | `c("Superior", "Michigan", "Huron", "Erie", "Ontario")` — W to E |
| `lake_colors` | Named vector, one color per lake |
| `marker_colors` | Named vector: Folmer, Leray, 18S, Morphology |
| `tax_ranks` | `c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")` |

---

## Common Patterns

### Aggregate to a rank, then subset

```r
ps_genus <- agg_rank(ps_folmer, "Genus")
ps_arthropods <- subset_taxa_custom(ps_genus, list(rank="Phylum", include="Arthropoda"))
```

### Get a filtered list of ps objects ready to plot

```r
ps_filt <- filter_ps_list(ps_all_methods,
                          markers = c("Folmer", "Leray"),
                          lakes = c("Erie", "Ontario"),
                          tsub = list(rank = "Phylum", include = "Arthropoda"))
```

### Build a long data frame for ggplot from multiple markers

```r
df <- build_long_df(ps_filt, rank = "Genus", relative = TRUE)
# returns: Sample, Marker, Lake, Mesh, Station_ID, taxon, Abundance
```

### Save a plot

```r
save_plot(p, "output/alpha/figures/alpha_boxplot_lake.pdf", width = 12, height = 8)
```

### Save stats as CSV + Word table

```r
save_stats(df, "output/alpha/stats/kruskal_results", caption = "Kruskal-Wallis results")
# creates: kruskal_results.csv and kruskal_results.docx
```

---

## Output Structure

Each analysis script writes to `output/<analysis>/stats/` and
`output/<analysis>/figures/`. Default output root is `output/` but
can be changed via `output_root` in `setup.R` or `run_all.R`.