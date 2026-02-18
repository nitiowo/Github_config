# Phyloseq Analysis Pipeline

## Overview

This directory contains a comprehensive phyloseq-based analysis pipeline for zooplankton biodiversity data from the Great Lakes. The pipeline processes multiple molecular markers (COI-Folmer, COI-Leray, 18S) alongside morphological identification data to compare community composition, diversity patterns, and taxonomic concordance across sampling methods and geographic locations.

## Main Analysis Files

- **`setup.R`** - Loads phyloseq objects, defines palettes, and sets global parameters
- **`alpha_Version2.R`** - Alpha diversity analysis (richness, diversity indices)
- **`beta_Version2.R`** - Beta diversity and ordination (NMDS, PCoA, PERMANOVA)
- **`composition_Version7.R`** - Community composition and taxonomic profiles
- **`differential_Version2.R`** - Differential abundance testing
- **`exploratory_Version7.R`** - Initial data exploration and summary statistics
- **`run_all_Version7.R`** - Master script to run all analyses sequentially
- **`zoop_functions_Version7.R`** - Custom functions for data manipulation and visualization

Additional analysis modules include geographic patterns, focal taxon analyses, heatmaps, rarefaction curves, overlap analyses, variance partitioning, and ground truth comparisons. See `cheatsheet.md` for detailed documentation on available objects, control options, and common usage patterns.
