# Phyloseq Analysis Pipeline

## Overview

This directory contains a comprehensive phyloseq-based analysis pipeline for zooplankton biodiversity data from the Great Lakes. The pipeline processes multiple molecular markers (COI-Folmer, COI-Leray, 18S) alongside morphological identification data to compare community composition, diversity patterns, and taxonomic concordance across sampling methods and geographic locations.

The analysis pipeline and its accompanying scripts are documented and organized to facilitate reproducibility and to allow traceability of analytical steps across runs and users.

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

## Notes on provenance and visibility

- This repository may be updated using manual edits, automation, or AI-assisted tools (for example, GitHub Copilot or other agents). When changes are pushed to a remote Git hosting service (such as GitHub), the visibility of AI involvement depends on the metadata included with the commit and any associated pull request (author/committer fields, commit message, and PR description).

- The purpose of this small README edit is to support an experiment comparing visibility of AI-assisted changes when made via different interfaces (e.g., GitHub Chat/Copilot Agent vs. local edits in VS Code). After pushing a test branch with this change, examine the remote commit and PR to determine what information about AI assistance is visible to collaborators or the public.

- Suggested next steps:
	- Commit this change to a new test branch and push to the remote.
	- Open a pull request and inspect the PR timeline, commit metadata, and any UI indicators that reference AI-assisted suggestions or automated tooling.
	- Optionally take screenshots or export the PR timeline for documentation.

If you want, I can commit this change here locally (with a branch and commit message) and help you push it to a remote test branch and open a PR so you can observe the visibility differences.
