# Copilot Context: Repositories, Workflow, and Specific Suggestions

This document captures the decisions, patterns, and repo-specific context from our discussion so you can reuse it in future chats. It includes practical workflow guidance for working across your local machine and the university HPC and concrete recommendations for the two repositories: NicloTemp and Metagenomics.

Use this file as shared context by referencing it at the start of a new chat: “Use docs/copilot-context.md as context.”

---

## Goals

- Keep analysis repos reproducible and professional for a specific study.
- Avoid committing ephemeral, environment-specific tweaks (HPC paths, emails, temporary outputs).
- Separate configuration from code using gitignored config files.
- Present repos clearly with strong READMEs, description/topics, licenses, and minimal hardcoded personal details.

---

## General Workflow for Local + HPC

- Edit and commit scripts primarily on your personal computer; push to GitHub.
- On HPC, avoid editing scripts except for legitimate bug fixes — instead, put run-specific values in gitignored config files (below).
- When you do fix a bug on HPC, commit and push, then pull locally.
- For occasional temporary changes: `git stash` → `git pull` → `git stash pop` (use sparingly).
- Keep data and outputs gitignored in both environments.

### Configuration Externalization Pattern

- Create and commit `config.example` files with placeholder variables; add the real `config` files to `.gitignore`.
  - Shell: `config.sh.example` → user copies to `config.sh` and edits for their environment (HPC vs local).
  - R: `config.R.example` → user copies to `config.R` and edits.
  - Nextflow/nf-core: `params.yml.example` → user copies to `params.yml` (path-only overrides).

- Scripts source these config files and refer to variables instead of hardcoding paths/emails.

- Benefits:
  - Identical scripts run everywhere.
  - No tiny commits for run-specific tweaks.
  - Clean pulls on HPC (no uncommitted changes blocking `git pull`).

---

## Repository 1: NicloTemp

RNA-seq analysis pipeline (QC, trimming, lane combining, mapping, alignment filtering, counting). References suggest Brugia (*Btru*) and *Shae* genomes.

### Current Structure (observed)

- `01_QC/QC_Scripts/`: trimming, FastQC/MultiQC, lane combining, sample sheet generation.
- `02_mapping/mapping_scripts/`: HISAT2 mapping, post-processing (primary-only BAMs), featureCounts, nf-core/rnaseq params and jobs.
- `.gitignore`: ignores data, outputs, logs, indexes (good).
- README.md: minimal (`# NicloTemp`).
- Visibility: private.

### Specific Suggestions

1. README (high impact)
   - Write a comprehensive README:
     - Title and one-line summary (e.g., “RNA-seq pipeline for niclosamide temperature-stress experiments in Brugia spp.”).
     - Scientific context: organism(s), experimental design, goal.
     - Pipeline overview: Raw reads → QC (FastQC/Trimmomatic) → lane combining → HISAT2 mapping → secondary alignment removal → featureCounts → downstream DE.
     - Directory structure and roles.
     - Software and versions (Trimmomatic, HISAT2, samtools, featureCounts, nf-core/rnaseq, conda envs).
     - How to run: logical steps; note cluster dependencies where relevant.
     - References to genome assemblies and literature.

2. Externalize configuration
   - Replace hardcoded absolute paths (e.g., `/temp180/mpfrende/nvincen2/...`) and email addresses in job scripts with variables from a gitignored `config.sh`.
   - Source `config.sh` at the top of job scripts and use `${PROJECT_DIR}`, `${FASTQ_DIR}`, etc.

3. Clean up duplicated/legacy files
   - Remove backup file `02_mapping/mapping_scripts/sam_secondaryRemove.job.save`.
   - Consolidate nf-core params:
     - Keep the one actually used; remove or clearly mark others as examples.
     - One param file references an unrelated organism (*A. baumannii*); remove or move it to an `examples/` folder.

4. Minor bug fix and hygiene
   - `01_QC/QC_Scripts/postTrim_multiqc.sh`: two `mv` commands incorrectly source from `fastqc52` and `trim_52` twice; adjust the second references to use `fastqc53` and `trim_53`.
   - `.gitignore`: typo — `hisat2-build_indeces/` → `hisat2-build_indices/`.

5. Repo metadata and visibility
   - Add a description and topics:
     - Description: “RNA-seq analysis pipeline for niclosamide temperature-stress experiments in Brugia spp.”
     - Topics: `rna-seq`, `bioinformatics`, `hisat2`, `trimmomatic`, `samtools`, `featurecounts`, `nf-core`, `nextflow`, `differential-expression`.
   - Add a LICENSE (MIT or GPL-3.0).
   - Make the repo public when ready.

### README Outline (starter)

- Title
- Abstract / Study context
- Pipeline overview diagram or bullet steps
- Directory structure
- Dependencies and environment setup (conda envs or modules)
- Running the pipeline (step-by-step)
- Configuration (explain `config.sh`)
- Outputs
- References (genomes, tools, papers)
- License and citation

---

## Repository 2: Metagenomics

Metabarcoding pipeline for Great Lakes zooplankton — COI (Folmer/Leray) and 18S markers — with database curation and rich phyloseq analyses.

### Current Structure (observed)

- `00_prep` → `01_QC` → `02_cutadapt` → `03_dadaASV` → `04_dadaTax` → `05_taxize` → `06_phyloseq`.
- `DB_scripts/`: BOLD database wrangling, clustering, taxonomy fixing (Python/R/Jupyter, shell).
- `nf-core_ampliseq/`: job scripts and params for Nextflow pipeline runs.
- `06_phyloseq/zoop_phyloseq.Rmd`: polished RMarkdown (master configuration).
- `06_phyloseq/plot_scripts/`: multiple R analysis scripts.
- `README.md`: minimal.
- Visibility: private.

### Specific Suggestions

1. README (high impact)
   - Write a comprehensive top-level README:
     - Title: “Great Lakes Zooplankton Metabarcoding Pipeline.”
     - Summary: markers used (COI/18S), sampling context, biological questions.
     - Pipeline overview: Raw reads → QC → demultiplexing → DADA2 ASVs → taxonomy assignment → taxize cleanup → phyloseq diversity/composition analyses.
     - Directory structure table with inputs/outputs per stage.
     - Software stack and versions (cutadapt, DADA2, nf-core/ampliseq, phyloseq, taxize, etc.).
     - Reproducibility: environment files, seeds, version pins.
     - How to reproduce: preferred “driver” (Nextflow or RMarkdown path), data locations (or sample/synthetic data).
     - Highlight `zoop_phyloseq.Rmd` as the main analysis notebook; optionally render and host via GitHub Pages.
     - References.

2. Externalize configuration
   - Many R scripts use `setwd('/Volumes/Samsung_1TB/Zooplankton/')` and absolute paths. Replace with `here::here()` and a gitignored `config.R` supplying base directories and key inputs.

3. Naming and cleanup
   - Rename `06_phyloseq/plot_scripts/Untitled.R` to something meaningful (e.g., `taxa_count_concordance.R`).
   - Fix typo: `beta_diversuty.R` → `beta_diversity.R`.
   - Reduce duplicated boilerplate across R scripts (library loads, data import):
     - Create `06_phyloseq/utils/setup.R` to load libraries and source config; `source(here::here("06_phyloseq/utils/setup.R"))` in each analysis script.

4. Subdirectory READMEs
   - Add READMEs for key stages:
     - `03_dadaASV/`: DADA2 parameters and outputs.
     - `04_dadaTax/`: taxonomy assignment approach.
     - `05_taxize/`: taxonomic cleanup strategy.
     - `DB_scripts/`: BOLD data processing overview, script purposes, expected inputs/outputs.
     - `06_phyloseq/`: analysis overview.

5. Repo metadata and visibility
   - Add a description and topics:
     - Description: “Metabarcoding analysis pipeline for Great Lakes zooplankton communities using COI and 18S markers.”
     - Topics: `metabarcoding`, `metagenomics`, `dada2`, `phyloseq`, `taxize`, `bioinformatics`, `nf-core`, `ampliseq`, `zooplankton`.
   - Add a LICENSE (MIT or GPL-3.0).
   - Make the repo public when ready.

### README Outline (starter)

- Title
- Abstract / Study context
- Pipeline overview (diagram or bullet steps)
- Directory structure
- Dependencies and environment setup
- How to run (end-to-end or stage-by-stage)
- Configuration (explain `config.R` and `params.yml`)
- Outputs (examples, key plots)
- Main analysis notebook (link to rendered HTML if available)
- References (methods, tools, databases)
- License and citation

---

## Documentation & Reproducibility Standards

- Environments:
  - Commit conda `environment.yml` (and `renv.lock` for R), Nextflow profiles; optionally Docker/Apptainer recipes.
  - Pin versions and record seeds for deterministic results.
- Releases:
  - Tag a release for publications and archive with Zenodo for DOI.
- Citation:
  - Add `CITATION.cff`.
- Data/outputs:
  - Keep large data, generated outputs, logs gitignored.
  - Optionally provide small synthetic/sample data for smoke tests.

---

## .gitignore Guidance

- Keep ignoring:
  - Large inputs (fastq/fq/fa/fna/faa/jsonl/zip), annotations (gff/gtf/gbff), outputs (pdf/csv/tsv/png/txt/html/log/json/bam/sam/ht2).
  - Tool outputs and caches (Nextflow `.nextflow/`, `work/`, `nfcore_output/`, `multiqc/`, job logs).
  - Environment-specific configs: `config.sh`, `config.R`, `params.yml`.
- Fix minor typos (e.g., `hisat2-build_indices/`).

---

## Example Config Templates (commit these as examples and gitignore the real files)

Refer to these examples from your scripts and R code. Place the real files alongside them and add to `.gitignore`.

```
- config.sh.example   → copy to config.sh and edit (HPC vs local)
- config.R.example    → copy to config.R and edit
- params.yml.example  → copy to params.yml and edit for Nextflow/nf-core
```

---

## How To Use This Context in Copilot Chats

- “Use docs/copilot-context.md as context.”
- “Based on the context file, propose a README for NicloTemp.”
- “Refactor job scripts to use config.sh variables instead of hardcoded paths.”
- “Create setup.R and update phyloseq plot scripts to source it.”
- “Suggest topics and descriptions for both repos based on this context.”

---
