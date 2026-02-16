# Metagenomics Repo Cleanup

## Repo structure (current)

```
.gitignore
README.md
Metagenomics.Rproj
Data/
  README.md
00_prep/
  prep_scripts/
    sampsheet.sh
01_QC/
  README.md
  QC_scripts/
    fastQC.sh
    multiqc.sh
    fullQC.job
02_cutadapt/
  README.md
  02a_reverse_demux/02a_reverse_demux_scripts/  (3 scripts)
  02b_forward_demux/...
  02c_rescue_demux/...
  demux_scripts/         (samplesheet_gen.sh, post_demux_qc.job)
  demux_total_scripts/   (demux_submit_array.job, TODO.md)
03_dadaASV/
  README.md  (EMPTY)
  dada_asv_scripts/
    dada2_asv.R
    dada2_asv_SE.R
    dada_asv.job
04_dadaTax/
  README.md  (EMPTY)
  dada_tax_scripts/
    dada_tax_assign.R
    dada_tax.job
  ref_DBs/   (empty or gitignored)
05_taxize/
  README.md  (EMPTY)
06_phyloseq/
  README.md  (EMPTY)
  setup.R
  (zoop_phyloseq.Rmd — on a different branch/commit)
```

---

## Part 1: Structural cleanup

### Commit 1 — Flatten unnecessary nested script dirs
**Why:** `01_QC/QC_scripts/`, `00_prep/prep_scripts/`, `03_dadaASV/dada_asv_scripts/`, `04_dadaTax/dada_tax_scripts/` are pointless nesting. Every step already has its own numbered folder. Having `03_dadaASV/dada_asv_scripts/` just means extra `../` everywhere. Move scripts up one level.

```
00_prep/sampsheet.sh           (was prep_scripts/sampsheet.sh)
01_QC/fastQC.sh                (was QC_scripts/fastQC.sh)
01_QC/multiqc.sh
01_QC/fullQC.job
03_dadaASV/dada2_asv.R         (was dada_asv_scripts/dada2_asv.R)
03_dadaASV/dada2_asv_SE.R
03_dadaASV/dada_asv.job
04_dadaTax/dada_tax_assign.R   (was dada_tax_scripts/dada_tax_assign.R)
04_dadaTax/dada_tax.job
```

Update any relative paths inside the scripts after moving (e.g. `../primers/` → `./primers/` or whatever is now correct).

**Commit message:** `flatten script subdirectories in 00-04`

### Commit 2 — Clean up 02_cutadapt structure
**Why:** You have `demux_scripts/`, `demux_total_scripts/`, AND `02a_reverse_demux/02a_reverse_demux_scripts/`. The naming is redundant (`02a_reverse_demux/02a_reverse_demux_scripts/`). Either:
- Keep the `02a/02b/02c` subdirs but rename the inner script dirs to just `scripts/`, OR
- Flatten entirely so `02_cutadapt/` has `reverse_demux.sh`, `forward_demux.sh`, `rescue.sh` etc.

Pick whichever matches how you actually ran things. The `demux_scripts/` and `demux_total_scripts/` folders look like drafts/alternatives — if they're superseded by the 02a/02b/02c structure, delete them.

Also **delete `demux_total_scripts/TODO.md`** — it's a to-do list for yourself, not analysis documentation.

**Commit message:** `clean up cutadapt directory structure, remove TODO`

### Commit 3 — Rename `Data/` to `data/`
Lowercase is standard convention. Capital `Data/` stands out.

**Commit message:** `rename Data/ to data/`

---

## Part 2: READMEs

### Commit 4 — Rewrite root README.md
Current content is just "Metagenomics project scripts and resources". Replace with something like:

```markdown
# Great Lakes Zooplankton Metabarcoding

Scripts for processing zooplankton metabarcoding samples from the Great Lakes.
Three primer sets were used (Folmer COI, Leray COI, 18S rRNA) alongside
morphological identification.

## Pipeline steps

| Dir | Step |
|-----|------|
| `00_prep` | Generate sample lists from raw fastq files |
| `01_QC` | FastQC + MultiQC on raw reads |
| `02_cutadapt` | Demultiplex by primer and trim adapters |
| `03_dadaASV` | DADA2 ASV inference |
| `04_dadaTax` | DADA2 taxonomic assignment |
| `05_taxize` | Taxonomic cleanup with taxize |
| `06_phyloseq` | Community analysis in R (phyloseq) |

Raw data and large outputs are gitignored. See `.gitignore` for details.
```

**Commit message:** `rewrite root README`

### Commit 5 — Fix existing READMEs, add missing ones
**`01_QC/README.md`** — Currently has busted markdown (heading levels are backwards: `##` then `#`). Also says "simplifies" (typo). Rewrite:

```markdown
# Quality control

Ran FastQC on all raw reads, then summarised with MultiQC.
`sampsheet.sh` (in `00_prep/`) generates the sample list first.

- `fastQC.sh` — runs fastqc on everything in the raw data dir
- `multiqc.sh` — aggregates fastqc output
- `fullQC.job` — SGE job that runs both
```

**`02_cutadapt/README.md`** — Has a typo ("Demulitplexing"). Rewrite briefly:

```markdown
# Primer demultiplexing and trimming

Cutadapt was used in 3 steps:
1. `02a_reverse_demux/` — separate by reverse primer (HCO2198 vs SSUR)
2. `02b_forward_demux/` — split HCO2198 reads by forward primer (LCO1490 vs mICOintF)
3. `02c_rescue_demux/` — check unmatched reads for anything missed in step 1
```

**`03_dadaASV/README.md`** — Currently empty. Add:

```markdown
# ASV inference (DADA2)

Ran DADA2 on each primer set separately.
- `dada2_asv.R` — paired-end pipeline
- `dada2_asv_SE.R` — single-end variant (used for LCO1490 which didn't merge well)
- `dada_asv.job` — SGE submission script
```

**`04_dadaTax/README.md`** — Empty. Add:

```markdown
# Taxonomic assignment (DADA2)

Used `assignTaxonomy` + `addSpecies` against a custom reference database.
- `dada_tax_assign.R` — takes a seqtab RDS and ref fastas, outputs taxa table
- `dada_tax.job` — SGE submission script
- `ref_DBs/` — reference databases (gitignored, too large)
```

**`05_taxize/README.md`** — Empty. Add whatever is actually in there, or:

```markdown
# Taxonomic cleanup

Used the `taxize` R package to standardise names from DADA2 output.
```

**`06_phyloseq/README.md`** — Empty. Add:

```markdown
# Community analysis

- `setup.R` — loads DADA2 outputs and metadata into phyloseq objects
- `zoop_phyloseq.Rmd` — main analysis notebook (diversity, ordination, etc.)
```

**`Data/README.md`** — Currently "Put raw data in this folder". Change to:

```markdown
# Data

Morphological count/biomass data and sample metadata.
Raw sequence data is not tracked (too large).
```

**`00_prep/`** — No README. Add:

```markdown
# Prep

- `sampsheet.sh` — generates `samples.txt` from filenames in the raw data directory
- `raw_dir.txt` — path to raw fastq location (not tracked)
```

**Commit message:** `add and fix READMEs for all pipeline steps`

---

## Part 3: AI tells in shell scripts

### Commit 6 — Clean up shell script comments

**`02_cutadapt/demux_scripts/samplesheet_gen.sh`** — The `abspath()` helper with three fallbacks (realpath → readlink → python3) is over-engineered for a personal analysis script. No grad student writes that. Replace with just `realpath "$1"` or even hardcoded paths. Also the comments read like documentation for a tool ("usage:", formal function docs). Simplify to just inline comments.

**`02_cutadapt/demux_total_scripts/demux_submit_array.job`** — The block comment at the top reads like generated docs:
```
# Submit hierarchical demux pipeline as an SGE array
# Usage: qsub submit_array.sh
#
# Notes:
#  - Make sure "samples.txt" is in ...
#  - Make sure reverse_demux.sh, forward_demux.sh, rescue.sh ...
```
Real lab scripts don't have structured "Notes:" blocks. Trim to 1-2 lines max.

**`01_QC/QC_scripts/fastQC.sh`** and **`multiqc.sh`** — Comments like `# Load fastqc (eg; HPC module or conda) if necessary` and `# Deactivate fastqc` / `# conda deactivate if necessary` are hand-holdy. Just leave the commented conda lines and drop the explanations.

**`00_prep/sampsheet.sh`** — `# This script generates a sample sheet from your raw data directory` + `# Usage: bash sampsheet.sh` is a minor flag. Change to just `# make sample list from fastq filenames` or similar.

**Across all scripts**: the phrase "if necessary" appears in several places — it's a pattern of hedging that reads like AI. Cut it.

**Commit message:** `trim comments in shell scripts`

---

## Part 4: AI tells in R scripts

### Commit 7 — Clean up dada2 R scripts

**`03_dadaASV/dada2_asv.R` and `dada2_asv_SE.R`:**

- `# This is inspired by https://benjjneb.github.io/dada2/tutorial_1_8.html` — fine, but the `# Usage:` line underneath reads tool-like. Change to just the link.
- The commented-out install block at top is fine (common pattern).
- `# Modify and add new parameter based on your data.` — vague and instructional. Replace with a comment about what you actually used, e.g. `# truncLen chosen based on quality profiles`.
- `# Change based on sample names` — same problem. Replace with what your naming convention actually is.
- Section dividers like `######################################################################################################` are fine but the length is excessive. Shorten to ~40 chars or use `# ---`.

**`04_dadaTax/dada_tax_assign.R`:**
- Clean, mostly fine. The `##### SCRIPT MODE #####` / `##### INTERACTIVE MODE #####` blocks are a legit pattern. Keep.
- `# The following line replaces any instances of NNNNNNNNNN with nothing` — fine but could just be `# strip Ns from justConcatenate joins`.

**Commit message:** `clean up R script comments in 03-04`

### Commit 8 — Fix setup.R

**`06_phyloseq/setup.R`:**
- `# This script loads dada2 output, standardizes metadata and phyloseq formats, and outputs phyloseq objects` — rewrite as `# build phyloseq objects from dada2 output + metadata`
- The hardcoded `setwd("/Volumes/Samsung_1TB/Zooplankton/Metagenomics/")` is fine (it's real), but note it won't work for anyone else. Add a comment like `# local path — adjust if needed`
- The script seems unfinished — it loads `morph_otu` but never creates `ps_morph` or saves anything. Either finish it or add a `# TODO: finish morph phyloseq + saveRDS` at the bottom.

**Commit message:** `tidy setup.R comments, note unfinished morph section`

---

## Part 5: AI tells in zoop_phyloseq.Rmd (the big one)

### Commit 9 — Strip "Master Configuration" framing

The entire config block screams AI:
- `★★★  MASTER CONFIGURATION  ★★★` — no human writes this
- `Change these values to instantly re-run the entire document with different settings.`
- `Each parameter is used throughout the notebook.`
- `**Change parameters here to control every downstream analysis and plot.**`
- Markdown bold + the star emoji

Replace the section header with just `## Parameters` and the code comment with:
```r
# ---- analysis parameters ----
```
Delete the markdown paragraph above the chunk entirely.

Kill all the `# Options: "ASV", or any element of tax_ranks`-style option menus in the comments. Those are for tool users. Just state what you used and why.

**Commit message:** `replace config block framing in phyloseq Rmd`

### Commit 10 — Rewrite helper function comments

Current style:
```r
# ─────────────────────────────────────────────────
#  Utility functions used throughout the document
# ─────────────────────────────────────────────────
```

Replace with `# --- helpers ---` or just delete. The box-drawing characters are a tell.

Other specific flags:
- `# Set Lake as ordered factor (West → East) + Mesh as factor` — the arrow character `→` is unusual. Use `->` or just write `W to E`.
- `# Agglomerate to a taxonomic rank; "ASV" = no agglomeration` — the semicolon usage is a GPT pattern. Just `# aggregate to rank (skip if ASV)`.
- `# Subset taxa by the global \`taxa_subset\` config (or a supplied list)` — backtick-quoting a variable name in a comment is a docs pattern. Just `# subset taxa if taxa_subset is set`.
- `# One-step "prepare" wrapper` — scare-quoting in comments is a tell. Just `# prep wrapper`.
- `# Morph-available sample names` — compound-adjective hyphenation in comments is unusual for code. `# samples that have morph data`.

**Commit message:** `simplify helper function comments`

### Commit 11 — Strip the verbose prose blocks

Throughout the Rmd, there are long narrative paragraphs between chunks like:

> "To evaluate whether the markers agree on which sites are more or less diverse, I computed Spearman rank correlations of site-level alpha diversity between every pair of methods."

This is way too polished for an analysis notebook. Either:
- Delete them entirely (the code + section headers are self-explanatory), OR
- Reduce each to 1 short sentence max

Also kill any `> **Note:**` blockquotes — that's a documentation pattern.

**Commit message:** `trim prose in phyloseq Rmd`

### Commit 12 — Fix comment style inside code chunks

Patterns to find-and-replace:
- `# ── text ──` (em-dash box style) → `# text` or `# --- text ---`
- Unicode arrows `→` → `->` or just words
- Comments that explain what the *next line does* in full English sentences — e.g. `# Identify top‑N taxa by total relative abundance` before a dplyr chain. Trim to `# top N taxa` or delete.
- `# 95% confidence ellipses when ≥ 2 groups with ≥ 3 points each` — the `≥` symbols are a tell. Use `>=`.
- Star emoji `★` anywhere

**Commit message:** `normalise comment formatting in Rmd`

### Commit 13 — Remove the parametric tool patterns

The Rmd has patterns designed for reusability that make no sense for a one-off analysis:
- The `taxa_subset` system with `list(rank = "Phylum", include = "Rotifera")` examples
- `active_markers` as a configurable vector
- `analysis_rank` as a global toggle
- `comparison_var` as a string that gets injected into formulas
- Functions that accept `rank`, `group_var`, `distance`, etc. as arguments

You don't need to rip all of this out — some of it is genuinely useful. But:
1. Remove the commented-out option menus (the `# Options:` lines)
2. Remove `taxa_subset` and the `subset_taxa_custom()` function if you never actually used subsetting
3. Hardcode values where you only ever use one setting (e.g. if you always use `"Lake"`, don't pass `comparison_var` — just write `"Lake"`)

**Commit message:** `remove unused configurability from Rmd`

---

## Part 6: Config files

### Commit 14 — Add config.sh to 02_cutadapt

The cutadapt scripts have primer paths and directory paths scattered across multiple files. Add a `02_cutadapt/config.sh`:

```bash
# paths used across demux steps
RAW_DIR=$(cat ../00_prep/raw_dir.txt)
PRIMER_REV=./primers/reverse_primers.fa
PRIMER_FWD_HCO=./primers/forward_primers_HCO.fa
OUTBASE=./demux_out
```

Then source it from the individual scripts: `source ../config.sh` (adjust relative path).

**Don't** do this for `03_dadaASV` or `04_dadaTax` — those already take paths as CLI arguments via argparser, which is fine.

**Commit message:** `add config.sh for cutadapt paths`

---

## Part 7: Gitignore and loose ends

### Commit 15 — Clean up .gitignore

The current gitignore is thorough but has some issues:
- `**/*output*` will match any file/dir with "output" in the name — this is probably too aggressive. You might want `*_output/` instead.
- `**/*txt` ignores ALL .txt files including `samples.txt` and `raw_dir.txt` which you probably want tracked. Either use `!00_prep/samples.txt` exceptions or be more specific with the patterns.
- `**/*csv` and `**/*tsv` ignore metadata CSVs you probably want tracked (like `zoop96_metadata.csv`). Add exceptions or restructure.
- `**/*messaround` and `**/*helper` and `**/*random` — these are dev patterns, fine to keep but they look a bit odd.
- Remove the duplicate log section (you have both `**/*log` under "Output files" and `**/*log*` under "Log files").

**Commit message:** `clean up gitignore, add exceptions for tracked data files`

### Commit 16 — Remove hardcoded absolute paths from job files

`03_dadaASV/dada_asv.job` and `04_dadaTax/dada_tax.job` have full paths like `/temp180/mpfrende/nvincen2/Metagenomics/...`. These are fine as a record of what you ran, but add a comment like `# ran on CRC cluster` at the top so it's clear these are specific to your setup.

Also `fullQC.job` has your email `nvincen2@nd.edu` — that's fine and actually makes it look more real. Keep it.

**Commit message:** `annotate job scripts with cluster info`

---

## Summary: commit order

| # | Message | What |
|---|---------|------|
| 1 | `flatten script subdirectories in 00-04` | Move scripts up, delete empty dirs |
| 2 | `clean up cutadapt directory structure, remove TODO` | Delete drafts + TODO.md |
| 3 | `rename Data/ to data/` | Lowercase |
| 4 | `rewrite root README` | Actual project description |
| 5 | `add and fix READMEs for all pipeline steps` | All subdirectory READMEs |
| 6 | `trim comments in shell scripts` | Remove tool-like/hedging language |
| 7 | `clean up R script comments in 03-04` | Section dividers, instructional comments |
| 8 | `tidy setup.R comments, note unfinished morph section` | setup.R specifically |
| 9 | `replace config block framing in phyloseq Rmd` | Kill stars, master config prose |
| 10 | `simplify helper function comments` | Box chars, GPT punctuation patterns |
| 11 | `trim prose in phyloseq Rmd` | Delete/shorten narrative paragraphs |
| 12 | `normalise comment formatting in Rmd` | Unicode, em-dashes, sentence comments |
| 13 | `remove unused configurability from Rmd` | Kill option menus, dead params |
| 14 | `add config.sh for cutadapt paths` | Centralise paths |
| 15 | `clean up gitignore, add exceptions for tracked data files` | Fix overly broad patterns |
| 16 | `annotate job scripts with cluster info` | Add context to hardcoded paths |

---

## Quick-reference: AI tells cheat sheet

| Tell | Where | Fix |
|------|-------|-----|
| `★` emoji, `★★★ MASTER CONFIGURATION ★★★` | Rmd config block | Delete |
| `# ── text ──` em-dash box comments | Rmd everywhere | Use `# --- text` or plain `#` |
| `─────────` box-drawing lines | Rmd helpers section | Delete or use `# ---` |
| Unicode `→` `≥` `‑` in comments | Rmd throughout | Use `->` `>=` `-` |
| `# Options: "X", "Y", or "Z"` menus | Rmd config block | Delete |
| `(eg; HPC module or conda) if necessary` | Shell QC scripts | Delete |
| `abspath()` with 3 fallbacks | samplesheet_gen.sh | Use just `realpath` |
| Structured `# Notes:` blocks | demux_submit_array.job | Trim to 1 line |
| `# This script generates a sample sheet from your raw data directory` | sampsheet.sh | Shorten |
| `# Modify and add new parameter based on your data.` | dada2_asv.R | State what you used |
| `# Change based on sample names` | dada2_asv.R | State your convention |
| Formal `> **Note:**` blockquotes | Rmd | Delete |
| Multi-sentence prose paragraphs between chunks | Rmd | Delete or 1 line max |
| Backtick-quoting variable names in `#` comments | Rmd helpers | Just name it |
| Compound-adjective hyphenation (`Morph-available`) | Rmd | Write normally |
| Semicolons in comments (`rank; "ASV" = no agglomeration`) | Rmd | Rephrase |