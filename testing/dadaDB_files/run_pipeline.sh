#!/usr/bin/env bash
#
# Master pipeline: BOLD → DADA2 reference database
#
# Usage:
#   bash run_pipeline.sh                          # run all steps
#   bash run_pipeline.sh 1,2,4,5                  # run only steps 1, 2, 4, and 5
#   bash run_pipeline.sh 1-3,5                    # run steps 1 through 3, plus step 5
#   bash run_pipeline.sh 1,2,4,5 --skiptaxonomy   # skip step 3; step 4 uses raw headers
#
# Steps:
#   1 - Filter relevant taxa
#   2 - Add outgroup sequences
#   3 - Resolve taxonomy (taxize/R)  [can be skipped with --skiptaxonomy]
#   4 - Format to dada2 FASTA
#   5 - Cluster & addSpecies format
#
set -euo pipefail

SCRIPTS="scripts"
SKIP_TAXONOMY=0

# Separate the step-selection argument from flags
STEP_ARG=""
for arg in "$@"; do
    if [[ "$arg" == "--skiptaxonomy" ]]; then
        SKIP_TAXONOMY=1
    else
        STEP_ARG="$arg"
    fi
done

# Parse step selection into an array of enabled step numbers
declare -A ENABLED_STEPS

if [ -z "$STEP_ARG" ]; then
    # No step argument: run all steps
    for s in 1 2 3 4 5; do ENABLED_STEPS[$s]=1; done
else
    # Parse comma-separated steps and ranges (e.g. "1,3-5")
    IFS=',' read -ra PARTS <<< "$STEP_ARG"
    for part in "${PARTS[@]}"; do
        if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            for s in $(seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"); do
                ENABLED_STEPS[$s]=1
            done
        elif [[ "$part" =~ ^[0-9]+$ ]]; then
            ENABLED_STEPS[$part]=1
        else
            echo "ERROR: unrecognized step specifier: '$part'" >&2
            exit 1
        fi
    done
fi

if [ "$SKIP_TAXONOMY" -eq 1 ]; then
    echo "  --skiptaxonomy: step 4 will use FASTA headers directly"
fi

echo ""
echo "  Steps to run: ${!ENABLED_STEPS[*]}"

# Required input files for each step — checked before the step runs
# Step 4 requirements change depending on --skiptaxonomy
declare -A STEP_INPUTS
STEP_INPUTS[1]="config/params.yaml"
STEP_INPUTS[2]="01_filtered/target_sequences.fasta"
STEP_INPUTS[3]="02_outgroups/combined_with_outgroups.fasta"
if [ "$SKIP_TAXONOMY" -eq 1 ]; then
    STEP_INPUTS[4]="02_outgroups/combined_with_outgroups.fasta"
else
    STEP_INPUTS[4]="02_outgroups/combined_with_outgroups.fasta 03_taxonomy_fixed/resolved_taxonomy.csv"
fi
STEP_INPUTS[5]="04_dada2_format/dada2_with_accession.fasta"

check_inputs() {
    local step=$1
    local desc=$2
    local missing=0
    for f in ${STEP_INPUTS[$step]}; do
        if [ ! -f "$f" ]; then
            echo "  ERROR: Step ${step} (${desc}) requires '${f}' but it does not exist." >&2
            missing=1
        fi
    done
    if [ "$missing" -eq 1 ]; then
        echo "  Ensure all prerequisite steps have been run first." >&2
        exit 1
    fi
}

run_step() {
    local step=$1
    local desc=$2
    local cmd=$3
    if [ -n "${ENABLED_STEPS[$step]+x}" ]; then
        check_inputs "${step}" "${desc}"
        echo ""
        echo "╔══════════════════════════════════════════════════╗"
        echo "║  Step ${step}: ${desc}"
        echo "╚══════════════════════════════════════════════════╝"
        echo ""
        eval "${cmd}"
    else
        echo "  Skipping step ${step}: ${desc}"
    fi
}

STEP4_FLAG=""
[ "$SKIP_TAXONOMY" -eq 1 ] && STEP4_FLAG="--skiptaxonomy"

run_step 1 "Filter relevant taxa"            "python3 ${SCRIPTS}/01_filter_taxa.py"
run_step 2 "Add outgroup sequences"           "python3 ${SCRIPTS}/02_add_outgroups.py"
run_step 3 "Resolve taxonomy (taxize/R)"      "Rscript ${SCRIPTS}/03_resolve_taxonomy.R"
run_step 4 "Format to dada2 FASTA"            "python3 ${SCRIPTS}/04_format_dada2.py ${STEP4_FLAG}"
run_step 5 "Cluster & addSpecies format"      "bash    ${SCRIPTS}/05_cluster_addspecies.sh"

echo ""
echo "════════════════════════════════════════════════════"
echo "  Pipeline complete!"
echo ""
echo "  Output files:"
echo "    03_taxonomy_fixed/"
echo "      1. all_taxa_with_ids.csv        - Raw parsed taxonomy"
echo "      2. unique_taxa.csv              - Unique taxa with IDs"
echo "      3. resolved_taxonomy.csv        - Corrected 7-rank taxonomy"
echo "      4. resolution_log.csv           - Per-attempt fate log"
echo "    04_dada2_format/"
echo "      5. dada2_with_accession.fasta   - Corrected + accession IDs"
echo "      6. dada2_assignTaxonomy.fasta   - assignTaxonomy format"
echo "      7. dada2_assignTaxonomy_genus.fasta - Genus-level only"
echo "    05_clustered/"
echo "      8. dada2_addSpecies.fasta       - addSpecies format"
echo "════════���═══════════════════════════════════════════"