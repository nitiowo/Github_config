#!/usr/bin/env bash
#
# Master pipeline: SILVA → DADA2 reference database
#
# Usage:
#   bash run_pipeline.sh [step_number]
#   bash run_pipeline.sh          # run all steps
#   bash run_pipeline.sh 3        # run from step 3 onward
#
set -euo pipefail

START_STEP=${1:-1}
SCRIPTS="scripts"

run_step() {
    local step=$1
    local desc=$2
    local cmd=$3
    if [ "${step}" -ge "${START_STEP}" ]; then
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

run_step 1 "Filter relevant taxa"            "python3 ${SCRIPTS}/01_filter_taxa.py"
run_step 2 "Add outgroup sequences"           "python3 ${SCRIPTS}/02_add_outgroups.py"
run_step 3 "Resolve taxonomy (taxize/R)"      "Rscript ${SCRIPTS}/03_resolve_taxonomy.R"
run_step 4 "Format to dada2 FASTA"            "python3 ${SCRIPTS}/04_format_dada2.py"
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