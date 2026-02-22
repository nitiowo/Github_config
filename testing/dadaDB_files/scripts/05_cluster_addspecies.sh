#!/usr/bin/env bash
#
# Step 5: Cluster at 100% identity WITHIN each unique taxonomy,
# then format for dada2::addSpecies as ">AccessionID Genus species".
#
# Uses Output 5 (dada2_with_accession.fasta) from Step 4.
# Format: >ACCESSION Kingdom;Phylum;Class;Order;Family;Genus;Species
#
# Clustering is done per-taxonomy-string so identical sequences with
# different taxonomy are NEVER collapsed.
#
# Output 8: dada2_addSpecies.fasta
#     >ACCESSION Genus species
#     ACGT...
#
# Usage:
#   bash 05_cluster_addspecies.sh

set -euo pipefail

INPUT="04_dada2_format/dada2_with_accession.fasta"
OUTDIR="05_clustered"
ADDSPECIES="${OUTDIR}/dada2_addSpecies.fasta"
LOG="${OUTDIR}/cluster_log.txt"
THREADS=8

mkdir -p "${OUTDIR}"

echo "=== Counting input sequences ==="
INPUT_COUNT=$(grep -c "^>" "${INPUT}")
echo "Input sequences: ${INPUT_COUNT}"

# ── Split by taxonomy, cluster per group, merge ──────────────
echo "=== Clustering at 100% identity per taxonomy ==="

# Run Python script ONCE and capture all output
PYTHON_OUTPUT=$(python3 - "${INPUT}" "${OUTDIR}" "${ADDSPECIES}" "${THREADS}" 2>&1 <<'PYEOF'
import sys
import os
import subprocess
from collections import defaultdict
from pathlib import Path

input_fasta = sys.argv[1]
work_dir    = sys.argv[2]
output_file = sys.argv[3]
threads     = sys.argv[4]

tmp_dir = os.path.join(work_dir, "tmp_per_taxon")
os.makedirs(tmp_dir, exist_ok=True)

# ── 1. Parse input and group by full taxonomy string ─────────
# Header format: >ACCESSION Kingdom;Phylum;Class;Order;Family;Genus;Species
groups = defaultdict(list)  # taxonomy_string → [(accession, sequence), ...]

with open(input_fasta) as f:
    acc = None
    tax = None
    seq_lines = []
    
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            # Save previous record
            if acc is not None and tax is not None:
                seq = "".join(seq_lines)
                groups[tax].append((acc, seq))
            
            # Parse new header
            header = line[1:]  # Remove >
            parts = header.split(None, 1)  # Split on first whitespace
            acc = parts[0]
            tax = parts[1] if len(parts) > 1 else ""
            seq_lines = []
        else:
            seq_lines.append(line)
    
    # Don't forget the last record
    if acc is not None and tax is not None:
        seq = "".join(seq_lines)
        groups[tax].append((acc, seq))

total_input_seqs = sum(len(v) for v in groups.values())
print(f"Unique taxonomy strings: {len(groups)}")
print(f"Total sequences parsed: {total_input_seqs}")

# ── 2. Cluster each taxonomy group independently ────────────
all_centroids = []  # list of (accession, genus, species, sequence)
clustering_stats = {
    "total_before_clustering": 0,
    "total_after_clustering": 0,
    "sequences_removed_by_clustering": 0
}

for tax_str, seqs in groups.items():
    ranks = tax_str.split(";")
    
    # Extract genus and species
    genus   = ranks[5].strip() if len(ranks) > 5 else "NA"
    species = ranks[6].strip() if len(ranks) > 6 else "NA"
    
    clustering_stats["total_before_clustering"] += len(seqs)
    
    # If only one sequence, no need to cluster
    if len(seqs) == 1:
        acc, seq = seqs[0]
        all_centroids.append((acc, genus, species, seq))
        clustering_stats["total_after_clustering"] += 1
        continue
    
    # Write temp FASTA for this taxonomy group
    safe_name = tax_str.replace(";", "_").replace(" ", "_").replace("/", "_")[:100]
    tmp_in  = os.path.join(tmp_dir, f"{safe_name}_in.fasta")
    tmp_out = os.path.join(tmp_dir, f"{safe_name}_out.fasta")
    
    with open(tmp_in, "w") as fh:
        for a, s in seqs:
            fh.write(f">{a}\n{s}\n")
    
    # Run vsearch 100% clustering
    cmd = [
        "vsearch",
        "--cluster_fast", tmp_in,
        "--id", "1.0",
        "--centroids", tmp_out,
        "--threads", threads,
        "--quiet",
    ]
    subprocess.run(cmd, check=True)
    
    # Read centroids
    with open(tmp_out) as fh:
        c_acc = None
        c_seq_lines = []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if c_acc is not None:
                    all_centroids.append((c_acc, genus, species, "".join(c_seq_lines)))
                    clustering_stats["total_after_clustering"] += 1
                c_acc = line[1:]
                c_seq_lines = []
            else:
                c_seq_lines.append(line)
        if c_acc is not None:
            all_centroids.append((c_acc, genus, species, "".join(c_seq_lines)))
            clustering_stats["total_after_clustering"] += 1
    
    # Clean up
    os.remove(tmp_in)
    os.remove(tmp_out)

# Calculate clustering reduction
clustering_stats["sequences_removed_by_clustering"] = (
    clustering_stats["total_before_clustering"] - 
    clustering_stats["total_after_clustering"]
)

print(f"\nClustering complete:")
print(f"  Sequences before clustering: {clustering_stats['total_before_clustering']}")
print(f"  Sequences after clustering:  {clustering_stats['total_after_clustering']}")
print(f"  Removed by clustering:       {clustering_stats['sequences_removed_by_clustering']}")

# ── 3. Write addSpecies FASTA (Output 8) ────────────────────
# Format: >AccessionID Genus species
written = 0
skipped_na = 0

with open(output_file, "w") as out:
    for acc, genus, species, seq in all_centroids:
        # Skip if genus or species is NA/empty
        if genus in ("NA", "") or species in ("NA", ""):
            skipped_na += 1
            continue
        
        # Build binomial
        if species.startswith(genus):
            binomial = species
        else:
            binomial = f"{genus} {species}"
        
        out.write(f">{acc} {binomial}\n{seq}\n")
        written += 1

print(f"\nFinal output:")
print(f"  addSpecies sequences written: {written}")
print(f"  Skipped (NA genus/species):   {skipped_na}")

# Write stats for bash script to capture (MUST BE LAST LINE)
print(f"STATS_FOR_LOG:{clustering_stats['total_before_clustering']},"
      f"{clustering_stats['total_after_clustering']},"
      f"{clustering_stats['sequences_removed_by_clustering']},"
      f"{written},{skipped_na}")
PYEOF
)

# Display Python output to console
echo "${PYTHON_OUTPUT}"

# ── Cleanup temp directory ───────────────────────────────────
rmdir "${OUTDIR}/tmp_per_taxon" 2>/dev/null || true

# ── Extract stats from Python output ─────────────────────────
# Python prints: STATS_FOR_LOG:before,after,clustered_removed,written,na_removed
STATS_LINE=$(echo "${PYTHON_OUTPUT}" | grep "^STATS_FOR_LOG:" | cut -d: -f2)

if [ -n "${STATS_LINE}" ]; then
    IFS=',' read -r BEFORE_CLUSTER AFTER_CLUSTER REMOVED_CLUSTER WRITTEN REMOVED_NA <<< "${STATS_LINE}"
else
    # Fallback if parsing fails
    OUTPUT_COUNT=$(grep -c "^>" "${ADDSPECIES}" || echo 0)
    BEFORE_CLUSTER="${INPUT_COUNT}"
    AFTER_CLUSTER="${OUTPUT_COUNT}"
    REMOVED_CLUSTER="?"
    WRITTEN="${OUTPUT_COUNT}"
    REMOVED_NA="?"
fi

# ── Generate log file ────────────────────────────────────────
{
  echo "═════════════════════���══════════════════════════════════════"
  echo "  Clustering & addSpecies Formatting Summary"
  echo "════════════════════════════════════════════════════════════"
  echo ""
  echo "Input sequences:                    ${INPUT_COUNT}"
  echo ""
  echo "After 100% identity clustering:"
  echo "  Sequences remaining:              ${AFTER_CLUSTER}"
  echo "  Sequences removed (clustering):   ${REMOVED_CLUSTER}"
  echo ""
  echo "After filtering for valid Genus/Species:"
  echo "  Final addSpecies sequences:       ${WRITTEN}"
  echo "  Removed (NA genus/species):       ${REMOVED_NA}"
  echo ""
  echo "Total reduction:                    $(( INPUT_COUNT - WRITTEN )) sequences"
  echo ""
  echo "════════════════════════════════════════════════════════════"
  echo ""
  echo "Breakdown of reduction:"
  echo "  1. Clustering (100% identical):   ${REMOVED_CLUSTER} sequences"
  echo "  2. Missing genus/species:         ${REMOVED_NA} sequences"
  echo ""
  echo "Output file:"
  echo "  8. addSpecies: ${ADDSPECIES}"
  echo "════════════════════════════════════════════════════════════"
} > "${LOG}"

cat "${LOG}"
echo ""
echo "=== Done ==="