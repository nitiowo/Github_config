#!/usr/bin/env python3
"""
Step 4: Rebuild corrected taxonomy FASTA files from resolved taxonomy.

Uses output file 3 (resolved_taxonomy.csv) from Step 3 to look up
corrected taxonomy for each accession, then writes:

  Output 5: dada2_with_accession.fasta
      >ACCESSION Kingdom;Phylum;Class;Order;Family;Genus;Species
      ACGT...

  Output 6: dada2_assignTaxonomy.fasta
      >Kingdom;Phylum;Class;Order;Family;Genus;Species;
      ACGT...

  Output 7: dada2_assignTaxonomy_genus.fasta
      >Kingdom;Phylum;Class;Order;Family;Genus;
      ACGT...

Usage:
    python 04_format_dada2.py
"""

import csv
import sys
from pathlib import Path
from Bio import SeqIO

# ── Paths ─────────────────────────────────────────────────────
input_fasta       = "02_outgroups/combined_with_outgroups.fasta"
resolved_csv      = "03_taxonomy_fixed/resolved_taxonomy.csv"
output_dir        = Path("04_dada2_format")
output_dir.mkdir(exist_ok=True)

file_with_acc     = output_dir / "dada2_with_accession.fasta"       # Output 5
file_assign_tax   = output_dir / "dada2_assignTaxonomy.fasta"       # Output 6
file_assign_genus = output_dir / "dada2_assignTaxonomy_genus.fasta" # Output 7
log_file          = output_dir / "format_log.txt"

RANK_COLS  = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
GENUS_COLS = RANK_COLS[:6]  # Up to Genus

# ═══════════════════════════════════════════════════════════════
# 1. BUILD ACCESSION → CORRECTED TAXONOMY LOOKUP
# ═══════════════════════════════════════════════════════════════

print("Loading resolved taxonomy...")

# Read resolved_taxonomy.csv and expand accession lists
# Each row has: taxa_id, Kingdom..Species, source_db, accessions (pipe-delimited)
acc_to_ranks = {}  # accession → dict of rank values

with open(resolved_csv, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Split pipe-delimited accession list
        accessions = row["accessions"].split("|")
        
        # Build rank dict for this taxa_id
        ranks = {}
        for r in RANK_COLS:
            val = row.get(r, "").strip()
            ranks[r] = val if val and val != "NA" else "NA"
        
        # Map every accession to these ranks
        for acc in accessions:
            acc = acc.strip()
            if acc:
                acc_to_ranks[acc] = ranks

print(f"Loaded taxonomy for {len(acc_to_ranks):,} accessions")

# ═══════════════════════════════════════════════════════════════
# 2. REBUILD FASTA FILES
# ═══════════════════════════════════════════════════════════════

print("Writing output FASTA files...")

written  = 0
skipped  = 0

with open(file_with_acc, "w") as out5, \
     open(file_assign_tax, "w") as out6, \
     open(file_assign_genus, "w") as out7:

    for record in SeqIO.parse(input_fasta, "fasta"):
        acc = record.id

        if acc not in acc_to_ranks:
            skipped += 1
            continue

        ranks = acc_to_ranks[acc]

        # Clean sequence: remove gaps and dots, uppercase
        seq_clean = str(record.seq).replace("-", "").replace(".", "").upper()

        # Build taxonomy strings
        full_tax  = ";".join(ranks.get(r, "NA") for r in RANK_COLS)
        genus_tax = ";".join(ranks.get(r, "NA") for r in GENUS_COLS)

        # Output 5: >ACCESSION Kingdom;Phylum;Class;Order;Family;Genus;Species
        out5.write(f">{acc} {full_tax}\n{seq_clean}\n")

        # Output 6: >Kingdom;Phylum;Class;Order;Family;Genus;Species;
        out6.write(f">{full_tax};\n{seq_clean}\n")

        # Output 7: >Kingdom;Phylum;Class;Order;Family;Genus;
        out7.write(f">{genus_tax};\n{seq_clean}\n")

        written += 1

# ═══════════════════════════════════════════════════════════════
# 3. LOG
# ═══════════════════════════════════════════════════════════════

summary = (
    f"Sequences written:              {written:,}\n"
    f"Skipped (not in resolved table): {skipped:,}\n"
    f"\nOutput files:\n"
    f"  5. With accession:    {file_with_acc}\n"
    f"  6. assignTaxonomy:    {file_assign_tax}\n"
    f"  7. Genus-level:       {file_assign_genus}\n"
)
print(summary)
with open(log_file, "w") as f:
    f.write(summary)