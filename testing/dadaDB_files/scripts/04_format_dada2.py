#!/usr/bin/env python3
"""
Step 4: Rebuild corrected taxonomy FASTA files from resolved taxonomy.

By default, uses output file 3 (resolved_taxonomy.csv) from Step 3 to look up
corrected taxonomy for each accession.

With --skiptaxonomy, skips Step 3 entirely and parses taxonomy directly from
the FASTA headers produced by Step 2 (expected format:
  >ACCESSION Kingdom;Phylum;Class;Order;Family;Genus;Species).

Outputs:
  Output 5: dada2_with_accession.fasta
      >ACCESSION Kingdom;Phylum;Class;Order;Family;Genus;Species

  Output 6: dada2_assignTaxonomy.fasta
      >Kingdom;Phylum;Class;Order;Family;Genus;Species;

  Output 7: dada2_assignTaxonomy_genus.fasta
      >Kingdom;Phylum;Class;Order;Family;Genus;

Usage:
    python 04_format_dada2.py
    python 04_format_dada2.py --skiptaxonomy
"""

import argparse
import csv
import sys
from pathlib import Path
from Bio import SeqIO

# ── Arguments ────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--skiptaxonomy", action="store_true",
                    help="Parse taxonomy from FASTA headers instead of resolved_taxonomy.csv")
args = parser.parse_args()

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
# 1. BUILD ACCESSION → TAXONOMY LOOKUP
# ═══════════════════════════════════════════════════════════════

acc_to_ranks = {}  # accession → dict of rank values

if args.skiptaxonomy:
    # Parse taxonomy directly from FASTA headers
    # Expected header format: >ACCESSION Kingdom;Phylum;Class;Order;Family;Genus;Species
    print("Loading taxonomy from FASTA headers (--skiptaxonomy mode)...")
    for record in SeqIO.parse(input_fasta, "fasta"):
        acc = record.id
        # Header description is everything after the first space
        desc = record.description[len(acc):].strip()
        parts = [p.strip() for p in desc.split(";")]
        # Pad or trim to exactly 7 ranks
        parts = parts[:7]
        while len(parts) < 7:
            parts.append("NA")
        ranks = {r: (v if v else "NA") for r, v in zip(RANK_COLS, parts)}
        acc_to_ranks[acc] = ranks
else:
    # Load corrected taxonomy from resolved_taxonomy.csv
    print("Loading resolved taxonomy from CSV...")
    with open(resolved_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            accessions = row["accessions"].split("|")
            ranks = {}
            for r in RANK_COLS:
                val = row.get(r, "").strip()
                ranks[r] = val if val and val != "NA" else "NA"
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

tax_source = "FASTA headers (--skiptaxonomy)" if args.skiptaxonomy else "resolved_taxonomy.csv"
summary = (
    f"Taxonomy source:                 {tax_source}\n"
    f"Sequences written:               {written:,}\n"
    f"Skipped (no taxonomy found):     {skipped:,}\n"
    f"\nOutput files:\n"
    f"  5. With accession:    {file_with_acc}\n"
    f"  6. assignTaxonomy:    {file_assign_tax}\n"
    f"  7. Genus-level:       {file_assign_genus}\n"
)
print(summary)
with open(log_file, "w") as f:
    f.write(summary)