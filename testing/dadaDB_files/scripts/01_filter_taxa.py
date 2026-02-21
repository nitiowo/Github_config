#!/usr/bin/env python3
"""
Step 1: Filter target taxa from reference database (SILVA or BOLD format).

Reads the input FASTA in whatever format is specified in config, extracts
target taxa and outgroup pools, and writes ALL outputs in SILVA format:
    >ACCESSION Rank1;Rank2;Rank3;...;RankN

This normalization means scripts 02-05 only ever see SILVA-formatted headers.

Usage:
    python 01_filter_taxa.py
"""

import random
import re
import sys
import yaml
from pathlib import Path
from Bio import SeqIO

# ═══════════════════════════════════════════════════════════════
# 0. CONFIGURATION
# ═══════════════════════════════════════════════════════════════

with open("config/params.yaml") as f:
    cfg = yaml.safe_load(f)

input_primary   = cfg["input_primary"]
input_secondary = cfg.get("input_secondary", "")
target_file     = cfg["target_classes_file"]
db_format       = cfg.get("database_format", "silva").lower()

out_dir = Path("01_filtered")
out_dir.mkdir(exist_ok=True)

file_targets     = out_dir / "target_sequences.fasta"
file_outgroups_1 = out_dir / "primary_outgroups.fasta"
file_outgroups_2 = out_dir / "secondary_outgroups.fasta"
file_log         = out_dir / "filter_log.txt"

# ═══════════════════════════════════════════════════════════════
# 1. HEADER PARSING → SILVA FORMAT
# ═══════════════════════════════════════════════════════════════

def parse_silva_header(description):
    """Parse SILVA header: >ACCESSION Rank1;Rank2;...;RankN"""
    parts = description.split(None, 1)
    accession = parts[0]
    taxonomy = parts[1] if len(parts) > 1 else ""
    ranks = [r.strip() for r in taxonomy.rstrip(";").split(";") if r.strip()]
    return accession, ranks


def parse_silva_notax_header(description, _counter=[0]):
    """Parse taxonomy-only SILVA header: >Rank1;Rank2;...;RankN (no accession).
    Synthesises a sequential placeholder accession (seq_00001, seq_00002, ...)."""
    _counter[0] += 1
    accession = f"seq_{_counter[0]:05d}"
    ranks = [r.strip() for r in description.rstrip(";").split(";") if r.strip()]
    return accession, ranks


def parse_bold_header(description):
    """Parse BOLD header: >AANIC174-10|COI-5P|Australia|Animalia,Arthropoda,..."""
    parts = description.split("|")
    accession = parts[0].strip()

    if len(parts) < 4:
        return accession, []

    taxonomy = parts[3].strip()
    ranks = [r.strip() for r in taxonomy.split(",") if r.strip() and r.strip() != "None"]
    return accession, ranks


def parse_header(description, fmt="silva"):
    """Universal parser. Returns (accession, [rank1, rank2, ...]).
    Supported formats: 'silva', 'silva_notax', 'bold'."""
    if fmt == "bold":
        return parse_bold_header(description)
    elif fmt == "silva_notax":
        return parse_silva_notax_header(description)
    else:
        return parse_silva_header(description)


def to_silva_header(accession, ranks):
    """Convert accession + ranks to SILVA format string (no leading >)."""
    tax_str = ";".join(ranks)
    return f"{accession} {tax_str}"


# ═══════════════════════════════════════════════════════════════
# 2. LOAD TARGET CLASSES
# ═══════════════════════════════════════════════════════════════

with open(target_file) as f:
    target_classes = set(line.strip() for line in f if line.strip() and not line.startswith("#"))

print(f"Target classes ({len(target_classes)}): {', '.join(sorted(target_classes))}")
print(f"Database format: {db_format.upper()}")

# ═══════════════════════════════════════════════════════════════
# 3. PROCESS PRIMARY FASTA
# ═══════════════════════════════════════════════════════════════

def process_fasta(fasta_path, fmt):
    """Parse a FASTA file, return lists of (record, accession, ranks) for targets and outgroups."""
    targets = []
    outgroups = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        accession, ranks = parse_header(record.description, fmt)

        # Check if any rank matches a target class
        is_target = any(r in target_classes for r in ranks)

        # Rewrite record header to SILVA format
        record.id = accession
        record.name = accession
        record.description = to_silva_header(accession, ranks)

        if is_target:
            targets.append(record)
        else:
            outgroups.append(record)

    return targets, outgroups


print(f"\nProcessing primary file: {input_primary}")
primary_targets, primary_outgroups = process_fasta(input_primary, db_format)
print(f"  Target sequences:   {len(primary_targets)}")
print(f"  Outgroup sequences: {len(primary_outgroups)}")

# ═══════════════════════════════════════════════════════════════
# 4. PROCESS SECONDARY FASTA (if provided)
# ═══════════════════════════════════════════════════════════════

secondary_targets = []
secondary_outgroups = []
has_secondary = bool(input_secondary and input_secondary.strip() and Path(input_secondary).exists())

if has_secondary:
    print(f"\nProcessing secondary file: {input_secondary}")
    secondary_targets, secondary_outgroups = process_fasta(input_secondary, db_format)
    print(f"  Target sequences:   {len(secondary_targets)}")
    print(f"  Outgroup sequences: {len(secondary_outgroups)}")
else:
    print("\nNo secondary file provided. All outgroups will come from primary file.")

# ═══════════════════════════════════════════════════════════════
# 5. WRITE OUTPUTS (all in SILVA format)
# ═══════════════════════════════════════════════════════════════

# Combine targets from both files
all_targets = primary_targets + secondary_targets
SeqIO.write(all_targets, str(file_targets), "fasta")
print(f"\nTargets written: {file_targets} ({len(all_targets)} sequences)")

# Write primary outgroups
SeqIO.write(primary_outgroups, str(file_outgroups_1), "fasta")
print(f"Primary outgroups written: {file_outgroups_1} ({len(primary_outgroups)} sequences)")

# Write secondary outgroups (if any)
if has_secondary:
    SeqIO.write(secondary_outgroups, str(file_outgroups_2), "fasta")
    print(f"Secondary outgroups written: {file_outgroups_2} ({len(secondary_outgroups)} sequences)")

# ═══════════════════════════════════════════════════════════════
# 6. LOG
# ═══════════════════════════════════════════════════════════════

log_lines = [
    f"Database format: {db_format.upper()}",
    f"Target classes ({len(target_classes)}): {', '.join(sorted(target_classes))}",
    "",
    f"Primary file: {input_primary}",
    f"  Total sequences:    {len(primary_targets) + len(primary_outgroups)}",
    f"  Target sequences:   {len(primary_targets)}",
    f"  Outgroup pool:      {len(primary_outgroups)}",
]

if has_secondary:
    log_lines += [
        "",
        f"Secondary file: {input_secondary}",
        f"  Total sequences:    {len(secondary_targets) + len(secondary_outgroups)}",
        f"  Target sequences:   {len(secondary_targets)}",
        f"  Outgroup pool:      {len(secondary_outgroups)}",
    ]
else:
    log_lines += ["", "Secondary file: not provided"]

log_lines += [
    "",
    f"Combined targets: {len(all_targets)}",
    "",
    "NOTE: All output files are in SILVA format regardless of input format.",
    f"  Header format: >ACCESSION Rank1;Rank2;...;RankN",
    f"  (For 'silva_notax' inputs, ACCESSION is a synthetic placeholder: seq_00001, seq_00002, ...)",
]

summary = "\n".join(log_lines)
print(f"\n{summary}")
with open(file_log, "w") as f:
    f.write(summary + "\n")