#!/usr/bin/env python3
"""
Step 2: Add outgroup sequences to target sequences.

Randomly samples outgroups from the pools created in Step 1.
If a secondary outgroup file exists, splits sampling across both pools.
If no secondary file exists, all outgroups come from the primary pool.

All files are already in SILVA format from Step 1.

Usage:
    python 02_add_outgroups.py
"""

import random
import yaml
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

# ═══════════════════════════════════════════════════════════════
# 0. CONFIGURATION
# ═══════════════════════════════════════════════════════════════

with open("config/params.yaml") as f:
    cfg = yaml.safe_load(f)

n_total    = cfg.get("n_outgroups_total", 200)
seed       = cfg.get("random_seed", 42)

target_fasta      = "01_filtered/target_sequences.fasta"
outgroup_primary   = "01_filtered/primary_outgroups.fasta"
outgroup_secondary = "01_filtered/secondary_outgroups.fasta"

out_dir = Path("02_outgroups")
out_dir.mkdir(exist_ok=True)

output_fasta = out_dir / "combined_with_outgroups.fasta"
log_file     = out_dir / "outgroup_log.txt"

random.seed(seed)

# ═══════════════════════════════════════════════════════════════
# 1. LOAD TARGETS
# ═══════════════════════════════════════════════════════════════

print("Loading target sequences...")
targets = list(SeqIO.parse(target_fasta, "fasta"))
print(f"  Target sequences: {len(targets)}")

# ═══════════════════════════════════════════════════════════════
# 2. LOAD OUTGROUP POOLS
# ═══════════════════════════════════════════════════════════════

has_secondary = Path(outgroup_secondary).exists() and Path(outgroup_secondary).stat().st_size > 0

print(f"\nLoading primary outgroup pool: {outgroup_primary}")
primary_pool = list(SeqIO.parse(outgroup_primary, "fasta"))
print(f"  Primary pool size: {len(primary_pool)}")

secondary_pool = []
if has_secondary:
    print(f"Loading secondary outgroup pool: {outgroup_secondary}")
    secondary_pool = list(SeqIO.parse(outgroup_secondary, "fasta"))
    print(f"  Secondary pool size: {len(secondary_pool)}")
else:
    print("No secondary outgroup pool found. All outgroups from primary.")

# ═══════════════════════════════════════════════════════════════
# 3. STRATIFIED SAMPLING BY CLASS
# ═══════════════════════════════════════════════════════════════

def extract_class_from_header(description):
    """Extract a rough class identifier for stratification.
    Uses the 3rd rank (typically Class in SILVA taxonomy)."""
    parts = description.split(None, 1)
    if len(parts) < 2:
        return "Unknown"
    ranks = [r.strip() for r in parts[1].split(";") if r.strip()]
    # Use 3rd rank if available (typically Class), else 2nd, else 1st
    if len(ranks) >= 3:
        return ranks[2]
    elif len(ranks) >= 2:
        return ranks[1]
    elif len(ranks) >= 1:
        return ranks[0]
    return "Unknown"


def stratified_sample(pool, n_sample):
    """Sample n sequences from pool, stratified by class."""
    if n_sample <= 0:
        return []
    if len(pool) <= n_sample:
        return list(pool)

    # Group by class
    by_class = defaultdict(list)
    for record in pool:
        cls = extract_class_from_header(record.description)
        by_class[cls].append(record)

    # Calculate proportional allocation
    sampled = []
    remaining = n_sample
    classes = sorted(by_class.keys())

    for i, cls in enumerate(classes):
        members = by_class[cls]
        if i == len(classes) - 1:
            # Last class gets whatever is remaining
            n_from_class = remaining
        else:
            # Proportional allocation
            n_from_class = max(1, round(len(members) / len(pool) * n_sample))
            n_from_class = min(n_from_class, remaining, len(members))

        sampled.extend(random.sample(members, min(n_from_class, len(members))))
        remaining = n_sample - len(sampled)

        if remaining <= 0:
            break

    return sampled[:n_sample]


# Determine how many to sample from each pool
if has_secondary:
    # Split evenly between primary and secondary
    n_primary   = n_total // 2
    n_secondary = n_total - n_primary
    print(f"\nSampling {n_primary} from primary, {n_secondary} from secondary")
else:
    # All from primary
    n_primary   = n_total
    n_secondary = 0
    print(f"\nSampling all {n_primary} outgroups from primary pool")

sampled_primary   = stratified_sample(primary_pool, n_primary)
sampled_secondary = stratified_sample(secondary_pool, n_secondary)

all_outgroups = sampled_primary + sampled_secondary

print(f"  Sampled from primary:   {len(sampled_primary)}")
if has_secondary:
    print(f"  Sampled from secondary: {len(sampled_secondary)}")
print(f"  Total outgroups:        {len(all_outgroups)}")

# ═══════════════════════════════════════════════════════════════
# 4. COMBINE AND WRITE
# ═══════════════════════════════════════════════════════════════

combined = targets + all_outgroups
SeqIO.write(combined, str(output_fasta), "fasta")
print(f"\nCombined output: {output_fasta} ({len(combined)} sequences)")

# ═══════════════════════════════════════════════════════════════
# 5. LOG
# ═══════════════════════════════════════════════════════════════

# Class breakdown for primary
primary_classes = defaultdict(int)
for r in sampled_primary:
    primary_classes[extract_class_from_header(r.description)] += 1

log_lines = [
    "═══════════════════════════════════════════════════════════",
    "  Outgroup Sampling Summary",
    "═══════════════════════════════════════════════════════════",
    "",
    f"Random seed: {seed}",
    f"Total outgroups requested: {n_total}",
    "",
    f"Target sequences: {len(targets)}",
    "",
    f"Primary outgroup pool: {len(primary_pool)} sequences",
    f"  Sampled: {len(sampled_primary)}",
]

for cls in sorted(primary_classes.keys()):
    log_lines.append(f"    {cls}: {primary_classes[cls]}")

if has_secondary:
    secondary_classes = defaultdict(int)
    for r in sampled_secondary:
        secondary_classes[extract_class_from_header(r.description)] += 1

    log_lines += [
        "",
        f"Secondary outgroup pool: {len(secondary_pool)} sequences",
        f"  Sampled: {len(sampled_secondary)}",
    ]
    for cls in sorted(secondary_classes.keys()):
        log_lines.append(f"    {cls}: {secondary_classes[cls]}")
else:
    log_lines += ["", "Secondary outgroup pool: not provided"]

log_lines += [
    "",
    f"Total outgroups sampled: {len(all_outgroups)}",
    f"Total combined sequences: {len(combined)}",
    "",
    f"Output: {output_fasta}",
    "═══════════════════════════════════════════════════════════",
]

summary = "\n".join(log_lines)
print(f"\n{summary}")
with open(log_file, "w") as f:
    f.write(summary + "\n")