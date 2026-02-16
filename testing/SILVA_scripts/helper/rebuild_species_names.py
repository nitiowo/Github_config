#!/usr/bin/env python3
"""
Standalone script to rebuild dada2 FASTA files with corrected species names.

Uses already-resolved taxonomy from Step 3 but re-applies species cleaning rules.

Inputs:
  - 03_taxonomy_fixed/resolved_taxonomy.csv  (has resolved higher ranks + accessions)
  - 03_taxonomy_fixed/resolution_log.csv     (has raw species names)
  - 04_dada2_format/dada2_with_accession.fasta (has sequences + accessions)

Outputs:
  - rebuilt/dada2_with_accession.fasta       (Output 5)
  - rebuilt/dada2_assignTaxonomy.fasta       (Output 6)
  - rebuilt/dada2_assignTaxonomy_genus.fasta (Output 7)

Usage:
    python rebuild_species_names.py
"""

import csv
import re
from pathlib import Path
from Bio import SeqIO

# ════════════════════════════════════════════════════���══════════
# 0. SETUP
# ═══════════════════════════════════════════════════════════════

resolved_csv = "03_taxonomy_fixed/resolved_taxonomy.csv"
log_csv      = "03_taxonomy_fixed/resolution_log.csv"
input_fasta  = "04_dada2_format/dada2_with_accession.fasta"

out_dir = Path("rebuilt")
out_dir.mkdir(exist_ok=True)

output_5 = out_dir / "dada2_with_accession.fasta"
output_6 = out_dir / "dada2_assignTaxonomy.fasta"
output_7 = out_dir / "dada2_assignTaxonomy_genus.fasta"

RANK_COLS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
GENUS_COLS = RANK_COLS[:6]

# ═══════════════════════════════════════════════════════════════
# 1. SPECIES NAME CLEANING FUNCTION
# ═══════════════════════════════════════════════════════════════

def clean_species(species_str):
    """
    Apply exception rules to species string.
    Returns cleaned species string or None.
    
    Rules:
    1. First letter must be capitalized
    2. No uncultured/metagenome/environmental terms
    3. Strip parenthetical common names
    4. sp./cf./aff. - keep only if something follows the period
    5. Strip _X suffixes
    """
    if not species_str or not species_str.strip():
        return None
    
    cleaned = species_str.strip()
    
    # RULE 5: Strip SILVA _X suffixes
    cleaned = re.sub(r'_X{1,3}$', '', cleaned)
    cleaned = re.sub(r'_X{1,3}\b', '', cleaned)  # Also mid-string
    
    # RULE 3: Strip parenthetical common names
    cleaned = re.sub(r'\s*\([^)]*\)', '', cleaned)
    cleaned = cleaned.strip()
    
    if not cleaned:
        return None
    
    # RULE 1: First letter must be capitalized
    if not cleaned[0].isupper():
        return None
    
    # RULE 2: No problematic terms
    problematic = r'\b(uncultured|unidentified|unknown|metagenome|metagenomic|environmental)\b'
    if re.search(problematic, cleaned, re.IGNORECASE):
        return None
    
    # RULE 4: Handle sp./cf./aff. qualifiers
    # Check if there's content after the qualifier
    qualifier_pattern = r'\b(sp\.|cf\.|aff\.)\s*(.*)$'
    match = re.search(qualifier_pattern, cleaned)
    if match:
        after_qualifier = match.group(2).strip()
        if not after_qualifier:
            # Nothing after qualifier - return None
            return None
        # Otherwise keep the full string
        return cleaned
    
    return cleaned

# ═══════════════════════════════════════════════════════════════
# 2. LOAD RESOLVED TAXONOMY
# ═══════════════════════════════════════════════════════════════

print("Loading resolved taxonomy...")

# Build taxa_id → higher ranks mapping
taxa_to_ranks = {}  # taxa_id → {Kingdom: ..., Phylum: ..., ..., Genus: ...}
taxa_to_accessions = {}  # taxa_id → [acc1, acc2, ...]

with open(resolved_csv, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        tid = row['taxa_id']
        
        # Get Kingdom through Genus (not Species yet)
        ranks = {}
        for rank in GENUS_COLS:
            val = row.get(rank, '').strip()
            ranks[rank] = val if val and val != 'NA' else None
        
        taxa_to_ranks[tid] = ranks
        
        # Parse accession list
        accs = [a.strip() for a in row['accessions'].split('|') if a.strip()]
        taxa_to_accessions[tid] = accs

print(f"  Loaded {len(taxa_to_ranks)} taxa IDs")

# ═══════════════════════════════════════════════════════════════
# 3. LOAD RAW SPECIES NAMES FROM LOG
# ═══════════════════════════════════════════════════════════════

print("Loading raw species names from resolution log...")

# Build taxa_id → raw species name
# Use the FIRST "resolved" entry for each taxa_id (when resolution succeeded)
taxa_to_raw_species = {}

with open(log_csv, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        tid = row['taxa_id']
        outcome = row.get('final_outcome', '')
        
        # Only look at successful resolutions
        if outcome == 'resolved' and tid not in taxa_to_raw_species:
            raw_name = row.get('raw_name', '').strip()
            if raw_name:
                taxa_to_raw_species[tid] = raw_name

print(f"  Found raw species names for {len(taxa_to_raw_species)} taxa IDs")

# ═══════════════════════════════════════════════════════════════
# 4. BUILD ACCESSION → FULL TAXONOMY MAPPING
# ═══════════════════════════════════════════════════════════════

print("Building accession → taxonomy mapping...")

acc_to_taxonomy = {}

for tid, accessions in taxa_to_accessions.items():
    # Get higher ranks
    higher_ranks = taxa_to_ranks.get(tid, {})
    
    # Get and clean species name
    raw_species = taxa_to_raw_species.get(tid, None)
    cleaned_species = clean_species(raw_species) if raw_species else None
    
    # Build full 7-rank taxonomy
    full_taxonomy = {}
    for rank in GENUS_COLS:
        full_taxonomy[rank] = higher_ranks.get(rank, None)
    full_taxonomy['Species'] = cleaned_species
    
    # Map all accessions to this taxonomy
    for acc in accessions:
        acc_to_taxonomy[acc] = full_taxonomy

print(f"  Mapped {len(acc_to_taxonomy)} accessions")

# ═══════════════════════════════════════════════════════════════
# 5. REBUILD FASTA FILES
# ═══════════════════════════════════════════════════════════════

print("\nRebuilding FASTA files...")

written_5 = 0
written_6 = 0
written_7 = 0
skipped = 0

with open(output_5, 'w') as out5, \
     open(output_6, 'w') as out6, \
     open(output_7, 'w') as out7:
    
    for record in SeqIO.parse(input_fasta, 'fasta'):
        acc = record.id
        
        if acc not in acc_to_taxonomy:
            skipped += 1
            continue
        
        taxonomy = acc_to_taxonomy[acc]
        
        # Clean sequence
        seq = str(record.seq).replace('-', '').replace('.', '').upper()
        
        # Build taxonomy strings (replace None with "NA")
        def rank_val(r):
            val = taxonomy.get(r)
            return val if val else "NA"
        
        full_tax = ";".join(rank_val(r) for r in RANK_COLS)
        genus_tax = ";".join(rank_val(r) for r in GENUS_COLS)
        
        # Output 5: >ACCESSION Kingdom;Phylum;...;Species
        out5.write(f">{acc} {full_tax}\n{seq}\n")
        written_5 += 1
        
        # Output 6: >Kingdom;Phylum;...;Species;
        out6.write(f">{full_tax};\n{seq}\n")
        written_6 += 1
        
        # Output 7: >Kingdom;Phylum;...;Genus;
        out7.write(f">{genus_tax};\n{seq}\n")
        written_7 += 1

# ═══════════════════════════════════════════════════════════════
# 6. REPORT
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*60)
print("  Rebuild Complete")
print("="*60)
print(f"Sequences written:              {written_5}")
print(f"Skipped (not in taxonomy):      {skipped}")
print()
print("Output files:")
print(f"  5. With accession:    {output_5}")
print(f"  6. assignTaxonomy:    {output_6}")
print(f"  7. Genus-level:       {output_7}")
print("="*60)

# ═══════════════════════════════════════════════════════════════
# 7. SPECIES CLEANING SUMMARY
# ═══════════════════════════════════════════════════════════════

print("\nSpecies name cleaning summary:")

total_taxa = len(taxa_to_raw_species)
cleaned_species = sum(1 for tax in acc_to_taxonomy.values() if tax['Species'])
na_species = total_taxa - cleaned_species

print(f"  Total taxa with raw species:  {total_taxa}")
print(f"  Species kept after cleaning:  {cleaned_species}")
print(f"  Species set to NA:            {na_species}")

# Show some examples
print("\nExample transformations:")
example_count = 0
for tid, raw in list(taxa_to_raw_species.items())[:10]:
    cleaned = clean_species(raw)
    if raw != cleaned:
        status = cleaned if cleaned else "NA"
        print(f"  '{raw}' → '{status}'")
        example_count += 1
        if example_count >= 5:
            break

print("\n✓ Done! You can now use the rebuilt files for steps 5+ of the pipeline.")
print("  Or manually run step 5 clustering on rebuilt/dada2_with_accession.fasta")