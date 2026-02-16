"""Shared utilities for FASTA parsing and SILVA header handling."""

import re
from typing import Dict, Tuple, List, Optional
from Bio import SeqIO


def parse_silva_header(header: str) -> Dict[str, str]:
    """
    Parse a SILVA-formatted FASTA header.
    
    SILVA headers typically look like:
    >ACCESSION.start.stop Taxonomy;Path;Down;To;Species
    
    Returns dict with 'accession' and 'taxonomy' (the full semicolon-delimited string).
    """
    parts = header.split(maxsplit=1)
    accession = parts[0]
    taxonomy = parts[1].strip() if len(parts) > 1 else ""
    return {"accession": accession, "taxonomy": taxonomy}


def get_taxonomy_ranks(taxonomy_str: str) -> List[str]:
    """Split SILVA semicolon-delimited taxonomy into a list of rank names."""
    return [t.strip() for t in taxonomy_str.rstrip(";").split(";") if t.strip()]


def get_lowest_level_name(taxonomy_str: str) -> str:
    """Return the lowest (most specific) taxonomic name from a SILVA taxonomy string."""
    ranks = get_taxonomy_ranks(taxonomy_str)
    return ranks[-1] if ranks else ""


def taxonomy_contains_any(taxonomy_str: str, target_names: set) -> bool:
    """Check if any rank in the taxonomy matches any name in target_names."""
    ranks = get_taxonomy_ranks(taxonomy_str)
    return bool(set(ranks) & target_names)


def load_target_classes(filepath: str) -> set:
    """Load target class names from a text file (one per line)."""
    with open(filepath) as f:
        return {line.strip() for line in f if line.strip() and not line.startswith("#")}


def iter_fasta_nonempty(filepath: str):
    """Iterate over FASTA records, skipping entries with empty sequences."""
    for record in SeqIO.parse(filepath, "fasta"):
        seq_str = str(record.seq).replace("-", "").replace(".", "").strip()
        if len(seq_str) > 0:
            yield record