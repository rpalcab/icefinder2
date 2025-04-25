#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
from Bio import SeqIO
import shutil

class SequenceProcessingError(Exception):
    """Custom exception for sequence processing errors."""
    pass

def check_fasta(input: Path) -> None:
    try:
        _ = next(SeqIO.parse(input, 'fasta'))
    except StopIteration:
        raise SequenceProcessingError("The input file is empty or not in FASTA format")
    except Exception as e:
        raise SequenceProcessingError(f"Error parsing FASTA: {e}")
    
def check_gff(input: Path) -> None:
    if not input.exists():
        raise FileNotFoundError(f"GFF file '{input}' does not exist.")
    if input.suffix not in [".gff", ".gff3"]:
        raise SequenceProcessingError(f"Invalid GFF file extension: '{input.suffix}'")
    if input.stat().st_size == 0:
        raise SequenceProcessingError("The GFF file is empty.")

def check_tools() -> None:
    for tool in ['aragorn', 'hmmscan2', 'mkvtree', 'vmatch']:
        tool_file = Path(__file__).parents[1] / 'tool' / tool
        if not tool_file.exists():
            raise FileNotFoundError(f"Necessary tool '{tool}' not found. Please, check {tool_file}")

def check_dbs() -> None:
    for db in ['ICEscan.hmm.*']:
        db_files = list((Path(__file__).parents[1] / 'data').glob(db))
        if not db_files:
            raise FileNotFoundError(f"Necessary database matching '{db}' not found. Please, check the data directory.")

def copy_files(src: Path, dst: Path) -> None:
    """
    Recursively copy a file or directory from `src` to `dst`.
    - If `src` is a file, it’s copied (with metadata) to `dst`.
    - If `src` is a directory, its entire tree is copied into `dst`.
    """
    if src.is_file():
        # Ensure parent dirs exist, then copy file
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    elif src.is_dir():
        # Copy entire directory tree; create dst if needed
        dst.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, dst, dirs_exist_ok=True)
    else:
        raise FileNotFoundError(f"Source not found: {src}")
