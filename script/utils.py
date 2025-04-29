#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
from Bio import SeqIO
import shutil
import re
import gffutils
from typing import Dict, List, Tuple, Optional

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
    for db in ['ICEscan.hmm.*', 'degradation.*', 'metal.*', 'oriT_db.*', 'resfinder.*', 'symbiosis.*', 'transposase.*', 'virulence.*']:
        db_files = list((Path(__file__).parents[1] / 'data').glob(db))
        if not db_files:
            raise FileNotFoundError(f"Necessary database matching '{db}' not found. Please, check the data directory -> {Path(__file__).parents[1] / 'data'}")

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

def zill(header: str, num: int, width: int = 5) -> str:
    return f"{header}_{str(num).zfill(width)}"

def get_gff(gff_path: str) -> Tuple[Dict[str, str], Dict[str, List[str]], str, int]:
    """
    Parse a GFF3 file with gffutils and return:
      - trna_dict: IDs → product for tRNA/tmRNA only
      - pos_dict:  IDs → [start, end, strand, product]
      - header:    last ID’s prefix before '_' 
      - totalnum:  last ID’s numeric suffix as int

    Builds (or loads) a SQLite database via gffutils.create_db().
    """
    # Load GFF file
    db = gffutils.create_db(
        str(gff_path),
        dbfn=':memory:',
        force=True,            # overwrite existing DB
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True
    )

    trna_dict = {}
    pos_dict = {}
    header = ''
    totalnum = 0
    # Iterate all features in genomic order
    for feature in db.all_features(order_by='start'):
        id_ = feature.id
        # Attributes are lists; take first product if present
        product = feature.attributes.get('product', [''])[0]
        pos_dict[id_] = [
            str(feature.start),
            str(feature.end),
            feature.strand,
            product
        ]
        # Only keep tRNA / tmRNA entries
        if feature.featuretype in {'tRNA', 'tmRNA'}:
            trna_dict[id_] = product

    header = id_.split('_')[0]
    totalnum = len(list(db.all_features()))

    return trna_dict, pos_dict, header, totalnum

def get_num(id_str: str) -> int:
    try:
        return int(id_str.split('_', 1)[1].lstrip('0') or '0')
    except (IndexError, ValueError):
        raise ValueError(f"Invalid ID format: '{id_str}'")

def find_max_distance(numbers: List[int]) -> Optional[Tuple[int, int]]:
    """
    Return the consecutive pair from `numbers` with the maximal absolute difference.
    If fewer than two numbers, return None.
    """
    if len(numbers) < 2:
        return None
    # Evaluate pairs and pick the one with max gap
    return max(zip(numbers, numbers[1:]), key=lambda pair: abs(pair[1] - pair[0]))

def candidate_setup(seq_record, outdir: Path) -> Tuple[Path, Path]:
	contig_folder = outdir / 'tmp' / f'{seq_record.id}'
	contig_folder.mkdir(exist_ok=True)
	contig_fasta = contig_folder / f'{seq_record.id}.fasta'
	with contig_fasta.open("w") as oh:
		SeqIO.write(seq_record, oh, "fasta")
	return contig_folder, contig_fasta

