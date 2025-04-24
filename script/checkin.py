#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import sys
from pathlib import Path
from typing import Tuple
from Bio import SeqIO

class SequenceProcessingError(Exception):
    """Custom exception for sequence processing errors."""
    pass

def is_fagb(filepath: Path) -> str:

    for fmt, label in (('fasta', 'fa'), ('gb', 'gb')):
        try:
            records = list(SeqIO.parse(str(filepath), fmt))
            if records:
                return label
        except Exception:
            continue
    return ''

def remove_folders_with_run_id(root_dir: Path, run_id: str) -> None:
    for path in root_dir.rglob(f'*{run_id}*'):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

def process_fasta(input_file: Path, run_id: str, fa_dir: Path, intype: str) -> Tuple[Path, str]:
    out = fa_dir / f'{run_id}.fa'
    shutil.copy(str(input_file), str(out))
    records = list(SeqIO.parse(str(input_file), 'fasta'))
    if len(records) == 1:
        return out, 'fa'
    if intype.lower() == 'metagenome':
        return out, 'multifa'
    raise SequenceProcessingError("FASTA input must contain exactly one record unless Metagenome")

def process_gbk(infile: Path, run_id: str, fa_dir: Path, gb_dir: Path) -> Tuple[Path, str]:
    """
    Process a GenBank input: ensure single record, valid features, and convert to FASTA.
    """
    records = list(SeqIO.parse(str(infile), 'gb'))
    if len(records) != 1:
        raise SequenceProcessingError("GenBank input must contain exactly one record")
    record = records[0]
    # Check for missing locus_tag qualifiers in CDS features
    missing_lt = sum(
        1 for f in record.features if f.type == 'CDS' and 'locus_tag' not in f.qualifiers
    )
    if missing_lt > 10:
        raise SequenceProcessingError(
            "Too many CDS without locus_tag in GenBank input"
        )
    # Check for degenerate sequences (e.g., a sequence of identical bases)
    if len(set(str(record.seq))) == 1:
        raise SequenceProcessingError(
            "Invalid GenBank sequence: appears degenerate"
        )
    # Assign new ID and write outputs
    record.id = run_id
    gb_out = gb_dir / f'{run_id}.gbk'
    fa_out = fa_dir / f'{run_id}.fa'
    SeqIO.write(record, str(gb_out), 'gb')
    SeqIO.write(record, str(fa_out), 'fasta')
    return fa_out, 'gb'

def get_fagb(run_id: str, input_file: Path, intype: str, tmp_dir: Path, fa_dir: Path, gb_dir: Path) -> None:

	remove_folders_with_run_id(gb_dir,run_id)
	remove_folders_with_run_id(tmp_dir,run_id)

	filetype = is_fagb(input_file)
	if not filetype:
		raise SequenceProcessingError(
               "The input file is not a standard FASTA/GenBank format"
        )

	if filetype == 'fa':
		return process_fasta(input_file, run_id, fa_dir, intype)
	return process_gbk(input_file, run_id, fa_dir, gb_dir)
