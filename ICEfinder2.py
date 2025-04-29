#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
from pathlib import Path
from typing import Tuple

from script.utils import check_fasta, check_gff, check_tools, check_dbs
from script.metaICE import _meta


logging.basicConfig(
        level=logging.INFO,
        datefmt='%m/%d/%Y %I:%M:%S %p',
        handlers=[logging.StreamHandler()]
    )

class SequenceProcessingError(Exception):
    """Custom exception for sequence processing errors."""
    pass

def get_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        prog="ICEfinder2.py",
        description="ICEfinder - Identify ICE elements in genomic data."
    )

    parser.add_argument(
        "-v", "--version", action="version", version="2.0", help="Show ICEfinder version"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        type=Path,
        help="Input file in FASTA/Genbank format (Genbank only for single genome)",
    )
    parser.add_argument(
        "-o", "--outdir", required=True,
        type=Path,
        help="Output directory for results",
    )
    parser.add_argument(
        "-g", "--gff",
        type=Path,
        help="Gene annotation file in gff format",
    )
    parser.add_argument(
        "-k", "--kraken_db",
        type=Path, default=None,
        help="Path to Kraken2 database. If not provided, taxonomical identification is skipped",
    )

    return parser.parse_args()


def setup_directories(outdir: Path) -> Tuple[Path, Path, Path]:
    """Creates and returns the required directory structure."""
    outdir.mkdir(parents=True, exist_ok=True)
    tmp_dir = outdir / "tmp"
    fa_dir = tmp_dir / "fasta"
    gb_dir = tmp_dir / "gbk"

    tmp_dir.mkdir(exist_ok=True)
    fa_dir.mkdir(exist_ok=True)
    gb_dir.mkdir(exist_ok=True)

    return tmp_dir, fa_dir, gb_dir


def main() -> None:
    logging.info("\n###################\nStarting ICEfinder2\n###################")
    args = get_args()

    # Initial parameters
    run_id = args.input.stem
    logging.info(
        f"\nInput data ->\n"
        f"\tFile: {args.input}\n"
        f"\tRun ID: {run_id}\n"
        f"\tOutput dir: {args.outdir}\n"
    )

    # Checks
    check_fasta(args.input)     # Fasta integrity

    ##TODO: Check GFF content is valid (look for python packages like seqio.bio but for gff)
    if args.gff:                # GFF integrity (if provided)
        check_gff(args.gff)

    if args.kraken_db:          # Kraken DB exists (if provided)
        if not args.kraken_db.exists() or not args.kraken_db.is_dir():
            raise FileNotFoundError(f"Kraken2 database directory '{args.kraken_db}' does not exist or is not a directory.")
    
    check_dbs()

    check_tools()
    
    # Output setup
    tmp_dir,fa_dir,gb_dir = setup_directories(args.outdir)
    
    # 
    logging.info("Executing meta mode")
    _meta(run_id, args.input, args.outdir, args.kraken_db, args.gff)

    logging.info(f"Analysis completed for {run_id}!")


if __name__ == "__main__":
    main()
