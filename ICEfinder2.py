#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
from pathlib import Path
from typing import Tuple

from script.checkin import get_fagb
from script.config import get_param
from script.metaICE import _meta
from script.single import _single


logging.basicConfig(
        level=logging.INFO,
        datefmt='%m/%d/%Y %I:%M:%S',
        handlers=[logging.StreamHandler()]
    )

def add_arguments_to_parser(parser: argparse.ArgumentParser) -> None:
    """Adds command-line arguments to the parser."""
    parser.add_argument(
        "-v", "--version", action="version", version="2.0", help="Show ICEfinder version"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Input file in FASTA/Genbank format (Genbank only for single genome)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        required=True,
        help="Output directory for results",
    )
    parser.add_argument(
        "-t",
        "--type",
        type=str,
        required=True,
        choices=["Single", "Metagenome"],
        help="Analysis type: Single genome or Metagenome",
    )


def setup_directories(outdir: Path) -> Tuple[Path, Path, Path]:
    """Creates and returns the required directory structure."""
    tmp_dir = outdir / "tmp"
    fa_dir = tmp_dir / "fasta"
    gb_dir = tmp_dir / "gbk"

    tmp_dir.mkdir(parents=True, exist_ok=True)
    fa_dir.mkdir(exist_ok=True)
    gb_dir.mkdir(exist_ok=True)

    return tmp_dir, fa_dir, gb_dir


def main() -> None:
    logging.info("\n###################\nStarting ICEfinder2\n###################")

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="ICEfinder - Identify ICE elements in genomic data",
        usage="python ICEfinder.py -i input_file -o output_dir -t {Single,Metagenome}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments_to_parser(parser)
    args = parser.parse_args()

    # Validate and create output directories
    args.outdir.mkdir(parents=True, exist_ok=True)
    tmp_dir, fa_dir, gb_dir = setup_directories(args.outdir)

    # Get parameters from configuration
    logging.debug("Loaded configuration parameters")

    # Set up analysis parameters
    run_id = args.input.stem
    logging.info(
        f"\nInput data -->\n"
        f"\tFile: {args.input}\n"
        f"\tType: {args.type}\n"
        f"\tRun ID: {run_id}\n"
        f"\tOutput dir: {args.outdir}\n"
    )

    # Process input file
    infile, filetype = get_fagb(run_id, args.input, args.type, tmp_dir, fa_dir, gb_dir)

    # Execute analysis pipeline
    if args.type == "Single":
        logging.info("Executing Single genome mode")
        _single(run_id, infile, filetype)
    else:
        logging.info("Executing Metagenome mode")
        _meta(run_id, infile)

    logging.info(f"Analysis completed for {run_id}!")


if __name__ == "__main__":
    main()
