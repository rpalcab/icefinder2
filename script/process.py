#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,time
import random,json
import string,shutil
import logging
from pathlib import Path
from typing import List, Union, Optional, Tuple, Dict
import subprocess
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline

from Bio import SeqIO
from ete3 import NCBITaxa
from Bio.SeqUtils import GC
from functools import cmp_to_key
from script.function import getblast
from script.config import get_param
from script.utils import copy_files

logging.basicConfig(
        level=logging.INFO,
        datefmt='%m/%d/%Y %I:%M:%S %p',
        handlers=[logging.StreamHandler()]
    )


# Output processing
def get_time() -> str:
	return time.asctime( time.localtime(time.time()) )

def process_hmmscan(scan_file: Path):
	df = pd.read_table(scan_file, sep=r"\s+", comment='#', header=None, usecols=range(5), 
					   names=['target_name', 'accession', 'query_name', 'accession_fs', 'E-value_fs'])
	df.drop_duplicates(subset='query_name', inplace=True)
	df['key'] = df['query_name'].str.split('_').str[:-1].str.join('_')

	chosen = []
	for k, group in df.groupby('key'):
		if scanf(group['target_name'].tolist()):
			chosen.append(k)
	return chosen

def scanf(hmmlist: List) -> bool:

	ice_count = []
	for line in hmmlist:
		if 'MOB' in line:
			ice_count.append('MOB')
		elif 't4cp' in line or 'tcpA' in line:
			ice_count.append('t4cp')
		elif 'FA' in line:
			ice_count.append('T4SS')
		elif line in [
					 'Phage_integrase', 'UPF0236',
					 'Recombinase', 'rve', 'TIGR02224',
					 'TIGR02249', 'TIGR02225', 'PB001819'
					 ]:
			ice_count.append('Int')
		else:
			ice_count.append('T4SS')
	if ice_count.count('MOB') and ice_count.count('t4cp') and ice_count.count('Int') and ice_count.count('T4SS') >= 5:
		return True
	else:
		return False

def process_kraken2(output: Path) -> Dict:
	spdict = {}
	with output.open() as fh:
		for line in fh:
			_, seq_id, taxid_str, *rest = line.rstrip("\n").split("\t")
			if taxid_str == "0":
				spdict[seq_id] = "-"
			else:
				spdict[seq_id], _ = get_ranks(int(taxid_str), NCBITaxa())
	return spdict

def get_ranks(taxid: int, _ncbi_db) -> Tuple[str, str]:
    """
    Given a taxonomic ID, return a tuple (species_name, strain_name).
    If no species or strain can be resolved, returns '-' or '' respectively.
    """
    # 1. Retrieve full lineage and ranks
    lineage = _ncbi_db.get_lineage(taxid)
    ranks = _ncbi_db.get_rank(lineage)               # {taxid: rank, ...}
    names = _ncbi_db.get_taxid_translator(lineage)   # {taxid: name, ...}

    # 2. Invert to map rank → taxid
    rank2tid = {r: tid for tid, r in ranks.items()}

    # 3. Pick the “best” taxid for species-level naming:
    #    species → genus → phylum
    for primary_rank in ("species", "genus", "phylum"):
        sp_tid = rank2tid.get(primary_rank)
        if sp_tid is not None:
            break
    else:
        sp_tid = None

    # 4. Optionally get the strain taxid (only valid if species was found)
    st_tid = rank2tid.get("strain") if rank2tid.get("species") else None

    # 5. Lookup names, with fallbacks
    sp_name = names.get(sp_tid, "-") if sp_tid else "-"
    st_name = names.get(st_tid, "")  if st_tid else ""

    return sp_name, st_name

def process_seqkit(result: Path, run_id: str) -> Dict:
	for raw in result.stdout.splitlines():
		line = raw.strip()
		if not line.startswith("#") and not line.lower().startswith("file"):
			cols = line.split()
			length = cols[4]
			contig_count = cols[3]
			n50 = cols[12]

	base_dict = {'JobID': run_id,
		    	 'Submission date': get_time(),
		    	 'Total length': f'{length} bp',
   		    	 'Contig number': contig_count,
		    	 'Sequence N50': f'{n50} bp'
		    }
	
	return base_dict

def parse_icescan(ice_res: Path) -> Tuple[Dict, Dict]:
	df_systems = pd.read_table(ice_res, comment='#')
	df_systems = df_systems[df_systems.model_fqn.str.contains('Chromosome')]
	df_systems = df_systems[~df_systems.model_fqn.str.contains('UserReplicon_IME')]

	ice_dict = {}
	info_dict = {}
	for _, row in df_systems.iterrows():
		gbk_name = row['hit_id']
		tags = get_feat(row.gene_name)
		mpf_type = row.model_fqn.split('/')[-1].split('_')[1] if 'IME' not in row.model_fqn else ''
		mob_type = tags.split('@')[1] if 'Relaxase@' in tags else ''
		ice_tag = f"ICE{row['sys_id'].split('_')[-1]}"

		ice_dict.setdefault(ice_tag, {})[gbk_name] = tags
		info_dict.setdefault(ice_tag, {'mob': [], 'mpf': []})
		if mob_type and mob_type not in info_dict[ice_tag]['mob']:
			info_dict[ice_tag]['mob'].append(mob_type)
		if mpf_type not in info_dict[ice_tag]['mpf']:
			info_dict[ice_tag]['mpf'].append(mpf_type)

	return ice_dict, info_dict

def get_feat(feat):

	integrases_list = ['Phage_integrase', 'UPF0236', 'Recombinase', 'rve', 
					   'TIGR02224', 'TIGR02249', 'TIGR02225', 'PB001819']

	if feat in integrases_list:
		return f'Integrase@{feat}'

	elif 'T4SS_MOB' in feat:
		tag = feat.split('_')[1]
		return f'Relaxase@{tag}'

	elif 't4cp' in feat or 'tcpA' in feat:
		tag = feat.split('_')[1]
		return f'T4CP@{tag}'

	else:
		return 'T4SS@'+feat.replace('T4SS_','')

def process_vmatch(result) -> pd.DataFrame:
    rows = []
    for raw in result.stdout.splitlines():
        line = raw.strip()
        if not line.startswith("#"):
            cols = line.split()
            start1 = int(cols[2])
            end1 = start1 + int(cols[0])
            start2 = int(cols[6])
            end2 = start2 + int(cols[4])
            rows.append([start1, end1, start2, end2])
    
    df = pd.DataFrame(rows, columns=["start1", "end1", "start2", "end2"])
    return df

def parse_defensefinder(infile: Path) -> Dict:
    results: Dict[str, str] = {}
    with infile.open() as fh:
        for line in fh:
            if line.startswith('replicon'):
                continue
            fields = line.rstrip("\n").split('\t')
            locus, gene = fields[1], fields[2]
            results[locus] = gene.replace('__', ',')
    return results

# Commands
def kraken2(run_id: str, infile: Path, outdir: Path, db: Path) -> Tuple[Dict, Path]:
	logging.info("Running Kraken2")

	report = outdir / 'tmp' / f'{run_id}_kraken.report'
	output = outdir / 'tmp' / f'{run_id}_kraken.output'

	cmd = [
		'kraken2',
		'--db', db,
		'--report', report,
		'--output', output,
		infile
		]
	
	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
		spdict = process_kraken2(output)
	except subprocess.CalledProcessError as e:
		logging.error(f"Kraken2 classification failed: {e}")

	return spdict,report

def prodigal(run_id: str, infile: Path, outdir: Path) -> Path:
	logging.info("Running Prodigal")

	anno_fa = outdir / 'tmp' / f'{run_id}.faa'
	anno_gff = outdir / 'tmp' / f'{run_id}.gff'
	cmd = [
		'prodigal',
		'-c',
		'-m',
		'-q',
		'-p', 'meta',
		'-f', 'gff',
		'-i', infile,
		'-a', anno_fa,
		'-o', anno_gff
		]
	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	except subprocess.CalledProcessError as e:
		logging.error(f"Prodigal annotation failed: {e}")

	return anno_fa

def seqkit(run_id: str, infile: Path, outdir: Path) -> Path:
	logging.info("Running Seqkit stats")

	cmd = [
		'seqkit',
		'stats',
		'-a', infile
		]
	try:
		result = subprocess.run(cmd, check=True, capture_output=True, text=True)
		base_file = outdir / 'tmp' / f'{run_id}_info.json'
		base_dict = process_seqkit(result, run_id)
		with base_file.open("w") as fh:
			json.dump(base_dict, fh, indent=4)
	except subprocess.CalledProcessError as e:
		logging.error(f"Seqkit stats failed: {e}")

	return base_file

def hmmscan(run_id: str, infile: Path, outdir: Path, anno_fa: Path) -> Path:
	logging.info("Running Hmmscan2")

	scan_file = outdir / 'tmp' / f'{run_id}_prescan.tsv'
	icescan_hmm = Path(__file__).parents[1] / 'data' / 'ICEscan.hmm'
	hmmscan_tool = Path(__file__).parents[1] / 'tool' / 'hmmscan2'
	cmd = [
		hmmscan_tool,
		'-E', str(0.00001),
		'--tblout', scan_file,
		icescan_hmm,
		anno_fa
		]
	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	except subprocess.CalledProcessError as e:
		logging.error(f"hmmscan search failed: {e}")

	return scan_file

def prokka(run_id: str, infile: Path, outdir: Path) -> Path:
	logging.info(f"Running Prokka on contig {run_id}")
	outdir = outdir / 'prokka'
	cmd = [
		'prokka',
		infile,
		'--force',
		'--fast',
		'--quiet', 
		'--cdsrnaolap',
		'--cpus', str(8),
		'--outdir', outdir,
		'--prefix', run_id
		]

	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	except subprocess.CalledProcessError as e:
		logging.error(f"Prokka search failed: {e}")

	return outdir

def ICEscan(contig_id: str, prokka_dir: Path, outdir: Path) -> Tuple[Dict, Dict]:
	logging.info(f"Running Macsyfinder with ICEscan model on contig {contig_id}")
	anno_fa = prokka_dir / f'{contig_id}.faa'
	outdir = outdir / 'mcsy_icescan'
	cmd = [
		'macsyfinder',
		'--db-type', 'ordered_replicon',
		'--models-dir', Path(__file__).parents[1] / 'data' / 'macsydata',
		'--models', 'ICEscan', 'all',
		'--replicon-topology', 'linear',
		'--coverage-profile', '0.3',
		'--sequence-db', anno_fa,
		'--force',
		'-o', outdir
	]

	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
		ice_dict, info_dict = parse_icescan(outdir / 'all_systems.tsv')
	except subprocess.CalledProcessError as e:
		logging.error(f"Macsyfinder execution failed: {e}")

	return ice_dict, info_dict

def get_dr(run_id: str, infile: Path, output: Path):
	logging.info(f"Running Mkctree and Vmatch on contig {run_id}")
	dr_index = output / f'{run_id}_DR'

	mkctree_cmd = [
		Path(__file__).parents[1] / 'tool' / 'mkvtree',
		'-db', infile,
		'-indexname', dr_index,
		'-dna',
		'-pl',
		'-lcp',
		'-suf',
		'-tis',
		'-ois',
		'-bwt',
		'-bck',
		'-sti1'
	]

	vmatch_cmd = [
		Path(__file__).parents[1] / 'tool' / 'vmatch',
		'-l', str(15),
		dr_index
	]

	try:
		subprocess.run(mkctree_cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
		result = subprocess.run(vmatch_cmd, check=True, capture_output=True, text=True)
		dr_list = process_vmatch(result)
	except subprocess.CalledProcessError as e:
		logging.error(f"DR search failed: {e}")

	return dr_list

def blastp(faa_file, db, outfile):
	logging.info(f"Running Blastp for {db.stem} detection")
	blastp_cline = NcbiblastpCommandline(cmd='blastp', query=faa_file, db=db, \
                       evalue=0.0001, num_threads=20, max_hsps=1, max_target_seqs=1, \
                       outfmt="6 std slen stitle", out=outfile)
	blastp_cline()

def blastn(fna_file, db, outfile):
	logging.info(f"Running Blastn for {db.stem} detection")
	blastp_cline = NcbiblastnCommandline(cmd='blastn', query=fna_file, db=db, \
                       evalue=0.0001, num_threads=20, max_hsps=1, max_target_seqs=1, \
                       outfmt="6 std slen stitle", out=outfile)
	blastp_cline()

def havalue(value,out):

	blast_filter = {}
	for line in open(out,'r').readlines():
		lines = line.strip().split('\t')
		havalue = (int(lines[3])/int(lines[12]))*float(lines[2])/100
		if havalue >= float(value):
			blast_filter[lines[0]]=lines[1].split('|')[1]
	return blast_filter

def blastp_all(infile, outdir):
	db_dict = {'transposase': 0.64, 'virulence': 0.64, 'degradation': 0.64, 'metal': 0.64, 'symbiosis': 0.64}
	results_dict = {}
	for db, ha in db_dict.items():
		outfile = outdir / f'{db}.out'
		blastp(infile, Path(__file__).parents[1] / 'data' / db, outfile)
		results_dict[f'{db}_dict'] = havalue(ha, outfile)

	return results_dict

def blastn_all(infile, outdir):
	outfile = outdir / 'resfinder.out'
	blastn(infile, Path(__file__).parents[1] / 'data' / 'resfinder', outfile)
	results_dict = {'resfinder': havalue('0.81',outfile)}

	return results_dict

def defense_finder(run_id: str, infile: Path, outdir: Path) -> Dict:
	logging.info("Running Defense_finder")
	outdir = outdir / 'defense_finder'
	cmd = [
		'defense-finder', 'run',
		'-w', str(8),
		'--models-dir', Path(__file__).parents[1] / 'data' / 'macsydata',
		'-o', outdir,
		infile
		]
	
	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
		df_dict = parse_defensefinder(outdir / f'{run_id}_defense_finder_genes.tsv')
	except subprocess.CalledProcessError as e:
		logging.error(f"Macsyfinder execution failed: {e}")

	return df_dict

