#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,time
import random,json
import string,shutil
import logging
from pathlib import Path
from typing import List, Union, Optional, Tuple
import subprocess
import pandas as pd

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

param = get_param()
workdir = param[0]
soft_dir = param[1]
blastn = param[5]
macsyfinder = param[9]
hmmsearch = param[10]

tmp_dir = os.path.join(workdir,'tmp') 
in_dir = os.path.join(tmp_dir,'fasta')
gb_dir = os.path.join(tmp_dir,'gbk') # AHORA ES EL OUTPUT DE PROKKA


# Refactored
def get_time() -> str:

	return time.asctime( time.localtime(time.time()) )

def taxonomy(run_id: str, input: Path, outdir: Path, db: Path) -> Tuple[Path, dict, Path]:
	logging.info("Running Kraken2")

	report = outdir / 'tmp' / f'{run_id}_kraken.report'
	output = outdir / 'tmp' / f'{run_id}_kraken.output'
	drawout = outdir / 'tmp' / 'kraken.html'

	cmd = [
		'kraken2',
		'--db', db,
		'--report', report,
		'--output', output,
		input
		]
	
	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	except subprocess.CalledProcessError as e:
		logging.error(f"Kraken2 classification failed: {e}")

	spdict = {}
	with output.open() as fh:
		for line in fh:
			_, seq_id, taxid_str, *rest = line.rstrip("\n").split("\t")
			if taxid_str == "0":
				spdict[seq_id] = "-"
			else:
				spdict[seq_id], _ = get_ranks(int(taxid_str), NCBITaxa())

	return drawout,spdict,report

def prodigal(run_id: str, input: Path, outdir: Path) -> Path:
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
		'-i', input,
		'-a', anno_fa,
		'-o', anno_gff
		]
	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	except subprocess.CalledProcessError as e:
		logging.error(f"Prodigal annotation failed: {e}")

	return anno_fa

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

def getbase(run_id: str, input: Path, outdir: Path) -> Path:
	logging.info("Getting basic sample stats")

	cmd = [
		'seqkit',
		'stats',
		'-a', input
		]
	try:
		result = subprocess.run(cmd, check=True, capture_output=True, text=True)
	except subprocess.CalledProcessError as e:
		logging.error(f"Seqkit stats failed: {e}")

	for raw in result.stdout.splitlines():
		line = raw.strip()
		if not line.startswith("#") and not line.lower().startswith("file"):
			cols = line.split()
			length = cols[4]
			contig_count = cols[3]
			n50 = cols[12]
			

	base_file = outdir / 'tmp' / f'{run_id}_info.json'
	base_dict = {'JobID': run_id,
		    	 'Submission date': get_time(),
		    	 'Total length': f'{length} bp',
   		    	 'Contig number': contig_count,
		    	 'Sequence N50': f'{n50} bp'
		    }
	with base_file.open("w") as fh:
		json.dump(base_dict, fh, indent=4)

	return base_file

def hmmscan(run_id: str, input: Path, outdir: Path, anno_fa: Path) -> Path:
	logging.info("Running hmmscan2")

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

def prescan(run_id: str, input: Path, outdir: Path) -> List:
	anno_fa = prodigal(run_id, input, outdir)
	scan_file = hmmscan(run_id, input, outdir, anno_fa)

	df = pd.read_table(scan_file, sep=r"\s+", comment='#', header=None, usecols=range(5), 
					   names=['target_name', 'accession', 'query_name', 'accession_fs', 'E-value_fs'])
	df.drop_duplicates(subset='query_name', inplace=True)
	df['key'] = df['query_name'].str.split('_').str[:-1].str.join('_')

	chosen = []
	for k, group in df.groupby('key'):
		if scanf(group['target_name'].tolist()):
			chosen.append(k)

	return chosen

def prokka(run_id: str, input: Path) -> Path:
	logging.info("Running Prokka")
	outdir = input.parent / 'prokka'
	cmd = [
		'prokka',
		input,
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

def ICEscan(contig_id, prokka_dir, ice_dir) -> Path:
	logging.info("Running Macsyfinder on ICEscan model")
	anno_fa = prokka_dir / f'{contig_id}.faa'
	ice_res = ice_dir / 'all_systems.tsv'

	cmd = [
		'macsyfinder',
		'--db-type', 'ordered_replicon',
		'--models-dir', Path(__file__).parents[1] / 'data' / 'macsydata',
		'--models', 'ICEscan', 'all',
		'--replicon-topology', 'linear',
		'--coverage-profile', '0.3',
		'--sequence-db', anno_fa,
		'-o', ice_res
	]

	try:
		subprocess.run(cmd, check=True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	except subprocess.CalledProcessError as e:
		logging.error(f"Macsyfinder execution failed: {e}")

	return ice_res



# Unsure
def rename(run_id, input):
	run_dir = os.path.join(tmp_dir,run_id)
	if not os.path.exists(run_dir):
		os.makedirs(run_dir)
	if not os.path.exists(gb_dir):
		os.mkdir(gb_dir)

	filename = os.path.basename(input)
	resultf = filename.rsplit('.', 1)[0]
	input1 = os.path.join(os.path.dirname(input), resultf)

	i = 1
	id_dict = {}
	newIDfa = os.path.join(run_dir,run_id+'_newID.fa')
	newfa = open(newIDfa,'w')

	for seq_record in SeqIO.parse(input, "fasta"):
		if len(seq_record.id) > 15:
			realID = seq_record.id[:15]+'...'
		else:
			realID = seq_record.id
		contigID = 'contig_' + str(i)
		seq_record.id = contigID
		seqfa = str(seq_record.seq)
		newfa.write(">%s\n%s\n" % (contigID,seqfa))
		id_dict[contigID] = realID
		i += 1
	return id_dict


def gc(fasta_file,start,end):

	record = SeqIO.read(fasta_file, "fasta")
	sequence = record.seq[start-1:end]
	gcs = str("%.2f"%GC(sequence))

	return gcs

def calculate_gc(fasta_file, start, end, window_size, step_size):

	record = SeqIO.read(fasta_file, "fasta")
	if start == 0:
		start = 1
	sequence = record.seq[start-1:end]

	windows = []
	gc_contents = []
	pos = []
	j = start/1000 + 0.025
	for i in range(0, len(sequence) - window_size + 1, step_size):
		window = sequence[i:i+window_size]
		gc_content = GC(window)
		gc_contents.append(gc_content)
		pos.append(round(j, 4))
		j += 0.05

	gcdict = {
		        'xData':pos,
		        'datasets':[{
		            'name':'',
		            'data':gc_contents,
		            'unit':'%',
		            'type':'line',
		            "valueDecimals": 1
		        }]
	    }

	return gcdict


def getgff(run_id):

	gffile = os.path.join(gb_dir, run_id + '.gff')
	with open(gffile,'r') as gffin:
		trnadict = {}
		posdict = {}
		for line in gffin.readlines():
			if 'ID=' in line:
				lines = line.strip().split('\t')
				ids = lines[8].split(';')[0].split('=')[1]
				header = ids.split('_')[0]
				product = lines[8].split('product=')[1]
				pos = [lines[3],lines[4],lines[6],product]
				if lines[2] == 'tRNA' or lines[2] == 'tmRNA':
					trnadict[ids] = product
				posdict[ids] = pos
				totalnum = getnum(ids)

	return trnadict,posdict,header,totalnum

def getnum(ID):

	return int(ID.split('_')[1].lstrip("0"))

def find_max_distance(numbers):
	max_distance = -1
	max_distance_index = -1

	for i in range(len(numbers) - 1):
		distance = abs(numbers[i] - numbers[i+1])
		if distance > max_distance:
			max_distance = distance
			max_distance_index = i

	if max_distance_index == -1:
		return None

	return numbers[max_distance_index], numbers[max_distance_index + 1]

def pos_tag(pos,posdict,ICE,final,totalnum,dirtag):

	tICE = ICE
	tfinal = final
	for k,v in posdict.items():
		vstart,vend = int(v[0]),int(v[1])
		if int(pos) <= vend:
			if dirtag == 's':
				tICE = getnum(k)
				tfinal = max(1, tICE - 5)
			else:
				if vstart > int(pos):
					tICE = getnum(k) - 1 
				else:
					tICE = getnum(k)
				tfinal = min(totalnum, tICE + 5)
			break
	return tICE, tfinal

def merge_tRNA(run_id,ICEdict,DRlist):

	trnadict,posdict,header,totalnum = getgff(run_id)
	fICE = getnum(next(iter(ICEdict)))
	eICE = getnum(list(ICEdict.keys())[-1])

	nfICEnum = max(1, fICE - 5)
	neICEnum = min(totalnum, eICE + 5)

	ICEtagnum = [nfICEnum,neICEnum]
	trnalist = []
	for key,value in trnadict.items():
		if nfICEnum <= getnum(key) <= neICEnum:
			ICEtagnum.append(getnum(key))
			trnalist.append(value)

	ICEtagnum.sort()
	finalstart,finalend = find_max_distance(ICEtagnum)

	myDR1 = posdict[zill(header,fICE)][0]
	myDR2 = ''
	myDR3 = ''
	myDR4 = posdict[zill(header,eICE)][1]

	if trnalist:
		if finalstart == nfICEnum:
			eICE = finalend
			finalend = min(totalnum, finalend + 5)
			myDR4 = posdict[zill(header,eICE)][1]					
			for line in DRlist:
				DRs = line.split('|')
				if int(DRs[3]) - int(DRs[0]) > 500000:
					continue
				if int(DRs[3]) - int(DRs[0]) < 5000:
					continue					
				if int(posdict[zill(header,eICE)][0]) < int(DRs[3]) < int(posdict[zill(header,eICE)][1]):
					checktrna = 0
					for key,value in trnadict.items():
						if int(DRs[0]) <= value[0] <= int(DRs[3]) and int(DRs[0]) <= value[1] <= int(DRs[3]):
							checktrna += 1
					if checktrna >= 2:
						break

					fICE,finalstart = pos_tag(DRs[0],posdict,fICE,finalstart,totalnum,'s')
					myDR1 = DRs[0]
					myDR2 = DRs[1]
					myDR3 = DRs[2]
					myDR4 = DRs[3]
					break

		elif finalend == neICEnum:
			fICE = finalstart
			finalstart =  max(1, finalstart - 5)
			myDR1 = posdict[zill(header,fICE)][0]
			for line in DRlist:
				DRs = line.split('|')
				if int(DRs[3]) - int(DRs[0]) > 500000:
					continue	
				if int(DRs[3]) - int(DRs[0]) < 5000:
					continue									
				if int(posdict[zill(header,fICE)][0]) < int(DRs[0]) < int(posdict[zill(header,fICE)][1]):
					checktrna = 0
					for key,value in trnadict.items():
						if int(DRs[0]) <= value[0] <= int(DRs[3]) and int(DRs[0]) <= value[1] <= int(DRs[3]):
							checktrna += 1
					if checktrna >= 2:
						break
					eICE,finalend = pos_tag(DRs[3],posdict,eICE,finalend,totalnum,'e')
					myDR1 = DRs[0]
					myDR2 = DRs[1]
					myDR3 = DRs[2]
					myDR4 = DRs[3]									
					break

	return myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,posdict,header,trnalist

def get_DR(run_id,input):

	DRindex = os.path.join(tmp_dir, run_id, run_id+'_DR')
	DRout = os.path.join(tmp_dir, run_id, run_id+'_DRout')
	maktree_cmd = './tool/mkvtree -db ' + input +' -indexname '+ DRindex +' -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1'
	vmatch_cmd = './tool/vmatch -l 15 '+ DRindex + ' > '+ DRout
	os.system(maktree_cmd)
	os.system(vmatch_cmd)

	DRlist = []
	with open(DRout,'r') as DRin:
		for line in DRin.readlines():
			lines = line.strip().split()
			if not line.startswith('#'):
				DR = [lines[2],str(int(lines[2])+int(lines[0])),lines[6],str(int(lines[6])+int(lines[4]))]
				DRlist.append('|'.join(DR))
	return DRlist

def get_ICE(contig_id, input, ice_dir, prokka_dir):

	ice_res = ICEscan(contig_id, prokka_dir, ice_dir)

	with open(ice_res,'r') as ICEin:
		ICEdict = {}
		infodict = {}
		for line in ICEin.readlines():
			if 'Chromosome' in line:
				lines = line.strip().split('\t')
				if 'UserReplicon_IME' not in lines[5]:
					gbname = lines[1]
					tags = get_feat(lines[2])
					mpf = lines[4].split('/')[-1].split('_')[1]

					if 'Relaxase@' in tags:
						mob = tags.split('@')[1]
					else:
						mob = ''
					ICEtag = 'ICE'+lines[5].split('_')[-1]

					ICEdict.setdefault(ICEtag,{})[gbname]=tags

					if ICEtag not in infodict:
						infodict[ICEtag] = {'mob': [], 'mpf': []}
					if mob not in infodict[ICEtag]['mob']:
						if mob:
							infodict[ICEtag]['mob'].append(mob)
					if mpf not in infodict[ICEtag]['mpf']:
						infodict[ICEtag]['mpf'].append(mpf)	

	dictICE = {}
	posdict = {}
	trnalist = []
	header = ''
	DRlist = get_DR(contig_id,input)
	for key,value in ICEdict.items():
		myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,posdict,header,trnalist = merge_tRNA(contig_id,value,DRlist)
		dictICE[key] = [myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend]

	return dictICE,ICEdict,posdict,header,trnalist,infodict

def args(run_id):

	return getblast(run_id)

def zill(header,num):

	return header+ '_' +str(num).zfill(5)

def get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product):

	feature = [feature]
	product = [product]

	if gene in argdict:
		feature.append('AR')
		product.append(argdict[gene])
	if gene in vfdict:
		feature.append('VF')
		product.append(vfdict[gene])
	if gene in isdict:
		feature.append('IS')
		product.append(isdict[gene])
	if gene in dfdict:
		feature.append('Defense')
		product.append(dfdict[gene])

	if gene in metaldict:
		feature.append('Metal')
		product.append(metaldict[gene])
	if gene in popdict:
		feature.append('Degradation')
		product.append(popdict[gene])
	if gene in symdict:
		feature.append('Symbiosis')
		product.append(symdict[gene])

	feature = '; '.join(list(filter(None, feature)))
	product = '; '.join(list(filter(None, product)))

	return feature,product

def oritseq(run_id, regi, input, start, end):

	oritseq = '-'
	fafile = os.path.join(tmp_dir,run_id,regi+'_fororit.fa')	
	with open(fafile,'w') as orif:
		seq = getfa(input,start,end)
		orif.write('>fororit\n')
		orif.write(seq)

	oriT_Database = os.path.join(workdir,'data','oriT_db')
	blastn_out = os.path.join(tmp_dir,run_id,regi+'_oriTout')
	blast_cmd = [blastn, "-db", oriT_Database, "-query", fafile, "-evalue 0.01 -word_size 11 -outfmt '6 std qlen slen' -num_alignments 1 -out", blastn_out,">/dev/null"]
	os.system(' '.join(blast_cmd))

	with open(blastn_out,'r') as oritout:
		for line in oritout.readlines():
			lines = line.strip().split()
			if lines[0]:
				matchl = int(lines[3])
				slen = int(lines[13])
				ident = float(lines[2])
				hvalue = (matchl/slen)*ident	
				if hvalue > 0.49:
					oritseq = getfa(fafile,str(int(lines[6])-1),lines[7])
					break
	return oritseq

def get_feat(feat):

	featuredict = {
		'Phage_integrase':'Integrase','UPF0236':'Integrase',
        'Recombinase':'Integrase','rve':'Integrase',
        'TIGR02224':'Integrase','TIGR02249':'Integrase',
        'TIGR02225':'Integrase','PB001819':'Integrase'
	}
	if feat in featuredict:
		return 'Integrase@'+feat
	elif 'T4SS_MOB' in feat:
		tag = feat.split('_')[1]
		return 'Relaxase@'+tag
	elif 't4cp' in feat:
		tag = feat.split('_')[1]
		return 'T4CP@' + tag
	elif 'tcpA' in feat:
		tag = feat.split('_')[1]
		return 'T4CP@' + tag
	else:
		return 'T4SS@'+feat.replace('T4SS_','')

def getcolor(feature,product):

	coldict = {'DR':'black','Gene':'#C0C0C0',
		   'Hyp':'#DCDCDC','Integrase':'blue',
		   'Transposase':'yellow','T4SS':'lightpink','T4CP':'orange',
		   'Relaxase':'brown','AR':'red','tRNA':'black',
		   'Flank':'gray','VF':'#ba8448','Defense':'#00B050',
		   'Metal':'#03A89E','Degradation':'#640B0F','Symbiosis':'#FFFFCD'
	}

	namedict = {'Hyp':'Hypothetical protein','Gene':'Other gene',
		    'AR':'Antibiotic resistance gene',
		    'VF':'Virulence factor','Metal':'Metal resistance',
		    'Flank':'Flank region','Defense':'Defense system',
		    'Transposase':'Transposase','Relaxase':'Relaxase',
		    'T4CP':'T4CP','T4SS':'T4SS','Integrase':'Integrase',
		    'Degradation':'Degradation','Symbiosis':'Symbiosis'
	}

	if 'Integrase' in feature:
		feature = 'Integrase'
	elif 'T4SS' in feature:
		feature = 'T4SS'
	elif 'T4CP' in feature:
		feature = 'T4CP'
	elif 'Relaxase' in feature:
		feature = 'Relaxase'
	elif 'IS' in feature:
		feature = 'Transposase'
	elif 'VF' in feature:
		feature = 'VF'
	elif 'AR' in feature:
		feature = 'AR'
	elif 'Defense' in feature:
		feature = 'Defense'
	elif 'Metal' in feature:
		feature = 'Metal'
	elif 'Degradation' in feature:
		feature = 'Degradation'
	elif 'Symbiosis' in feature:
		feature = 'Symbiosis'

	elif feature == 'Flank':
		feature == 'Flank'
	elif feature == '':
		if product == 'hypothetical protein':
			feature = 'Hyp'
		else:
			feature = 'Gene'
	else:
		feature = 'Gene'

	return coldict[feature], namedict[feature]

def gstrand(instra):

	strands = {'+' : 1, '-' : -1}
	return strands[instra]

def getfa(input,s,e):

	seq_record = SeqIO.read(input, "fasta")
	sequence = seq_record.seq[int(s):int(e)]
	return str(sequence)

def get_map(contig_id: str, id_dict, fasta_file: Path, outdir: Path, prokka_dir: Path):
	print(outdir, contig_id)
	js_dir = outdir / 'js'
	js_dir.mkdir(exist_ok=True)
	gc_map = outdir / 'script' / 'js' / 'gcmap.js'
	view_file = outdir / 'script' / 'js' / 'view.html'
	ice_dir = outdir / f'{contig_id}_ICE'
	# Aquí da fallo, ¿no está borrando?
	if ice_dir.exists() and ice_dir.is_dir():
		print(f"removing {ice_dir}")
		shutil.rmtree(ice_dir)
	ice_dir.mkdir()

	dictICE,ICEdict,posdict,header,trnalist,infodict = get_ICE(contig_id, fasta_file, ice_dir, prokka_dir)
	argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict = args(contig_id)

	ICEss = {}
	for key,value in dictICE.items():
		genelist = []
		regi = contig_id+'_'+key
		regijs = 'contig_'+contig_id.split("_contig_", 1)[-1] +'_'+key
		genefile = os.path.join(outdir,regi+'_gene.json')
		infofile = os.path.join(outdir,regi+'_info.json')
		gcjson = os.path.join(js_dir,regijs+'_gc.js')
		mapfile = os.path.join(js_dir,regijs+'.js')
		htmlfile = os.path.join(outdir,regi+'.html')
		[myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend] = value

		start = finalstart
		while start < fICE:
			gene = zill(header,start)
			s,e,strand,pro = posdict[gene]
			pos = s+'..'+e+' ['+strand+'], '+str(int(e)-int(s)+1)

			feature = 'Flank'
			product = pro
			feature,product = get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product)
			if 'hypothetical protein;' in product:
				product = product.replace('hypothetical protein;','')

			start += 1
			content = {
					'gene':gene,
					'pos':pos,
					'prod': product,
					'featu': feature
			}
			genelist.append(content)

		mov = fICE
		while mov <= eICE:
			gene = zill(header,mov)
			s,e,strand,pro = posdict[gene]
			pos = s+'..'+e+' ['+strand+'], '+str(int(e)-int(s)+1)

			if gene in ICEdict[key]:
				[feature,pro11] = ICEdict[key][gene].split('@')
			else:
				feature,pro11 = '',''

			if pro11:
				if pro == 'hypothetical protein':
					product = pro11
				else:
					product = pro+', '+ pro11
			else:
				product = pro

			feature,product = get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product)
			mov += 1
			content = {
					'gene':gene,
					'pos':pos,
					'prod': product,
					'featu': feature
			}
			genelist.append(content)				

		while mov <= finalend:
			gene = zill(header,mov)
			s,e,strand,pro = posdict[gene]
			pos = s+'..'+e+' ['+strand+'], '+str(int(e)-int(s)+1)

			feature = 'Flank'
			product = pro
			feature,product = get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product)
			if 'hypothetical protein;' in product:
				product = product.replace('hypothetical protein;','')

			mov += 1
			content = {
					'gene':gene,
					'pos':pos,
					'prod': product,
					'featu': feature
			}
			genelist.append(content)

		with open(genefile,'w') as gene_file:
			json.dump(genelist, gene_file, indent=4)

		contigID = contig_id.split('_', 1)[1]

		sgene = zill(header,fICE)
		egene = zill(header,eICE)
		s1,e1,strand1,pro1 = posdict[sgene]
		s2,e2,strand2,pro2 = posdict[egene]
		if myDR1 == '0':
			myDR1 = '1'

		ICEss[regi] = '|'.join([myDR1,myDR4,str(fICE),str(eICE)])

		# host = spdict[contigID]
		host = "-"
		gcc = gc(fasta_file,int(myDR1),int(myDR4))
		source = id_dict[contigID]

		if myDR2:
			DR1 = getfa(fasta_file,myDR1,myDR2)
			DR2 = getfa(fasta_file,myDR3,myDR4)
			DRw = 'attL:'+myDR1+'..'+myDR2+'('+DR1+')  '+'attR:'+myDR3+'..'+myDR4+'('+DR2+')'
		else:
			DRw = '-'

		oritseqs = oritseq(contig_id, regi, fasta_file, myDR1,myDR4)
#		oritdesc = "<br>".join([oritseqs[i:i+63] for i in range(0, len(oritseqs), 63)])

		ICEinfo = {
			'Contig source':source,
			'Host Strain':host,
			'GC Content (%)':gcc,
			'Length (bp)':str(int(e2)-int(s1)+1),
			'oriT seq':oritseqs,
			'DRs':DRw,
			'Relaxase Type': ','.join(infodict[key]['mob']),
			'Mating pair formation systems':','.join(infodict[key]['mpf']),
			'Close to tRNA':','.join(trnalist)
		}
		with open(infofile,'w') as info_file:
			json.dump(ICEinfo, info_file, indent=4)

		i = 1
		mapzlist = []
		mapflist = []
		for gene in genelist:
			color, name = getcolor(gene['featu'],gene['prod'])
			start = gene['pos'].split(' ')[0].split('..')[0]
			end = gene['pos'].split(' ')[0].split('..')[1]
			strand = gstrand(gene['pos'].split('[')[1].split(']')[0])
			product = gene['prod']

			if product == '':
				product = 'hypothetical protein'

			anno = {
					'start' : start,
					'end' : end,
					'strand' : strand,
					'locus_tag' : 'M'+ str(i),
					'type' : 'others',
					'color' : color,
					'description' : 'Location: '+gene['pos'].split(' ')[0]+' ('+gene['pos'].split(' ')[2]+' bp)<br>Type: ' +name +'<br>Detail: '+ product
				}
			if strand == 1:
				mapzlist.append(anno)
			else:
				mapflist.append(anno)
			i += 1

		head = 'var borders = [];\nvar tta_codons = [];\nvar orfs ='
		s = genelist[0]['pos'].split(' ')[0].split('..')[0]
		e = genelist[-1]['pos'].split(' ')[0].split('..')[1]

		gcdict = calculate_gc(fasta_file, int(s), int(e), 500, 50)
		with open(gc_map, 'r') as original_file:
			original_content = original_file.read()
		with open(gcjson,'w') as gein2:
			gein2t = 'var jsonData = ' + str(gcdict)+';'
			gein2.write(gein2t)
			gein2.write(original_content) 

#		with open(gcjson,'w') as gein2:
#			json.dump(gcdict, gein2, indent=4)

		maps = str(mapzlist)+';\nvar orfs2 ='+str(mapflist)+';\nvar clusterf2 = { start: '+s+', end: '+ \
				  e+', idx: 1, orfs: orfs, borders: borders, tta_codons:tta_codons,\
				  label: \'\', unordered: true };\nvar clusterr2 = { start: '+ s+', end: '+ \
				  e+', idx: 2, orfs: orfs2, borders: borders, tta_codons:tta_codons,\
				  label: \'\', unordered: true };\nsvgene.drawClusters("'+regijs+'", [clusterf2, clusterr2], 50, 920);'
		with open(mapfile,'w') as map_file:
			map_file.write(head + maps)

		with open(view_file, 'r') as file:
			file_content = file.read()
		new_content = file_content.replace('XXXX', regijs)
		with open(htmlfile, 'w') as file:
			file.write(new_content)

	return ICEss

def delete_folders_starting_with_keyword(dir, keyword):
	for dirpath, dirnames, filenames in os.walk(dir, topdown=False):
		for dirname in dirnames:
			if dirname.startswith(keyword):
				folder_to_remove = os.path.join(dirpath, dirname)
				shutil.rmtree(folder_to_remove)

def getfasta(run_id,resultdir,id_dict,key,s,e,stag,etag):

	fafile = os.path.join(tmp_dir, run_id, run_id+'.fa')	
	faafile = os.path.join(tmp_dir, 'gbk', run_id+'.faa')

	outfa = os.path.join(resultdir,key+'.fa')
	outfaa = os.path.join(resultdir,key+'.faa')

	seq_record = SeqIO.read(fafile, "fasta")
	with open(outfa, "w") as output_handle1:
		sequence = seq_record.seq[int(s)-1:int(e)]
		ID = '_'.join(seq_record.id.split('_')[-2:])
		seq_record.description = ''
		seq_record.seq = sequence
		seq_record.id = id_dict[ID] + ' ' + s +'-'+e
		SeqIO.write(seq_record, output_handle1, "fasta")

	faa_records = SeqIO.parse(faafile, "fasta")
	with open(outfaa, "w") as output_handle2:
		for faa_record in faa_records:
			seq_id = getnum(faa_record.id)
			if int(stag) <= seq_id <= int(etag):
				SeqIO.write(faa_record, output_handle2, "fasta")

def _meta(run_id: str, input: Path, outdir: Path, kraken_db: Union[Path, None], gff: Union[Path, None]) -> None:

	jsback = outdir / 'script' / 'js'

	# En cuarentena, quizá recomendable borrar?
	id_dict = rename(run_id,input)

	# ICE markers annotations
	chosen_contigs = prescan(run_id, input, outdir)

	# Taxonomy assignation (optional)
	drawout = None
	spdict = {}
	report = None
	if kraken_db:
		drawout,spdict,report = taxonomy(run_id, input, outdir, kraken_db)
		copy_files(report, outdir)
		copy_files(drawout, outdir)

	# Get basic information from input
	basefile = getbase(run_id, input, outdir)
	copy_files(basefile, outdir)

	#############################################################
	##################### Por aquí me quedé #####################
	#############################################################
	ICEsum = outdir / tmp_dir / f'{run_id}_ICEsum.json'
	i = 1 
	ICEsumlist = []
	for seq_record in SeqIO.parse(input, "fasta"):
		if seq_record.id in chosen_contigs:

			# Extract contig of interest in fasta format
			contig_folder = outdir / tmp_dir / f'{seq_record.id}'
			contig_folder.mkdir(exist_ok=True)
			contig_fasta = contig_folder / f'{seq_record.id}.fasta'
			with contig_fasta.open("w") as oh:
				SeqIO.write(seq_record, oh, "fasta")

			# Annotate contig
			prokka_dir = prokka(seq_record.id, contig_fasta)
			ICEss = get_map(seq_record.id, id_dict, contig_fasta, contig_folder, prokka_dir)
			copy_files(contig_folder, outdir)

			if ICEss:
				for key,value in ICEss.items():
					[s,e,stag,etag] = value.split('|')
					lengt = int(e) - int(s) + 1
					ICEs = {
						'id' : str(i),
				        'seqid': id_dict[seq_record.id],
				        # 'species':spdict[seq_record.id],
						'species':'',
				        'location': s+'..'+e,
				        'length': lengt,
				        'detail': key
				    }
					ICEsumlist.append(ICEs)
					getfasta(seq_record.id,outdir,id_dict,key,s,e,stag,etag)
					i += 1 

	with open(ICEsum,'w') as ice_file:
		json.dump(ICEsumlist, ice_file, indent=4)

	copy_files(ICEsum, outdir)
	jsdir = os.path.join(outdir,'js')
	copy_files(jsback, jsdir)
	delete_folders_starting_with_keyword(tmp_dir, run_id)
