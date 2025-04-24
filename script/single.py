#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,time
import json
import shutil
import subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqUtils import GC
import logging
from typing import List
from pathlib import Path
from Bio.SeqFeature import CompoundLocation
from script.function import getblast

logging.basicConfig(
		level=logging.INFO,
		datefmt='%m/%d/%Y %I:%M:%S %p'
	)

workdir = './result'

tmp_dir = os.path.join(workdir,'tmp') 
in_dir = os.path.join(tmp_dir,'fasta')
gb_dir = os.path.join(tmp_dir,'gbk')



def get_time():

	return time.asctime( time.localtime(time.time()) )

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

def prokka(run_id: str, infile: Path, cpus: int = 8):
	logging.info("Running Prokka")
	cmd = [
		'prokka',
		infile,
		'--force',
		'--fast',
		'--quiet',
		'--cdsrnaolap',
		'--cpus', str(cpus),
		'--outdir', gb_dir,
		'--prefix', run_id
	]

	try:
		subprocess.run(cmd, check=True)
	except subprocess.CalledProcessError as e:
		logging.error(f"Prokka annotation failed: {e}")

def ICEscan(run_id = str, gb_dir = Path, ICE_dir = Path, cov_profile: float = 0.3):
	logging.info("Running Macsyfinder on ICEscan model")
	anno_fa = gb_dir / f'{run_id}.faa'

	cmd = [
		'macsyfinder',
		'--db-type', 'ordered_replicon',
		'--models-dir', './data/macsydata/',
		'--models', 'ICEscan', 'all',
		'--replicon-topology', 'linear',
		'--coverage-profile', str(cov_profile),
		'--sequence-db', anno_fa,
		'-o', ICE_dir
	]
	try:
		subprocess.run(cmd, check=True)
	except subprocess.CalledProcessError as e:
		logging.error(f"ICEscan annotation failed: {e}")

def getgff1(run_id):
	gffile = os.path.join(gb_dir, run_id + '.gff')
	locusdict = {}
	with open(gffile,'r') as gffin:
		trnadict = {}
		posdict = {}
		for line in gffin.readlines():
			if 'ID=' in line:
				lines = line.strip().split('\t')
				ids = lines[8].split(';')[0].split('=')[1]
				locusdict[ids] = ids		
				header = ids.split('_')[0]
				product = lines[8].split('product=')[1]
				pos = [lines[3],lines[4],lines[6],product]
				if lines[2] == 'tRNA' or lines[2] == 'tmRNA':
					trnadict[ids] = pos
				posdict[ids] = pos
				totalnum = getnum(ids)

	return [trnadict,posdict,header,totalnum,locusdict]

def getgff(run_id: str, gb_dir: Path) -> List:
	gb_file = gb_dir / f'{run_id}.gbk'
	faa_file = gb_dir / f'{run_id}.faa'
	ffn_file = gb_dir / f'{run_id}.ffn'

	locus_dict = {}
	pos_dict = {}
	trna_dict = {}

	prefix = "TMPID"
	count = 0

	with faa_file.open('w') as faa_out, ffn_file.open('w') as ffn_out:
		for record in SeqIO.parse(gb_file, "genbank"):
			for feature in record.features:
				quals = feature.qualifiers
				if "locus_tag" not in quals:
					continue
				orig_id = quals["locus_tag"][0]

				seq_nt = feature.location.extract(record).seq

				prod = (quals.get("product", ["tmRNA"])[0]
					if feature.type != "tmRNA"
					else "tmRNA")

				new_id = zill(prefix, count + 1)
				locus_dict[new_id] = orig_id
				start = int(feature.location.start)
				end = int(feature.location.end)
				strand = gstrand1(feature.location.strand)
				pos_dict[new_id] = [str(start), str(end), str(strand), prod]

				if feature.type in ("CDS", "rRNA"):
					if feature.type == "CDS":
						aa_seq = seq_nt.translate(to_stop=True)
						faa_out.write(f">{new_id} {prod}\n{aa_seq}\n")
					ffn_out.write(f">{new_id} {prod}\n{seq_nt}\n")

				if feature.type in ("tRNA", "tmRNA"):
					trna_dict[new_id] = pos_dict[new_id]

				count += 1
	return [trna_dict, pos_dict, prefix, count-1, locus_dict]

def getnum(ID):

	return int(ID.split('_')[1].lstrip("0"))

def zill(header,num):

	return header+ '_' +str(num).zfill(5)

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

def merge_tRNA(run_id,ICEdict,DRlist,listgff):

	[trnadict,posdict,header,totalnum,locusdict] = listgff

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
		if finalend == neICEnum:
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
						if int(DRs[0]) <= int(value[0]) <= int(DRs[3]) and int(DRs[0]) <= int(value[1]) <= int(DRs[3]):
							checktrna += 1
					if checktrna >= 2:
						break

					eICE,finalend = pos_tag(DRs[3],posdict,eICE,finalend,totalnum,'e')
					myDR1 = DRs[0]
					myDR2 = DRs[1]
					myDR3 = DRs[2]
					myDR4 = DRs[3]									
					break

		elif finalstart == nfICEnum:
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
						if int(DRs[0]) <= int(value[0]) <= int(DRs[3]) and int(DRs[0]) <= int(value[1]) <= int(DRs[3]):
							checktrna += 1
					if checktrna >= 2:
						break

					fICE,finalstart = pos_tag(DRs[0],posdict,fICE,finalstart,totalnum,'s')
					myDR1 = DRs[0]
					myDR2 = DRs[1]
					myDR3 = DRs[2]
					myDR4 = DRs[3]
					break
	return myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,posdict,header,trnalist,locusdict

def get_DR(run_id,infile):

	DRindex = os.path.join(tmp_dir, run_id, run_id+'_DR')
	DRout = os.path.join(tmp_dir, run_id, run_id+'_DRout')
	maktree_cmd = './tool/mkvtree -db ' + infile +' -indexname '+ DRindex +' -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1'
	vmatch_cmd = './tool/vmatch -l 15 '+ DRindex + ' > '+ DRout
	os.system(maktree_cmd)
	os.system(vmatch_cmd)

	DRlist = []
	with open(DRout,'r') as DRin:
		for line in DRin.readlines():
			lines = line.strip().split()
			if not line.startswith('#'):
				DR = [str(int(lines[2])+1),str(int(lines[2])+int(lines[0])),str(int(lines[6])+1),str(int(lines[6])+int(lines[4]))]
				DRlist.append('|'.join(DR))
	return DRlist

def oritseq(run_id, regi, infile, start, end):
	logging.info("Running Blast for OriT detection")
	oritseq = '-'
	fafile = os.path.join(tmp_dir,run_id,regi+'_fororit.fa')	
	with open(fafile,'w') as orif:
		seq = getfa(infile,start,end)
		orif.write('>fororit\n')
		orif.write(seq)

	oriT_Database = os.path.join(workdir,'data','oriT_db')
	blastn_out = os.path.join(tmp_dir,run_id,regi+'_oriTout')
	blast_cmd = ['blastn', "-db", oriT_Database, "-query", fafile, "-evalue 0.01 -word_size 11 -outfmt '6 std qlen slen' -num_alignments 1 -out", blastn_out,">/dev/null"]
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

def ICE_filter(ice_res):
	logging.info("Filtering ICEs")
	if not ice_res.is_file():
		raise FileNotFoundError(f"ICE result file not found: {ice_res}")

	# Initialize dictionaries and lists
	ice_dict = {'ICE': []}
	ime_dict = defaultdict(list)
	ime_genes = defaultdict(list)
	aice_dict = defaultdict(list)
	ice_tags = []
	aice_tags = []

	# Read and filter lines
	with ice_res.open('r') as infile:
		for line in infile:
			if 'Chromosome' not in line:
				continue
			cols = line.rstrip().split('\t')
			# Expect at least 8 columns
			if len(cols) < 8 or cols[7] != '1':
				continue

			id_tag = cols[3]
			ice_tag = cols[5]
			feature_source = cols[2]

			if 'IME' in ice_tag:
				ime_dict[ice_tag].append(id_tag)
				ime_genes[ice_tag].append(feature_source)

			elif 'AICE' in ice_tag:
				aice_dict[ice_tag].append(id_tag)
				if ice_tag not in aice_tags:
					aice_tags.append(ice_tag)

			else:
				ice_dict['ICE'].append(id_tag)
				if ice_tag not in ice_tags:
					ice_tags.append(ice_tag)

	# Define gene indicator mapping
	integrase_map = {
		'Phage_integrase': 'Integrase', 'UPF0236': 'Integrase',
		'Recombinase': 'Integrase', 'rve': 'Integrase',
		'TIGR02224': 'Integrase', 'TIGR02249': 'Integrase',
		'TIGR02225': 'Integrase', 'PB001819': 'Integrase'
	}

	# Identify IME tags with both MOB and Integrase
	valid_ime_tags: List[str] = []
	for tag, genes in ime_genes.items():
		gene_types = set()
		for gene in genes:
			if 'Relaxase_' in gene or 'T4SS_MOB' in gene:
				gene_types.add('MOB')
			elif gene in integrase_map:
				gene_types.add('Int')
		if gene_types == {'MOB', 'Int'}:
			valid_ime_tags.append(tag)

	# Filter IME tags where their IDs are not a subset of ICE IDs
	ime_tags: List[str] = []
	ice_ids = set(ice_dict['ICE'])
	for tag, ids in ime_dict.items():
		if tag in valid_ime_tags and not set(ids).issubset(ice_ids):
			ime_tags.append(tag)

	# Combine results preserving order: ICE, IME, AICE
	return ice_tags + ime_tags + aice_tags

def get_ICE(run_id = str, infile = Path, listgff = List, outdir = Path, tmp_dir = Path, gb_dir = Path):

	ICE_dir = tmp_dir / f'{run_id}_ICE'
	ICE_dir.mkdir(parents=True, exist_ok=True)
	ICE_res = ICE_dir / 'all_systems.tsv'

	shutil.rmtree(ICE_dir, ignore_errors=True)
	ICEscan(run_id, gb_dir, ICE_dir)
	ftag = ICE_filter(ICE_res)

	with ICE_res.open('r') as ICEin:
		ICEdict = {}
		infodict = {}
		for line in ICEin.readlines():
			if 'Chromosome' in line:
				lines = line.strip().split('\t')
				if lines[7] != '1':
					continue
				elif lines[5] not in ftag:
					continue
				else:									
					gbname = lines[1]
					tags = get_feat(lines[2])

					if 'T4SS' in lines[4]:
						mpf = lines[4].split('/')[-1].split('_')[1]
					else:
						mpf = ''
					if 'Relaxase@' in tags:
						mob = tags.split('@')[1]
					else:
						mob = ''
					ICEtag = 'ICE'+lines[5].split('_')[-1]

					if 'IME' in lines[5]:
						ICEtag = 'IME'+lines[5].split('_')[-1]
					elif 'AICE' in lines[5]:
						ICEtag = 'AICE'+lines[5].split('_')[-1]
					else:
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
	DRlist = get_DR(run_id,infile)

	for key,value in ICEdict.items():
		myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,posdict,header,trnalist,locusdict = merge_tRNA(run_id,value,DRlist,listgff)
		dictICE[key] = [myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,trnalist,locusdict]

	return dictICE,ICEdict,posdict,header,infodict					

def args(run_id):

	return getblast(run_id)

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

def get_feat(feat):
	logging.info(f"Extracting ICE features: {feat}")
	featuredict = {
		'Phage_integrase':'Integrase','UPF0236':'Integrase',
		'Recombinase':'Integrase','rve':'Integrase',
		'TIGR02224':'Integrase','TIGR02249':'Integrase',
		'TIGR02225':'Integrase','PB001819':'Integrase',
		'RepSAv2':'Rep','DUF3631':'Rep','Prim-Pol':'Rep',
		'FtsK_SpoIIIE':'Tra'}

	if feat in featuredict:
		return featuredict[feat]+'@'+feat
	elif 'T4SS_MOB' in feat:
		tag = feat.split('_')[1]
		return 'Relaxase@'+tag
	elif 'Relaxase_' in feat:
		tag = feat.split('_')[1:]
		return 'Relaxase@'+'_'.join(tag)
	elif 't4cp' in feat:
		tag = feat.split('_')[1]
		return 'T4CP@' + tag
	elif 'tcpA' in feat:
		tag = feat.split('_')[1]
		return 'T4CP@' + tag
	elif 'FATA_' in feat or 'FA_' in feat:
		tag = feat.split('_')[1]
		return 'T4SS@' + tag
	else:
		return 'T4SS@'+feat.replace('T4SS_','')

def getcolor(feature,product):

	coldict = {'DR':'black','Gene':'#C0C0C0',
		   'Hyp':'#DCDCDC','Integrase':'blue',
		   'Transposase':'yellow','T4SS':'lightpink','T4CP':'orange',
		   'Relaxase':'brown','AR':'red','tRNA':'black',
		   'Flank':'gray','VF':'#ba8448','Defense':'#00B050',
		   'Metal':'#03A89E','Degradation':'#640B0F','Symbiosis':'#FFFFCD',
		   'Rep':'black','Tra':'black'
	}

	namedict = {'Hyp':'Hypothetical protein','Gene':'Other gene',
			'AR':'Antibiotic resistance gene',
			'VF':'Virulence factor','Metal':'Metal resistance',
			'Flank':'Flank region','Defense':'Defense system',
			'Transposase':'Transposase','Relaxase':'Relaxase',
			'T4CP':'T4CP','T4SS':'T4SS','Integrase':'Integrase',
			'Degradation':'Degradation','Symbiosis':'Symbiosis',
			'Rep':'Rep','Tra':'Tra'
	}

	if 'Integrase' in feature:
		feature = 'Integrase'
	elif 'T4SS' in feature:
		feature = 'T4SS'
	elif 'T4CP' in feature:
		feature = 'T4CP'
	elif 'Relaxase' in feature:
		feature = 'Relaxase'
	elif 'Rep' in feature:
		feature = 'Rep'
	elif 'Tra' in feature:
		feature = 'Tra'
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

def gstrand1(instra):

	strands = {1 : '+', -1 : '-'}
	return strands[instra]

def getfa(infile,s,e):

	seq_record = SeqIO.read(infile, "fasta")
	sequence = seq_record.seq[int(s):int(e)]
	return str(sequence)

def get_map(run_id = str, infile = Path, listgff = List, outdir = Path, tmp_dir = Path, gb_dir = Path, js_dir = Path, jsback = Path):
	gcmap = jsback / 'gcmap.js'
	viewfile = jsback / 'view.html'
	dictICE,ICEdict,posdict,header,infodict = get_ICE(run_id, infile, listgff, outdir, tmp_dir, gb_dir)

	argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict = args(run_id)

	ICEss = {}
	for key,value in dictICE.items():
		genelist = []
		regi = run_id+'_'+key
		regijs = key
		genefile = os.path.join(outdir,regi+'_gene.json')
		infofile = os.path.join(outdir,regi+'_info.json')
		gcjson = os.path.join(js_dir,regijs+'_gc.js')
		mapfile = os.path.join(js_dir,regijs+'.js')
		htmlfile = os.path.join(outdir,regi+'.html')
		[myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,trnalist,locusdict] = value

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
					'gene':locusdict[gene],
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
					'gene':locusdict[gene],
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
					'gene':locusdict[gene],
					'pos':pos,
					'prod': product,
					'featu': feature
			}
			genelist.append(content)

		with open(genefile,'w') as gene_file:
			json.dump(genelist, gene_file, indent=4)

		sgene = zill(header,fICE)
		egene = zill(header,eICE)
		s1,e1,strand1,pro1 = posdict[sgene]
		s2,e2,strand2,pro2 = posdict[egene]
		if myDR1 == '0':
			myDR1 = '1'
		gcc = gc(infile,int(myDR1),int(myDR4))
		ICEss[regi] = '|'.join([myDR1,myDR4,str(fICE),str(eICE),gcc])

		if myDR2:
			DR1 = getfa(infile,myDR1,myDR2)
			DR2 = getfa(infile,myDR3,myDR4)
			DRw = 'attL:'+myDR1+'..'+myDR2+'('+DR1+')  '+'attR:'+myDR3+'..'+myDR4+'('+DR2+')'
		else:
			DRw = '-'
		if trnalist:
			trnaout = trnalist[0][3]+' ('+trnalist[0][0]+'..'+trnalist[0][1]+') ['+trnalist[0][2]+']'
		else:
			trnaout = '-'

		oritseqs = oritseq(run_id, regi, infile, myDR1,myDR4)
#		oritdesc = "<br>".join([oritseqs[i:i+63] for i in range(0, len(oritseqs), 63)])

		if 'IME' in regi:
			typeIE = 'IME'
		elif 'AICE' in regi:
			typeIE = 'AICE'
		else:
			typeIE = 'T4SS-type ICE'

		ICEinfo = {
			'Type':typeIE,
			'Location (nt)':myDR1+'..'+myDR4,
			'Length (bp)':str(int(myDR4)-int(myDR1)+1),
			'GC Content (%)':gcc,
			'oriT seq':oritseqs,
			'DRs':DRw,
			'Relaxase Type': ','.join(infodict[key]['mob']),
			'Mating pair formation systems':','.join(infodict[key]['mpf']),
			'Close to tRNA':trnaout
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

		gcdict = calculate_gc(infile, int(s), int(e), 500, 50)

		with open(gcmap, 'r') as original_file:
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

		with open(viewfile, 'r') as file:
			file_content = file.read()
		new_content = file_content.replace('XXXX', regijs)
		with open(htmlfile, 'w') as file:
			file.write(new_content)

	return ICEss

def get_color(region):

		coldict = {'T4SS-type ICE':'fill:rgba(0, 128, 164,0.9)',
						'IME':'fill:rgba(41,76,166,0.9)',
						'AICE':'fill:rgba(255, 84,0,0.9)'
		}
		return coldict[region]

def copy_files(source_dir, destination_dir):

	if os.path.isfile(source_dir):  # 如果源路径是文件
		shutil.copy(source_dir, destination_dir)
	elif os.path.isdir(source_dir):  # 如果源路径是文件夹
		files = os.listdir(source_dir)
		for file in files:
			source_file = os.path.join(source_dir, file)
			if os.path.isfile(source_file):
				shutil.copy(source_file, destination_dir)
			elif os.path.isdir(source_file):
				destination_subdir = os.path.join(destination_dir, file)
				copy_files(source_file, destination_subdir)

def getfasta(run_id,infile,key,s,e,stag,etag,locusdict):

	faafile = os.path.join(tmp_dir, 'gbk', run_id+'.faa')
	outfa = os.path.join(workdir,'result',run_id,key+'.fa')
	outfaa = os.path.join(workdir,'result',run_id,key+'.faa')

	seq_record = SeqIO.read(infile, "fasta")
	with open(outfa, "w") as output_handle1:
		sequence = seq_record.seq[int(s)-1:int(e)]
		ID = '_'.join(seq_record.id.split('_')[-2:])
		seq_record.description = ''
		seq_record.seq = sequence
		seq_record.id = ID + ' ' + s +'-'+e
		SeqIO.write(seq_record, output_handle1, "fasta")

	faa_records = SeqIO.parse(faafile, "fasta")
	with open(outfaa, "w") as output_handle2:
		for faa_record in faa_records:
			old_id = faa_record.id
			seq_id = getnum(old_id)
			if int(stag) <= seq_id <= int(etag):
				new_id = locusdict.get(old_id, old_id)
				faa_record.id = new_id
				SeqIO.write(faa_record, output_handle2, "fasta")

def getbase(run_id,filetype,homelist,final_dir):

	infile = os.path.join(in_dir,run_id)
	final_dir = os.path.join(workdir,'result',run_id)
	basefile = os.path.join(final_dir, run_id+'_info.json')

	if filetype == 'gb':
		filet = 'genbank'
	else:
		filet = 'fasta'
	for seq_record in SeqIO.parse(infile, filet):
		realID = seq_record.id
		realID = realID[:15]
		desc = seq_record.description
		lengt = len(seq_record.seq)
		gcs = "%.2f"%GC(seq_record.seq)

	basedict = {'JobID':run_id,
			'Submission date':get_time(),
			'Sequence name':realID,
			'Length': str(lengt)+' bp',
   			'GC Content': str(gcs) + ' %'
			}
	with open(basefile,'w') as gein2:
		json.dump(basedict, gein2, indent=4)

def _single(run_id: str,infile: Path,filetype: str, outdir: Path.resolve, tmp_dir: Path, fa_dir: Path, gb_dir: Path):

	js_dir = outdir / 'js'
	js_dir.mkdir(exist_ok=True)
	jsback = outdir / 'script' / 'js'
	if  filetype == 'fa':
		prokka(run_id,infile)
		listgff = getgff1(run_id)
	else:
		listgff = getgff(run_id, gb_dir)

	ICEss = get_map(run_id, infile, listgff, outdir, tmp_dir, gb_dir, js_dir, jsback)

	i = 1 
	ICEsumlist = []
	homelist = []
	if ICEss:
		for key,value in ICEss.items():
			[s,e,stag,etag,gc] = value.split('|')
			lengt = int(e) - int(s) + 1

			if 'IME' in key:
				typeIE = 'IME'
			elif 'AICE' in key:
				typeIE = 'AICE'
			else:
				typeIE = 'T4SS-type ICE'

			ICEs = {
				'id' : str(i),
				'region':'Region'+ str(i),
				'location': s+'..'+e,
				'length': str(lengt),
				'gc':gc,
				'type': typeIE,
				'detail': key
			}

			homedict = {
			  	'start' : int(s),
			   	'end' : int(e),
			   	'color' : get_color(typeIE),
			   	'info' : '_'.join([str(lengt),str(gc),typeIE]),
			   	'text' : 'Region'+ str(i)
			   }
			homelist.append(homedict)
			ICEsumlist.append(ICEs)
			getfasta(run_id,infile,key,s,e,stag,etag,locusdict)
			i += 1

	ICEsum = os.path.join(outdir, run_id+'_ICEsum.json')
	with open(ICEsum,'w') as ice_file:
		json.dump(ICEsumlist, ice_file, indent=4)

	copy_files(jsback, js_dir)
	getbase(run_id,filetype,homelist,outdir)

	tmpfile = os.path.join(tmp_dir,run_id)
	shutil.rmtree(tmpfile)

