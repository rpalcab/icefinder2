#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,time
import json
import shutil
import logging
from pathlib import Path
from typing import List, Union, Tuple, Dict, Any, Sequence
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils import GC
from script.config import get_param
from script.utils import copy_files, zill, get_gff, get_num, find_max_distance, candidate_setup
from script.process import kraken2, prodigal, seqkit, hmmscan, process_hmmscan, prokka, ICEscan, get_dr, blastn_all, blastp_all, defense_finder

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

def prescan(run_id: str, input: Path, outdir: Path) -> List:
	anno_fa = prodigal(run_id, input, outdir)
	scan_file = hmmscan(run_id, input, outdir, anno_fa)
	candidates = process_hmmscan(scan_file)
	return candidates

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

def pos_tag(position: int, pos_dict: Dict[str, Tuple[str, str]], current_index: int, current_final: int, total_genes: int, direction: str) -> Tuple[int, int]:
	"""
	Adjust an ICE gene index and its boundary based on a genomic position.

	- `direction` 's' updates start side; 'e' updates end side.
	"""
	# Iterate through sorted entries to find the first covering interval
	for key, (start_s, end_s) in pos_dict.items():
		start_i, end_i = int(start_s), int(end_s)
		if position <= end_i:
			if direction == 's':
				new_index = get_num(key)
				new_final = max(1, new_index - 5)
			else:
				# Move index back one if the DR start is before this interval
				idx = get_num(key)
				new_index = idx - 1 if start_i > position else idx
				new_final = min(total_genes, new_index + 5)
			return new_index, new_final

	# If not found, return original
	return current_index, current_final

# def merge_tRNA(run_id: str, ice_dict: dict, dr_df: pd.DataFrame, prokka_dir: Path) -> Tuple[List, Dict, str, List]:
# 	# Load annotations
# 	gff_path = prokka_dir / f"{run_id}.gff"
# 	trna_dict, pos_dict, header, total_genes = get_gff(gff_path)

# 	fICE = get_num(next(iter(ice_dict)))
# 	eICE = get_num(list(ice_dict.keys())[-1])

# 	nfICEnum = max(1, fICE - 5)
# 	neICEnum = min(total_genes, eICE + 5)

# 	ICEtagnum = [nfICEnum,neICEnum]
# 	trnalist = []
# 	for key,value in trna_dict.items():
# 		if nfICEnum <= get_num(key) <= neICEnum:
# 			ICEtagnum.append(get_num(key))
# 			trnalist.append(value)

# 	ICEtagnum.sort()
# 	finalstart,finalend = find_max_distance(ICEtagnum)

# 	myDR1 = pos_dict[zill(header,fICE)][0]
# 	myDR2 = ''
# 	myDR3 = ''
# 	myDR4 = pos_dict[zill(header,eICE)][1]

# 	if trnalist:
# 		if finalstart == nfICEnum:
# 			eICE = finalend
# 			finalend = min(total_genes, finalend + 5)
# 			myDR4 = pos_dict[zill(header,eICE)][1]					
# 			for line in dr_list:
# 				DRs = line.split('|')
# 				if int(DRs[3]) - int(DRs[0]) > 500000:
# 					continue
# 				if int(DRs[3]) - int(DRs[0]) < 5000:
# 					continue					
# 				if int(pos_dict[zill(header,eICE)][0]) < int(DRs[3]) < int(pos_dict[zill(header,eICE)][1]):
# 					checktrna = 0
# 					for key,value in trna_dict.items():
# 						if int(DRs[0]) <= value[0] <= int(DRs[3]) and int(DRs[0]) <= value[1] <= int(DRs[3]):
# 							checktrna += 1
# 					if checktrna >= 2:
# 						break

# 					fICE,finalstart = pos_tag(DRs[0],pos_dict,fICE,finalstart,total_genes,'s')
# 					myDR1 = DRs[0]
# 					myDR2 = DRs[1]
# 					myDR3 = DRs[2]
# 					myDR4 = DRs[3]
# 					break

# 		elif finalend == neICEnum:
# 			fICE = finalstart
# 			finalstart =  max(1, finalstart - 5)
# 			myDR1 = pos_dict[zill(header,fICE)][0]
# 			for line in dr_list:
# 				DRs = line.split('|')
# 				if int(DRs[3]) - int(DRs[0]) > 500000:
# 					continue	
# 				if int(DRs[3]) - int(DRs[0]) < 5000:
# 					continue									
# 				if int(pos_dict[zill(header,fICE)][0]) < int(DRs[0]) < int(pos_dict[zill(header,fICE)][1]):
# 					checktrna = 0
# 					for key,value in trna_dict.items():
# 						if int(DRs[0]) <= value[0] <= int(DRs[3]) and int(DRs[0]) <= value[1] <= int(DRs[3]):
# 							checktrna += 1
# 					if checktrna >= 2:
# 						break
# 					eICE,finalend = pos_tag(DRs[3],pos_dict,eICE,finalend,total_genes,'e')
# 					myDR1 = DRs[0]
# 					myDR2 = DRs[1]
# 					myDR3 = DRs[2]
# 					myDR4 = DRs[3]									
# 					break

# 	return [myDR1, myDR2, myDR3, myDR4, fICE, eICE, finalstart, finalend], pos_dict, header, trnalist

def load_and_collect(
    run_id: str,
    ice_dict: Dict[str, Any],
    prokka_dir: Path
) -> Tuple[Dict[str, Sequence[str]], str, int, int, int, int, List[Tuple[int,int]], List[int]]:
    """
    Load GFF, extract ICE gene window, and gather tRNA positions.

    Returns:
      - trna_data, pos_data, header, total_genes,
      - start_win, end_win,
      - tRNA_ranges, ice_positions
    """
    gff_path = prokka_dir / f"{run_id}.gff"
    trna_data, pos_data, header, total_genes = get_gff(gff_path)

    ice_nums = sorted(get_num(k) for k in ice_dict)
    first_ice, last_ice = ice_nums[0], ice_nums[-1]
    start_win = max(1, first_ice - 5)
    end_win = min(total_genes, last_ice + 5)

    tRNA_ranges: List[Tuple[int,int]] = []
    ice_positions = [start_win, end_win]
    for key, coords in trna_data.items():
        idx = get_num(key)
        if start_win <= idx <= end_win:
            ice_positions.append(idx)
            try:
                s, e = int(coords[0]), int(coords[1])
                tRNA_ranges.append((s, e))
            except Exception:
                continue
    ice_positions.sort()
    return trna_data, pos_data, header, total_genes, start_win, end_win, tRNA_ranges, ice_positions


def adjust_window(
    start_win: int,
    end_win: int,
    ice_positions: List[int]
) -> Tuple[int, int]:
    """
    Adjust the ICE gene window by the largest gap in tRNA positions.
    """
    gap = find_max_distance(ice_positions)
    if gap:
        return gap
    return start_win, end_win


def refine_dr_boundaries(
    dr_df: pd.DataFrame,
    pos_data: Dict[str, Sequence[str]],
    header: str,
    total_genes: int,
    tRNA_ranges: List[Tuple[int,int]],
    start_win: int,
    end_win: int,
    final_start: int,
    final_end: int
) -> Tuple[str, str, str, str, int, int, int, int]:
    """
    Find and refine direct repeat (DR) boundaries.
    """
    dr1 = pos_data[zill(header, start_win)][0]
    dr4 = pos_data[zill(header, end_win)][1]
    dr2 = dr3 = ''

    extend_start = (final_start == start_win)
    extend_end = (final_end == end_win)
    if not (extend_start or extend_end):
        return dr1, dr2, dr3, dr4, start_win, end_win, final_start, final_end

    if extend_start:
        end_win = min(total_genes, final_end + 5)
        gene_key = zill(header, end_win)
    else:
        start_win = max(1, final_start - 5)
        gene_key = zill(header, start_win)

    gene_s, gene_e = map(int, pos_data[gene_key][:2])

    for _, row in dr_df.iterrows():
        a, b, c, d = map(int, (row['start1'], row['end1'], row['start2'], row['end2']))
        length = d - a
        if not 5000 <= length <= 500000:
            continue
        coord = d if extend_start else a
        if not (gene_s < coord < gene_e):
            continue
        inside = sum(1 for ts, te in tRNA_ranges if a <= ts <= d and a <= te <= d)
        if inside >= 2:
            continue
        if extend_start:
            start_win, final_start = pos_tag(a, pos_data, start_win, final_start, total_genes, 's')
        else:
            end_win, final_end = pos_tag(d, pos_data, end_win, final_end, total_genes, 'e')
        dr1, dr2, dr3, dr4 = map(str, (a, b, c, d))
        break

    return [dr1, dr2, dr3, dr4, start_win, end_win, final_start, final_end]


def merge_trna(
    run_id: str,
    ice_dict: Dict[str, Any],
    dr_df: pd.DataFrame,
    prokka_dir: Path
) -> Tuple[
    str, str, str, str,
    int, int, int, int,
    Dict[str, Sequence[str]], str,
    List[Tuple[int, int]]
]:
    trna_data, pos_data, header, total_genes, start_win, end_win, tRNA_ranges, ice_positions = load_and_collect(run_id, ice_dict, prokka_dir)

    final_start, final_end = adjust_window(start_win, end_win, ice_positions)

    dr_data = refine_dr_boundaries(
        dr_df, pos_data, header, total_genes,
        tRNA_ranges, start_win, end_win, final_start, final_end
    )
    return (dr_data, pos_data, header, tRNA_ranges)


def ice_markers(contig_id: str, prokka_dir: Path, ice_dict: dict, dr_df: pd.DataFrame):

	coords_dict = {}
	for key, value in ice_dict.items():
		coords_dict[key], pos_dict, header, trnalist = merge_trna(contig_id, value, dr_df, prokka_dir)
	return coords_dict, ice_dict, pos_dict, header, trnalist

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

def charact_contig(contig_id: str, id_dict, fasta_file: Path, outdir: Path):

	# Output paths and files
	# ice_dir = outdir / f'{contig_id}_ICE'
	# shutil.rmtree(ice_dir,ignore_errors=True)
	# ice_dir.mkdir()

	# Contig annotation
	prokka_dir = prokka(contig_id, fasta_file, outdir)
	ice_dict, info_dict = ICEscan(contig_id, prokka_dir, outdir)
	dr_df = get_dr(contig_id, fasta_file, outdir)

	########################################
	################# AQUI #################
	########################################
	# Calculate ICE limits
	coords_dict, ann_dict, pos_dict, header, trnalist = ice_markers(contig_id, prokka_dir, ice_dict, dr_df)

	# Detailed genes of interest search
	outblast = outdir / 'blast'
	outblast.mkdir(parents=True, exist_ok=True)
	argdict = blastn_all(prokka_dir / f'{contig_id}.ffn', outdir / 'blast')
	isdict, vfdict, metaldict, popdict, symdict = blastp_all(prokka_dir / f'{contig_id}.faa', outdir / 'blast')
	dfdict = defense_finder(contig_id, prokka_dir / f'{contig_id}.faa', outdir / 'defense_finder')

	ICEss = {}
	for key,value in coords_dict.items():
		genelist = []
		regi = contig_id+'_'+key
		regijs = 'contig_'+contig_id.split("_contig_", 1)[-1] +'_'+key
		genefile = os.path.join(outdir,regi+'_gene.json')
		infofile = os.path.join(outdir,regi+'_info.json')
		js_dir = outdir / 'js'
		js_dir.mkdir(exist_ok=True)
		gc_map = outdir / 'script' / 'js' / 'gcmap.js'
		view_file = outdir / 'script' / 'js' / 'view.html'
		gcjson = os.path.join(js_dir,regijs+'_gc.js')
		mapfile = os.path.join(js_dir,regijs+'.js')
		htmlfile = os.path.join(outdir,regi+'.html')
		[myDR1,myDR2,myDR3,myDR4,first_gene,last_gene,finalstart,finalend] = value

		start = finalstart
		while start < first_gene:
			gene = zill(header,start)
			s,e,strand,pro = pos_dict[gene]
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

		mov = first_gene
		while mov <= last_gene:
			gene = zill(header,mov)
			s,e,strand,pro = pos_dict[gene]
			pos = s+'..'+e+' ['+strand+'], '+str(int(e)-int(s)+1)

			if gene in ann_dict[key]:
				[feature,pro11] = ann_dict[key][gene].split('@')
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
			s,e,strand,pro = pos_dict[gene]
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

		sgene = zill(header,first_gene)
		egene = zill(header,last_gene)
		s1,e1,strand1,pro1 = pos_dict[sgene]
		s2,e2,strand2,pro2 = pos_dict[egene]
		if myDR1 == '0':
			myDR1 = '1'

		ICEss[regi] = '|'.join([myDR1,myDR4,str(first_gene),str(last_gene)])

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
			'Relaxase Type': ','.join(info_dict[key]['mob']),
			'Mating pair formation systems':','.join(info_dict[key]['mpf']),
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
			seq_id = get_num(faa_record.id)
			if int(stag) <= seq_id <= int(etag):
				SeqIO.write(faa_record, output_handle2, "fasta")

def _meta(run_id: str, input: Path, outdir: Path, kraken_db: Union[Path, None], gff: Union[Path, None]) -> None:

	jsback = outdir / 'script' / 'js'

	# En cuarentena, quizá borrar?
	id_dict = rename(run_id,input)

	# ICE markers annotations
	candidates = prescan(run_id, input, outdir)

	# Taxonomy assignation (optional)
	spdict = {}
	report = None
	if kraken_db:
		spdict,report = kraken2(run_id, input, outdir, kraken_db)
		# copy_files(report, outdir)

	# Get basic information from input
	basefile = seqkit(run_id, input, outdir)
	# copy_files(basefile, outdir)
	
	i = 1 
	ice_summary = []
	for seq_record in SeqIO.parse(input, "fasta"):
		if seq_record.id in candidates:
			logging.info(f"Analyzing contig {seq_record.id}")

			# Extract contig of interest in fasta format
			contig_folder, contig_fasta = candidate_setup(seq_record, outdir)

			# Characterize contig
			ICEss = charact_contig(seq_record.id, id_dict, contig_fasta, contig_folder)
			# copy_files(contig_folder, outdir)

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
					ice_summary.append(ICEs)
					getfasta(seq_record.id,outdir,id_dict,key,s,e,stag,etag)
					i += 1 

	ICEsum = outdir / tmp_dir / f'{run_id}_ICEsum.json'
	with open(ICEsum,'w') as ice_file:
		json.dump(ice_summary, ice_file, indent=4)

	copy_files(ICEsum, outdir)
	jsdir = os.path.join(outdir,'js')
	copy_files(jsback, jsdir)
	delete_folders_starting_with_keyword(tmp_dir, run_id)
