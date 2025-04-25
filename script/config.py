#!/usr/bin/env python
# -*- coding: utf-8 -*-

import configparser
import os

conf = configparser.ConfigParser()
conf.read("./config.ini")

def get_param():
	options = conf.options("Param")
	workdir = conf.get("Param", "workdir")
	kraken = conf.get("Param", "kraken")
	krakenDB = conf.get("Param", "krakenDB")
#	mmseqs = conf.get("Param", "mmseqs")
#	mmseqDB = conf.get("Param", "mmseqDB")
	defensefinder = conf.get("Param", "defensefinder")
	blastp = conf.get("Param", "blastp")
	blastn = conf.get("Param", "blastn")
	seqkit = conf.get("Param", "seqkit")
	prodigal = conf.get("Param", "prodigal")
	prokka = conf.get("Param", "prokka")
	macsyfinder = conf.get("Param", "macsyfinder")
	hmmsearch = conf.get("Param", "hmmsearch")	

	return workdir,kraken,krakenDB,defensefinder,blastp,blastn,seqkit,prodigal,prokka,macsyfinder,hmmsearch
