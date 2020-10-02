#!/usr/bin/env python

# inputs: gene name, mutation, assembly (hg19, hg38, mm9, ...), gtf file
# outputs:
#	donor DNA (~120 bp)
#	gRNA
#	diagnostic restriction enzyme

import gzip
import os
import sys
import argparse
import pandas as pd
from Bio import Entrez, SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gene", help="name of the gene to be mutated")
parser.add_argument("-m", "--mutation", help="missense mutation to be introduced in the gene")
parser.add_argument("-G", "--genome", help="genome assembly to be used (genomes available: hg38)")
parser.add_argument("-e", "--email", help="email needed from fetching sequences from Entrez database")

args = parser.parse_args()

gene_name = args.gene
mutation = args.mutation
genome = args.genome
email = args.email

print(gene_name, mutation, genome)

class mutation:
	def __init__(self, mut_name):
		self.name = mut_name
		self.wtaa = mut_name[0]
		self.aapos = mut_name[1:-1]
		self.mutaa = mut_name[-1]
		
def file_to_list(path):
	with gzip.open(path, "r") as f:
		text = f.readlines()
		text = [line.decode("utf-8") for line in text]
		text = [line.strip() for line in text]
		text = [line.split('\t') for line in text]
	return text

def build_gene_id_table(genome):
	parent_path = os.path.dirname(os.path.realpath(__file__))
	exon_path = parent_path + "/knownGene-" + genome + ".txt.gz"
	ref_path = parent_path + "/kgXref-" + genome + ".txt.gz"
	
	exon_text = file_to_list(exon_path)
	exon_table = pd.DataFrame(exon_text, columns = ["geneid", "chr", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "proteinID", "alignID"])
	exon_table = exon_table.set_index("geneid")
	
	name_text = file_to_list(ref_path)
	name_table = pd.DataFrame(name_text, columns = ["kgID", "mRNA_ID", "spID", "spDisplayID", "geneSymbol", "refseq", "protAcc", "description"])
	name_table = name_table.set_index("kgID")
	
	ref_table = exon_table.join(name_table)
	
	return ref_table

def build_genome_chrom_dict(genome):
	parent_path = os.path.dirname(os.path.realpath(__file__))
	chrom_id_path = parent_path + "/data/chr_ids/" + genome + ".txt"
	
	with open(chrom_id_path, 'r') as f:
		text = f.readlines()
		text = [line.strip() for line in text]
		text = [line.split('\t') for line in text]
		
	chrom_dict = {}
	for line in text:
		chrom_dict[line[0]] = line[1]
	return chrom_dict

def getDNA(email, genome, chr, start, stop):
	chrom_dict = build_genome_chrom_dict(genome)
	Entrez.email = email
	handle = Entrez.efetch(db = "nucleotide",
					   id = chrom_dict[chr],
					   rettype = "fasta",
					   strand = 1,
					   seq_start = start + 1,
					   seq_stop = stop)
	record = SeqIO.read(handle, "fasta")
	handle.close()
	return record.seq

def main():
	# take inputs (gene, mutation, genome)
	# extract DNA from genome (exons)
		# load in table
		# get exons
		# concat exons
		# revcomp if on minus strand
		# translate and check wt AA is correct
		# check if mutation occurs in exon-exon junction
	# create mutant
		# get mutant codon
		# replace in wt sequence
		# split back into exons, insert introns
	# get sgRNA
		# search for closet GG
		# return sgRNA sequence
	# insert silent mutations to block sgRNA binding
	return True

# exon_table = build_gene_id_table(genome)
# print(exon_table.head(10))

# ref_table = build_gene_ref_table(genome)
# print(ref_table.head(10))

# x = exon_table.join(ref_table)
# print(x.head(10))

ref_table = build_gene_id_table(genome)

idx = ref_table.index[ref_table["geneSymbol"] == gene_name]
print(idx)

# for i in idx:
	# print(i)
	
x = idx[0]
print(x)
	
def parse_exons(tran_id, ref_table):
	start_stop_pairs = []
	starts = ref_table.at[tran_id, "exonStarts"]
	starts = starts.split(",")
	if starts[-1] == '':
		starts.pop()
		
	stops = ref_table.at[tran_id, "exonEnds"]
	stops = stops.split(",")
	if stops[-1] == '':
		stops.pop()
	if len(starts) == len(stops):
		for i in range(len(starts)):
			pair = (starts[i], stops[i])
			start_stop_pairs.append(pair)
	print(tran_id)
	print(starts)
	print(stops)
	print(start_stop_pairs)
	return start_stop_pairs

exons = parse_exons(x, ref_table)

def build_cDNA(tran_id, ref_table, exons, email, genome):
	cDNA = ''
	chr = ref_table.at[x, "chr"]
	for exon in exons:
		start = int(exon[0])
		stop = int(exon[1])
		dna = getDNA(email, genome, chr, start, stop)
		cDNA += dna
	start_idx = int(ref_table.at[tran_id, "cdsStart"])
	cds_start = start_idx - int(exons[0][0])
	stop_idx = int(ref_table.at[tran_id, "cdsEnd"])
	cds_stop = stop_idx - int(exons[-1][-1])
	cDNA = cDNA[cds_start:cds_stop]
	return cDNA

seq = build_cDNA(x, ref_table, exons, email, genome)
# start_idx = int(ref_table.at[x, "cdsStart"])
# cds_start = start_idx - int(exons[0][0])
# stop_idx = int(ref_table.at[x, "cdsEnd"])
# cds_stop = stop_idx - int(exons[-1][-1])
# seq = seq[cds_start:cds_stop]


# print(start_idx)
# print(cds_start)
# print(stop_idx)
# print(cds_stop)
print(seq)
print(len(seq))