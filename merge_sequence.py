#!/usr/bin/env python3.8

## This script is used to select subset of sequences from fasta file
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description = 'Add sequences to data frame that contains gene ID, start and stop codon positions')
parser.add_argument('-f','--file', help = 'file with Gene ID, Transcript ID and start and stop position')
args = parser.parse_args()

Genes_seq = dict()
for record in SeqIO.parse("config/Mmusculus.GRCm38.100.cdna.ensembl.fa", "fasta") :
    IDs = record.id 
    seq = record.seq
    Genes_seq[IDs]=seq

with open(args.file) as pos_file :
    for line in pos_file :
        fields = line.split(" ")
        fileID = fields[0]
        start = int(fields[1])
        stop = int(fields[2])
        count = fields[3]
        if fileID in Genes_seq :
            geneSeq = Genes_seq[fileID]
            uORFseq = geneSeq[start+15:stop+15] 
            print(fileID,"\t",start,"\t",stop,"\t",uORFseq,"\t",count, end="")
