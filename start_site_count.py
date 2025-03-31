#!/usr/bin/env python3.8
import argparse

parser =argparse.ArgumentParser(description = 'calculates the proportion of uORFs starting with AUG, CUG, GUG & UUG')
parser.add_argument('file', help='input file contining the uORFs sequences')
args = parser.parse_args()

ATG = 0
CTG = 0
GTG = 0
TTG = 0
with open(args.file) as file :
    for line in file :
        fields = line.split("\t")
        geneID = fields[0]
        seq = str(fields[3])
        start = seq[1:4]
        if start == "ATG":
            ATG = ATG+1
        elif start == "CTG":
            CTG = CTG+1
        elif start == "GTG":
            GTG = GTG+1
        elif start == "TTG":
            TTG = TTG+1
    print("AUG", ATG)
    print("CUG", CTG)
    print("GUG", GTG)
    print("UUG", TTG)
    total = ATG+CTG+GTG+TTG
    ATG_p = ATG/float(total)
    CTG_p = CTG/float(total) 
    GTG_p = GTG/float(total) 
    TTG_p = TTG/float(total)
    print("AUG %", ATG_p)
    print("CUG %", CTG_p)
    print("GUG %", GTG_p)
    print("UUG %", TTG_p)
