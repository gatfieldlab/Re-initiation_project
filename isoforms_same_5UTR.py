#!/usr/bin/env python3.8

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description = 'select genes for which all isoforms have the same 5UTRend position')
parser.add_argument('file', help = 'File with Gene ID, transcript ID and 5UTR end position, seperated by \t')
args = parser.parse_args()

with open(args.file) as file :
    previousID = 1
    previousEnd = 1
    previousTranscriptID = 1
    for lines in file :
        field = lines.split("\t")
        GeneID = field[0]
        TranscriptID = field[1]
        End = int(field[2])
        if GeneID == previousID :
            if 0 <= End - previousEnd <= 5 or 0 >= End - previousEnd >= -5 :
                if previousEnd >=1 End :
                    if geneID == previousPreviousID && previousPreviousEnd >= End :
                        print(previousPreviousID,"\t",previousPreviousTranscriptID,"\t",previousPreviousEnd)
                    else :
                        print(GeneID,"\t",TranscriptID,"\t",End)
                else :
                    print(GeneID,"\t",previousTranscriptID,"\t",previousEnd)
                previousTranscriptID = TranscriptID
                previousEnd = End
                previousID = GeneID
        else :
            previousID = GeneID
            previousEnd = End
            previousTranscriptID = TranscriptID
