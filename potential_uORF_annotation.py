#!/usr/bin/env python3.8
import argparse
import re

parser =argparse.ArgumentParser(description = 'create file with gene ID,start and in frame stop codon positions of potential uORFs (no overlapping uORFs !)')
parser.add_argument('file', help='input file')
#parser.add_argument("-o", "--output")
args = parser.parse_args()

start = ("ATG","CTG","GTG","TTG")
stop = ("TAG","TGA","TAA")
with open(args.file) as file :
    for line in file :
        fields = line.split("\t")
        geneID = fields[0]
        GeneID = geneID.rstrip()
        transcriptID = fields[1]
        transcriptID1 = transcriptID.rstrip()
        TranscriptID = transcriptID1.lstrip()
        seq = fields[2]
        IDs = GeneID+'|'+TranscriptID
        for startcodons in start:
            for match1 in re.finditer(startcodons,seq):
                startPos = match1.start()
                afterStart = seq[startPos:]
                beforeStart = seq[:startPos]
                stops = dict()
                for stopcodons in stop:
                    for match2 in re.finditer(stopcodons,afterStart):
                        stopPos = match2.start() + len(beforeStart)
                        inframe = (stopPos - startPos) % 3
                        if inframe == 0:
                            stops[stopcodons]=stopPos
                            break
                if len(stops) != 0:
                    stopsorted = dict(sorted(stops.items(), key=lambda items:items[1]))
                    minstopPos = list(stopsorted.values())[0]
                    sequence = seq[startPos:minstopPos+3]
                    if startPos > 15 :
                        startPos2 = startPos - 16
                        minstopPos2 = minstopPos - 16
                        print (IDs,"\t",startPos2,"\t",minstopPos2,"\t","protein_coding|protein_coding",sep='')
#        output_file.write()

