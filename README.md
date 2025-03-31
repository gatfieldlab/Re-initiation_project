This repository contains the scripts used in the publication "MCTS2 and distinct eIF2D roles in uORF-dependent translation regulation revealed by in vitro re-initiation assays" (https://doi.org/10.1038/s44318-024-00347-3) from the Gatfield lab at the University of Lausanne.

1. Read counts
Rpf-counts :
This script was used to count the number of CDS or 5’ UTR reads per gene using bed files and an annotation file containing the ID, start and stop position of the features.

2. Counts normalisation and translation efficiency calculation
EdgeR_TMM_normalisation_TE.R :
R script to normalise RNA and RPF read counts using the EdgeR TMM normalisation method and calculate the CDS translation efficiency.

3. 5’ UTR relative occupency
5UTR_ribosome_relative_occupancy.R :
Allows to calculate and plot the relative ribosomal occupency of the 5’ UTR of transcripts with sufficient 5’ UTR coverage.

4. Differential translation analysis
Delta_TE_analysis.R :
This script is based on the delta TE method (Chothani et al., 2019) used to calculate the differential translation efficiency of transcripts. Takes files with gene counts and a design file summarizing the different conditions and batches of the experiment.

5. uORFs annotation
prepare_cds.sh :
A BAM file is used to create a CSV file containing the list of expressed transcripts in a sample. This information was then used to annotate uORFs of representative expressed transcripts.
isoforms_same_5UTR.py :
From a list of transcripts and their 5’UTR end position, creates a list of genes for which all isoforms have the same 5’ UTR end position.
merge_sequence.py :
This script merges a file with gene ID, start and stop positions of an RNA feature (e.g. uORF or 5’ UTR) with a file that contains the gene IDs and the transcript sequences.
potential_uORF_annotation.py :
Requires a file that countains the 5’ UTR sequence of transcripts. Annotates all potential non-overlapping uORFs with a start codon and an in frame stop codon.
uORF_annotations.R :
This code was used to annotate uORFs and plot different features of DENR- and eIF2D-dependent transcripts with uORFs.
start_site_count.py :
This script measures the proportion uORFs starting with of each type of start codons (AUG, CUG, GUG & UUG).
