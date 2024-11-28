#!/bin/bash

# check commands and perhaps versions?
# samtools, stringtie

get_version () {
  local case_command=$1
  local case_version=$2
  if command -v ${case_command} >/dev/null 2>&1; then
    version_str="$( ${case_command} ${case_version} 2>&1 )"
    if [ -z "${version_str}" ]; then
      version_str="UNDEF"
    fi
  fi
  echo "${version_str}"
}


# Get project name
if [ "$#" -ne 2 ]; then
  if [[ -z "${PROJECT_NAME}" ]]; then
  	echo "Project name cannot be deduced. Suply as following:"
  	echo "make_densities INFILE PROJECT_NAME" && exit 1
  fi
else
 	PROJECT_NAME="$2"
fi

# Set project home-folder
if [[ -z "${PROJECT_HOME}" ]]; then
  PROJECT_HOME="$(pwd)"
fi

echo "Using '${PROJECT_NAME}' as project's name"
echo "Using '${PROJECT_HOME}' as project's root path"

# Check for config file and create if it doesn't exist
CONFIG_PATH="${PROJECT_HOME}/config"
CONFIG_FILE="${CONFIG_PATH}/prepare_cds.conf"
if [ ! -f "${CONFIG_FILE}" ]; then
  echo "Creating ${CONFIG_FILE} ... Check and run again."
  mkdir -p ${CONFIG_PATH}
  cat >${CONFIG_FILE} <<EOL
# Project specs
export CDS_DIR="cds_preparation"
export MAP_DIR="mapping_data/STAR_genome"
export MAP_FILE="Aligned.sortedByCoord.out.bam"
export GTF_FILE="config/Mmusculus.GRCm38.97.gtf"
export FA_FILE="config/Mmusculus.GRCm38.97.cdna.ensembl.fa"
export PROT_FILE="config/Mmusculus.GRCm38.97_prot_coding_genes.tsv"
# stringtie specs
export STRING_DIR="stringtie_data"
export STRING_PROCESS="4"
export STRING_CPU="8"
EOL
  exit 0
fi

# Check stringtie
stringtie_ver=$( get_version "stringtie" "--version" )
[ -z "${stringtie_ver}" ] && echo "Need 'stringtie'. Exiting" && exit 1

# Check samtools
samtools_ver=$( get_version "samtools" "--version" )
[ -z "${samtools_ver}" ] && echo "Need 'samtools'. Exiting" && exit 1

echo "Found 'stringtie': ${stringtie_ver}"
echo "Found 'samtools': ${samtools_ver}"

# Set variables
source ${CONFIG_FILE}
export CDS_DIR="${PROJECT_HOME}/${CDS_DIR}"
export MAP_DIR="${PROJECT_HOME}/${MAP_DIR}"
export GTF_FILE="${PROJECT_HOME}/${GTF_FILE}"
export STRING_DIR="${PROJECT_HOME}/${STRING_DIR}"
export FILTERED="${CDS_DIR}/${PROJECT_NAME}_filtered_genes_transcripts.tsv"
LIMIT=${STRING_PROCESS}
PROT_CODING="${PROJECT_HOME}/${PROT_FILE}"
FASTA="${PROJECT_HOME}/${FA_FILE}"
TR_FILE="${CDS_DIR}/${PROJECT_NAME}_filtered_transcripts.tsv"
GENE_FILE="${CDS_DIR}/${PROJECT_NAME}_filtered_genes.tsv"
PROT_FILE="${CDS_DIR}/${PROJECT_NAME}_filtered_prot_coding_genes.tsv"
PREPARED_CDS="${CDS_DIR}/${PROJECT_NAME}_prepared_cds.tsv"
PREPARED_LOG="${CDS_DIR}/${PROJECT_NAME}_prepared_cds_info.txt"
ANNOT_FILE="${CDS_DIR}/${PROJECT_NAME}_gene_annotation.tsv"
P_GTF="${CDS_DIR}/${PROJECT_NAME}_filtered.gtf"
IGV_GTF="${CDS_DIR}/${PROJECT_NAME}_igv.gtf"
IGV_FASTA="${CDS_DIR}/${PROJECT_NAME}_igv.fa"
LOCI="${CDS_DIR}/${PROJECT_NAME}_loci_IDs.txt"

mkdir -p ${STRING_DIR}
mkdir -p ${CDS_DIR}
touch ${FILTERED}


echo "Transcript analysis with StringTie"
cat $1 | xargs -n 1 -P $LIMIT -I INFILE bash -c \
'echo "... processing INFILE"; '\
'stringtie ${MAP_DIR}/INFILE/${MAP_FILE} '\
'-o ${STRING_DIR}/INFILE/out.gtf -p ${STRING_PROCESS} -G ${GTF_FILE} '\
'-A ${STRING_DIR}/gene_abund.tab -C ${STRING_DIR}/cov_refs.gtf -B -e '\
'&& filter_tdata ${STRING_DIR}/INFILE/t_data.ctab >> ${FILTERED} '\
'&& echo "... finished processing of INFILE"'

echo "Sorting and removing redundancy ..."
sort ${FILTERED} | uniq --repeated > _temp_sort_file && mv _temp_sort_file ${FILTERED}

echo "Splitting Gene and Transcript IDs ..."
cut -f2 ${FILTERED} | sort > ${TR_FILE}
cut -f1 ${FILTERED} | sort | uniq > ${GENE_FILE}

echo "Filtering protein coding genes ..."
filter_by_ids --not-skip-firstline ${TR_FILE} ${PROT_CODING} > ${PROT_FILE}

echo "Preparing CDS models for filtered protein coding transcipts"
prepare_cds_models --list ${PROT_FILE} ${GTF_FILE} >${PREPARED_CDS} 2>${PREPARED_LOG}

echo "Preparing project GTF"
filter_gtf_by_ids --gene-id --tr-id --filter-file ${FILTERED} ${GTF_FILE} > ${P_GTF}

echo "Preparing project annotations"
filter_summarize_gtf --filter-file ${FILTERED} --gene-id --tr-id \
  -o gene_id,gene_name,gene_biotype,transcript_id,transcript_name,transcript_biotype \
  ${P_GTF} | sort | uniq > ${ANNOT_FILE}

echo "Preparing GTF and FASTA for IGV"
make_gtf_for_igv ${ANNOT_FILE} ${PREPARED_CDS} > ${IGV_GTF}
awk 'BEGIN {OFS="|"} {print $1,$2}' ${FILTERED} > ${LOCI}
filter_fasta_by_ids --filter-file ${LOCI} ${FASTA} > ${IGV_FASTA}
samtools faidx ${IGV_FASTA}

echo "Script completed."
exit 0
