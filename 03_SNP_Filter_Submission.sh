#!/bin/bash

#usage: sh 04_SNP_Filter_Submission.sh ${SPECIES} ${REF} ${PROJECTFOLDER} ${CPUs}

#catch variables
SPECIES=$1 #Species_name
REF=$2 #reference genome
PROJECTFOLDER=$3 #working directory
CPUs=$4

INPUTFOLDER="${PROJECTFOLDER}/snps.raw/raw"
OUTPUTFOLDER="${PROJECTFOLDER}/snps.filtered"
REGIONS="${PROJECTFOLDER}/regions.lst"

mkdir -p ${OUTPUTFOLDER}

cd ${OUTPUTFOLDER}

qsub -pe smp ${CPUs} -q small.q,medium.q,large.q -v "REF=${REF}" -v "INPUTFOLDER=${INPUTFOLDER}" -v "OUTPUTFOLDER=${OUTPUTFOLDER} -v "REGIONS=${REGIONS}" 04_SNP_Filter.sh
