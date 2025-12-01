#!/bin/bash

#usage: sh 03_SNP_Call_Submission.sh $REF $INDIR $BAMDIR $CPUs

set -x

#load modules
module load freebayes/1.3.2

#catch variables
REF=$1 #reference genome
INDIR=$2 #working directory
BAMDIR=$3 #directory with bam files
CPUs=$4

cd ${INDIR}

OUTPUTFOLDER="${INDIR}/snps.raw"
mkdir -p ${OUTPUTFOLDER}
mkdir -p ${OUTPUTFOLDER}/raw
ls -1 ${BAMDIR}/*.markdup.sorted.bam > ${BAMDIR}/bam.lst
BAMLIST="${BAMDIR}/bam.lst"

fasta_generate_regions.py ${REF}.fai 1000000 > ${INDIR}/regions.lst
cat "${INDIR}/regions.lst" | wc -l

REGIONS="${INDIR}/regions.lst"

qsub -pe smp ${CPUs} -q large.q,medium.q,small.q -v "BAMLIST=${BAMLIST}" -v "REGIONS=${REGIONS}" -v "REF=${REF}" -v "OUTPUTFOLDER=$OUTPUTFOLDER" 03_SNP_Call.sh
