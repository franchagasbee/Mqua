#!/bin/bash

#usage: sh 02_Read_Mapping_Submission.sh $SPECIES $REF $READFOLDER $OUTPUTDIR $SAMPLELIST $CPUs.

set -x

SPECIES=$1 #Species_name
REF=$2 #Reference genome
READFOLDER=$3 #Directory with filtered paired end reads
OUTPUTDIR=$4 #Directory where output should be written
SAMPLELIST=$5 #List of samples, one per line
CPUs=$6

cd ${OUTPUTDIR}

N=$(wc -l ${SAMPLELIST} | cut -d' ' -f 1) && echo $N
for (( i = 1 ; i < $N+1; i++))
    do
        SAMPLE=$(cat ${SAMPLELIST} | cut -f1 | sed -n $i'p')
        echo ${SAMPLE}
        FWD=$(echo ${READFOLDER}/${SAMPLE}_1.dedupe.filtered.fq.gz)
        ls $FWD
        REV=$(echo ${READFOLDER}/${SAMPLE}_2.dedupe.filtered.fq.gz)
        ls $REV
        echo "submitting mapping for ${SAMPLE}"
        qsub -pe smp $CPUs -q large.q,medium.q,small.q -v "SAMPLE=${SAMPLE}" -v "FWD=${FWD}" -v "REV=${REV}" -v "REF=${REF}" -v "OUTPUTDIR=${OUTPUTDIR}" 02_Read_Mapping.sh
    done
