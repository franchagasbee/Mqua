#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M yourname@email.com
#$ -m ae
#$ -N mapping

#this script maps filtered reads to the reference genome, marks duplicates
##usage: qsub -pe smp $CPUs -q small.q,medium.q,large,q -v "SAMPLE=$SAMPLE" -v "FWD=$FWD" -v "REV=$REV" -v "OUTPUTDIR=$OUTPUTDIR" -v "REF=$REF" 02_Read_Mapping.sh

set -x

#load modules
module load samtools/1.19.2
module load bwa-mem2/2.2.1
modula load java/jdk-11.0.7

#catch variables
CPUs=${NSLOTS}
SAMPLE=${SAMPLE} #sample ID
FWD=${FWD} #FWD reads
REV=${REV} #REV reads
OUTPUTDIR=${OUTPUTDIR} #Directory where output should be written
REF=${REF} #Reference genome

#READgroups, based on Input Sample ID
SAMPLEID=${SAMPLE}
ID=${SAMPLE}
CLADE=${SAMPLE}
LIBRARY="NEBnext"
PLATFORM="IlluminaNovaseq6k"
TIMEOFSEQ="2025"
READGROUPHEADER="@RG\tID:$ID\tSM:$SAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY\tPU:$CLADE\tDT:$TIMEOFSEQ"

#map reads
time bwa2 mem -t ${CPUs} -S ${REF} ${FWD} ${REV} -R ${READGROUPHEADER} \
| samtools view --threads ${CPUs} --reference ${REF} -b -u - \
| samtools sort -l 5 -m 1024M -T ${HOME}/tmp/samtools.${SAMPLE}.tmp --threads ${CPUs} --output-fmt BAM -o ${OUTPUTDIR}/${SAMPLE}.bam -
samtools index ${CPUs} ${OUTPUTDIR}/${SAMPLE}.bam
samtools flagstats -@ ${CPUs} -O tsv ${OUTPUTDIR}/${SAMPLE}.bam > ${OUTPUTDIR}/${SAMPLE}.bam.flagstat
echo "${OUTPUTDIR}/${SAMPLE}.bam finished"

time samtools sort -n -l 0 -m 1024M -T ${HOME}/tmp/samtools.${SAMPLE}.tmp --threads ${CPUs} --output-fmt BAM ${OUTPUTDIR}/${SAMPLE}.bam \
| samtools fixmate -m --threads ${CPUs} --reference ${REF} --output-fmt BAM - - \
| smtools sort -l 8 -m 1024M --threads ${CPUs} --output-fmt BAM -T ${HOME}/tmp/samtools.2.${SAMPLE}.tmp - \
| samtools markdup -S -s --mode t --output-fmt BAM --reference ${REF} --threads ${CPUs} --write-index -f ${OUTPUTDIR}/${SAMPLE}.markdup.bam.stats - ${OUTPUTDIR}/${SAMPLE}.markdup.bam
samtools index -@ ${CPUs} ${OUTPUTDIR}/${SAMPLE}.markdup.bam
samtools flagstat -@ ${CPUs} -O tsv ${OUTPUTDIR}/${SAMPLE}.markdup.bam > ${OUTPUTDIR}/${SAMPLE}.markdup.bam.flagstat
echo "${OUTPUTDIR}/${SAMPLE}.markdup.bam finished"

unset DISPLAY #to avoid java issues
time qualimap bamqc -bam ${OUTPUTDIR}/${SAMPLE}.markdup.bam -nw 4000 -nr 10000 -os -hm 6 -nt ${CPUs} -outformat PDF:HTML -outdir ${OUTPUTDIR}/${SAMPLE}.markdup.bam.bamqc -outfile ${SAMPLE}.markdup.bam.qualimap.pdf
echo "qualimap finished for ${SAMPLE}"
