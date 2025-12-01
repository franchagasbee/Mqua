#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M e.karalashvili@leibniz-lib.de
#$ -m e
#$ -N SNP_filterrepeats

#this script filters the SNP list by repeatitive regions using the merged repeat file
##usage: qsub -pe smp $CPUs -q small.q,medium.q,large.q -v "VCF=$VCF" -v "REPEATS=$REPEATS" -v "VCFOUT=$VCFOUT" ~/shells/mapSNP/2025_07_09_SNP_repeat_filter.sh

module load bedtools/2.29.2
module load htslib/1.19.1

## catch variables
VCF=$VCF
REPEATS=$REPEATS
VCFOUT=$VCFOUT
CPUs=$NSLOTS

echo "filtering SSRs and Ns"
bedtools intersect -a ${VCF} -b ${REPEATS} -v -header > ${VCFOUT}
vt peek ${VCFOUT} 2> ${VCFOUT}.stats
bgzip -f -@ $CPUs ${VCFOUT}
tabix -p vcf ${VCFOUT}.gz
