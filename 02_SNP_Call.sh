#!/bin/bash
#
#$ -S /bin/bash
#$ -j n
#$ -M yourname@email.com
#$ -m e
#$ -N SNP_call

#this script performs a SNP call using freebayes, parallelizes over regions
#usage: qsub -pe smp $CPUs -q large.q,medium.q,small.q -v "REGIONS=${REGIONS}" -v "BAMLIST=${BAMLIST}" -v "REF=${REF}" -v "OUTPUTFOLDER=${OUTPUTFOLDER} 03_SNP_call.sh

set -x

#load modules
module load htslib/1.19.1
module load samtools/1.19.2
module load bcftools/1.19
module load freebayes/1.3.2
module load vcftools/0.1.16
module load mambaforge/23.1.0
module load parallel/20201222

#catch variables
CPUs=${NSLOTS}
REGIONS=${REGIONS} #regions of reference genome to paralellize over
BAMLIST=${BAMLIST} #list of bam files to jointly call SNPs
REF=${REF} #reference genome
OUTPUTFOLDER=$OUTPUTFOLDER #directory where output should be written

MINALTFRAC=0.35 #threshold ALT allele (fraction)
MINALTN=4 #threshold ALT reads
MINCOV=6 #threshold coverage

#run SNP call
cat ${REGIONS} | parallel -k -j $CPUs "echo {} && freebayes --region {} \
--fasta-reference ${REF} \
--ploidy 2 \
--report-genotype-likelihood-max \
--use-mapping-quality \
--genotype-qualities \
--use-best-n-alleles 3 \
--haplotype-length 1 \
--min-mapping-quality 30 \
--min-base-quality 30 \
--min-alternate-fraction ${MINALTFRAC} \
--min-alternate-total ${MINALTN} \
--min-coverage ${MINCOV} \
--use-reference-allele \
--bam-list ${BAMLIST} |\
bgzip -f -@ 1 -c /dev/stdin > ${OUTPUTFOLDER}/raw/{}.vcf.gz \
&& tabix -fp vcf ${OUTPUTFOLDER}/raw/{}.vcf.gz"

echo "SNP call finished"
