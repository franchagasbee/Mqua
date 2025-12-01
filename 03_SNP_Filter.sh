#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -m e
#$ -M yourname@email.com
#$ -N SNP_Filter

#this script filters SNPs by coverage/allele frequencies 
#usage: qsub -pe smp ${CPUs} -q large.q,medium.q,small.q -v "REF=${REF}" -v "INPUTFOLDER=${INPUTFOLDER}" -v "OUTPUTFOLDER=${OUTPUTFOLDER}" -v "REGIONS=${REGIONS} 04_SNP_Filter.sh

set -x

#load modules
module load htslib/1.19.1
module load samtools/1.19.2
module load bcftools/1.19
module load freebayes/1.3.2
module load vcftools/0.1.16
module load mambaforge/23.1.0
module load parallel/20201222
module load vcflib/1.0.3

REF=${REF} #reference genome
INPUTFOLDER=${INPUTFOLDER} #directory with vcf files
OUTPUTFOLDER=${OUTPUTFOLDER} #directory where output should be written
REGIONS=${REGIONS} #regions of the reference genome to parallelize over
CPUs=$NSLOTS

MINALTFRAC=0.35
MINALTN=4
MINCOV=6

cat ${REGIONS} | parallel -k -j ${CPUs} "echo {} && \
bcftools filter -e 'INFO/DP<10 | QUAL<30' ${INPUTFOLDER}/{}.vcf.gz | bcftools norm -m -any - -o - -Oz | \
bcftools norm --rm-dup all --fasta-ref ${REF} --old-rec-tag PRENORMALIZE - | \
bcftools view -m2 -M2 -v snps -i'AC=2 | INFO/AF > 0.01 | FORMAT/DP>2' | \
vcffixup - | vcfnulldotslashdot | vcfstreamsort -a | sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.#\t./.#g" | \
bcftools view -v snps --write-index --output-type z5 --output ${OUTPUTFOLDER}/{}.vcf.gz"

echo "SNP filter done"
