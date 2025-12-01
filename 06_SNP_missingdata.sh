#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M e.karalashvili@leibniz-lib.de
#$ -m e
#$ -N SNPmissingdata

#this script filters the SNP list by missing data
## qsub -pe smp $CPUs -q small.q,medium.q,large.q,fast.q -v "VCF=$VCF" -v "INDIR=$INDIR" -v "OUTPUTDIR=$OUTPUTDIR" $HOME/scripts/cluster-SNP-missing-data-filter.sh

set -x

## require local vt installation

### load modules
module load vcftools/0.1.16
module load vcflib/1.0.3
module load htslib/1.19.1

CPUs=$NSLOTS
VCF=$VCF
INDIR=$INDIR
OUTPUTDIR=$OUTPUTDIR

echo "filtering ${VCF}"
VCFPATH=$(echo ${VCF} | rev | cut -d "." -f3- | rev)
VCFFILE=$(echo ${VCFPATH} | rev | cut -d"/" -f1 | rev)

NSAMPLES=$(vcfsamplenames ${INDIR}/${VCFFILE}.vcf.gz | wc -l) && echo ${NSAMPLES}
NMISSING25=$(echo -n ${NSAMPLES} | awk '{ $1=sprintf("%.0f",$1*0.25*2)} {print $1;}') && echo ${NMISSING25}
NMISSING=$((echo ${NSAMPLES}-${NMISSING25}/2 | bc -l) | awk '{ $1=sprintf("%.0f",$1)} {print $1;}') && echo ${NMISSING}

zcat ${INDIR}/${VCFFILE}.vcf.gz | vcftools --recode --recode-INFO-all -c --vcf - --max-missing-count ${NMISSING} |\
bgzip -f -@ ${CPUs} -c /dev/stdin > ${OUTPUTDIR}/${VCFFILE}.MaxMissing25.vcf.gz
tabix -fp vcf ${OUTPUTDIR}/${VCFFILE}.MaxMissing25.vcf.gz

vt peek ${OUTPUTDIR}/${VCFFILE}.MaxMissing25.vcf.gz 2> ${OUTPUTDIR}/${VCFFILE}.MaxMissing25.vcf.stats
cat ${OUTPUTDIR}/${VCFFILE}.MaxMissing25.vcf.stats
