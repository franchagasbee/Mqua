#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -m e
#$ -M e.karalashvili@leibniz-lib.de
#$ -N VCFstats

## usage: qsub -pe smp $CPUs -q small.q,medium.q,large.q -v "OUTPUTFOLDER=$OUTPUTFOLDER" -v "SPECIES=$SPECIES" -v "VCF=$VCF" -v "REF=$REF" ~/shells/mapSNP/2025_08_28_VCF_stats.sh

set -x

## load modules
module load vcftools/0.1.16
module load vcflib/1.0.3
module load htslib/1.19.1
module load bcftools/1.19
module load beagle/22Jul22
module load pigz/2.7
module load parallel/20201222

SPECIES=$SPECIES
OUTPUTDIR=$OUTPUTDIR
REF=$REF
VCF=$VCF
PHASINGDIR=$PHASINGDIR
CPUs=$NSLOTS

echo "filtering $VCF"
VCFPATH=$(echo $VCF | rev | cut -d"." -f3- | rev)
VCFFILE=$(echo ${VCFPATH} | rev | cut -d"/" -f1 | rev)
pigz -dkf ${VCF}

# SNPs per scaffold, filter low information scaffolds
bcftools index -s $VCF > ${OUTPUTDIR}/${VCFFILE}.indexstats
cat ${OUTPUTDIR}/${VCFFILE}.indexstats | awk ' $3 < 10 ' > ${OUTPUTDIR}/${VCFFILE}.scf.lessthan10.SNPs.tsv
cat ${OUTPUTDIR}/${VCFFILE}.indexstats | awk ' $3 < 10 ' |cut -f1 > ${OUTPUTDIR}/${VCFFILE}.scf.lessthan10.SNPs.lst
cat ${OUTPUTDIR}/${VCFFILE}.indexstats | awk ' $3 >= 10 ' > ${OUTPUTDIR}/${VCFFILE}.scf.morethan10.SNPs.tsv
cat ${OUTPUTDIR}/${VCFFILE}.indexstats | awk ' $3 >= 10 ' | cut -f1 > ${OUTPUTDIR}/${VCFFILE}.scf.morethan10.SNPs.lst
cat ${REF}.fai | cut -f1-2 | awk ' $2 >= 10000 ' | cut -f1 > ${REF}.10kb.lst
cat ${REF}.fai | sort -k2,2nr | head -n30 | cut -f1 > ${REF}.30largest.lst
cat ${OUTPUTDIR}/${VCFFILE}.scf.morethan10.SNPs.tsv | grep -f ${REF}.10kb.lst > ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.tsv
cat ${REF}.fai | wc -l > ${REF}.fai.n
cat ${REF}.10kb.lst | wc -l > ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.tsv.n
cat ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.tsv | wc -l > ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.tsv.n

cat ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.tsv | bioawk -t '{print $1,"1",$2}' > ${REF}.10kb10SNPs.regions
tabix -h --regions ${REF}.10kb10SNPs.regions $VCF | tee ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf | sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.\|.#g" | sed -r "s#\t\.:\.:\.:\.,\.:\.:\.:\.:\.:\.#\t.|.#g" | sed 's#0/0#0\|0#g;s#1/1#1\|1#g' > ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.vcf
vt peek ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.vcf 2> ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.vcf.stats

# phasing using beagle
mkdir ${PHASINGDIR}
time java -Xms1512m -Xmx40g -jar $BEAGLE/beagle.22Jul22.46e.jar gt=${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.vcf out=${PHASINGDIR}/${SPECIES}.phased impute=false ne=1000 nthreads=$CPUs window=10
tabix -fp vcf ${PHASINGDIR}/${SPECIES}.phased.vcf.gz
vt peek ${PHASINGDIR}/${SPECIES}.phased.vcf.gz 2> ${PHASINGDIR}/${SPECIES}.phased.vcf.gz.stats
cat ${PHASINGDIR}/${SPECIES}.phased.vcf.gz.stats

# vcflib vcfhetcount
vcfhetcount ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf > ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf.hetcount

################ windowed popStats ################3
WINDOW=10000

# Tajima's D per window
vcftools --vcf ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf --TajimaD ${WINDOW} --out ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.windowedTajimaD

# nucleotide diversity per window (pi)
vcftools --vcf ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf --window-pi ${WINDOW} --out ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.windowedPi

# heterozygocity per window
vcftools --vcf ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf --windowed-het ${WINDOW} --out ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.windowedHet

# SNP density per window
vcftools --vcf ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf --SNPdensity ${WINDOW} --out ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.SNPdensity

####################################
# vcflib popStats
SAMPLENR=$(echo $(tabix -H ${VCF} | tail -n1 | tr '\t' '\n' | wc -l) - 10 | bc -l) && echo ${SAMPLENR}
SAMPLESEQ=$(seq 0 ${SAMPLENR} | tr '\n' ',' | sed -e "s/,$//g") && echo ${SAMPLESEQ}
printf id"\t"pos"\t"AF"\t"He"\t"Ho"\t"nHets"\t"nHomRef"\t"nHomAlt"\t"FIS"\n" > ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf.popStats
popStats --type GL $(printf ${SAMPLESEQ}) --file ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf >> ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf.popStats

# vcflib genotypeSummary
SAMPLENR=$(echo $(tabix -H ${VCF} | tail -n1 | tr '\t' '\n' | wc -l) - 10 | bc -l) && echo ${SAMPLENR}
SAMPLESEQ=$(seq 0 ${SAMPLENR} | tr '\n' ',' | sed -e "s/,$//g") && echo ${SAMPLESEQ}
time genotypeSummary --type GL --target $(printf ${SAMPLESEQ}) --file ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf --snp ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.originalformat.vcf.genotypesummary

# vcflib sequenceDiversity
SAMPLENR=$(echo $(tabix -H ${VCF} | tail -n1 | tr '\t' '\n' | wc -l) - 10 | bc -l) && echo ${SAMPLENR}
SAMPLESEQ=$(seq 0 ${SAMPLENR} | tr '\n' ',' | sed -e "s/,$//g") && echo ${SAMPLESEQ}
printf id"\t"start"\t"end"\t"pi"\t"eHH"\n" > ${PHASINGDIR}/${SPECIES}.phased.vcf.sequenceDiversity
sequenceDiversity --type GT --target $(printf ${SAMPLESEQ}) --file ${PHASINGDIR}/${SPECIES}.phased.vcf.gz >> ${PHASINGDIR}/${SPECIES}.phased.vcf.sequenceDiversity

# vcflib iHS
printf id"\t"pos"\t"AF"\t"EHHref"\t"iHS"\t"x"\t"x"\n" > ${PHASINGDIR}/${SPECIES}.phased.vcf.header.iHS
SAMPLENR=$(echo $(tabix -H ${VCF} | tail -n1 | tr '\t' '\n' | wc -l) - 10 | bc -l) && echo ${SAMPLENR}
SAMPLESEQ=$(seq 0 ${SAMPLENR} | tr '\n' ',' | sed -e "s/,$//g") && echo ${SAMPLESEQ}
cat ${REF}.30largest.lst | parallel -j 10 "iHS --region {} --type GT --target $(printf ${SAMPLESEQ}) --file ${PHASINGDIR}/${SPECIES}.phased.vcf.gz > ${PHASINGDIR}/${SPECIES}.phased.vcf.{}.iHS"

#ROH
bcftools roh ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.vcf --threads $CPUs --AF-dflt 0.4 --GTs-only 30 --skip-indels --output ${OUTPUTDIR}/${SPECIES}.0.4.RHO
~/shells/mapSNP/2025_08_28_ROH_viz --min-length 1e6 -i ${OUTPUTDIR}/${SPECIES}.0.4.RHO -v ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.vcf -o ${OUTPUTDIR}/${SPECIES}.0.4.RHO.html
~/shells/mapSNP/2025_08_28_ROH_viz --min-length 1e6 -l 1 -i ${OUTPUTDIR}/${SPECIES}.0.4.RHO -v ${OUTPUTDIR}/${VCFFILE}.10kb10SNPs.vcf -o ${OUTPUTDIR}/${SPECIES}.0.4.l1.RHO.html

# ngsrelate
vcfsamplenames ${VCF} > ${VCF}.samples
ngsRelate -h ${VCF} -T GT -c 1 -A AF -z ${VCF}.samples -p ${CPUs} -O ${OUTPUTDIR}/${SPECIES}.ngsRelate

echo "done"
