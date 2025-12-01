#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M e.karalashvili@leibniz-zfmk.de
#$ -m n
#$ -N SNPmerge

#this script merges the per-region vcf files into a single SNP masterlist
## qsub -pe smp $CPUs -q fast.q -v "REGIONS=$REGIONS" -v "INPUTFOLDER=$INPUTFOLDER" -v "OUTPUTFOLDER=$OUTPUTFOLDER" $HOME/scripts/cluster-SNP-merge.sh

### load modules
module load htslib/1.19.1
module load samtools/1.19.2
module load bcftools/1.19
module load freebayes/1.3.2
module load vcftools/0.1.16
module load mambaforge/23.1.0
module load parallel/20201222
module load vcflib/1.0.3

## FYI install VT locally
## VT
#mkdir -p ~/bin
#mkdir -p ~/progz
#cd ~/progz
#git clone https://github.com/atks/vt.git  
#cd vt
#git submodule update --init --recursive 
#make -j 10 
#make test
#ln -s -f $PWD/vt ~/bin

REGIONS=$REGIONS
INPUTFOLDER=$INPUTFOLDER
OUTPUTFOLDER=$OUTPUTFOLDER

CPUs=$NSLOTS

FIRSTSCF=$(head -n1 $REGIONS)
tabix -H $INPUTFOLDER/$FIRSTSCF.vcf.gz > $OUTPUTFOLDER/header
cat $OUTPUTFOLDER/header > $OUTPUTFOLDER/$SPECIES.merged.vcf
zcat $INPUTFOLDER/*.vcf.gz | grep -v '#' >> $OUTPUTFOLDER/$SPECIES.merged.vcf
cat $OUTPUTFOLDER/$SPECIES.merged.vcf | vcfstreamsort -a > $OUTPUTFOLDER/$SPECIES.merged.sorted.vcf
bgzip -f -@ $CPUs $OUTPUTFOLDER/$SPECIES.merged.sorted.vcf
tabix -fp vcf $OUTPUTFOLDER/$SPECIES.merged.sorted.vcf.gz
vt peek $OUTPUTFOLDER/$SPECIES.merged.sorted.vcf.gz 2> $OUTPUTFOLDER/$SPECIES.merged.sorted.vcf.stats
cat $OUTPUTFOLDER/$SPECIES.merged.sorted.vcf.stats

#rm -f $OUTPUTFOLDER/$SPECIES.merged.vcf

echo "SNP merging done"
