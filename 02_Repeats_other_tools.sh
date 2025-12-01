#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M e.karalashvili@leibniz-zfmk.de
#$ -m n
#$ -N repeats2

#this script uses various tools to find repeat elements in the genome
# qsub -pe smp $CPUs -q large.q,small.q,medium.q -v "OUTPUTFOLDER=$OUTPUTFOLDER" -v "INPUTFASTA=$INPUTFASTA" -v "SPECIES=$SPECIES" ~/scripts/submission.cluster-repeats2-OTHER.sh

# catch and process inputs
CPUs=$NSLOTS
INPUT="$INPUT"
SPECIESNAME="$SPECIESNAME"
OUT="$OUT"
mkdir -p $OUT

## test if input exists
if [ -f "$INPUT" ]; then
    echo "INPUT exists"

else
    echo "INPUT does not exist: $INPUT" && exit 1
fi

## check if input is gzipped, unpack if needed, et File and Foldernames
if [[ $INPUT == *.gz ]]
then 
 echo "gzipped input, unpacking"
 pigz -dfk $INPUT
 INPUTFILE=$(echo $INPUT | rev | cut -d"/" -f1 | cut -d"." -f2- | rev) && echo $INPUTFILE
 INPUTFILENAME=$(echo -n $INPUTFILE | rev | cut -d"." -f2- | rev ) && echo $INPUTFILENAME
 INPUTFOLDER=$(echo $INPUT | rev | cut -d"/" -f2- | rev) && echo $INPUTFOLDER
else
 echo "not gzipped input"
  INPUTFILE=$(echo $INPUT | rev | cut -d"/" -f1 | cut -d"." -f1- | rev) && echo $INPUTFILE
  INPUTFILENAME=$(echo -n $INPUTFILE | rev | cut -d"." -f2- | rev ) && echo $INPUTFILENAME
  INPUTFOLDER=$(echo $INPUT | rev | cut -d"/" -f2- | rev) && echo $INPUTFOLDER
fi

CURRENTFOLDER=$PWD
cd $OUT
OUTPUT="$OUT/$INPUTFILENAME.repeats"
mkdir -p $OUTPUT
VARIABLESFILE="$OUTPUT/repeat-analyses.variables2.txt"

echo $INPUT > $VARIABLESFILE
echo $CPUs >> $VARIABLESFILE
echo $SPECIESNAME >> $VARIABLESFILE
echo $CURRENTFOLDER >> $VARIABLESFILE
echo $INPUTFILE >> $VARIABLESFILE
echo $INPUTFILENAME >> $VARIABLESFILE
echo $INPUTFOLDER >> $VARIABLESFILE
echo $OUT >> $VARIABLESFILE
echo $(date) >> $VARIABLESFILE


# load bashrc or $HOME/bin and conda_env of earlygrey
source ~/.bashrc
#export PATH=$PATH:$HOME/bin
module load mambaforge/23.1.0
mamba activate $HOME/progz/conda_envs/earlgrey

## re-asssign variables (after conda activation)
INPUT=$(cat $VARIABLESFILE | sed -n '1p')
CPUs=$(cat $VARIABLESFILE | sed -n '2p')
SPECIESNAME=$(cat $VARIABLESFILE | sed '3q;d')
CURRENTFOLDER=$(cat $VARIABLESFILE | sed '4q;d')
INPUTFILE=$(cat $VARIABLESFILE | sed '5q;d')
INPUTFILENAME=$(cat $VARIABLESFILE | sed '6q;d')
INPUTFOLDER=$(cat $VARIABLESFILE | sed '7q;d')
OUT=$(cat $VARIABLESFILE | sed '8q;d')
STARTDATE=$(cat $VARIABLESFILE | sed '9q;d')

## make outputfolder
SPECIES="$SPECIESNAME"
OUTPUT="$OUT/$INPUTFILENAME.repeats"
mkdir -p $OUTPUT
REF="$INPUTFOLDER/$INPUTFILE"
echo "$OUTPUT" >> $VARIABLESFILE

########################
########################
########################
## other repeat analyses

# Fasta Index
samtools faidx $REF
cat $REF.fai | cut -f1,2 > $REF.length

# occurences of N
printf "ctg\tstart\tend\t0\n" > $OUTPUT/$INPUTFILENAME.NNN.bed
seqtk cutN -g -n 1 $REF | bioawk -t '{print $0,$3-$2}' >> $OUTPUT/$INPUTFILENAME.NNN.bed
cat $OUTPUT/$INPUTFILENAME.NNN.bed | bioawk '{sum+=$4;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.NNN.bed.sum
printf "N\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.NNN.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt

# long homopolymer runs
seqtk hrun $REF | sed "s/ /\t/g" | bioawk -t '{print $0,$3-$2}' > $OUTPUT/$INPUTFILENAME.hrun.bed
cat $OUTPUT/$INPUTFILENAME.hrun.bed | awk '$5>=10' > $OUTPUT/$INPUTFILENAME.hrun.10.bed
cat $OUTPUT/$INPUTFILENAME.hrun.bed | bioawk '{sum+=$5;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.hrun.bed.sum
cat $OUTPUT/$INPUTFILENAME.hrun.10.bed | bioawk '{sum+=$5;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.hrun.10.bed.sum
printf "hrun\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.hrun.10.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt

# sdust
sdust -w 64 -t 20 $REF | bioawk -t '{print $0,$3-$2}' > $OUTPUT/$INPUTFILENAME.sdust.bed
cat $OUTPUT/$INPUTFILENAME.sdust.bed | awk '$4>=10' > $OUTPUT/$INPUTFILENAME.sdust.10.bed
cat $OUTPUT/$INPUTFILENAME.sdust.bed | bioawk '{sum+=$4;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.sdust.bed.sum
cat $OUTPUT/$INPUTFILENAME.sdust.10.bed | bioawk '{sum+=$4;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.sdust.10.bed.sum
printf "sdust\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.sdust.10.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt

## SSRs
#eTRF
etrf -m 200 -l 13 $REF | bioawk -t '{print $0,$3-$2}' > $OUTPUT/$INPUTFILENAME.SSRs.etrf.bed
cat $OUTPUT/$INPUTFILENAME.SSRs.etrf.bed | bioawk '{sum+=$7;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.SSRs.etrf.bed.sum
printf "SSRs.etrf\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.SSRs.etrf.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt

#SA-SSR
sa-ssr --min-repeats 5 --min-nucs 10 --max-ssr-len 12 --max-seq-len 25000000 -t $CPUs $REF $OUTPUT/$INPUTFILENAME.SSRs.sa-ssr.tsv
cat $OUTPUT/$INPUTFILENAME.SSRs.sa-ssr.tsv | tail -n+2 | sed 's/\t/____/' | sed 's/\s.*____/\t/' | sed 's/____/\t/' | bioawk -t '{$5=length($2)}{$6=$5*$3}{$7=$4+$6}NF' | tabtk cut -r -f1,4,7,2,5,3,6 > $OUTPUT/$INPUTFILENAME.SSRs.sa-ssr.bed
cat $OUTPUT/$INPUTFILENAME.SSRs.sa-ssr.bed | bioawk '{sum+=$7;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.SSRs.sa-ssr.bed.sum
printf "SSRs.sa-ssr\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.SSRs.sa-ssr.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.SSRs.sa-ssr.bed | cut -f1-3 | bioawk -t '{print $0,$3-$2}' | bedtools sort > $OUTPUT/$INPUTFILENAME.SSRs.sa-ssr.formatted.bed

export PATH=$HOME/progz/gcc12/bin:$PATH
export LD_LIBRARY_PATH=$HOME/progz/gcc12/lib64:$LD_LIBRARY_PATH
export CC=$HOME/progz/gcc12/bin/gcc
export CXX=$HOME/progz/gcc12/bin/g++
export FC=$HOME/progz/gcc12/bin/gfortran
#module load gcc/11.2.0  #newer gcc needs to be loaded for tantan to run (despite otherwise it compiled)

#Tantan
tantan -f4 $REF | bioawk -t '{print $0,$3-$2}' > $OUTPUT/$INPUTFILENAME.SSRs.tantan.bed
cat $OUTPUT/$INPUTFILENAME.SSRs.tantan.bed | awk '$5>=5' | bedtools sort > $OUTPUT/$INPUTFILENAME.SSRs.tantan.5.bed
cat $OUTPUT/$INPUTFILENAME.SSRs.tantan.bed | bioawk '{sum+=$8;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.SSRs.tantan.bed.sum
cat $OUTPUT/$INPUTFILENAME.SSRs.tantan.5.bed | bioawk '{sum+=$8;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.SSRs.tantan.5.bed.sum
printf "SSRs.tantan\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.SSRs.tantan.5.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt

#ULTRA --minunit 4 --minlen 14
ultra --min_unit 4 --min_length 10 --period 300 --threads $CPUs $REF | grep -v "^\"" | grep -v "^\{" | grep -v "^\}" > $OUTPUT/$INPUTFILENAME.SSRs.ultra.bed
cat $OUTPUT/$INPUTFILENAME.SSRs.ultra.bed | bioawk '{sum+=$9;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.SSRs.ultra.bed.sum
cat $OUTPUT/$INPUTFILENAME.SSRs.ultra.bed | bioawk -t '{print $0,$3-$2}' | bioawk '$11 >100 ' > $OUTPUT/$INPUTFILENAME.SSRs.ultra.100bp.bed
cat $OUTPUT/$INPUTFILENAME.SSRs.ultra.100bp.bed | bioawk '{sum+=$9;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.SSRs.ultra.100bp.bed.sum
printf "SSRs.ultra\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.SSRs.ultra.100bp.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.SSRs.ultra.bed | cut -f1-3 | sed -E 's/(.*)_/\1./' | bioawk -t '{print $0,$3-$2}' | bedtools sort > $OUTPUT/$INPUTFILENAME.SSRs.ultra.formatted.bed
cat $OUTPUT/$INPUTFILENAME.SSRs.ultra.100bp.bed | cut -f1-3 | sed -E 's/(.*)_/\1./' | bioawk -t '{print $0,$3-$2}' | bedtools sort > $OUTPUT/$INPUTFILENAME.SSRs.ultra.100bp.formatted.bed

#srf
mkdir -p $OUTPUT/$INPUTFILENAME.srf_tmp_dir
kmc -fm -k171 -t$CPUs -ci20 -cs100000 $REF $OUTPUT/$INPUTFILENAME.count.kmc $OUTPUT/$INPUTFILENAME.srf_tmp_dir
kmc_dump $OUTPUT/$INPUTFILENAME.count.kmc $OUTPUT/$INPUTFILENAME.count.txt
srf -p $REF $OUTPUT/$INPUTFILENAME.count.txt > $OUTPUT/$INPUTFILENAME.srf.fa
srfutils.js enlong $OUTPUT/$INPUTFILENAME.srf.fa > $OUTPUT/$INPUTFILENAME.srf.enlong.fa
minimap2 -c -N1000000 -f1000 -r100,100 -t$CPUs $OUTPUT/$INPUTFILENAME.srf.enlong.fa $REF > $OUTPUT/$INPUTFILENAME.srf.paf
touch $OUTPUT/$INPUTFILENAME.srf.bed
srfutils.js paf2bed $OUTPUT/$INPUTFILENAME.srf.paf > $OUTPUT/$INPUTFILENAME.srf.bed
srfutils.js bed2abun $OUTPUT/$INPUTFILENAME.srf.bed > $OUTPUT/$INPUTFILENAME.srf.bed.abundance
trf-mod $OUTPUT/$INPUTFILENAME.srf.fa > $OUTPUT/$INPUTFILENAME.srf.decomposed.fa
# TRF-mod: decompose HORs to monomers; #many similar seq are mapped in same repeat? rerun and filter customly: minimap2 -c -N1000 <(./srfutils.js enlong -d srf.fa) srf.fa
cat $OUTPUT/$INPUTFILENAME.srf.bed | bioawk -t '{print $0,$3-$2}' > $OUTPUT/$INPUTFILENAME.srf.with.size.bed
cat $OUTPUT/$INPUTFILENAME.srf.with.size.bed | bioawk '{sum+=$9;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.srf.with.size.bed.sum
printf "Sats.srf\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.srf.with.size.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
rm -rf $OUTPUT/$INPUTFILENAME.srf_tmp_dir $OUTPUT/$INPUTFILENAME.count.*

# mappability
# FYI mappability=1: k-mer occurs once with up to e errors, low mappability: k-mer belongs to a repetitive region
genmap index --fasta-file $REF --index $OUTPUT/$INPUTFILENAME.genmap.index
genmap map -K 35 -E 2 --threads $CPUs --bedgraph --wig --txt --output $OUTPUT/$INPUTFILENAME.genmap -I $OUTPUT/$INPUTFILENAME.genmap.index
cat $OUTPUT/$INPUTFILENAME.genmap.bedgraph | awk '$4<=0.05' | bedtools merge | bioawk -t '{print $0,$3-$2}' > $OUTPUT/$INPUTFILENAME.genmap.0.05.bed
cat $OUTPUT/$INPUTFILENAME.genmap.0.05.bed | bioawk '{sum+=$4;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.genmap.0.05.bed.sum
cat $OUTPUT/$INPUTFILENAME.genmap.0.05.bed | bioawk -t '$4 >= 100' > $OUTPUT/$INPUTFILENAME.genmap.0.05.100bp.bed
cat $OUTPUT/$INPUTFILENAME.genmap.0.05.100bp.bed | bioawk '{sum+=$4;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.genmap.0.05.100bp.bed.sum
cat $OUTPUT/$INPUTFILENAME.genmap.0.05.bed | bioawk -t '$4 >= 1000' > $OUTPUT/$INPUTFILENAME.genmap.0.05.1000bp.bed
cat $OUTPUT/$INPUTFILENAME.genmap.0.05.1000bp.bed | bioawk '{sum+=$4;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.genmap.0.05.1000bp.bed.sum
printf "genmap.0.05\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.genmap.0.05.100bp.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
printf "genmap.0.05-1kb\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.genmap.0.05.1000bp.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
rm -rf $OUTPUT/$INPUTFILENAME.genmap.index $OUTPUT/$INPUTFILENAME.genmap.bedgraph $OUTPUT/$INPUTFILENAME.genmap.wig $OUTPUT/$INPUTFILENAME.genmap.txt

ENDTIME=$(date)
echo "finished repeat analysis"
echo "started: "$STARTTIME
echo "finished at: "$ENDTIME
cat $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cd $CURRENTFOLDER
mamba deactivate
echo "repeats2 done"
