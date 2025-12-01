#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M e.karalashvili@leibniz-zfmk.de
#$ -m n
#$ -N repeats3

#this script runs trf-mod for repeat element discovery
# qsub -pe smp $CPUs -q large.q,small.q,medium.q -v "OUTPUTFOLDER=$OUTPUTFOLDER" -v "INPUTFASTA=$INPUTFASTA" -v "SPECIES=$SPECIES" ~/scripts/submission.cluster-repeats3-TRFmod.sh

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
VARIABLESFILE="$OUTPUT/repeat-analyses.variables3a.txt"

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

## 
SPECIES="$SPECIESNAME"
REF="$INPUTFOLDER/$INPUTFILE"
OUTPUT="$OUT/$INPUTFILENAME.repeats"
mkdir -p $OUTPUT
echo "$OUTPUT" >> $VARIABLESFILE

########################
########################
########################
## trf-mod

## ATTENTION. extremely long runtime of TRF on long centromere regions, TRF will get stuck
## solution: set a higher value for -l
## "-l to >100 for chm13-T2T helped in my case" , memory usage can get pretty high
## -l INT     maximum TR length expected (in millions) [2]

# run trf-mod - everything is done in case this hangs (then you need to kill the job)
## trf-mod takes extrmeel long sometimes (eg some scaffolds in Bombus terrestris)

time trf-mod -p 300 -l 25 $REF | bioawk -t '{print $0,$3-$2}' | bedtools sort > $OUTPUT/$INPUTFILENAME.SSRs.trfmod.bed
cat $OUTPUT/$INPUTFILENAME.SSRs.trfmod.bed | bioawk '{sum+=$4;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.SSRs.trfmod.bed.sum
printf "SSRs.trfmod\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $OUTPUT/$INPUTFILENAME.SSRs.trfmod.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt

cd $CURRENTFOLDER
mamba deactivate
echo "TRFmod done"
