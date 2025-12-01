#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M e.karalashvili@leibniz-zfmk.de
#$ -m e
#$ -N repeats4-merge

#this script merges all the outputs of the repeat annotation scripts into a single file
set -x

CPUs=$NSLOTS
INPUT="$INPUT"
SPECIESNAME="$SPECIESNAME"
OUT="$OUT"

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
OUTPUT="$OUT"
VARIABLESFILE="$OUTPUT/repeat-analyses.variables4.txt"

# load bashrc or $HOME/bin and conda_env of earlygrey
source ~/.bashrc
module load mambaforge/23.1.0
mamba activate $HOME/progz/conda_envs/earlgrey

## 
SPECIES="$SPECIESNAME"
REF="$INPUTFOLDER/$INPUTFILE"
OUTPUT="$OUT"
echo "$OUTPUT" >> $VARIABLESFILE

cat $OUTPUT/$INPUTFILE.earlgrey.filteredRepeats.withSize.bed \
$OUTPUT/$INPUTFILE.NNN.bed \
$OUTPUT/$INPUTFILE.sdust.10.bed \
$OUTPUT/$INPUTFILE.SSRs.ultra.formatted.bed \
$OUTPUT/$INPUTFILE.SSRs.etrf.bed \
$OUTPUT/$INPUTFILE.SSRs.tantan.5.bed \
$OUTPUT/$INPUTFILE.SSRs.sa-ssr.bed \
$OUTPUT/$INPUTFILE.srf.with.size.bed \
$OUTPUT/$INPUTFILE.SSRs.trfmod.bed \
$OUTPUT/$INPUTFILE.hrun.10.bed \
  | cut -f1-3 | sort -k1,1 -k2,2n \
  | bedtools merge -d 2 | bedtools sort \
  | bioawk -t '{print $0,$3-$2}' \
  | bgzip -f -@ $CPUs -c /dev/stdin > $OUTPUT/$INPUTFILENAME.repetitive.bed.gz
tabix -fp bed $OUTPUT/$INPUTFILENAME.repetitive.bed.gz
zcat $OUTPUT/$INPUTFILENAME.repetitive.bed.gz | bioawk '{sum+=$4;} END{print sum;}' > $OUTPUT/$INPUTFILENAME.repetitive.bed.sum

printf "Merged.repetitive\t" >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt
cat $REF.repetitive.bed.sum >> $OUTPUT/$INPUTFILENAME.repeats.stats.txt

cd $CURRENTFOLDER
mamba deactivate
echo "repeat merging done"
