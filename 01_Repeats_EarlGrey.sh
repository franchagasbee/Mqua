#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M e.karalashvili@leibniz-lib.de
#$ -m n
#$ -N earlgrey

#this script runs the EarlGrey pipeline for repeat discovery and classification
## usage: qsub -pe smp $CPUs -q small.q,medium.q,large.q -v "INPUTFASTA=$INPUTFASTA" -v "SPECIES=$SPECIES" ~/shells/repeats/01_Repeats_EarlGrey.sh

set -x
CPUs=$NSLOTS

## check if input is gzipped. unpack if needed.

if [[ $INPUTFASTA == *.gz ]]
then
    echo "gzipped input, unpacking"
    pigz -dfk $INPUTFASTA
else
    echo "input not gzipped, proceeding..."
fi

#STARTDATE=$(date)

## activate conda environment
source ~/.bashrc
module load mambaforge/23.1.0
mamba activate $HOME/progz/conda_envs/earlgrey

OUTPUT="${INPUTFASTA}.earlgrey"

echo "Input fasta: ${INPUTFASTA}"
echo "Species: ${SPECIES}"
echo "Output: ${OUTPUT}"
echo "CPUs: ${CPUs}"

## run EarlGrey pipeline
time earlGrey -g $INPUTFASTA -s $SPECIES -o $OUTPUT -t $CPUs -c yes -m yes -d yes

ENDDATE=$(date)

## process some files
cat $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'EarlGrey'.log | grep -v "% completed" | grep -v "\[0m" | grep -v "########" | grep -v 'running on ctg' > $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'EarlGrey'.short.log
zcat $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.bed.gz | bioawk -t '{print $0,$3-$2}' > $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed
cat $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed | bioawk '{sum+=$7;} END{print sum;}' > $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed.sum

## copy summary file(s))
cp $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed ./$INPUTFASTA.earlgrey.filteredRepeats.withSize.bed
cp $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed.sum ./$INPUTFASTA.earlgrey.filteredRepeats.withSize.bed.sum
cat $INPUTFASTA.earlgrey.filteredRepeats.withSize.bed | grep -vP 'Simple_repeat|Low_complexity|Satellite|Unknown' > $INPUTFASTA.earlgrey.filteredRepeats.withSize.TEonly.bed 
cat $INPUTFASTA.earlgrey.filteredRepeats.withSize.TEonly.bed | bioawk '{sum+=$7;} END{print sum;}' > $INPUTFASTA.earlgrey.filteredRepeats.withSize.TEonly.bed.sum

## write some stats
printf "EarlGrey\t" > repeats.stats.txt
cat $INPUTFASTA.earlgrey.filteredRepeats.withSize.bed.sum >> repeats.stats.txt
printf "EarlGreyTEonly\t" >> repeats.stats.txt
cat $INPUTFASTA.earlgrey.filteredRepeats.withSize.TEonly.bed.sum >> repeats.stats.txt
printf "Starting time: ${STARTDATE}" >> repeats.stats.txt
printf "Ending time: ${ENDDATE}" >> repeats.stats.txt
