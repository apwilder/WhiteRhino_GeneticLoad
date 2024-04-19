#!/bin/bash

REFERENCE=$1 # path to reference fasta file, e.g /hive/users/awilder/rhinoref/GCA_000283155.1_CerSimSim1.0_genomic.fna.gz
REFBASENAME=$2 # reference name to add to output files, e.g. CerSimSim1.0_genomic
SAMPLELIST=$3 # path to table fastq files with info on fastq files, e.g. /hive/users/awilder/SampleLists/Run1SampleTable.txt
FASTQDIR=$4 # path to directory of fastq files, e.g. /hive/users/awilder/RhinoFastqs/
OUTPUTDIR=$5 # path to directory where output files are to be written, e.g. /hive/users/awilder/Mapping/

cd $OUTPUTDIR
##### RUN EACH SAMPLE THROUGH PIPELINE #######
# Loop over each sample
for ID in `cat $SAMPLELIST | cut -f1`; do

#### MAPPING THE GENOMIC READS TO THE REFERENCE ####
# Extract relevant values from a table of sample ID, library ID, and platform unit (here in columns 2, 3, and 4, respectively) for each sequenced library
SAMPLENUM=`grep "$ID" $SAMPLELIST | cut -f2`
PU=`grep "$ID" $SAMPLELIST | cut -f3`
LIBRARY_ID=`grep "$ID" $SAMPLELIST | cut -f5`
NAME=`grep "$ID" $SAMPLELIST | cut -f6`
#SAMPLE1=$FASTQDIR'P9109_'$ID'_S'$SAMPLENUM'_L00'$LANE'_AdptTrim_P1.fastq.gz' #for Love fastq files
#SAMPLE2=$FASTQDIR'P9109_'$ID'_S'$SAMPLENUM'_L00'$LANE'_AdptTrim_P2.fastq.gz'
SAMPLE1=$FASTQDIR$ID'_1_sequence.txt.gz' #for Run1 fastq files
SAMPLE2=$FASTQDIR$ID'_2_sequence.txt.gz'

# Map the paired-end reads: use --very-sensitive-local for reference transcriptome and --very-sensitive for reference genome
bowtie2 -q --phred33 --very-sensitive -p 12 -I 0 -X 2000 --fr --rg-id $NAME --rg SM:$NAME \
--rg LB:$LIBRARY_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLE1 -2 $SAMPLE2 \
| samtools view -bS -F 4 > $NAME'_AdptTrimPaired_bt2_'$REFBASENAME'.bam'

done
