#!/bin/bash

SAMPLELIST=$1 # Filename and path to list of fastq files
OUTPUTDIR=$2 # path to directory where output files are to be written

##### RUN EACH SAMPLE THROUGH PIPELINE #######
# Loop over each sample
for ID in `cat $SAMPLELIST | cut -f1`; do
SAMPLENUM=`grep "$ID" $SAMPLELIST | cut -f2`
LANE=`grep "$ID" $SAMPLELIST | cut -f3`
SAMPLE=$OUTPUTDIR'P9109_'$ID'_S'$SAMPLENUM'_L00'$LANE

#### CLEANING THE READS ####
# Remove adapter sequence with Trimmomatic. 
trimmomatic PE -threads 6 -phred33 -trimlog \
$SAMPLE'_trimlog.out' $SAMPLE'_R1_001.fastq.gz' $SAMPLE'_R2_001.fastq.gz' \
$SAMPLE'_AdptTrim_P1.fastq.gz' $SAMPLE'_AdptTrim_U1.fastq.gz' \
$SAMPLE'_AdptTrim_P2.fastq.gz' $SAMPLE'_AdptTrim_U2.fastq.gz' \
ILLUMINACLIP:/hive/users/awilder/bin/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:4:true
done
