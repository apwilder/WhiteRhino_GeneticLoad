#!/bin/bash

export JAVA_HOME=/hive/users/awilder/bin/jre1.8.0_152
export PATH="$JAVA_HOME/bin:$PATH"

REFBASENAME=$1 # reference name to add to output files, e.g. CerSimSim1.0_genomic
SAMPLELIST=$2 # Filename and path to list of fastq files, e.g. ../SampleLists/LoveSampleTable.txt
OUTPUTDIR=$3 # path to directory where output files are to be written

cd $OUTPUTDIR
mkdir tmp
##### RUN EACH SAMPLE THROUGH PIPELINE #######
# Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST | cut -f6`; do

SAMPLE=$OUTPUTDIR$SAMPLEFILE

samtools sort -t 8 -o $SAMPLE'_AdptTrimPaired_bt2_'$REFBASENAME'_sort.bam' \
$SAMPLE'_AdptTrimPaired_bt2_'$REFBASENAME'.bam'

samtools index $SAMPLE'_AdptTrimPaired_bt2_'$REFBASENAME'_sort.bam'

java -Xmx60g -Djava.io.tmpdir=`pwd`/tmp -XX:ParallelGCThreads=8 -jar picard.jar MarkDuplicates \
I=$SAMPLE'_AdptTrimPaired_bt2_'$REFBASENAME'_sort.bam' O=$SAMPLE'_AdptTrimPaired_bt2_'$REFBASENAME'_dedup.bam' \
M=$SAMPLE'_AdptTrimPaired_bt2_'$REFBASENAME'_dupstat.txt' \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true TMP_DIR=`pwd`/tmp

done
