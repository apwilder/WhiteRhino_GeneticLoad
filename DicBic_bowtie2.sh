#!/bin/bash

bowtie2 -q --phred33 --very-sensitive -p 16 -I 0 -X 2000 --fr \
--rg-id DicBicMic --rg SM:DicBicMic --rg LB:Broad --rg PU:Broad \
--rg PL:ILLUMINA -x CerSimSim1.0_genomic -1 /public/groups/cgl/nwr/BlackRhino/DicBicMic_1_R1.fastq.gz \
-2 /public/groups/cgl/nwr/BlackRhino/DicBicMic_1_R2.fastq.gz | \
/public/groups/cgl/nwr/bin/samtools-1.8/samtools view -bS -F 4 \
> DicBicMic_AdptTrimPaired_bt2_CerSim1.0_genomic.bam

bowtie2 -q --phred33 --very-sensitive -p 16 -I 0 -X 2000 --fr \
--rg-id DicBicMic --rg SM:DicBicMic --rg LB:Broad --rg PU:Broad \
--rg PL:ILLUMINA -x CerSimSim1.0_genomic -1 DicBicMic_1_R1.fastq.gz \
-2 DicBicMic_1_R2.fastq.gz | samtools view -bS -F 4 \
> DicBicMic_AdptTrimPaired_bt2_CerSim1.0_genomic.bam

bowtie2 \
-q --phred33 --very-sensitive -p 12 \
-I 0 -X 2000 --fr --rg-id DicBicMin --rg SM:DicBicMin \
--rg LB:Broad --rg PU:Broad --rg PL:ILLUMINA -x CerSimSim1.0_genomic \
-1 SRR12010331_AdaptTrim_P1.fastq.gz -2 SRR12010331_AdaptTrim_P2.fastq.gz \
| samtools view -bS -F 4 > DicBicMin_AdptTrimPaired_bt2_CerSim1.0_genomic.bam

