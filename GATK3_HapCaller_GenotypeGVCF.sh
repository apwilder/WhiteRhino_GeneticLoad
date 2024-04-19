#!/bin/bash

cd /hive/users/awilder/Mapping/

for BAMFILE in `ls *_dedup.bam`; do
SAMPLE=`echo $BAMFILE | sed s/_AdptTrimPaired_bt2_CerSimSim1.0_genomic_dedup.bam//g`
java -jar GenomeAnalysisTK.jar \
-T HaplotypeCaller -ERC GVCF \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-I $BAMFILE \
-o $SAMPLE'_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz' \
-nct 12 >& $SAMPLE'_gvcf.out'
done

cd /hive/users/awilder/GATK/

mkdir -p tmp

java -Djava.io.tmpdir=`pwd`/tmp -jar GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V NWR_KB3731_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V NWR_KB5763_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V NWR_KB5764_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V NWR_KB5766_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V NWR_KB6571_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V NWR_KB8174_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V NWR_KB8175_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V NWR_KB9939_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V NWR_KB9947_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_102_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_103_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_104_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_105_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_106_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_107_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_108_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_109_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_110_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_KB13306_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_KB5892_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_KB6974_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-V SWR_KB7062_AdptTrimPaired_bt2_CerSimSim1.0_raw.g.vcf.gz \
-o AllWR_CerSimSim1.0_raw.vcf.gz -nt 8 >& AllWR_raw_vcf.out


