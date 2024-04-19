#!/bin/bash

java -jar GenomeAnalysisTK.jar \
-T HaplotypeCaller -ERC GVCF \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-I DicBicMic_AdptTrimPaired_bt2_CerSim1.0_genomic_sort.bam \
-I DicBicMin_AdptTrimPaired_bt2_CerSim1.0_genomic_sort.bam \
-L AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2.vcf.gz \
--output_mode EMIT_ALL_SITES \
-o DicBic_AdptTrimPaired_bt2_CerSim1.0_raw.g.vcf.gz -nct 12 >& DicBic_gvcf.out

java -jar GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna \
--includeNonVariantSites \
--variant DicBic_AdptTrimPaired_bt2_CerSim1.0_raw.g.vcf.gz \
-o DicBic_AdptTrimPaired_bt2_CerSim1.0_raw.vcf.gz -nt 12 >& DicBic.out

vcftools --gzvcf DicBic_AdptTrimPaired_bt2_CerSim1.0_raw.vcf.gz \
--extract-FORMAT-info GT --out DicBic_AdptTrimPaired_bt2_CerSim1.0_raw

#select and tabulate raw SNPs:
java -jar GenomeAnalysisTK.jar \
-T SelectVariants -R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_raw.vcf.gz -selectType SNP -selectType NO_VARIATION \
-o DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_snps.vcf.gz #done

java -jar GenomeAnalysisTK.jar \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna -T VariantsToTable \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_snps.vcf.gz -F CHROM -F POS -F AN -F BaseQRankSum \
-F DP -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR \
--out DicBicMic_AdptTrimPaired_bt2_CerSim1.0_raw_snps.tab #done

#select and tabulate raw indels:
java -jar GenomeAnalysisTK.jar \
-T SelectVariants -R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_raw.vcf.gz -selectType INDEL \
-o DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_indels.vcf.gz >& SubsetBRindels.nohup

java -jar GenomeAnalysisTK.jar \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna -T VariantsToTable \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_indels.vcf.gz -F CHROM -F POS -F AN -F BaseQRankSum \
-F DP -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR \
--out DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_indels.tab >& TabulateBRindels.nohup

