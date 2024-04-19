#!/bin/bash

#select and tabulate raw SNPs:
java -jar GenomeAnalysisTK.jar \
-T SelectVariants -R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V AllWR_CerSimSim1.0_raw.vcf.gz -selectType SNP -selectType NO_VARIATION \
-o AllWR_CerSimSim1.0_raw_snps.vcf.gz #done

java -jar GenomeAnalysisTK.jar \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna -T VariantsToTable \
-V AllWR_CerSimSim1.0_raw_snps.vcf.gz -F CHROM -F POS -F AN -F BaseQRankSum \
-F DP -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR \
--out AllWR_CerSimSim1.0_raw_snps.tab #done

#select and tabulate raw indels:
java -jar GenomeAnalysisTK.jar \
-T SelectVariants -R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V AllWR_CerSimSim1.0_raw.vcf.gz -selectType INDEL \
-o AllWR_CerSimSim1.0_raw_indels.vcf.gz >& SubsetBRindels.nohup &

java -jar GenomeAnalysisTK.jar \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna -T VariantsToTable \
-V AllWR_CerSimSim1.0_raw_indels.vcf.gz -F CHROM -F POS -F AN -F BaseQRankSum \
-F DP -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR \
--out AllWR_CerSimSim1.0_raw_indels.tab >& TabulateBRindels.nohup

vcftools --gzvcf AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2.vcf.gz \
--out AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3 \
--minGQ 20 --max-missing 0.3 --recode --recode-INFO-all \
>& vcftools_recode_GQ20mm3.nohup

vcftools --vcf AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3.recode.vcf \
--extract-FORMAT-info GT --out AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3 \
>& WR_GT_GQ20_mm3_extraction.nohup

