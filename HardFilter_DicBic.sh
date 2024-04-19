#!/bin/bash

java -jar GenomeAnalysisTK.jar \
-T VariantFiltration \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_snps.vcf.gz \
--filterExpression "QD < 22.0" --filterName "QD" --filterExpression "FS > 32.0" \
--filterName FS --filterExpression "MQ < 32.0" --filterName MQ \
--filterExpression "MQRankSum < -5.0" --filterName MQRankSum \
--filterExpression "ReadPosRankSum < -2.0" --filterName ReadPosRankSum \
--filterExpression "SOR > 3.0" --filterName "SOR" \
-o DicBic_AdptTrimPaired_bt2_CerSim1.0_HardFilt_snps.vcf.gz >& BR_FiltSNPS.nohup

java -jar GenomeAnalysisTK.jar \
-T SelectVariants \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_HardFilt_snps.vcf.gz \
-ef -o DicBic_AdptTrimPaired_bt2_CerSim1.0_QD22_FS32_MQ32_MQRS5_RPRS2_SOR3_snps.vcf.gz \
>& GATK_rmHardFilt_snps.nohup

java -jar GenomeAnalysisTK.jar \
-T VariantFiltration \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_indels.vcf.gz \
--filterExpression "QD < 22.0" --filterName "QD" --filterExpression "FS > 32.0" --filterName FS \
--filterExpression "MQRankSum < -5.0" --filterName MQRankSum \
--filterExpression "ReadPosRankSum < -2.0" --filterName ReadPosRankSum \
--filterExpression "SOR > 3.0" --filterName "SOR" \
-o DicBic_AdptTrimPaired_bt2_CerSim1.0_HardFilt_indels.vcf.gz \
>& GATK_HardFilt_BRindels.nohup

java -jar GenomeAnalysisTK.jar \
-T VariantFiltration \
-R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_HardFilt_indels.vcf.gz \
-ef -o DicBic_AdptTrimPaired_bt2_CerSim1.0_QD22_FS32_MQRS5_RPRS2_SOR3_indels.vcf.gz \
>& GATK_rmHardFilt_BRindels.nohup

#combine filtered variants
java -jar GenomeAnalysisTK.jar \
-T CombineVariants -R GCA_000283155.1_CerSimSim1.0_genomic.fna \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_QD22_FS32_MQRS5_RPRS2_SOR3_indels.vcf.gz \
-V DicBic_AdptTrimPaired_bt2_CerSim1.0_QD22_FS32_MQ32_MQRS5_RPRS2_SOR3_snps.vcf.gz \
-assumeIdenticalSamples --genotypemergeoption UNSORTED \
-o DicBic_AdptTrimPaired_bt2_CerSim1.0_QD22_FS32_MQ32_MQRS5_RPRS2_SOR3.vcf.gz \
>& GATK_combineFiltVariants.nohup

#get filtered genotypes
vcftools --gzvcf DicBic_AdptTrimPaired_bt2_CerSim1.0_QD22_FS32_MQ32_MQRS5_RPRS2_SOR3.vcf.gz \
--extract-FORMAT-info GT --out DicBic_AdptTrimPaired_bt2_CerSim1.0_QD22_FS32_MQ32_MQRS5_RPRS2_SOR3 \
>& BR_GTextraction.nohup

