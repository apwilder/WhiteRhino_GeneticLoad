#!/bin/bash

vcftools --gzvcf AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3.vcf.gz \
--out AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink --plink

plink2 --threads 8 --memory 2500 --bfile AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink \
--freq --keep SWR_samplelist_plink.txt --out SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink \
>& SWR_AF_plink.nohup

plink2 --threads 8 --memory 2500 --bfile  AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink \
--freq --keep NWR_samplelist_plink.txt --out NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink \
>& NWR_AF_plink.nohup

sed 's/:/\t/g' SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink.afreq \
> SWR__CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink_AF.txt
sed 's/:/\t/g' NWR__CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink.afreq \
> NWR__CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink_AF.txt

