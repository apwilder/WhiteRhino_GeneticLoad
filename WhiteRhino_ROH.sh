#!/bin/bash

bcftools-1.8/bcftools roh -e SWR_samplelist.txt -o SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh.txt \
AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2.vcf.gz >& SWRfilt_roh.nohup
bcftools-1.8/bcftools roh -e NWR_samplelist.txt -o NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh.txt \
AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2.vcf.gz >& NWRfilt_roh.nohup

awk '$2 ~ /NWR_/' NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh.txt > \
NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_NWRfilt.txt
awk '$2 ~ /SWR_/' SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh.txt > \
SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_SWRfilt.txt

awk '$1=="RG"' SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_SWRfilt.txt > \
SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_SWRfilt_RG.txt
awk '$1=="ST"' SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_SWRfilt.txt > \
SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_SWRfilt_ST.txt

awk '$1=="RG"' NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_NWRfilt.txt > \
NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_NWRfilt_RG.txt
awk '$1=="ST"' NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_NWRfilt.txt > \
NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_NWRfilt_ST.txt
