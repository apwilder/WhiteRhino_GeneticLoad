#!/bin/bash

#build rhino database in snpEff
#java -jar snpEff.jar build -gtf22 -v CerSim1.0
java -Xmx4g -jar snpEff.jar build -gff3 -v CerSim1.0 >& BuildDB_noCGP.nohup

#run snpEff
java -Xmx4g -jar snpEff.jar CerSim1.0 \
AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2.vcf.gz -o vcf \
-s AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_snpEffnoCGP.summary.html \
-lof | awk -F"\t|ANN=" -v OFS="\t" '{print $1, $2, $4, $5, $9}' | grep -v "#" > \
AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_snpEffnoCGP_ann.txt
