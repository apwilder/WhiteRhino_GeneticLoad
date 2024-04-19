White rhino genomic data and analysis overview

from

Genetic load and viability of a future restored northern white rhino population 

Aryn P. Wilder, Cynthia C. Steiner, Sarah Hendricks, Benjamin C. Haller, Chang Kim, Marisa L. Korody, Oliver A. Ryder

https://doi.org/10.1111/eva.13683


1.	Read processing and calling genotypes for white rhinos
	a.	Trim reads with trimmomatic (WhiteRhino_Trimmomatic.sh)

	b.	Map with Bowtie2 in very-sensitive mode (WhiteRhino_Bowtie2.sh)

	c.	Deduplicate (WhiteRhino_MarkDups.sh)

	d.	Call gVCFs for each sample and genotype all samples together (GATK3_HapCaller_GenotypeGVCF.sh)

3.	Map black rhino short reads to SWR genome with bowtie2
	a.	Map: DicBic_bowtie2.sh

	b.	Genotype at WR SNPs: Genotype_DicBic.sh

5.	Choose hard filters, filter and extract variants
	a.	ChooseHardFilters_WhiteRhino.R

	b.	HardFilter_WhiteRhino.sh

	c.	ChooseHardFilters_DicBic.R

	d.	HardFilter_DicBic.sh

7.	Estimate WR allele frequencies in plink (WhiteRhino_Plink_AF.sh)
8.	Get GERP scores for all variants
	a.	Download GERPs from sidow lab site (from alignment w/out WR)

	b.	Liftover CerSim to Hg18 to Hg19 and output GERPs for each scaffold (GetGERPs.R)

	c.	Unzip and cat chr_*_allFilt_GERPscores.txt > AllChr_allFilt_GERPscores.txt

	d.	Merge GERP scores with variant table and polarize ancestral/derived alleles with black rhino genome (WhiteRhino_MergeGERP.R)

10.	Annotate protein-coding variants with snpEff (WhiteRhino_snpEff.sh)
11.	Estimate ROH in genome with bcftools (WhiteRhino_ROH.sh)
12.	Some downstream analyses of genetic load in SWR and NWR (DownstreamGenomicAnalysis.R)
13.	Simulate future dynamics of genetic load and fitness in SLiM
a.	Create input genomes to start simulations (CreateSLiMinput.R)
b.	Simulate 10 generations with and without reintroduction of founders (Rhino_supp.slim and Rhino_nosupp.slim)
c.	Summarize SLiM output files (ReadSLiMoutfiles4server_slimv4.R)
d.	Plot results (PlotSimLoad.R)


