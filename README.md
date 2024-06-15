# High-Runner selection experiment genomic comparisons between High-Runner and Control lines generations at 61

This code is based those used in:
	Hillis D. A., L. Yadgary, G. M. Weinstock, F. Pardo-Manuel de Villena, D. Pomp, et al., 2020 Genetic 
	basis of aerobically supported voluntary exercise: results from a selection experiment with house mice. 
	Genetics 216: 781–804. https://doi.org/10.1534/genetics.120.303668


Data can be found either in the at figshare: https://doi.org/10.25386/ genetics.12436649

For the replication of the publication's complete analyses, the Rmarkdown and SAS code need to be run 
interchangeably. The progression is described in the Rmarkdown file:

Workflow for WGS SNP analyses [R code] [SAS Code]
  # 1 Subset Allelic Data for Simple Analyses (only for 100 loci sample)
  # 2 Prep Allelic Data for Mixed Model Analyses
  # 3 Create Template Files
  WGS SAS Code
  # 4 Merge Results
  # 5 AIC Model Selection
  # 6 Identify P-Value
  # 7 No Within-Line Variance
  # 8 Write Results
  # 9 Group SNPs
  # 10 Select Maxima from Groups

Haplotype
  # H11 Recode Alleles
  # 3 Create Template Files (Not necessary if files still remain from SNP analyses)
  Haplotype SAS Code
  # H12 Functions for Combining Haplotype Results
  # H13 Read Haplotype Results
  # H14 Run Haplotype Functions
  # H15 Combine into One File
 
Running the code in its entirety will address only a subset of loci. However, the first R chunk can be 
skipped and "File_S2_sample.csv" replaced with "File S2.csv" in the second chunk to analyze the whole data
set.

In short, this code takes the alleleic data and organizes it (in R) for mixed-model ANOVA analyses for four 
different variance structures (in SAS). Then (in R) the best variance structure is selected for each locus 
and the results are aggregated and summarized. Haplotype analyses are similar to WGS analyses but due to 
the use of three major haplotypes of many loci, two analyses are run for each of these loci and p-values 
are combined. 

The both SAS and R are used to minimize processing time. 
R for data wrangling and SAS for mixed-model analyses.