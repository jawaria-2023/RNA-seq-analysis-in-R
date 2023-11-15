# RNA-seq-analysis-in-R

Differential expression analysis
30 November 2020
Authors: Belinda Phipson, Anna Trigos, Matt Ritchie, Shian Su, Maria Doyle, Harriet Dashnow, Charity Law

Data files and Resources
Data files are available from: https://figshare.com/s/1d788fd384d33e913a2a. You should download the files listed below and place them into a folder called data in your working directory.

Data files:

GSE60450_Lactation-GenewiseCounts.txt
SampleInfo.txt
SampleInfo_Corrected.txt
Data files were originally obtained from:
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz
http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.all.v7.1.entrez.rds
http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.h.all.v7.1.entrez.rds

This material has been inspired by the following resources:
http://www.statsci.org/smyth/pubs/QLedgeRPreprint.pdf (Lun, Chen, and Smyth 2016)
http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html

____________Bioconductor/R Packages_____________________
Packages used:

limma
edgeR
Glimma
org.Mm.eg.db
gplots
RColorBrewer
