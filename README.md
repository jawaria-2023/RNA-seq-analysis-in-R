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

_____________Bioconductor/R Packages_____________________
Packages used:

1. **limma**
2. **edgeR**
3. **Glimma**
4. **org.Mm.eg.db**
5. **gplots**
6. **RColorBrewer**

_________________Overview________________________________
1. **Reading in table of counts**
2. **Adding annotation**
3. **Filtering lowly expressed genes**
4. **Quality control**
5. **Normalization for composition bias**
6. **Differential expression analysis**
7. **Testing relative to a threshold**
8. **Visualization**
9. **Gene set testing**

Introduction and data import: 
Over the past two decades, the widespread practice of measuring gene expression on a genome-wide scale has evolved, primarily utilizing microarrays before 2008. The landscape shifted with the introduction of next-generation sequencing technology in 2008, leading to a growing preference among scientists for utilizing this technology to explore and comprehend changes in gene expression within complex biological systems. The decreasing costs of sequencing have made RNA-Seq an increasingly accessible method, allowing researchers to efficiently measure the expression of tens of thousands of genes across multiple samples. Consequently, the focus has shifted from the generation of data to the challenges of storing and analyzing the vast datasets produced.

The analysis of an RNA-Seq experiment involves several key steps. The process commences with the sequencing of reads, which are then aligned to a reference genome. Subsequently, the counts of reads mapped to each gene are determined, resulting in a table of counts. This count data serves as the foundation for statistical analyses conducted in R. While mapping and counting are crucial tasks in the analysis pipeline, our focus today is on delving into the analysis phase starting directly from the count data.











