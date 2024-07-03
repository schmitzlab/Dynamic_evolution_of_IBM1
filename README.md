# Introduction
Here are the scripts used for data analysis in the research study titled 'Dynamic evolution of the heterochromatin sensing histone demethylase IBM1'.

# Input dataset
The input files required to run the scripts mentioned below can be accessed here:

# Scripts


## Analysis for arabidopsis population

`extract.1kb.intronseq.pl`   This script extract intron sequence based on given location information.

`scan1kb_get.methylation.pl`  This script scans for mCHG-enriched introns longer than 1kb and quantifying methylation level on introns.

`IBM1_isoform_mutant.r`  This script visualizes IBM1 isoform quantification results for mutant lines across various public datasets.

`IBM1_Methylation_and_Expression_Analysis.R`  This script performs an integrated analysis of IBM1 gene methylation and expression across multiple datasets. It reads, preprocesses, and merges methylation data with gene expression levels, subsequently analyzing correlations and generating visualizations to uncover correlation between IBM1 intron methylation and isoform expression. The results are meticulously plotted and saved.

## Pipeline for generating IBM1 gene tree and intron structure

`check.match.pl` This script verifies the alignment of amino acid sequences with their corresponding coding sequences (CDS) from input files. It checks if the amino acids match the reading frames of the CDS, filters out sequences with mismatches or errors, and generates statistics about the sequence verification process. 

## IBM1 intron sequence annotation
 
