# Introduction
Here are the scripts used for data analysis in the research study titled 'Dynamic evolution of the heterochromatin sensing histone demethylase IBM1'.

# Input dataset
The input files required to run the scripts mentioned below can be accessed here: [https://osf.io/gk6tn/files/osfstorage](https://osf.io/gk6tn)

# Scripts


## Analysis for arabidopsis population

`extract.1kb.intronseq.pl`   This script extract intron sequence based on given location information.

`scan1kb_get.methylation.pl`  This script scans for mCHG-enriched introns longer than 1kb and quantifying methylation level on introns.

`IBM1_isoform_mutant.r`  This script visualizes IBM1 isoform quantification results for mutant lines across various public datasets.

`IBM1_Methylation_and_Expression_Analysis.R`  This script performs an integrated analysis of IBM1 gene methylation and expression across multiple datasets. It reads, preprocesses, and merges methylation data with gene expression levels, subsequently analyzing correlations and generating visualizations to uncover correlation between IBM1 intron methylation and isoform expression. The results are meticulously plotted and saved.

## Pipeline for generating IBM1 gene tree and intron structure

`check.match.pl` This script verifies the alignment of amino acid sequences with their corresponding coding sequences (CDS) from input files. It checks if the amino acids match the reading frames of the CDS, filters out sequences with mismatches or errors, and generates statistics about the sequence verification process.   

`convertto.phylip.pl` This script converts sequence alignment files from Gblock to Phylip format, facilitating the construction of maximum likelihood trees using RAxML.

`frompep.to.codon.pl` This script converts peptide alignments into codon alignments.

`extract.intron.seq.species.pl` This script extracts specific intron sequences from a given genome assembly using gene IDs and location information from GFF annotations.

`covert.intron.location.pl` This script identifies the relative coordinates of targeted introns, exons, and protein domains of interest based on provided GFF files.

`methylation_structurev2.0.r` This script generates a gene tree for IBM1 and plots intron methylation data.

## IBM1 intron sequence annotation
 
`extract.snp.seq.longintron.pl` This script extracts SNPs from specific intron regions using gene and CDS annotations in GFF files.  

`GC.calculation.r` This script calculates the GC content of specified intron sequences.  

`snp.diversity.r` This script calculates SNP density and Tajima's D for specified intron regions.  

`species.ex.heatmap.r` This R script is designed to process and analyze gene expression data across multiple species. It reads in gene expression data and gene body methylation (gbM) counts from text files for two main visualizations: a heatmap that displays gene expression levels and a scatter plot that shows the ratio of gbM genes.   

`species.ex.methyl.r` This R script analyzes the relationship between gene expression and intron methylation status across various species. The script reads and merges relevant datasets, then performs statistical comparisons and visualizations  

`summaryTE.v2.r` This R script performs an integrated analysis combining gene expression data, intron methylation status, and transposable element (TE) annotations across multiple species to explore potential correlations between these genomic features
