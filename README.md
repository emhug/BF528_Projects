# BF528_Projects
Example code scripted as part of the BF528: Applications of Translational Bioinformatics (Spring 2020, Dr. Adam Labadorf) course.
These projects were completed as part of a group, where each member was assigned a different role in the analysis. The roles were: Data Curator, Programmer, Analyst, and Biologist. 
  Data curator: identify, download, and describe relevant datasets and literature
  Programmer: write code to transform the downloaded data into an interpretable form
  Analyst: examine and visualize processed data to aid in interpretation
  Biologist: connect the processed data into meaningful biological interpretation

## Project 3 files - Analyst role
This project aimed to replicate Wang, Charles, Binsheng Gong, Pierre R. Bushel, Jean Thierry-Mieg, Danielle Thierry-Mieg, Joshua Xu, Hong Fang, et al. 2014. “A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data” Nature Biotechnology 32 (9): 926–32. PMID: 4243706.

### run_limma.R
This file uses microarray results and the tox-group information as input to perform a differential expression analysis using limma. The output is three CSV files that contain information about the differential expression of each probe.

### make_histograms_significant_genes.R
This file uses the three CSV files from the run_limma.R script and a CSV file of the probes mapping to genes as input. This code generates histograms of the log fold change values of the significant genes for each chemical. This code also generates volcano plots of the significant genes for each chemical.

### compare_limma_deseq.R
This file uses the same input files as make_histograms_significant_genes.R and the DESeq files for each chemical. This code maps the probes in the limma output to the refseqids in the DESeq output and calculates the concordance between the two files. This concordance is then plotted as a function of the number of differentially expressed genes in each analysis. The concordance values for above- and below-median expression genes was also calculated. These results were plotted on a bar graph that displayed the overall, above- and below-median concordance values for each chemical.

### Output files
- DEGs_scatterplot.png
- chemical_histograms.png
- concordance_barplot.png

## Project 5 files - Biologist role
This project aimed to replicate O’Meara et al. Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration. Circ Res. Feb 2015. PMID: 25477501. The purpose of the biologist role was to analyze the FPKM (fragments per kilobase of exon model per million reads mapped) matrices to biologically interpret the experiments. 

### biologist_project5.R
This file uses the FPKM matrix information to generate a heatmap and retrieve the top 500 differentially expressed genes.

### Input files
- Ad_1.fpkm_tracking
- Ad_2.fpkm_tracking
- P0_2.fpkm_tracking
- P4_1.fpkm_tracking
- P4_2.fpkm_tracking
- P7_1.fpkm_tracking
- P7_2.fpkm_tracking
- gene_exp.diff

### Output files
- heatmap.png
- fpkm_through_time.png
