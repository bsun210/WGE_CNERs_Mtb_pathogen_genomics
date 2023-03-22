# A hybridization capture approach for pathogen genomics

This repository contains scripts used to analyze data and to make figures reported in Sundararaman et al, mBio (2023b).

We developed a cost-effective method named Circular Nucleic acid Enrichment Reagent synthesis (CNERs) to generate whole-genome enrichment probes. We demonstrated the method by producing probes for Mycobacterium tuberculosis which we used to enrich M. tuberculosis DNA that had been spiked at concentrations as low as 0.01% and 100 genome copies against human DNA background to 1225-fold and 4636-fold. 

## Data analyses
* mtb_fq2counts.sh:

This bash scripts processes paired-end sequencing data and creates a summary file. Takes raw fastq files, trims adapters using Cutadapt (https://cutadapt.readthedocs.io/en/stable/). Maps reads as PE mode using bwa mem (http://bio-bwa.sourceforge.net/bwa.shtml) to M. tuberculosis H37Rv reference genome (NC_000962.3). Bam files merged, duplicates removed, sorted and indexed using samtools (http://www.htslib.org/doc/samtools.html). Mapped read are counted using samtools stats and genome coverage deteremined using bedtools (https://bedtools.readthedocs.io/en/latest/index.html).

*mtb_100bp_cov_multiSamp_corr.py:

This python script is used to make Figures S3, S9-11. Takes 100bp raw coverage tsvs for multiple samples as comma separated list to make correlation matrix of normalized coverage across 100bp bins of the H37Rv reference genome. 

*mtb_genomCov.py:
This python script is used to make Figures 2A, B, S8. Takes a tsv with genome coverage at each base for multiple samples (10000, 1000, 100 copy capture experiments) and plots histogram of cumulative genome coverage.  

