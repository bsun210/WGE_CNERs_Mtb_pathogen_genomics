# A hybridization capture approach for pathogen genomics

This repository contains scripts used to analyze data and to make figures reported in Sundararaman et al, mBio (2023b).

We developed a cost-effective method named Circular Nucleic acid Enrichment Reagent synthesis (CNERs) to generate whole-genome enrichment probes. We demonstrated the method by producing probes for Mycobacterium tuberculosis which we used to enrich M. tuberculosis DNA that had been spiked at concentrations as low as 0.01% and 100 genome copies against human DNA background to 1225-fold and 4636-fold. 

## Data analyses
* mtb_fq2counts.sh:
  This bash scripts processes paired-end sequencing data and creates a summary file. Takes raw fastq files, trims adapters using Cutadapt (https://cutadapt.readthedocs.io/en/stable/). Removes low complexity reads from merged and unmerged reads using prinseq (https://github.com/uwb-linux/prinseq). Maps merged reads as SE and unmereged read-pairs as PE mode using bwa aln (http://bio-bwa.sourceforge.net/bwa.shtml). Bam files merged, sorted and indexed using samtools (http://www.htslib.org/doc/samtools.html). Duplicates marked and removed using Picard MarkDuplicates (https://broadinstitute.github.io/picard/). Individual SNP coverage depth and coverage around the SNP region are collected using bedtools (https://bedtools.readthedocs.io/en/latest/index.html).

M. tuberculosis H37Rv reference genome (NC_000962.3) 
