# ESSIS (Effector Somatic SNV and INDEL Seeker)
In silico pipeline for predicting functionally significant somatic single-nucleotide variants and indels in multiple regions of sample genomes.
<img src="images/ESSIS_logo.png" alt="ESSIS workflow" width="400" height="400"/>
## Table of Contents

● [Introduction](#introduction)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;○ [What is ESSIS?](#what-is-essis)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;○ [How does it work in summary?](#how-does-it-work-in-summary)

● [Dependencies and Installation](#dependencies-and-installation)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;○ [Install dependencies via conda](#install-dependencies-via-conda)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;○ [Install auxiliary dependencies](#install-auxiliary-dependencies)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;○ [Install ESSIS](#install-essis)

● [Input](#input)  
● [Usage](#usage)  
● [Output](#output)

## Introduction  
### What is ESSIS?
ESSIS is a bioinformatic pipeline that takes somatic mutations stored as a GRanges object in R (in Rds format) and predicts their functional significance (effectors) based on the scores of parameters representing known, biologically defined mechanisms—such as truncations of proteins caused by nonsense mutations or disruptions of transcription factor-binding motif sites induced by promoter mutations. ESSIS can process mutations from any type of high-throughput next-generation sequencing platform, including **whole-genome (WGS), whole-exome (WXS), targeted** and **bulk RNA sequencing**. It can also annotate mutations within both coding and non-coding genomic regions and, if a reference GTF file is provided, intergenic regions as well. Moreover, it can operate on your local computer or on a computational cluster in both multiprocessing and multithreaded modes.

### How does it work in summary?  
Pipeline ESSIS is split into six main steps:  
I.  Parse genomic elements from the mandatory file (in Rds format) provided by the user (see format in `examples/`).  
II. Identify regions with an excess number of mutations by running fishHook[^1] on each tumor-cohort + mutation-class combination.  
III. Rank all mutations in each hypermutated region and flag the highest-ranked ones as effectors based on the excess mutation count and specified biological criteria.  
IV. Merge all non-flagged mutations across cohorts and rerun fishHook separately for each mutation class to catch excess mutations in regions with low statistical power.  
V. Re-rank and flag effectors in these newly identified regions.  
VI. Apply a rule-based filter to flag any remaining mutations satisfying specific biological criteria.  

For a detailed description of biological criteria applied to each element type, please see [documentation](#).






