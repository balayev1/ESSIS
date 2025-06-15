# ESSIS (Effector Somatic SNV and INDEL Seeker)

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
ESSIS is a bioinformatic pipeline that takes somatic mutations stored as a GRanges object in R (in Rds format) and predicts their functional significance (effectors) based on the scores of parameters representing known, biologically defined mechanisms—such as protein truncations caused by nonsense mutations or disruptions of transcription factor-binding motif sites induced by promoter mutations. ESSIS can process mutations from any type of high-throughput next-generation sequencing platform, including **whole-genome (WGS), whole-exome (WXS), targeted** and **bulk RNA sequencing**. It can also annotate mutations within both coding and non-coding genomic regions and, if a reference GTF file is provided, intergenic regions as well. Moreover, it can operate on your local computer or on a computational cluster in both multiprocessing and multithreaded modes.

### How does it work in summary?  
Pipeline ESSIS is split into six main steps:  
I.  Parse genomic elements from the mandatory file (in Rds format) provided by the user (see format in `data/`).  
II. Identify regions with an excess number of mutations by running [fishHook](https://github.com/mskilab-org/fishHook) on each tumor-cohort + mutation-class combination.  
III. Rank all mutations in each hypermutated region and flag the highest-ranked ones as effectors based on the excess mutation count and specified biological criteria.  
IV. Merge all non-flagged mutations across cohorts and rerun fishHook separately for each mutation class to catch excess mutations in regions with low statistical power.  
V. Re-rank and flag effectors in these newly identified regions.  
VI. Apply a rule-based filter to flag any remaining mutations satisfying specific biological criteria.  

For a detailed description of biological criteria applied to each element type, please see [documentation](#).

## Dependencies and Installation
Pipeline ESSIS was implemented on Linux operating system. You’ll need both a conda environment and several R packages installed via `devtools::install_github()`.

### Install dependencies via conda
We provide `environment.yml` for reproducibility. First, make sure conda is installed. From your terminal, run:

```bash
if ! command -v conda >/dev/null 2>&1; then
  echo "Error: Conda is not installed. Please install [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) or load conda environment on a cluster"
  exit 1
fi
```

To create and activate the environment, run:
```bash
conda env create -f environment.yml    # create env with all required packages
conda activate essis_env               # activate the ESSIS environment
```

### Install auxiliary dependencies
Open an R session (or RStudio) and run:
```r
install.packages("devtools")
devtools::install_github("mskilab/fishHook")
devtools::install_github("mskilab/gUtils")
devtools::install_github("mskilab/gTrack")
devtools::install_github("mskilab/gChain")
devtools::install_github("mskilab/skitools")
```

### Install ESSIS
```bash
git clone https://github.com/balayev1/ESSIS.git  # clone the repo
cd ESSIS # change dir to ESSIS/
```

**Warning:** Follow each installation step sequentially! Always run ESSIS from within the `ESSIS/` directory.

## Input

Required input files:
```bash
tumor	class	mut_path	ave_read_depth
PAAD	SNV	data/PAAD.sSNV.ann.Rds	data/Ave_RD_by1kbwindows_paad.Rds
PAAD	INDEL	data/PAAD.sINDEL.ann.Rds	data/Ave_RD_by1kbwindows_paad.Rds
```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;● A GRanges object in Rds format with genomic regions and required column `type`. These are first lines of the file with genomic element regions:
```bash
GRanges object with 889848 ranges and 3 metadata columns:
           seqnames            ranges strand |              gene_name
              <Rle>         <IRanges>  <Rle> |            <character>
       [1]     chr1       63419-65418      + |                  OR4F5
       [2]     chr1       65419-65433      + |                  OR4F5
       [3]     chr1       65434-65519      + |                  OR4F5
       [4]     chr1       65520-65564      + |                  OR4F5
       [5]     chr1       65565-65573      + |                  OR4F5
       ...      ...               ...    ... .                    ...
  [889844]     chrY 21820455-21820734      * | OFD1P16Y,ENSG0000028..
  [889845]     chrY 26315142-26315389      * |             PPP1R12BP1
  [889846]     chrY 26441509-26441720      * |             intergenic
  [889847]     chrY 26563051-26563261      * |                TPTE2P4
  [889848]     chrY 56685452-56685708      * |             intergenic
                        gene_type        type
                      <character> <character>
       [1]         protein_coding    PROMOTER
       [2]         protein_coding        5UTR
       [3]         protein_coding    PROMOTER
       [4]         protein_coding        5UTR
       [5]         protein_coding         CDS
       ...                    ...         ...
  [889844] unprocessed_pseudoge..        CTCF
  [889845] unprocessed_pseudogene        CTCF
  [889846]                     NA        CTCF
  [889847] unprocessed_pseudogene        CTCF
  [889848]                     NA        CTCF
```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;●Output directory to store results
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;●Column for sample ID
**Warning:** Input mutation GRanges object must contain columns:
● REF: denotes reference allele
● ALT: denotes alternative allele
● CONSEQUENCE: comma-separated list of ALT consequences (see `data/PAAD.sSNV.ann.Rds`)
● NUM_TOOLS: number of ALT-finding tools 
● LC: linguistic complexity within 80 base-pair vicinity of ALT
● CADD_scaled: scaled CADD score for ALT

## Usage
```bash
Usage: 
Rscript essis.R --ge-file GE-FILE --mut-file MUT-FILE --output-dir OUTPUT-DIR  --id-col IDCOL [options]
GE-FILE: full path to GRanges Rds obj with genomic element regions and types
MUT-FILE: full path to tab-delimited file with sample info on tumor, mutation class and GRanges Rds obj with annotated mutations
OUTPUT-DIR: directory to store results
ID-COL: sample ID column, used as fishHook covariate

Options:
--gtf-path=GTF-PATH
       full path to GTF file, considers all regions outside of types in GE-FILE as 'intergenic'
--repli=REPLI
       full path to file with replication times per genomic region, used as fishHook covariate
--dinuc=DINUC
       full path to file with dinucleotide percentages per genomic region, used as fishHook covariate
--trinuc=TRINUC
       full path to file with trinucleotide percentages per genomic region, used as fishHook covariate
--id-cap=ID-CAP
       upper mutation limit to consider per genomic region per sample, used as fishHook covariate
--cores=CORES
       number of utilized threads, used as fishHook covariate
--tile-size=TILE-SIZE
       tile size, used as fishHook covariate
--tasks=TASKS
       number of processes, used fo fishHook run (not recommended to set > 1)
--min-cadd=MIN-CADD
      minimum scaled CADD score to consider mutation for effector status
--ranking-types=RANKING-TYPES
      comma-separated list of element types to rank
--coding-types=CODING-TYPES
      comma-separated list of coding element types; requires binary column LOFTEE (HC/NC) in mutation file
--regulatory-types=REGULATORY-TYPES
      comma-separated list of regulatory element types; requires columns motif_gain_status and motif_break_status in mutation file with factors 0 and 1.
--mirna_binding-types=MIRNA_BINDING-TYPES
      comma-separated list of miRNA-binding element types; requires column mirBS_status in mutation file with factors 0 and 1.
--mirna-type=MIRNA-TYPE
      miRNA label
--ctcf-type=CTCF-TYPE
      CTCF-binding motif label; requires columns motif_gain_status and motif_break_status with factors 0 and 1, and column motif_name_break with motif name in mutation file
--intron-types=INTRON-TYPES
      comma-separated list of intronic element types to identify splice sites; requires column spliceSiteMod_status in mutation file with factors 0 and 1.
```

## Output
ESSIS writes parsed genomic element Rds files to `OUTPUT-DIR/parsed_elements`, fishHook results to `OUTPUT-DIR/fishHook` and mutations with effector flags to `OUTPUT-DIR/final_results`.
Final mutation files contain same columns as input Rds files with additional columns `effector_flag_byrank`, `effector_flag_byrule` and `effector_flag` that combines effectors from both rank-based and rule-based steps.


