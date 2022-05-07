<div align="left">
<a href=""><img src="https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink&label=UNMASC" height="100" /></a>
</div>

<!-- badges: start -->
![C++](https://img.shields.io/badge/C++-%2300599C.svg?style=square&logo=c%2B%2B&logoColor=gold)
![R](https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink)
![CRAN status](https://www.r-pkg.org/badges/version/UNMASC)
[![DOI](https://zenodo.org/badge/DOI/10.1093/narcan/zcab040.svg)](https://doi.org/10.1093/narcan/zcab040)
<!-- badges: end -->

## What is this for?

One goal of cancer genomics is to identify DNA variants specific to the cancer tissue within an individual. Perhaps a researcher would like to identify mutated genes and design a cancer treatment or therapy specific to that individual's cancer. These cancer variants are considered **somatic** or variants that cannot be inherited. Our normal tissue harbors inherited DNA variants called  **germline** variants that are present and identical across all normal tissue. 

If one sequences an individual's matched normal DNA (e.g. from blood or adjacent tissue) and tumor DNA, one can identify both **germline** and **somatic** mutations and more importantly, distinguish between them. However, without the matched normal DNA serving as a control, the performance of somatic mutation callers (MuTect2, Seurat, Indelocator, Varscan, Strelka, Strelka2, etc.) drops off in terms of recall (sensitivity) and precision (positive predictive value). Perhaps the tumor sample:

* lacks an available matched normal (e.g. patient is unavailable, has leukemia)
* sample contamination, poorly sequenced normal, 
* insufficient budget to sequence both samples per patient

A third set of detected and unavoidable variants are false positives or **artifacts** that can arise from several sources including poor sequencing, sample storage, read misalignment to the reference genome, etc. UNMASC attempts to identify somatic variants from tumor samples without an adequate matched normal.

<p align="center">
<img src="images/workflow.PNG" width="70%" />
<p align="center"><em>UNMASC workflow for a single tumor sample against Z unmatched normal controls. SB = strand bias, SEG = segmentation, OXOG = oxoG artifacts, FFPE = paraffin artifacts.</em></p>
</p>

## Description

This package is designed to filter and annotate tumor-only variant calls through the integration of public database annotations, clustering, and segmentation to provide the user with a clear characterization of each variant when called against a set of unmatched normal controls.

## Citation
Little, P., Jo, H., [Hoyle, A.](https://github.com/alanhoyle), [Mazul, A.](https://surgery.wustl.edu/people/angela-mazul/), [Zhao, X.](https://github.com/xiaobeizhao), Salazar, A.H., [Farquhar, D.](https://www.med.unc.edu/ent/directory/douglas-farquhar-md-mph/), [Sheth, S.](https://www.med.unc.edu/medicine/directory/siddharth-sheth-md/), [Masood, M.](https://www.linkedin.com/in/maheer-masood-03661685), Hayward, M.C., Parker, J.S., [Hoadley, K.A.](https://unclineberger.org/directory/katherine-hoadley/), [Zevallos, J.](https://www.linkedin.com/in/jpzevallosmd) and [Hayes, D.N.](https://hayeslab.lab.uthsc.edu/) (2021). UNMASC: tumor-only variant calling with unmatched normal controls. *NAR Cancer*, 3(4), zcab040. [[HTML](https://academic.oup.com/narcancer/article/3/4/zcab040/6382329), [PDF](https://academic.oup.com/narcancer/article-pdf/3/4/zcab040/40514892/zcab040.pdf), [Supplement](https://academic.oup.com/narcancer/article/3/4/zcab040/6382329#supplementary-data)]

## Installation

R/RStudio code to check, install, and load libraries.

```R
all_packs = as.character(installed.packages()[,1])
pandoc = Sys.getenv("RSTUDIO_PANDOC")
build_vign = !is.null(pandoc) && file.exists(pandoc)

if( !("smarter" %in% all_packs) ){
	stop("Check https://github.com/pllittle/smarter for installation")
}

library(smarter)
cran_packs = c("devtools","Rcpp","RcppArmadillo","emdbook",
	"scales","BiocManager","parallel","doParallel",
	"data.table","grDevices","foreach")
bioc_packs = c("seqTools","Rsamtools","GenomicRanges",
	"IRanges")

smarter::smart_packDeps(
	cran_packs = cran_packs,
	bioc_packs = bioc_packs,
	github_packs = c("pllittle/UNMASC"),
	pandoc = pandoc,
	build_vign = build_vign)

```

## UNMASC inputs

* annotated variant calls (e.g. Strelka/Strelka2 + VEP)
* target capture bed file: contains contig, start position, end position columns
* centromere start/end bed file
* dict_chrom file: Run `samtools view -H tumor.bam` and save the output.
* tumor bam filename

## Workflow template code for UNMASC inputs

UNMASC's benchmark samples were run with Strelka. Assuming 

* [Strelka](https://sites.google.com/site/strelkasomaticvariantcaller/)
* [Strelka2](https://github.com/Illumina/strelka), 
* [GATK](https://github.com/broadinstitute/gatk), and
* [VEP](https://uswest.ensembl.org/info/docs/tools/vep/index.html) 

are installed along with corresponding dependencies (Perl, HTSlib, etc.), 
Linux commands are provided below to run these software for variant calling 
and annotation. Running our customized VEP annotation requires downloading 
a COSMIC database VCF. For example, CosmicCodingMuts.vcf.gz for GRCh37 with 
the latest release can be found at 
[here](https://cancer.sanger.ac.uk/cosmic/download?genome=37). 
We have instructed VEP to annotate variants with 1000 Genomes population 
allele frequencies, ExAC/gnomAD population allele frequencies, variant transcripts, 
impacts/consequences, and COSMIC counts with stable and legacy IDs.

Refer to our documentation for steps.

## Prepare UNMASC's annotated VCF

Code below will import merged Strelka2 vcfs, VEP annotation vcf, and target bed file 
to prepare UNMASC's main `vcf` input R data.frame.

```R
vcf = prep_UNMASC_VCF(outdir,DAT,FILTER,
	target_fn,anno_fn,ncores)
```

## Run UNMASC

Template code to run UNMASC

```R
# Example Inputs
tumorID = "tumor01"
outdir 	= file.path(".",tumorID)
vcf     = <a data.frame of required column names, refer to documentation>
tBAM_fn = "path/to/tumor/bam"
bed_centromere_fn = "path/to/centromere/start/end/bed/file"
dict_chrom_fn = "path/to/chromosome/length/file"
qscore_thres = 30 # Qscore threshold
exac_thres = 5e-3 # ExAC population allele frequency threshold for germline filtering
ad_thres = 5
rd_thres = 10
cut_BAF = 5e-2
minBQ = 13
minMQ = 40
eps_thres = 0.5
psi_thres = 0.02
hg = "19"
binom = TRUE
gender = NA
ncores = 1 # number of threads/cores available

# Package main function
UNMASC::run_UNMASC(tumorID = tumorID,outdir = outdir,
	vcf = vcf,tBAM_fn = tBAM_fn,bed_centromere_fn = bed_centromere_fn,
	dict_chrom_fn = dict_chrom_fn,qscore_thres = qscore_thres,
	exac_thres = exac_thres,ad_thres = ad_thres,rd_thres = rd_thres,
	cut_BAF = cut_BAF,minBQ = minBQ,minMQ = minMQ,eps_thres = eps_thres,
	psi_thres = psi_thres,hg = "19",binom = binom,gender = gender,
	ncores = ncores)
```

<!---
## Output

The column definitions of the sample output `tumorOnly_VCs.tsv` are described using the code below.

```
library(UNMASC)
print(UNMASC::readme_VC,right = FALSE)
```

## Vignette
To open the package vignette, run ```vignette("UNMASC")```.

-->

## Future directions

* Workflow containers
* Shell script functions for dependency installations to variant calling to UNMASC
* Add MuTect/MuTect2 sample code
* Applying UNMASC toward circulating plasma tumor cell DNA
* Identifying somatic mutations missed by tumors with matched normals

## FAQs



