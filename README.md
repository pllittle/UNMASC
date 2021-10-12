# UNMASC - Unmatched Normals and Mutant Allele Status Characterization
This package is designed to filter and annotate tumor-only variant calls through the integration of public database annotations, clustering, and segmentation to provide the user with a clear characterization of each variant when called against a set of unmatched normal controls.

## Citation
Little, L., Jo, H., Hoyle, A., Mazul, A., Zhao, X., Salazar, A.H., Farquhar, D., Sheth, S., Masood, M., Hayward, M.C., Parker, J.S., Hoadley, K.A., Zevallos, J. and Hayes, D.N. (2021). UNMASC: tumor-only variant calling with unmatched normal controls. *NAR Cancer*, 3(4), zcab040. [[HTML](https://academic.oup.com/narcancer/article/3/4/zcab040/6382329), [PDF](https://academic.oup.com/narcancer/article-pdf/3/4/zcab040/40514892/zcab040.pdf)]

## Installation
```
library(devtools)

# Check dependencies
all_packs = as.character(installed.packages()[,1])
req_packs = c("Rcpp","RcppArmadillo","emdbook",
	"scales","seqTools","parallel","doParallel",
	"data.table","grDevices","Rsamtools","GenomicRanges",
	"IRanges","foreach")
if( !all(req_packs %in% all_packs) ){
	miss_packs = req_packs[!(req_packs %in% all_packs)]
	stop(sprintf("Missing package(s): %s",paste(miss_packs,collapse = ",")))
}

# UNMASC
if( !("UNMASC" %in% all_packs) )
	devtools::install_github("pllittle/UNMASC")
```

## UNMASC inputs

* Annotated variant calls (e.g. Strelka/Strelka2 + VEP)
* Target capture bed file
* Centromere start/end bed file
* Tumor bam filename

## Workflow template code for UNMASC inputs

UNMASC's benchmark samples were run with Strelka. Assuming Strelka/Strelka2, GATK, and VEP are installed along with corresponding dependencies, Linux commands are provided below to run Strelka, Strelka2, GATK and VEP for variant calling and annotation. Running our VEP annotation requires downloading a COSMIC database VCF. For example, CosmicCodingMuts.vcf.gz for GRCh37 with the latest release can be found at [here](https://cancer.sanger.ac.uk/cosmic/download?genome=37). We've instructed VEP to annotate variants with 1000 Genomes population allele frequencies, ExAC/gnomAD population allele frequencies, variant transcripts, impacts/consequences, and COSMIC counts with legacy IDs.

```
# ----------
# Directories and filenames
# ----------
stk_dir=<strelka directory>
out_dir=<output directory>
stk2_dir=<strelka2 directory>
gatk_dir=<GATK directory>
vep_dir=<VEP directory>

nbam=<control bam filename>
tbam=<tumor bam filename>
fasta=<reference fasta filename>
nthreads=<number of threads>
genome=GRCh37
cosmic_fn=<COSMIC coding variants VCF>

# ----------
# Strelka commands for a tumor and control pair.
# ----------
$stk_dir/bin/configureStrelkaWorkflow.pl \
	--normal=$nbam --tumor=$tbam --ref=$fasta \
	--config=$out_dir/config.ini \
	--output-dir=$out_dir

make -C $out_dir/ -j $nthreads

# ----------
# Strelka2 commands for a tumor and control pair.
# ----------
$stk2_dir/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam $nbam --tumorBam $tbam --referenceFasta $fasta \
	--disableEVS --exome --runDir $out_dir

$out_dir/runWorkflow.py -m local -j $nthreads

# ----------
# GATK command to merge SNV and INDEL VCFs.
# ----------
$gatk_dir/gatk MergeVcfs \
	--INPUT $out_dir/results/variants/somatic.snvs.vcf.gz \
	--INPUT $out_dir/results/variants/somatic.indels.vcf.gz \
	--OUTPUT $out_dir/somatic.vcf.gz

# ----------
# Variant Effect Predictor (VEP)
# ----------
vep_fields=IMPACT,Consequence,SYMBOL,HGVSc,HGVSp,AF
vep_fields="$vep_fields,gnomAD_AF,COSMIC,COSMIC_CNT"
vep_fields="$vep_fields,COSMIC_LEGACY_ID"

zcat $out_dir/somatic.vcf.gz > $out_dir/somatic.vcf

$vep_dir/vep --format vcf --species homo_sapiens \
	-i $out_dir/somatic.vcf -o $out_dir/somatic.vep.vcf \
	--fork $nthreads --cache --dir_cache $vep_dir \
	--cache_version 104 --assembly $genome --fasta $fasta \
	--force_overwrite --no_stats --domains --hgvs --af \
	--af_gnomad --vcf \
	--custom $cosmic_fn,COSMIC,vcf,exact,0,CNT,LEGACY_ID \
	--fields "$vep_fields"

```

## Simulated Example

## Output

The column definitions of the sample output `tumorOnly_VCs.tsv` are described using the code below.

```
library(UNMASC)
print(UNMASC::readme_VC,right = FALSE)
```

## Vignette
To open the package vignette, run ```vignette("UNMASC")```.

