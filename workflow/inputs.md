# Preparing UNMASC Inputs

## Outline

* [Shell Directories and Variables](https://github.com/pllittle/UNMASC/blob/main/workflow/inputs.md#setting-directories-and-variables)
* [Sourcing Workflow Scripts](https://github.com/pllittle/UNMASC/blob/main/workflow/inputs.md#source-scripts)
* [Custom Local Installations](https://github.com/pllittle/UNMASC/blob/main/workflow/inputs.md#custom-local-installations)
* [Get COSMIC vcf](https://github.com/pllittle/UNMASC/blob/main/workflow/inputs.md#get-cosmic-vcf)
* [TO_workflow](https://github.com/pllittle/UNMASC/blob/main/workflow/inputs.md#to_workflow)
* [Constructing UNMASC's `vcf`](https://github.com/pllittle/UNMASC/blob/main/workflow/inputs.md#preparing-unmascs-vcf)
* [Execution](https://github.com/pllittle/UNMASC/blob/main/workflow/inputs.md#execution)
* [Outputs](https://github.com/pllittle/UNMASC/blob/main/workflow/inputs.md#outputs)

## Setting Directories and Variables

<details>
<summary>Click to expand!</summary>

Below are fixed variables to specify.

```Shell
gatk_dir=; [ -z "$gatk_dir" ] \
	&& echo "Set gatk_dir, GATK directory" >&2 \
	&& return 1

git_dir=; [ -z "$git_dir" ] \
	&& echo "Set git_dir, location to store GitHub repos" >&2 \
	&& return 1
[ ! -d $git_dir ] && mkdir $git_dir

stk2_dir=; [ -z "$stk2_dir" ] \
	&& echo "Set stk2_dir, location of Strelka2 dir!" >&2 \
	&& return 1

vep_dir=; [ -z "$vep_dir" ] \
	&& echo "Set vep_dir, VEP directory" >&2 \
	&& return 1

vep_rel=; [ -z "$vep_rel" ] \
	&& echo "Set vep_rel, VEP release number, preferably 105 with GRCh37 and GRCh38 supported" >&2 \
	&& return 1

vep_cache=; [ -z "$vep_cache" ] \
	&& echo "Set vep_cache, VEP homo_sapiens cache, should be vep, refseq, or merged" >&2 \
	&& return 1

hts_dir=; [ -z "$hts_dir" ] \
	&& echo "Set hts_dir, the HTS directory" >&2 \
	&& return 1

cosm_dir=; [ -z "$cosm_dir" ] \
	&& echo "Set cosm_dir, directory containing COSMIC vcf" >&2 \
	&& return 1

cosm_ver; [ -z "$cosm_ver" ] \
	&& echo "Set cosm_ver, COSMIC version number, e.g. 95" >&2 \
	&& return 1

fasta_fn=; [ -z "$fasta_fn" ] \
	&& echo "Set fasta_fn, the reference FASTA" >&2 \
	&& return 1

nthreads=; [ -z "$nthreads" ] \
	&& echo "Set nthreads, number of threads or cores" >&2 \
	&& return 1

genome=; [ -z "$genome" ] \
	&& echo "Set genome, like GRCh37 or GRCh38" >&2 \
	&& return 1

```

Sample-specific Variables

```Shell
out_dir=; [ -z "$out_dir" ] \
	&& echo "Set out_dir, output directory" >&2 \
	&& return 1

nbams=; [ -z "$nbams" ] \
	&& echo "Set nbams, file listing all normal bam full paths" >&2 \
	&& return 1
echo -e "Detected $(cat $nbams | wc -l) normal controls." >&2

tbam=; [ -z "$tbam" ] \
	&& echo "Set tbam, tumor bam full path" >&2 \
	&& return 1

```

</details>

## Source scripts

<details>
<summary>Click to expand!</summary>

You are welcome to install your own programs and dependencies. I have
provided below the steps and scripts I use to automate the programs 
related to UNMASC.

```Shell
# Pull my functions

cd $git_dir
[ ! -d baSHic ] && git clone https://github.com/pllittle/baSHic.git >&2
[ -d baSHic ] && cd baSHic && git pull >&2

cd $git_dir
[ ! -d UNMASC ] && git clone https://github.com/pllittle/UNMASC.git >&2
[ -d UNMASC ] && cd UNMASC && git pull >&2

# Source functions

. $git_dir/baSHic/scripts/genomic.sh
[ ! $? -eq 0 ] && echo "Some error in sourcing genomic.sh" >&2 && return 1

. $git_dir/UNMASC/workflow/tumor_only.sh
[ ! $? -eq 0 ] && echo "Some error in sourcing tumor_only.sh" >&2 && return 1

```

</details>

## Custom local installations

<details>
<summary>Click to expand!</summary>

* gcc, libtool, perl, bzip2, xz, zlib, curl, expat, db,
* HTSlib, VEP, Strelka2

One can control where these programs are installed by adding `-a $apps_dir`
to the code below. `apps_dir` is just wherever you would like to install these programs.

```Shell
install_gcc
install_libtool
install_perl
install_bzip2
install_xz
install_zlib
install_curl
install_expat
install_db
install_htslib
install_VEP -r $vep_rel
install_strelka2

```

If you rely on these functions for installations, these will determine 
`hts_dir`, `stk2_dir`, and `vep_dir` definitions.

</details>

## Get COSMIC VCF

<details>
<summary>Click to expand!</summary>

Run the following code to obtain the appropriate COSMIC VCF.
The function below will prompt the user to input COSMIC website's login
and password. The generated file is needed for `TO_workflow()` below.

```Shell
get_COSMIC_canonical -g $genome \
	-v $cosm_ver -h $hts_dir \
	-c $cosm_dir

```

The function definition is located 
[here](https://github.com/pllittle/baSHic/blob/main/scripts/genomic.sh#L480).

</details>

## TO_workflow()

<details>
<summary>Click to expand!</summary>

Currently, `TO_workflow` is designed for `Strelka2` and `VEP`. If others have
alternate workflows, you are welcome to inspect this function's 
[definitions](https://github.com/pllittle/UNMASC/blob/main/workflow/tumor_only.sh#L10)
to extract specific steps to execute :smile:.

```Shell
TO_workflow -c $nthreads -f $fasta_fn -g $genome \
	-d $cosm_dir -e $cosm_ver -h $hts_dir -k $gatk_dir \
	-n $nbams -o $out_dir -s $stk2_dir -t $tbam \
	-v $vep_dir -r $vep_rel -a $vep_cache

```

</details>

## Preparing UNMASC's `vcf`

<details>
<summary>Click to expand!</summary>

Assuming `UNMASC` is successfully installed and the above `TO_workflow()` was
run successfully, the instructions below describes the inputs to construct 
UNMASC's main input `vcf`. Otherwise, the user needs to construct `vcf` ensuring
all required columns are available and formatted correctly (refer to 
`?UNMASC::run_UNMASC`).

* `outdir` String instructing where UNMASC outputs should be stored. Notice
this is different from the Shell variable `out_dir` from above.
* `DAT` A R data.frame containing `FILENAME` for full path to each 
control VCF and `STUDYNUMBER` to identify each normal control.
* `FILTER` R list used to specify some liberal pre-filtering of loci
based on total read depth and quality score. 
* `target_fn` String containing the full path to a target BED file
with tab-delimited columns with headers `Chr` (e.g. chr1), 
`Start` (start position), and `End` (end position).
* `anno_fn` String containing the full path to the annotated vcf file.
If `TO_workflow()` was used, the R string `anno_fn` is 
`$out_dir/allvar_ann.vcf.gz`.
* `nlines` Positive integer specifying how many lines into the `anno_fn`
to initially read in to determine if its formatting matches what 
`prep_UNMASC_VCF()` is expecting.
* `ncores` Positive integer specifying how many threads/cores are available.
This is used here to reduce the computational time to import each control VCF.
This argument may be more handy if samples underwent WGS or WES sequencing.

```R
vcf = UNMASC::prep_UNMASC_VCF(
	outdir = outdir,
	DAT = DAT,
	FILTER = NULL,
	target_fn = target_fn,
	anno_fn = anno_fn,
	nlines = 100,
	ncores = 1)

```

</details>

## Execution

<details>
<summary>Click to expand!</summary>

Arguments for `run_UNMASC()` with template inputs. 

* `tumorID = "tumor01"`, tumor sample ID
* `outdir = file.path(".",tumorID)`, output location for the tumor sample
* `vcf = vcf`, the data.frame generated by `prep_UNMASC_VCF()`
* `tBAM_fn = "path/to/tumor/bam"`
* `bed_centromere_fn = "path/to/centromere/start/end/bed/file"`, tab-delimited
file without headers with three columns containing contig, start position, 
and end position.
* `dict_chrom_fn = "path/to/chromosome/length/file"`, stored output from 
`samtools view -H tumor.bam`
* `qscore_thres = 30`, Qscore threshold
* `exac_thres = 5e-3`, gnomAD/ExAC population allele frequency threshold for 
germline filtering
* `ad_thres = 5`, alternate depth threshold
* `rd_thres = 10`, total depth threshold
* `cut_BAF = 5e-2`, cutoff for variants to exclude before running segmentation
* `minBQ = 13`, minimum base quality
* `minMQ = 40`, minimum mapping quality
* `eps_thres = 0.5`, noise mixture proportion threshold for determining a H2M 
segment
* `psi_thres = 0.02`, over-dispersion of beta-binomial threshold for determining 
a H2M segment
* `hg = "19"`, labeling for output figures
* `binom = TRUE`, set to `TRUE` to model read counts with binomial distribution. 
Set to `FALSE` to explore over both binomial and beta-binomial distributions
* `gender = NA`, set to `NA` if gender is unknown. Otherwise set to `"MALE"` 
or `"FEMALE"`
* `ncores = 1`, number of threads/cores available, aids with computational runtime 
in extracting strand-specific read counts per locus.

```R
# Package main function
UNMASC::run_UNMASC(
	tumorID = tumorID,
	outdir = outdir,
	vcf = vcf,
	tBAM_fn = tBAM_fn,
	bed_centromere_fn = bed_centromere_fn,
	dict_chrom_fn = dict_chrom_fn,
	qscore_thres = qscore_thres,
	exac_thres = exac_thres,
	ad_thres = ad_thres,
	rd_thres = rd_thres,
	cut_BAF = cut_BAF,
	minBQ = minBQ,
	minMQ = minMQ,
	eps_thres = eps_thres,
	psi_thres = psi_thres,
	hg = hg,
	binom = binom,
	gender = gender,
	ncores = ncores)

```

</details>

## Outputs

<details>
<summary>Click to expand!</summary>

Below are the main outputs from `UNMASC` to determine if
all steps ran smoothly and to assess if the package or 
input files requires debugging. If any of these directories 
or files are missing, check the corresponding Rout file 
for any error or warning messages or flags for low quality 
sample information.

* `image.rds`: Stores the comprehensive inputs for `UNMASC`.
Once this file is created, `run_UNMASC()` skips past the 
time consuming pre-processing of input files. This image 
can be useful for re-producing results and inspecting errors. 
To have a clean restart or reset, first remove this file.
* `nCLUST`: Figures of normal VAF clustering. If three clusters
of normal VAF are not present, the loci supplied to `UNMASC` may
have undergone pre-filtering of the normal VAF.
* `nSEG`: Figures of normal VAF segmentation. This is useful for 
visualizing hard-to-map (H2M) regions and whether or not read counts
should be modeled with a binomial or beta-binomial distribution.
* `tSEG`: Figures of tumor VAF segmentation. This is useful for
assessing the degree of sparsity when characterizing the local germline 
cluster behavior relative to each potential somatic locus. Also
these figures may prove useful for inspecting copy number aberrations
due to allelic imbalance (B allele frequencies deviating from 0.5).
* `tumor_genotype.tsv`: Experimental output. Attempting to 
reverse-genotype an individual based on their tumor genomics. Loci
are annotated with inferred genotypes, H2M status, strand bias, etc. 
to aid in isolating higher quality genotype calls. May be useful
for performing genotype PCA or mapping NGS-based sample data to 
microarray sample data.
* `tumorOnly_VCs.tsv`: `UNMASC`'s tumor-only variant calls with 
comprehensive annotation contained in the `LABEL` column for 
prioritizing variants. Additional columns contain metrics used
to construct the `LABEL` columns annotations such as strand bias,
oxoG artifact, paraffin artifact, H2M status, germline-like loci, 
etc.

</details>
