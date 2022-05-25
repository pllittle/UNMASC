# Preparing UNMASC Inputs

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

## Custom local installations

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

## Get COSMIC VCF

Run the following code to obtain the appropriate COSMIC VCF.
The function below will prompt the user to input COSMIC website's login
and password.

```Shell
get_COSMIC_canonical -g $genome \
	-v $cosm_ver -h $hts_dir \
	-c $cosm_dir

```

## My custom function for Strelka2 + VEP

Currently, `TO_workflow` is designed for `Strelka2` and `VEP`. If others have
alternate workflows, you are welcome to inspect this function's 
[definitions](https://github.com/pllittle/UNMASC/blob/main/workflow/tumor_only.sh)
to extract specific steps to execute :smile:.

```Shell
TO_workflow -c $nthreads -f $fasta_fn -g $genome \
	-d $cosm_dir -e $cosm_ver -h $hts_dir -k $gatk_dir \
	-n $nbams -o $out_dir -s $stk2_dir -t $tbam \
	-v $vep_dir -r $vep_rel -a $vep_cache

```