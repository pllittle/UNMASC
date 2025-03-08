#!/bin/sh

[ -z "$git_dir" ] && git_dir=$(cd $(dirname $BASH_SOURCE)/../..; pwd)

odir=$(pwd)

for repo in baSHic copythatdna somdna UNMASC; do
	
	repo_dir="$git_dir/$repo"
	tmp_url=https://github.com/pllittle/$repo.git
	
	if [ ! -d "$repo_dir" ]; then
		cd "$git_dir"
		git clone "$tmp_url" >&2
		[ $? -eq 0 ] && continue
	else
		cd "$repo_dir"
		git pull >&2
		[ $? -eq 0 ] && continue
	fi
	
	echo -e "Error clone/pull $repo" >&2
	return 1
	
done

cd "$odir"

for fn in base colors getEnv install linux_latex \
	linux_python linux_perl; do
	. "$git_dir/baSHic/scripts/$fn.sh"
	[ $? -eq 0 ] && continue
	echo -e "Error src-ing baSHic's $fn.sh" >&2
	return 1
done

for fn in getCounts; do
	. "$git_dir/copythatdna/scripts/$fn.sh"
	[ $? -eq 0 ] && continue
	echo -e "Error src-ing copythatdna's $fn.sh" >&2
	return 1
done

for fn in callingAllVariants; do
	. "$git_dir/somdna/scripts/$fn.sh"
	[ $? -eq 0 ] && continue
	echo -e "Error src-ing somdna's $fn.sh" >&2
	return 1
done

unset fn repo repo_dir

# Tumor-only Workflow, credit to Heejoon Jo
TO_workflow(){
	local genome cosm_dir cosm_ver stk_dir gatk_dir vep_dir vep_rel cache_type status
	local cache_dir nbam nbam2 nbams tbam fasta_fn out_dir cosmic_fn species
	local gnomad_dir gnomad_fn ncores work_dir hts_dir cosm_muts cnt tbam2 cmd
	
	while [ ! -z $1 ]; do
		case $1 in
			-c | --ncores )
				shift
				ncores="$1"
				;;
			-d | --cosm_dir )
				shift
				cosm_dir="$1"
				;;
			-e | --cosm_ver )
				shift
				cosm_ver="$1"
				;;
			-f | --fasta_fn )
				shift
				fasta_fn="$1"
				;;
			-g | --genome )
				shift
				genome="$1"
				;;
			-h | --hts_dir )
				shift
				hts_dir="$1"
				;;
			-k | --gatk )
				shift
				gatk_dir="$1"
				;;
			-n | --nbams )
				shift
				nbams="$1"
				;;
			-o | --out_dir )
				shift
				out_dir="$1"
				;;
			-s | --stk_dir )
				shift
				stk_dir="$1"
				;;
			-t | --tbam )
				shift
				tbam="$1"
				;;
			-v | --vep_dir )
				shift
				vep_dir="$1"
				;;
			-r | --vep_rel )
				shift
				vep_rel="$1"
				;;
			--cache_type )
				shift
				cache_type="$1"
				;;
			--cache_dir )
				shift
				cache_dir="$1"
				;;
			--gnomad_dir )
				shift
				gnomad_dir="$1"
				;;
			--species )
				shift
				species="$1"
				;;
			--cosm_muts )
				shift
				cosm_muts="$1"
				;;
		esac
		shift
	done
	
	# Check inputs
	[ -z $ncores ] 		&& echo "Add -c <number of threads, cores>" >&2 && return 1
	[ -z "$fasta_fn" ] 	&& echo "Add -f <fasta_fn>" >&2 && return 1
	[ -z $genome ] 		&& echo "Add -g <genome, e.g. GRCh37, GRCh38>" >&2 && return 1
	[ -z "$hts_dir" ] 	&& echo "Add -h <hts_dir, the ../../ directory of htsfile>" >&2 && return 1
	[ -z "$gatk_dir" ] 	&& echo "Add -k <GATK dir>" >&2 && return 1
	[ -z "$nbams" ] 		&& echo "Add -n <normal bams file list>" >&2 && return 1
	[ -z "$out_dir" ] 	&& echo "Add -o <tumor output dir>" >&2 && return 1
	[ -z "$stk_dir" ] 	&& echo "Add -s <Strelka2 dir>" >&2 && return 1
	[ -z "$tbam" ] 			&& echo "Add -t <tumor bam>" >&2 && return 1
	[ -z "$vep_dir" ] 	&& echo "Add -v <vep_dir>" >&2 && return 1
	[ -z "$vep_rel" ] 	&& echo "Add -r <vep_rel, VEP release>" >&2 && return 1
	[ -z "$cosm_dir" ] 	&& echo "Add -d <COSMIC dir>" >&2 && return 1
	[ -z $cosm_ver ] 	&& echo "Add -e <COSMIC version, e.g. 94, 95, 101>" >&2 && return 1
	[ -z "$cache_type" ] && echo "Add --cache_type <VEP cache, e.g. vep, refseq, merged>" >&2 && return 1
	[ -z "$cache_dir" ] && echo "Add --cache_dir <dir>" >&2 && return 1
	[ -z "$gnomad_dir" ] && echo "Add --gnomad_dir <dir>" >&2 && return 1
	[ -z "$species" ] && species=homo_sapiens
	[ -z "$cosm_muts" ] && echo "Add --cosm_muts <coding, noncoding, all>" >&2 && return 1
	
	check_array $cache_type vep refseq merged
	[ ! $? -eq 0 ] && echo "Error: VEP cache should be vep, refseq, or merged" >&2 && return 1
	
	check_array "$cosm_muts" coding noncoding all
	[ ! $? -eq 0 ] && echo "Error: cosm_muts value incorrect" >&2 && return 1
	
	new_mkdir "$out_dir"
	tbam2=`echo $tbam | sed 's|/|\n|g' | tail -n 1`
	cnt=1
	
	# Run variant calling with strelka
	for nbam in `cat $nbams`; do
		nbam2=`echo $nbam | sed 's|/|\n|g' | tail -n 1`
		echo -e "`date`: $cnt) tumor = $tbam2 vs normal = $nbam2" >&2
		
		work_dir="$out_dir/normal_$cnt"
		new_mkdir "$work_dir"
		
		# Run Strelka2
		run_strelka2_soma -s "$stk_dir" -g "$gatk_dir" \
			-t "$tbam" -n "$nbam" -r "$fasta_fn" -o "$work_dir" \
			-c $ncores
		
		# Append merged vcf to txt: retain chr1-chr22,chrX,chrY
		[ ! -s "$work_dir/somatic.vcf.gz" ] \
			&& echo -e "Missing/empty $work_dir/somatic.vcf.gz" >&2 \
			&& return 1
		zgrep -v ^# "$work_dir/somatic.vcf.gz" | cut -f1-5 \
			| awk '$1 ~ /^chr(1?[0-9]|2[0-2]|X|Y)$/' \
			> "$out_dir/somatic_$cnt.txt"
		
		let cnt=cnt+1
	done
	
	# Combine across normals
	new_rm "$out_dir/allvar.vcf" "$out_dir/allvar.vcf.gz"
	echo -e "`date`: Get unique loci calls for annotation" >&2
	cat $out_dir/somatic_*.txt | sort -k1,1V -k2,2n \
		| uniq > $out_dir/allvar.vcf
	rm $out_dir/somatic_*.txt
	
	# Get COSMIC if missing
	cosmic_fn="$cosm_dir"
	[ "$cosm_muts" == "coding" ] && cosmic_fn="$cosmic_fn/CosmicCodingMuts"
	[ "$cosm_muts" == "noncoding" ] && cosmic_fn="$cosmic_fn/CosmicNonCodingMuts"
	[ "$cosm_muts" == "all" ] && cosmic_fn="$cosmic_fn/CosmicMuts"
	cosmic_fn="${cosmic_fn}_${genome}_v${cosm_ver}_canonical.vcf.gz"
	
	if [ ! -f "$cosmic_fn" ]; then
		cmd="get_COSMIC_canonical -g $genome -v $cosm_ver"
		cmd="$cmd -c \"$cosm_dir\" -h \"$hts_dir\""
		echo -e "Run this command to get all the COSMIC vcfs:\n\n\t${purple}${cmd}${NC}\n" >&2
		return 1
	fi
	
	# Check gnomad file
	gnomad_fn="$gnomad_dir/gnomad.genomes.r2.0.1.sites.noVEP.$species.$genome.vcf.gz"
	[ ! -f "$gnomad_fn" ] && echo -e "$gnomad_fn missing" >&2 && return 1
	
	# Run VEP
	run_VEP -c "$cosmic_fn" \
		-a "$gnomad_fn" \
		-s "$species" \
		-f "$fasta_fn" \
		-g $genome \
		-i "$out_dir/allvar.vcf" \
		-n $ncores \
		-o "$out_dir/allvar_ann.vcf" \
		-v "$vep_dir" \
		-r "$vep_rel" \
		-h "$cache_dir" \
		-t "$cache_type"
	[ ! $? -eq 0 ] && echo "Error with VEP" >&2 && return 1
	new_rm "$out_dir/allvar.vcf"
	
	return 0
	
}


###

