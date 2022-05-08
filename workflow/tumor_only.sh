#!/bin/sh

[ -z "$git_dir" ] && git_dir=$(cd $(dirname $BASH_SOURCE)/../..; pwd)

for fn in base colors getEnv genomic; do
	. $git_dir/baSHic/scripts/$fn.sh
done

# Tumor-only Workflow, credit to Heejoon Jo
TO_workflow(){
	local genome cosm_dir cosm_ver stk_dir gatk_dir vep_dir vep_rel status
	local nbam nbam2 nbams tbam fasta_fn out_dir cosmic_fn vep_cache
	local ncores work_dir hts_dir cnt tbam2 cmd
	
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
			-a | --vep_cache )
				shift
				vep_cache="$1"
				;;
		esac
		shift
	done
	
	# Check inputs
	[ -z $ncores ] 		&& echo "Add -c <number of threads, cores>" >&2 && return 1
	[ -z $fasta_fn ] 	&& echo "Add -f <fasta_fn>" >&2 && return 1
	[ -z $genome ] 		&& echo "Add -g <genome, e.g. GRCh37, GRCh38>" >&2 && return 1
	[ -z $hts_dir ] 	&& echo "Add -h <hts_dir>" >&2 && return 1
	[ -z $gatk_dir ] 	&& echo "Add -k <GATK dir>" >&2 && return 1
	[ -z $nbams ] 		&& echo "Add -n <normal bams file list>" >&2 && return 1
	[ -z $out_dir ] 	&& echo "Add -o <tumor output dir>" >&2 && return 1
	[ -z $stk_dir ] 	&& echo "Add -s <Strelka2 dir>" >&2 && return 1
	[ -z $tbam ] 			&& echo "Add -t <tumor bam>" >&2 && return 1
	[ -z $vep_dir ] 	&& echo "Add -v <vep_dir>" >&2 && return 1
	[ -z $vep_rel ] 	&& echo "Add -r <vep_rel, VEP release>" >&2 && return 1
	[ -z $cosm_dir ] 	&& echo "Add -d <COSMIC dir>" >&2 && return 1
	[ -z $cosm_ver ] 	&& echo "Add -e <COSMIC version, e.g. 94, 95>" >&2 && return 1
	[ -z $vep_cache ] && echo "Add -a <VEP cache, e.g. vep, refseq, merged>" >&2 && return 1
	
	check_array $vep_cache vep refseq merged
	[ ! $? -eq 0 ] && echo "Error: VEP cache should be vep, refseq, or merged" >&2 && return 1
	
	new_mkdir $out_dir
	tbam2=`echo $tbam | sed 's|/|\n|g' | tail -n 1`
	cnt=1
	
	# Run variant calling with strelka
	for nbam in `cat $nbams`; do
		nbam2=`echo $nbam | sed 's|/|\n|g' | tail -n 1`
		echo -e "`date`: $cnt) tumor = $tbam2 vs normal = $nbam2" >&2
		
		work_dir=$out_dir/normal_$cnt
		new_mkdir $work_dir
		
		# Run Strelka2
		run_strelka2_soma -s $stk_dir -g $gatk_dir \
			-t $tbam -n $nbam -r $fasta_fn -o $work_dir \
			-c $ncores
		
		# Append merged vcf to txt: retain chr1-chr22,chrX,chrY
		[ ! -s $work_dir/somatic.vcf.gz ] \
			&& echo -e "Missing/empty $work_dir/somatic.vcf.gz" >&2 \
			&& return 1
		zgrep -v ^# $work_dir/somatic.vcf.gz | cut -f1-5 \
			| awk '$1 ~ /^chr(1?[0-9]|2[0-2]|X|Y)$/' \
			> $out_dir/somatic_$cnt.txt
		
		let cnt=cnt+1
	done
	
	# Combine across normals
	new_rm $out_dir/allvar.vcf $out_dir/allvar.vcf.gz
	echo -e "`date`: Get unique loci calls for annotation" >&2
	cat $out_dir/somatic_*.txt | sort -k1,1V -k2,2n \
		| uniq > $out_dir/allvar.vcf
	rm $out_dir/somatic_*.txt
	
	# Get COSMIC if missing
	cosmic_fn=$cosm_dir/CosmicCodingMuts_${genome}_v${cosm_ver}_canonical.vcf.gz
	if [ ! -f $cosmic_fn ]; then
		cmd="get_COSMIC_canonical -g $genome -v $cosm_ver"
		cmd="$cmd -c $cosm_dir -h $hts_dir"
		echo -e "Run this command to get the COSMIC vcf:\n\n\t${purple}${cmd}${NC}\n" >&2
		return 1
	fi
	
	# Run VEP
	run_VEP -c $cosmic_fn -f $fasta_fn -g $genome \
		-i $out_dir/allvar.vcf -n 1 -o $out_dir/allvar_ann.vcf \
		-v $vep_dir -r $vep_rel -a $vep_cache
	[ ! $? -eq 0 ] && echo "Error with VEP" >&2 && return 1
	new_rm $out_dir/allvar.vcf
	
	return 0
	
}


###

