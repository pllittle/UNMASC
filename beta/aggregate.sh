#!/bin/sh

agg_UNMASC(){
	local filename fn header out_dir tag
	local nsamps ii SAMP SAMP_PATH out_fn
	
	header=no
	
	# Parse cmdline arguments
	while [ ! -z "$1" ]; do
		case $1 in
			-f | --filename )
				shift
				filename="$1"
				;;
			-h | --header )
				header=yes
				;;
			-o | --out_dir )
				shift
				out_dir="$1"
				;;
			-t | --tag )
				shift
				tag="$1"
				;;
		esac
		shift
	done
	
	# Check inputs
	[ -z "$filename" ] 	&& echo "Add -f <filename with list of files>" >&2 && return 1
	[ -z "$out_dir" ] 	&& echo "Add -o <output directory>" >&2 && return 1
	[ -z "$tag" ] 			&& echo "Add -t <tag, for output filename with extension>" >&2 && return 1
	
	# Check stuff
	[ ! -f "$filename" ] 	&& echo -e "$filename is missing" >&2 && return 1
	[ ! -d "$out_dir" ] 	&& mkdir "$out_dir"
	
	# Process file
	fn=$(echo $filename | sed 's|/|\n|g' | tail -n 1)
	[ "$header" == "no" ] && cp "$filename" "$out_dir/$fn"
	[ "$header" == "yes" ] && tail -n +2 "$filename" > "$out_dir/$fn"
	filename="$out_dir/$fn"
	out_fn="$out_dir/$tag"
	
	# Print some info to user
	nsamps=$(cat "$filename" | wc -l)
	echo -e "Your file contains $nsamps samples to aggregate together" >&2
	
	# Parse file, aggregate
	for ii in $(seq $nsamps); do
		
		SAMP=$(sed -n ${ii}p "$filename" | cut -f 1)
		SAMP_PATH=$(sed -n ${ii}p "$filename" | cut -f 2)
		
		if [ $ii -eq 1 ]; then
			head -n 1 "$SAMP_PATH" | awk '{print "ID\t" $0}' > "$out_fn"
		fi
		
		tail -n +2 "$SAMP_PATH" | awk -v ID=$SAMP '{print ID "\t" $0}' >> "$out_fn"
		
	done
	
	echo -e "Output file located at $out_fn" >&2
	
	return 0
	
}




##



##