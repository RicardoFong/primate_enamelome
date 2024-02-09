#!/bin/bash

### submitted as job array in SLURM workload manager

######
### Arguments
######

inaln=$1
outdir=$2
spec=$3


######
### Prepare modules and directory
######


module load seqtk/1.3-GCC-11.2.0


mkdir -p "$outdir/tmp_$SLURM_ARRAY_TASK_ID"


######
### Randomly sample one indiv per species
######

if [ -f "$outdir/random_MSA_headers.txt" ]; then
rm "$outdir/random_MSA_headers.txt"
fi

touch "$outdir/tmp_$SLURM_ARRAY_TASK_ID/$random_MSA_headers.txt"


	while read species;do 

		# get file with all individuals per species	
		cat $inaln | grep $species > "$outdir/tmp_$SLURM_ARRAY_TASK_ID/tmp_spec.txt"
		
		# get number of species to set up  random number generator
		no_spec=$(wc -l "$outdir/tmp_$SLURM_ARRAY_TASK_ID/tmp_spec.txt" | cut -d" " -f1)

		# run random number generator, use species number as argument to define range of number from 1 - no_spec
			if (( $no_spec > 1 ))
			then	
				rand=$((1 + $RANDOM % $no_spec))
				echo "individual no. $rand out of $no_spec individuals of $species chosen."
			else
				rand=1
				echo "only one individual of $species."
			fi

		# Get random line	
		cat "$outdir/tmp_$SLURM_ARRAY_TASK_ID/tmp_spec.txt" | sed -n "$rand"p >> "$outdir/tmp_$SLURM_ARRAY_TASK_ID/random_MSA_headers.tmp"


	done < "$spec"

# get fasta entries of all headers
cat "$outdir/tmp_$SLURM_ARRAY_TASK_ID/random_MSA_headers.tmp" | cut -d">" -f2 > "$outdir/tmp_$SLURM_ARRAY_TASK_ID/random_MSA_headers.txt"

seqtk subseq "$inaln" "$outdir/tmp_$SLURM_ARRAY_TASK_ID/random_MSA_headers.txt" > "$outdir/random_MSA_$SLURM_ARRAY_TASK_ID.fasta"

rm -r "$outdir/tmp_$SLURM_ARRAY_TASK_ID"
