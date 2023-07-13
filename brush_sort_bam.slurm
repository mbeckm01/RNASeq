#!/bin/bash

#SBATCH --job-name=sort_test
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=8:00:00

module load samtools

in=/projects/mougeots_research/micaela_beckman/BAMs
out=/projects/mougeots_research/micaela_beckman/sorted_BAMs/


for f in $(ls $in/*out.bam)
do
	out_prefix=${f#"$in/"}
	out_prefix=${out_prefix%"_sorted.bam"}
	samtools sort -t  -o $out_prefix -O bam -T $out_prefix
done

