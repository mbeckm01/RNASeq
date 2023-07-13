#!/bin/bash
  
#SBATCH --job-name=BamSort
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=8:00:00

module load samtools

#input directory
in=/projects/mougeots_research/micaela_beckman/sams
#output directory
out=/projects/mougeots_research/micaela_beckman/BAMs

#change directories to output directory
cd $out
#run samtools on all sam files
for f in $(ls $in/*.sam)
do
	samtools view -bS -h $f >${f/%.sam/.bam} | samtools sort -t -o ${f/.out.bam/_sorted.bam} 
done
