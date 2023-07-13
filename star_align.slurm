#!/bin/bash

#SBATCH --job-name=align_brush_star
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=8:00:00
#SBATCH --mem=64GB

IFS=$'\n'
set -eu
umask 007
shopt -s nullglob

IN_DIR=/projects/mougeots_research/micaela_beckman/RNASeq_Brush_Analysis/brush_trim
GENOME_DIR=/projects/mougeots_research/micaela_beckman/RNASeq_Brush_Analysis/star_alignment/out
OUT_DIR=/projects/mougeots_research/micaela_beckman/RNASeq_Brush_Analysis/bam_trim
THREADS=32


module load star/2.7.9a

cd $OUT_DIR

# Aligning brush samples to human genome using star
printf "Aligning brush samples to human genome using star\n"

i=0
for file in $(ls $IN_DIR/*MB[0-9]*)
do
	printf "Aligning brush sample $file\n"
	out_prefix=${file#"$IN_DIR/"}
	out_prefix=${out_prefix%"fastq"}
	STAR --runThreadN $THREADS --genomeDir $GENOME_DIR --readFilesIn $file --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 33437379172 --outFileNamePrefix $out_prefix

	if [ $i == 0 ];
	then
		echo "file number_of_input_reads average_input_read_length #_mapped_uniquely %_mapped_uniquely #_mapped_to_multiple_loci %_mapped_to_multiple_loci #_mapped_to_too_many_loci %_mapped_to_too_many_loci #_unmapped_too_many_mismatches %_unmapped_too_many_mismatches #_unmapped_too_short %_unmapped_too_short #_unmapped_other %_unmapped_other #_chimeric_reads %_chimeric_reads" > $OUT_DIR/results.txt
	fi

	temp00=$(grep "Number of input reads" $out_prefix"Log.final.out" | awk -F ' ' '{print $6}')
	temp01=$(grep "Average input read length" $out_prefix"Log.final.out" | awk -F ' ' '{print $6}')
	temp02=$(grep "Uniquely mapped reads number" $out_prefix"Log.final.out" | awk -F ' ' '{print $6}')
	temp03=$(grep "Uniquely mapped reads %" $out_prefix"Log.final.out" | awk -F ' ' '{print $6}')
	temp04=$(grep "Number of reads mapped to multiple loci" $out_prefix"Log.final.out" | awk -F ' ' '{print $9}')
	temp05=$(grep "% of reads mapped to multiple loci" $out_prefix"Log.final.out" | awk -F ' ' '{print $9}')
	temp06=$(grep "Number of reads mapped to too many loci" $out_prefix"Log.final.out" | awk -F ' ' '{print $10}')
	temp07=$(grep "% of reads mapped to too many loci" $out_prefix"Log.final.out" | awk -F ' ' '{print $10}')
	temp08=$(grep "Number of reads unmapped: too many mismatches" $out_prefix"Log.final.out" | awk -F ' ' '{print $9}')
	temp09=$(grep "% of reads unmapped: too many mismatches" $out_prefix"Log.final.out" | awk -F ' ' '{print $9}')
	temp10=$(grep "Number of reads unmapped: too short" $out_prefix"Log.final.out" | awk -F ' ' '{print $8}')
	temp11=$(grep "% of reads unmapped: too short" $out_prefix"Log.final.out" | awk -F ' ' '{print $8}')
	temp12=$(grep "Number of reads unmapped: other" $out_prefix"Log.final.out" | awk -F ' ' '{print $7}')
	temp13=$(grep "% of reads unmapped: other" $out_prefix"Log.final.out" | awk -F ' ' '{print $7}')
	temp14=$(grep "Number of chimeric reads" $out_prefix"Log.final.out" | awk -F ' ' '{print $6}')
	temp15=$(grep "% of chimeric reads" $out_prefix"Log.final.out" | awk -F ' ' '{print $6}')

	echo ${out_prefix%"."} $temp00 $temp01 $temp02 $temp03 $temp04 $temp05 $temp06 $temp07 $temp08 $temp09 $temp10 $temp11 $temp12 $temp13 $temp14 $temp15 >> $OUT_DIR/results.txt

	i=$((i+1))

	printf "Done\n\n"
done

printf "DONE\n\n"
