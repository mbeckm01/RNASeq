#!/bin/bash
  
#SBATCH --job-name=align_2
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=32G

module load bowtie2
module load samtools

#input directory
in=/nobackup/mougeots_research/micaela_beckman/RNASeq_Plaque/kneaddata_data
#output directory
out=/nobackup/mougeots_research/micaela_beckman/RNASeq_Plaque/ncbi_align

cd $out

#bowtie2-build -f --large-index $out/ncbi_datasets.fna ncbi_index


#bowtie2 -q -a -x ncbi_index $in/A2MAMP3-11.IonXpressRNA_007_trimmed_kneaddata.fastq -S A2MAMP3-11_007_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS A2MAMP3-11_007_trimmed_kneaddata_bowtie2_ncbi.sam > A2MAMP3-11_007.bam
#samtools sort A2MAMP3-11_007.bam -o A2MAMP3-11_007_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/A1FAMP2-24.IonXpressRNA_009_trimmed_kneaddata.fastq -S A1FAMP2-24_009_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS A1FAMP2-24_009_trimmed_kneaddata_bowtie2_ncbi.sam > A1FAMP2-24_009.bam
#samtools sort A1FAMP2-24_009.bam -o A1FAMP2-24_009_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/A2MPMP3-11.IonXpressRNA_012_trimmed_kneaddata.fastq -S A2MPMP3-11_012_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS A2MPMP3-11_012_trimmed_kneaddata_bowtie2_ncbi.sam > A2MPMP3-11_012.bam
#samtools sort A2MPMP3-11_012.bam -o A2MPMP3-11_012_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/A3MAMP3-11.IonXpressRNA_009_trimmed_kneaddata.fastq -S A3MAMP3-11_009_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS A3MAMP3-11_009_trimmed_kneaddata_bowtie2_ncbi.sam > A3MAMP3-11_009.bam
#samtools sort A3MAMP3-11_009.bam -o A3MAMP3-11_009_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/A3MAMP3-11.IonXpressRNA_009_trimmed_kneaddata.fastq -S A3MAMP3-11_009_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS A3MAMP3-11_009_trimmed_kneaddata_bowtie2_ncbi.sam > A3MAMP3-11_009.bam
#samtools sort A3MAMP3-11_009.bam -o A3MAMP3-11_009_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/A3MPMP3-11.IonXpressRNA_013_trimmed_kneaddata.fastq -S A3MPMP3-11_013_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS A3MPMP3-11_013_trimmed_kneaddata_bowtie2_ncbi.sam > A3MPMP3-11_013.bam
#samtools sort A3MPMP3-11_013.bam -o A3MPMP3-11_013_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/A4FAMP3-11.IonXpressRNA_010_trimmed_kneaddata.fastq -S A4FAMP3-11_010_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS A4FAMP3-11_010_trimmed_kneaddata_bowtie2_ncbi.sam > A4FAMP3-11_010.bam
#samtools sort A4FAMP3-11_010.bam -o A4FAMP3-11_010_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/A4FPMP3-11.IonXpressRNA_014_trimmed_kneaddata.fastq -S A4FPMP3-11_014_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS A4FPMP3-11_014_trimmed_kneaddata_bowtie2_ncbi.sam > A4FPMP3-11_014.bam
#samtools sort A4FPMP3-11_014.bam -o A4FPMP3-11_014_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B1FAMP2-24.IonXpressRNA_001_trimmed_kneaddata.fastq -S B1FAMP2-24_001_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B1FAMP2-24_001_trimmed_kneaddata_bowtie2_ncbi.sam > B1FAMP2-24_001.bam
#samtools sort B1FAMP2-24_001.bam -o B1FAMP2-24_001_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B1FAMP3-11.IonXpressRNA_014_trimmed_kneaddata.fastq -S B1FAMP3-11_014_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B1FAMP3-11_014_trimmed_kneaddata_bowtie2_ncbi.sam > B1FAMP3-11_014.bam
#samtools sort B1FAMP3-11_014.bam -o B1FAMP3-11_014_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B1FPMP2-24.IonXpressRNA_007_trimmed_kneaddata.fastq -S B1FPMP2-24_007_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B1FPMP2-24_007_trimmed_kneaddata_bowtie2_ncbi.sam > B1FPMP2-24_007.bam
#samtools sort B1FPMP2-24_007.bam -o B1FPMP2-24_007_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B2FAMP2-24.IonXpressRNA_008_trimmed_kneaddata.fastq -S B2FAMP2-24_008_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B2FAMP2-24_008_trimmed_kneaddata_bowtie2_ncbi.sam > B2FAMP2-24_008.bam
#samtools sort B2FAMP2-24_008.bam -o B2FAMP2-24_008_sorted.bam

bowtie2 -q -a -x ncbi_index $in/B2FAMP3-11.IonXpressRNA_015-modified_trimmed_kneaddata.fastq -S B2FAMP3-11_015_trimmed_kneaddata_bowtie2_ncbi.sam
samtools view -bS B2FAMP3-11_015_trimmed_kneaddata_bowtie2_ncbi.sam > B2FAMP3-11_015.bam
samtools sort B2FAMP3-11_015.bam -o B2FAMP3-11_015_sorted.bam

bowtie2 -q -a -x ncbi_index $in/B3MAMP2-24.IonXpressRNA_002-modified_trimmed_kneaddata.fastq -S B3MAMP2-24_002_trimmed_kneaddata_bowtie2_ncbi.sam
samtools view -bS B3MAMP2-24_002_trimmed_kneaddata_bowtie2_ncbi.sam > B3MAMP2-24_002.bam
samtools sort B3MAMP2-24_002.bam -o B3MAMP2-24_002_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B3MAMP3-11.IonXpressRNA_016_trimmed_kneaddata.fastq -S B3MAMP3-11_016_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B3MAMP3-11_016_trimmed_kneaddata_bowtie2_ncbi.sam > B3MAMP3-11_016.bam
#samtools sort B3MAMP3-11_016.bam -o B3MAMP3-11_016_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B3MPMP2-24.IonXpressRNA_008_trimmed_kneaddata.fastq -S B3MPMP2-24_008_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B3MPMP2-24_008_trimmed_kneaddata_bowtie2_ncbi.sam > B3MPMP2-24_008.bam
#samtools sort B3MPMP2-24_008.bam -o B3MPMP2-24_008_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B4MAMP2-24.IonXpressRNA_003_trimmed_kneaddata.fastq -S B4MAMP2-24_003_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B4MAMP2-24_003_trimmed_kneaddata_bowtie2_ncbi.sam > B4MAMP2-24_003.bam
#samtools sort B4MAMP2-24_003.bam -o B4MAMP2-24_003_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B4MAMP3-11.IonXpressRNA_001_trimmed_kneaddata.fastq -S B4MAMP3-11_001_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B4MAMP3-11_001_trimmed_kneaddata_bowtie2_ncbi.sam > B4MAMP3-11_001.bam
#samtools sort B4MAMP3-11_001.bam -o B4MAMP3-11_001_sorted.bam

#bowtie2 -q -a -x ncbi_index $in/B4MPMP2-24.IonXpressRNA_008_trimmed_kneaddata.fastq -S B4MPMP2-24_008_trimmed_kneaddata_bowtie2_ncbi.sam
#samtools view -bS B4MPMP2-24_008_trimmed_kneaddata_bowtie2_ncbi.sam > B4MPMP2-24_008.bam
#samtools sort B4MPMP2-24_008.bam -o B4MPMP2-24_008_sorted.bam

