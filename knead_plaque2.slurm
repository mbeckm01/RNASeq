#!/bin/bash
#SBATCH --job-name=kneadpl2
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=8:00:00


module load kneaddata
module load bowtie2
module load samtools

#input directory
in=/nobackup/mougeots_research/micaela_beckman/RNASeq_Plaque/trimmed_plaque_fastqs/
#output directory
out=/nobackup/mougeots_research/micaela_beckman/RNASeq_Plaque/kneaddata_data/homo_sapiens_index

cd $out

#kneaddata_database --download human_genome bowtie2 $out

#test the kneaddata workflow

kneaddata --input $in/A1FAMP2-24.IonXpressRNA_009_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B1FAMP3-11.IonXpressRNA_014_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B1FPMP2-24.IonXpressRNA_007_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B2FAMP2-24.IonXpressRNA_008_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B2FAMP3-11.IonXpressRNA_015-modified_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B3MAMP2-24.IonXpressRNA_002-modified_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B3MAMP3-11.IonXpressRNA_016_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B3MPMP2-24.IonXpressRNA_008_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B4MAMP2-24.IonXpressRNA_003_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B4MAMP3-11.IonXpressRNA_001_trimmed.fq --reference-db $out/ --output $out
kneaddata --input $in/B4MPMP2-24.IonXpressRNA_008_trimmed.fq --reference-db $out/ --output $out
