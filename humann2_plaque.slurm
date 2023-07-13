#!/bin/bash

#SBATCH --job-name=humann2
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=8:00:00


module load humann2/0.11.2
module load samtools/1.10

#input directory
in=/nobackup/mougeots_research/micaela_beckman/RNASeq_Plaque/kneaddata_data
#output directory
in2=/nobackup/mougeots_research/micaela_beckman/RNASeq_Plaque
out=/nobackup/mougeots_research/micaela_beckman/RNASeq_Plaque/humann2

cd $out
#test humann2 dependencies
humann2 --version

#run humann2 on one sample as a test

#humann2 --input $in/A2MAMP3-11.IonXpressRNA_007_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/A2MPMP3-11.IonXpressRNA_012_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/A3MAMP3-11.IonXpressRNA_009_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/A3MPMP3-11.IonXpressRNA_013_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/A4FAMP3-11.IonXpressRNA_010_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/A4FPMP3-11.IonXpressRNA_014_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B1FAMP2-24.IonXpressRNA_001_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B1FAMP3-11.IonXpressRNA_014_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B1FPMP2-24.IonXpressRNA_007_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B2FAMP2-24.IonXpressRNA_008_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B2FAMP3-11.IonXpressRNA_015-modified_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B3MAMP2-24.IonXpressRNA_002-modified_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B3MAMP3-11.IonXpressRNA_016_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B3MPMP2-24.IonXpressRNA_008_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B4MAMP2-24.IonXpressRNA_003_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B4MAMP3-11.IonXpressRNA_001_trimmed_kneaddata.trimmed.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#humann2 --input $in/B4MPMP2-24.IonXpressRNA_008_trimmed_kneaddata.fastq --resume --threads 12 --output $out/ --nucleotide-database $in2/chocophlan --protein-database $in2/uniref

#join tables together
#humann2_join_tables  -i $out -o merged_plaque_genefamilies.tsv --file_name genefamilies

#renomalize the merged table
#humann2_renorm_table -i merged_plaque_genefamilies.tsv -o renorm_genefamilies_cpm.tsv --units cpm

#merge pathway abundances together
#humann2_join_tables -i $out -o merged_plaque_pathabundances.tsv --file_name pathabundance

#renormalize the merged table
#humann2_renorm_table -i merged_plaque_pathabundances.tsv -o renorm_pathabundance_cpm.tsv --units cpm

#associate
#humann2_associate --input STIME_renorm_pathabundance_cpm.tsv --last-metadatum STIME --focal-metadatum STIME --focal-type categorical --output stats.txt
#plot 
humann2_barplot --input STIME_renorm_pathabundance_cpm.tsv --sort sum --focal-feature PWY-6151 --focal-metadatum STIME --last-metadatum STIME --output plaque_pathabundance2_PWY6151_time.png 

#sum metadata 
humann2_barplot --input STIME_renorm_pathabundance_cpm.tsv --sort sum metadata --focal-feature PWY-6151 --focal-metadatum STIME --last-metadatum STIME --output plaque_pathabundance2_PWY6151_sumtime.png
#try barplot
#humann2_barplot --input $out/A1FAMP2-24.humann2/A1FAMP2-24.IonXpressRNA_009_trimmed_pathabundance.tsv --sort sum --output test_A1FAMP2-24.png
