# RNASeq
All data and accompanying files for "Human oral mucosa and oral microbiome interactions following supragingival plaque reconstitution: a diet-controlled balanced design proof-of-concept model for oral toxicities"


Files are as follows:

align_ncbi.slurm: slurm script to align kneaddata output (bacterial RNA) to large ncbi index, convert sam to bam and then sort bam files

bams.slurm: slurm script to convert sam to bam files and sort the bam files for human plaque samples

brush_sort_bam.slurm: slurm script to convert brush biopsies from sam to bam and sort

featurecountscommandtrimmed.txt: command used to run featurecounts on RNASeq data

deseq_commands.R: R script used to run DESeq2 on human brush/plaque samples, annotate dataset, and return differentially expressed genes, go ters, and kegg pathways

humann2_plaque.slurm: slurm script used to run HumanN2 on bacterial plaque samples

knead_plaque2.slurm: slurm script used to remove bacterial RNA data from human plaque samples

star_align.slurm: slurm script to align brush and human plaque samples
