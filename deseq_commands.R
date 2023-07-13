#load packages
library(gage)
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)

#read in a counts file and a metadata file
counts <- read_csv("~/BrushCountsTrim_ZerosRemoved.csv")
View(counts)
meta <- read_csv("~/BrushMetatrim.csv")
View(meta)

#Change the data to be dataframes
counts=as.data.frame(counts)
meta=as.data.frame(meta)

#check to see if counts and meta data have the same genes
names(counts)[-1]
meta$id

#set up DESeq NOTE: tidy must be true if the data is in the form of a dataframe
dds=DESeqDataSetFromMatrix(countData=counts,colData=meta,design=~dex,tidy=TRUE)
dds

#run DESeq
dds=DESeq(dds)

#save the results
res=results(dds,tidy=TRUE)

#change the results to a tibble
res2=as_tibble(res)

#filter out the genes for significance <0.05
res2 %>%
  filter(padj<0.05) %>%	
  #write results to a csv
  write_csv("sigresults.csv")

#plot strip plot of a gene
plotCounts(dds,gene='FAM161A',intgroup='dex')
plotCounts(dds,gene='FAM161A',intgroup='dex',returnData=TRUE)

#make a boxplot of the data
plotCounts(dds,gene='FAM161A',intgroup='dex',returnData=TRUE) %>%
  ggplot(aes(dex,count))+geom_boxplot(aes(fill=dex))+scale_y_log10() + ggtitle("FAM161A")

#mutate results to add a column called sig that evaluates to TRUE if padj<0.05 and false if not, and NA if padj is also NA
res3 = res2 %>% mutate(sig=padj<0.05)
res3 %>%
  group_by(sig)%>%
  summarize(n=n())

#variance stabilizing transformation which removes dependence of the variance on the mean
vsdata=vst(dds,blind=FALSE)

#plot PCA
plotPCA(vsdata,intgroup='dex')

#save results to different variable after ordering adjusted pvalues
res4=res2[order(res2$padj),]

#plot volcano
par(mfrow=c(1,1))
with(res4,plot(log2FoldChange,-log10(pvalue),pch=20,main='Volcano Plot',xlim=c(-3,3)))
with(subset(res4,padj<0.05),points(log2FoldChange,-log10(pvalue),pch=20,col='blue'))
with(subset(res4,padj<0.05 & abs(log2FoldChange)>2),points(log2FoldChange,-log10(pvalue),pch=20,col='red'))

#Load libraries for Pathway Analysis
library("AnnotationDbi")
library("org.Hs.eg.db")
library(gage)
library(gageData)
library(pathview)
#reorder by p-value
result=results(dds,contrast=c('dex','AM','PM'))
result=result[order(result$pvalue),]
summary(result)

#use the mapIds function to add columns to the results table; specify that keytype is ENSEMBL, column tells function what information we want, multiVals tells function what to do if there are multiple possible values for a single input value. Asing to give us the first oen that occurs in the database
result$name=mapIds(org.Hs.eg.db,key=row.names(result),column='GENENAME',keytype='SYMBOL',multiVals='first')
result$symbol=mapIds(org.Hs.eg.db,key=row.names(result),column='SYMBOL',keytype='SYMBOL',multiVals='first')
result$entrez=mapIds(org.Hs.eg.db,key=row.names(result),column='ENTREZID',keytype='SYMBOL',multiVals='first')

#viewheader
head(result,10)

#kegg.sets.hs is a named list of 229 elements. Each element is a character vector of member gene Entrez IDs for a single KEGG pathway. (See also go.sets.hs). sigmet.idx.hs is an index of numbers of sinaling and metabolic pathways in kegg.set.gs.
data(kegg.sets.hs)
data(sigmet.idx.hs)

#gives you the 'cleaner' gene sets of signaling and metabolic pathways only
kegg.sets.hs=kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)

#create a named vector of fold changes, where the names of the values are Entrez gene ids
foldchanges=result$log2FoldChange
names(foldchanges)=result$entrez
head(foldchanges)

#get the kegg results
keggresult=gage(foldchanges,gsets=kegg.sets.hs,same.dir=TRUE)

#look at both greater, down, and the statistics
lapply(keggresult,head)

#get the kegg pathways, set as a dataframe,filter and view the first 5
keggresultpaths=data.frame(id=rownames(keggresult$less),keggresult$less)%>%
  as.data.frame() %>%
  filter(row_number()<=5) %>%
  .$id %>%
  as.character()
keggresultpaths

#get the IDs
keggresultids=substr(keggresultpaths,start=1,stop=8)
keggresultids

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresultids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

#grab all the go terms
data(go.sets.hs)
data(go.subs.hs)

#specifically grab biological process terms
gobpsets=go.sets.hs[go.subs.hs$BP]

#kegg the GO bp results
gobpresult=gage(foldchanges,gsets=gobpsets,same.dir=TRUE)

#look at up and down regulated and the statistics
lapply(gobpresult,head)

#look at the results
View(gobpresult)