

# Import counts, sample table and annotation
A673.counts<-read.table("tablecounts_raw_RREB1_040822.txt", sep = "\t", header = TRUE)

colnames(A673.counts)[1]<-"gene.id"

sample.table<-read.table("sample description table RREB1_040822.txt", sep="\t", header=TRUE)

gene.names.annotation<-read.table("tableannot.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)





# Order both tables (sample table, gene counts) in ascending, matching order

sample.table2<-sample.table[order(sample.table$Sample.ID),]

sample.names<-as.character(sample.table2$Sample.ID)

A673.counts.ordered<-(A673.counts[sample.names])


#add row names
rownames(A673.counts.ordered)<-A673.counts$gene.id


#set sample table "Set" and "Condition" columns to factor

sample.table2$Set<-as.factor(sample.table2$Set)
sample.table2$Condition<-as.factor(sample.table2$Condition)


rownames(sample.table2)<-sample.table2$Sample.ID


sample.table3<-sample.table2[,c(2:4), drop=FALSE]



#check to see if the row names in sample table match with column names in counts table

all(rownames(sample.table3) == colnames(A673.counts.ordered))
#TRUE




#create the DESeq object

library("DESeq2")


dds<-DESeqDataSetFromMatrix(countData = A673.counts.ordered, colData = sample.table3, design = ~ Set + Condition)




########### For GSEA anaylsis, Create normalised counts table############# 


dds <- estimateSizeFactors(dds)

normalised_counts_a673<-counts(dds, normalized=TRUE)


# filter rows with <6 total counts (low/ non-expressed genes), for GSEA table


dim(normalised_counts_a673)
#[1] 57820     6

keep.norm<-rowSums(normalised_counts_a673) >= 6

normalised_counts_a673.filtered<-normalised_counts_a673[keep.norm,]

dim(normalised_counts_a673.filtered)

# 23100     6


all(row.names(sample.table3)==colnames(normalised_counts_a673.filtered))

colnames(normalised_counts_a673.filtered)<-sample.table3$Biological.name



# adding a value of 1 to all, to remove zeros (to avoid infinity scores with fold change calculations)

normalised_counts_a673.filtered<-normalised_counts_a673.filtered + 1

write.csv(normalised_counts_GSEA_a673.siCTvsiRREB1.si1_plus1, file="080822_a673_siCT_siRREB1.si1_DEseq_normalised_filtered_rowsums_6_plus_1.csv")









##### PCA visualisation ####

rld <- rlog(dds, blind = TRUE)

plotPCA(rld, intgroup=c("Condition"))+
theme_bw()





####Differential expression#####

# filter rows with <6 total counts (low/ non-expressed genes), for DESeq2 analysis

keep<-rowSums(counts(dds)) >= 6

dds<-dds[keep,]


#Build the model and carry out DE

#reorder conditions to set siCT as the reference level

dds$Condition<-relevel(dds$Condition, ref = "siCT")

dds<-DESeq(dds)

res<-results(dds)



#####plot RREB1 expression######


# Verify RREB1 knockdown 

library(ggplot2)


RREB1.plotdata<-plotCounts(dds, gene="ENSG00000124782", intgroup ="Condition", returnData = TRUE)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

ggplot(data=RREB1.plotdata, aes(x=Condition, y=count)) +
  ggtitle("RREB1 expression") +
  geom_point(position = position_jitter(w=0.12), size=1.8, colour="firebrick3")+
  theme_bw()+
  ylab("Normalised Read Count")+
  theme(plot.title = element_text(hjust = 0.5))



res_si1_vs_siCT<-results(dds, contrast = c("Condition", "si1", "siCT"))


#removing genes with NAs
res_si1_vs_siCT.2<-res_si1_vs_siCT[complete.cases(res_si1_vs_siCT),]




####Get DE genes for si1#########
## DE genes at p-value (padj)<0.05, and a fold change >1.5 

# Order genes based on log2fold change
si1.ordered<-res_si1_vs_siCT.2[order(res_si1_vs_siCT.2$log2FoldChange),]

si1.upregulated<-si1.ordered[si1.ordered$padj<0.05 & si1.ordered$log2FoldChange>0.585,]

dim(si1.upregulated)
#2350
 

si1.downregulated<-si1.ordered[si1.ordered$padj<0.05 & si1.ordered$log2FoldChange<(-0.585),]
 
dim(si1.downregulated)
 #1438


########annotate with gene names

#si1 upregulated
si1.upregulated$gene_id<-rownames(si1.upregulated)

si1.upregulated<-as.data.frame(si1.upregulated)

si1.upregulated<-merge(si1.upregulated, gene.names.annotation, by="gene_id")

si1.upregulated<-si1.upregulated[order(si1.upregulated$log2FoldChange, decreasing = TRUE),]

write.csv(si1.upregulated, "070822_A673_RREB1_si1_upregulated_genes_1.5fold_padj_0.05.csv")




#si1 downregulated
si1.downregulated$gene_id<-rownames(si1.downregulated)

si1.downregulated<-as.data.frame(si1.downregulated)

si1.downregulated<-merge(si1.downregulated, gene.names.annotation, by="gene_id")

si1.downregulated<-si1.downregulated[order(si1.downregulated$log2FoldChange),]

write.csv(si1.downregulated, "070822_A673_RREB1_si1_downregulated_genes_1.5fold_padj_0.05.csv")









#### Constructing a volcano plot #####

si1.res.datafr<-as.data.frame(res_si1_vs_siCT.2)


si1.res.datafr$Significance<-ifelse((si1.res.datafr$padj<0.05&si1.res.datafr$log2FoldChange<(-0.585))|(si1.res.datafr$padj<0.05&si1.res.datafr$log2FoldChange>(0.585)), "Significant", "Not Significant")


si1.volcano.plot<-ggplot(data=si1.res.datafr, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(color=Significance), show.legend=FALSE)+
  ggtitle("RREB1 knockdown \nsiRREB1#1 v/s siCT.1") +
  scale_color_manual(values=c("grey","dodgerblue3"))+
  xlim(-8, 8)+
  xlab("Log2 Fold change")+
  ylab("-Log10 FDR")


si1.volcano.plot



################## Heatmap of genes ############


library(gplots)
library(pheatmap)

normalised_counts_A673<-counts(dds, normalized=TRUE)

# Exract DE genes from normalised table for heatmap


si1_DEgenes<-normalised_counts_A673[c(si1.upregulated$gene_id,si1.downregulated$gene_id),]


col.namesA673<-as.character(sample.table3$Biological.name)

colnames(si1_DEgenes)<-col.namesA673

#reorder columns to group treatments

si1_DEgenes<-si1_DEgenes[,c(1,3,5,2,4,6)]


par(oma=c(5,2,2,2))

pheatmap(si1_DEgenes, scale = "row", color = colorRampPalette(c("blue3","white","red3"))(256), cluster_rows=FALSE, cluster_cols = FALSE, show_rownames = FALSE, cellwidth = 20)





################GO and KEGG analysis ##############

library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGG.db)



# GO analysis si DE genes 

si1.downreg.GObp<-enrichGO(gene = si1.downregulated$gene_id, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = "BP", 
                           pAdjustMethod =  "BH")

dotplot(si1.downreg.GObp, showCategory=10, title="GO:BP siRREB1 si1 downregulated", font.size=11)


si1.upreg.GObp<-enrichGO(gene = si1.upregulated$gene_id, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = "BP", 
                                   pAdjustMethod =  "BH")

dotplot(si1.upreg.GObp, showCategory=12, title="GO:BP siRREB1 si1 upregulated-\n relaxed DE threshold")



########### KEGG #########



si1.upregulated.ENTEREZ<-bitr(si1.upregulated$gene_id, fromType = 'ENSEMBL', 
                              toType = c('ENTREZID'), OrgDb = org.Hs.eg.db, drop = TRUE)


si1.upreg.KEGG<-enrichKEGG(gene = si1.upregulated.ENTEREZ$ENTREZID, organism = "hsa", keyType = 'ENTEREZID', use_internal_data = TRUE)


dotplot(si1.upreg.KEGG, showCategory=20, title="KEGG siRREB1 si1 upregulated")



si1.downregulated.ENTEREZ<-bitr(si1.downregulated$gene_id, fromType = 'ENSEMBL', 
                                      toType = c('ENTREZID'), OrgDb = org.Hs.eg.db, drop = TRUE)

si1.downregulated.KEGG<-enrichKEGG(gene = si1.downregulated.ENTEREZ$ENTREZID, organism = "hsa", keyType = 'ENTEREZID', use_internal_data = TRUE)


dotplot(si1.downregulated.KEGG, showCategory=20, title="KEGG siRREB1 si1 downregulated", font.size=11)





####### Overlap with EWSR1-FL1  DE genes ########

library(stringr)


EF1.DEgenes<-read.table(file = "RNAseq_se_paper_20190726_DE_celline_Yes_A673_High_Yes_A673_Low_p0.01_logfc1_exp10.txt", sep="\t", header = TRUE)

EF1.up<-EF1.DEgenes[EF1.DEgenes$DE=="Up",]

EF1.down<-EF1.DEgenes[EF1.DEgenes$DE=="Down",]

EF1.up<-EF1.up$Gene_id

EF1.down<-EF1.down$Gene_id

EF1.up<-gsub("\\..*$","", EF1.up)

EF1.down<-gsub("\\..*$","", EF1.down)


length(EF1.up)
#1836

length(EF1.down)
#1600

length(si1.upregulated$gene_id)
#2350

length(si1.downregulated$gene_id)
#1438


length(intersect(EF1.up, si1.upregulated$gene_id))
#479

length(intersect(EF1.down, si1.downregulated$gene_id))
#266



# Test statistical significance of overlap

library(GeneOverlap)

#57820 is number of total rows (genes) in the RNA -Seq counts table     
#21754 is number of rows after filtering low expression

gen.size<-23365

go.obj<-newGeneOverlap(EF1.up, si1.upregulated$gene_id, genome.size=gen.size)

go.obj<-testGeneOverlap(go.obj)

go.obj

# GeneOverlap object:
#   listA size=1836
# listB size=2350
# Intersection size=479
# Overlapping p-value=6.9e-95
# Jaccard Index=0.1



go.obj<-newGeneOverlap(EF1.down, si1.downregulated$gene_id, genome.size=gen.size)

go.obj<-testGeneOverlap(go.obj)

go.obj

# GeneOverlap object:
#   listA size=1600
# listB size=1438
# Intersection size=266
# Overlapping p-value=2e-53
# Jaccard Index=0.1



