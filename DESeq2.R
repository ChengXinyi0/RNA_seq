#=======================preprocess
library(DESeq2)
library(XML)
#data_set #<- read.csv('1220featureCounts_all.csv', header = TRUE, row.names = 1)
setwd("/Users/chengxinyi/Desktop/test")
data_set <- read.csv('gene_matrix_count.csv', header = TRUE, row.names = 1)
any(is.na(data_set))
colnames(data_set) <- c('HER21', 'Normal2', 'NonTNBC3', 'TNBC1', 'NonTNBC2', 'Normal3', 
                        'TNBC2', 'HER22', 'Normal1', 'HER23', 'TNBC3', 'NonTNBC1')
data_set <- data_set[ , order(names(data_set))]

#sample <- c('HER21', 'HER22', 'HER23', 'NonTNBC1', 'NonTNBC2', 'NonTNBC3', 'Normal1', 'Normal2', 'Normal3', 'TNBC1', 'TNBC2', 'TNBC3') 
#type <- c('HER2', 'HER2', 'HER2','NonTNBC', 'NonTNBC', 'NonTNBC', 'Normal', 'Normal', 'Normal', 'TNBC', 'TNBC', 'TNBC')
coldata #<- cbind(sample, type)

dds_pre <- DESeq2::DESeqDataSetFromMatrix(countData = data_set,  colData = coldata, design = ~ type)
dds <- DESeq2::DESeq(dds_pre)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
resultpca <- plotPCA(vsd, intgroup=c("sample", "type"), returnData=TRUE)


dds_t <- vst(dds, blind=TRUE)
DESeq2::rlog(dds, blind = TRUE)

#=======================heatmap===================
library('pheatmap')
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
ntd <- normTransform(dds)
df <- as.data.frame(colData(dds)[,c("sample","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
#vsd
hp <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
#rld
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
#Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
hp
#Principal component plot of the samples
PCA <- plotPCA(vsd, intgroup=c("type"))#change 'condition' to 'sample'
#merge together
layout(matrix(c(3,1), nrow = 1), widths = c(2,1), heights = c(5,3))
layout()
layout.show(n)
plot_grid(hp, PCA,ncol = 2, nrow = 1)

#=====================get DEGs up/down in different pairs==================================================================

#Contrasts -> calculate the log2FC and padj. e.g.HER2 vs Normal 
con_HER2vsNor <- results(dds, contrast=c("type","HER2","Normal"))
#filtering NA padj
con_HER2vsNor_df <- as.data.frame(con_HER2vsNor)#convert to data.frame
con_HER2vsNor_df_na <- na.omit(con_HER2vsNor_df)#del the gene ID containing NA
#get degs' ENSEMBL ID, preparation for GO analysis
HER2vsNor_deg <- row.names(con_HER2vsNor_df_na[con_HER2vsNor_df_na$padj < 0.05, ])
                                               #& abs(con_HER2vsNor_df_na$log2FoldChange) > 1, ])

dim.data.frame(HER2vsNor_deg)#count the number of DEGs
HER2vsNor_deg_up <- row.names(con_HER2vsNor_df_na[con_HER2vsNor_df_na$padj < 0.05
                                               & (con_HER2vsNor_df_na$log2FoldChange) > 0, ])
dim.data.frame(HER2vsNor_deg_up)
HER2vsNor_deg_down <- row.names(con_HER2vsNor_df_na[con_HER2vsNor_df_na$padj < 0.05
                                                  & (con_HER2vsNor_df_na$log2FoldChange) < 0, ])
dim.data.frame((HER2vsNor_deg_down))
#same for other pairs
#NonTNBC vs Normal
con_NTvsNor <- results(dds, contrast=c("type","NonTNBC","Normal"))
all_gene_IDs <- row.names(con_NTvsNor)
con_NTvsNor_df <- as.data.frame(con_NTvsNor)
con_NTvsNor_df_na <- na.omit(con_NTvsNor_df)
NTvsNor_deg <- row.names(con_HER2vsNor_df_na[con_HER2vsNor_df_na$padj < 0.05, ])
#                                             & abs(con_HER2vsNor_df_na$log2FoldChange) > 1, ])
dim.data.frame(NTvsNor_deg)
NTvsNor_deg_up <- row.names(con_HER2vsNor_df_na[(con_HER2vsNor_df_na$padj < 0.05 
                                                 & con_HER2vsNor_df_na$log2FoldChange) > 0, ])
dim.data.frame(NTvsNor_deg_up)
NTvsNor_deg_down <- row.names(con_HER2vsNor_df_na[(con_HER2vsNor_df_na$padj < 0.05 
                                                 & con_HER2vsNor_df_na$log2FoldChange) < 0, ])
dim.data.frame(NTvsNor_deg_down)

#TNBC vs Normal
con_TvsNor <- results(dds, contrast=c("type","TNBC","Normal"))
con_TvsNor_df <- as.data.frame(con_NTvsNor)
con_TvsNor_df_na <- na.omit(con_TvsNor_df)
TvsNor_deg <- row.names(con_TvsNor_df_na[con_TvsNor_df_na$padj < 0.05, ])
#                                         & abs(con_TvsNor_df_na$log2FoldChange) > 1, ])
dim.data.frame(TvsNor_deg)
TvsNor_deg_up <- row.names(con_TvsNor_df_na[con_TvsNor_df_na$padj < 0.05 
                                            & (con_TvsNor_df_na$log2FoldChange) > 0, ])

dim.data.frame(TvsNor_deg_up)
TvsNor_deg_down <- row.names(con_TvsNor_df_na[con_TvsNor_df_na$padj < 0.05 
                                            & (con_TvsNor_df_na$log2FoldChange) < 0, ])
dim.data.frame(TvsNor_deg_down)

#TNBC vs non TNBC
con_TvsNT <- results(dds, contrast=c("type","TNBC","NonTNBC"))
con_TvsNT_df <- as.data.frame(con_TvsNT)
con_TvsNT_df_na <- na.omit(con_TvsNT_df)
TvsNT_deg <- row.names(con_TvsNT_df_na[con_TvsNT_df_na$padj < 0.05, ])
#                                       & abs(con_TvsNT_df_na$log2FoldChange) > 1, ])
dim.data.frame(TvsNT_deg)
TvsNT_deg_up <- row.names(con_TvsNT_df_na[con_TvsNT_df_na$padj < 0.05
                                       & (con_TvsNT_df_na$log2FoldChange) > 0, ])
dim.data.frame(TvsNT_deg_up)
TvsNT_deg_down <- row.names(con_TvsNT_df_na[con_TvsNT_df_na$padj < 0.05
                                          & (con_TvsNT_df_na$log2FoldChange) < 0, ])
dim.data.frame(TvsNT_deg_down)
DESeq2::counts(TvsNT_deg)
#HER2 vs NonTNBC
con_HER2vsNT_df <- as.data.frame(con_HER2vsNT)
con_HER2vsNT_df_na <- na.omit(con_HER2vsNT_df)
HER2vsNT_deg <- row.names(con_HER2vsNT_df_na[con_HER2vsNT_df_na$padj < 0.05, ])
#                                          & abs(con_HER2vsNT_df_na$log2FoldChange) > 1, ])
dim.data.frame(HER2vsNT_deg)
HER2vsNT_deg_up <- row.names(con_HER2vsNT_df_na[con_HER2vsNT_df_na$padj < 0.05 
                                                & (con_HER2vsNT_df_na$log2FoldChange) > 0, ])
dim.data.frame(HER2vsNT_deg_up)
HER2vsNT_deg_down <- row.names(con_HER2vsNT_df_na[con_HER2vsNT_df_na$padj < 0.05 
                                                & (con_HER2vsNT_df_na$log2FoldChange) < 0, ])
dim.data.frame(HER2vsNT_deg_down)

#HER2 vs T
con_HER2vsT_df <- as.data.frame(con_HER2vsT)
con_HER2vsT_df_na <- na.omit(con_HER2vsT_df)
HER2vsT_deg <- row.names(con_HER2vsT_df_na[con_HER2vsT_df_na$padj < 0.05, ])
#                                          & abs(con_HER2vsT_df_na$log2FoldChange) > 1, ])
dim.data.frame(HER2vsT_deg)
HER2vsT_deg_up <- row.names(con_HER2vsT_df_na[con_HER2vsT_df_na$padj < 0.05 
                                              & (con_HER2vsT_df_na$log2FoldChange) > 0, ])
dim.data.frame(HER2vsT_deg_up)
HER2vsT_deg_down <- row.names(con_HER2vsT_df_na[con_HER2vsT_df_na$padj < 0.05 
                                              & (con_HER2vsT_df_na$log2FoldChange) < 0, ])
dim.data.frame(HER2vsT_deg_down)
#=====================volcano diagram=============================================================================
library(EnhancedVolcano)
#volcano diagram, e.g. HER2 vs Normal, but not pretty :(
plot(con_HER2vsNor_df_na$log2FoldChange, -log10(con_HER2vsNor_df_na$padj), pch = 16, main = 'HER2vsNormal',
     cex.lab=1.5,xlab = "log2FoldChange",ylab = "-log10Pvalue")

#pretty volcano diagram

H_Tv <- EnhancedVolcano(con_HER2vsT_df_na,
                        lab = rownames(con_HER2vsT_df_na),
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        title = "HER2 vs TNBC",
                        subtitle = 'DEGs: 4942, up: 1848, down: 3094',
                        colAlpha = 0.5,
                        legendPosition = 'bottom',
                        legendLabSize = 12,
                        legendIconSize = 4.0,
                        labSize = 1)
H_NTv <- EnhancedVolcano(con_HER2vsNT_df_na,
                         lab = rownames(con_HER2vsNT_df_na),
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         title = "HER2 vs Non TNBC",
                         subtitle = 'DEGs: 4721, up: 2368, down: 2353',
                         colAlpha = 0.5,
                         legendPosition = 'bottom',
                         legendLabSize = 12,
                         legendIconSize = 4.0,
                         labSize = 1)
H_Nv <- EnhancedVolcano(con_HER2vsNor_df_na,
                        lab = rownames(con_HER2vsNor_df_na),
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        title = "HER2 vs Normal",
                        subtitle = 'DEGs: 14397, up: 9436, down: 4961',
                        colAlpha = 0.5,
                        legendPosition = 'bottom',
                        legendLabSize = 12,
                        legendIconSize = 4.0,
                        labSize = 1)
NT_Nv <- EnhancedVolcano(con_NTvsNor_df_na,
                         lab = rownames(con_NTvsNor_df_na),
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         title = "Non TNBC vs Normal",
                         subtitle = 'DEGs: 14397, up: 14397, down: 0',
                         colAlpha = 0.5,
                         legendPosition = 'bottom',
                         legendLabSize = 12,
                         legendIconSize = 4.0,
                         labSize = 1)
T_NTv <- EnhancedVolcano(con_TvsNT_df_na,
                         lab = rownames(con_TvsNT_df_na),
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         title = "TNBC vs Non TNBC",
                         subtitle = 'DEGs: 1730, up: 807, down: 923',
                         colAlpha = 0.5,
                         legendPosition = 'bottom',
                         legendLabSize = 12,
                         legendIconSize = 4.0,
                         labSize = 1)
T_Nv <- EnhancedVolcano(con_TvsNor_df_na,
                        lab = rownames(con_TvsNor_df_na),
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        title = "TNBC vs Normal",
                        subtitle = 'DEGs: 12525, up: 8570, down: 3955',
                        colAlpha = 0.5,
                        legendPosition = 'bottom',
                        legendLabSize = 12,
                        legendIconSize = 4.0,
                        labSize = 1)
plot_grid(H_Nv, T_Nv, NT_Nv, H_NTv, H_Tv, T_NTv, ncol = 3, nrow = 2)

#example
EnhancedVolcano(con_HER2vsT_df_na,
                lab = rownames(con_HER2vsT_df_na),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "HER2 vs TNBC",
                subtitle = element_blank(),
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                pointSize = 1,
                labSize = 4)
#p < 0.05, logFC > 1
points(con_HER2vsNor_df_na[con_HER2vsNor_df_na$log2FoldChange > 1 & con_HER2vsNor_df_na$padj < 0.05,]$log2FoldChange,
       pch = 16, -log10(con_HER2vsNor_df_na[con_HER2vsNor_df_na$log2FoldChange > 1 & con_HER2vsNor_df_na$padj < 0.05,]$padj), col = 'red')



#p < 0.05, logFC < -1
points(con_HER2vsNor_df_na[con_HER2vsNor_df_na$log2FoldChange < -1 & con_HER2vsNor_df_na$padj < 0.05,]$log2FoldChange,
       pch = 16, -log10(con_HER2vsNor_df_na[con_HER2vsNor_df_na$log2FoldChange < -1 & con_HER2vsNor_df_na$padj < 0.05,]$padj), col = 'green')

abline(h= -log10(0.05),v=c(-1,1))
#still need some change
EnhancedVolcano(con_HER,
                lab = rownames(con_HER2vsNor_df_na),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "HER2 vs Normal",
                subtitle = element_blank(),
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                pointSize = 1,
                labSize = 4)


#=====================GO analysis================
GO_database <- 'org.Hs.eg.db'
library(clusterProfiler)
#GO analysis with ENSEMBL
GO_NTvsNor <- enrichGO( 
  gene = NTvsNor_deg,
  OrgDb = GO_database,
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  universe = all_gene_IDs,
  qvalueCutoff = 0.05,
  readable = T
  )
head(GO_NTvsNor)
#GO analysis with ENSEMBL TNBC vs Normal
GO_TvsNor <- enrichGO(    
  gene = TvsNor_deg,
  OrgDb = GO_database,
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = T
  )
head(GO_TvsNor)

GO_TvsNT <- enrichGO(    
  gene = TvsNT_deg,
  OrgDb = GO_database,
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = T
  )
head(GO_TvsNT)

GO_HER2vsNT <- enrichGO(
  gene = HER2vsNT_deg,
  OrgDb = GO_database,
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = T
  )
head(GO_HER2vsNT)

GO_HER2vsT <- enrichGO(
  gene = HER2vsT_deg,
  OrgDb = GO_database,
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = T
)
head(GO_HER2vsT)
GO_HER2vsNor <- enrichGO(
  gene = HER2vsNor_deg,
  OrgDb = GO_database,
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = T
)
head(GO_HER2vsNor)
#GO_result <- [GO_NTvsNor_2, GO_TvsNor, GO_TvsNT, GO_HER2vsNT, HER2vsT, GO_HER2vsNor]

#GO plot
library(ggplot2)
library(cowplot)
#barpolt
GO_NTvsNor_f <- barplot(GO_NTvsNor, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis barplot: NonTNBC vs Normal')#barplotGO
GO_TvsNor_f <- barplot(GO_TvsNor, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis barplot: TNBC vs Normal')
GO_TvsNT_f <- barplot(GO_TvsNT, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis barplot: TNBC vs NonTNBC')
GO_HER2vsNT_f <- barplot(GO_HER2vsNT, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis barplot: HER2 vs NonTNBC')
GO_HER2vsT_f <- barplot(GO_HER2vsT, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis barplot: HER2 vs TNBC')
GO_HER2vsNor_f <- barplot(GO_HER2vsNor, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis barplot: HER2 vs Normal')
plot_grid(GO_NTvsNor_f,GO_TvsNor_f,GO_TvsNT_f,GO_HER2vsNT_f,GO_HER2vsT_f,GO_HER2vsNor_f,ncol = 3, nrow = 2)


#dotplot
GO_NTvsNor_d <- dotplot(GO_NTvsNor, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis dotplot: NonTNBC vs Normal')
GO_TvsNor_d <- dotplot(GO_TvsNor, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis dotplot: TNBC vs Normal')
GO_TvsNT_d <- dotplot(GO_TvsNT, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis dotplot: TNBC vs NonTNBC')
GO_HER2vsNT_d <- dotplot(GO_HER2vsNT, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis dotplot: HER2 vs NonTNBC')
GO_HER2vsT_d <- dotplot(GO_HER2vsT, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis dotplot: HER2 vs TNBC')
GO_HER2vsNor_d <- dotplot(GO_HER2vsNor, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") + ggtitle('GO analysis dotplot: HER2 vs Normal')
plot_grid(GO_HER2vsT_f, GO_HER2vsT_d,ncol = 2, nrow = 1)

#just plot
GO_BP_NTvsNor <- enrichGO( NTvsNor_deg,#GO:BP
                 OrgDb = GO_database,
                 keyType = "ENSEMBL",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_BP_NTvsNor)#not good

#have no idea
selected_genst=tT[tT$P.Value<0.05&abs(tT$logFC)>1,]
dim(selected_genst)
cC<-read.csv(file.choose(),stringsAsFactors = F)
selected_gensc=cC[tT$P.Value<0.05&abs(tT$logFC)>1,]
aA<-read.csv(file.choose(),stringsAsFactors = F)
selected_gensa=aA[tT$P.Value<0.05&abs(tT$logFC)>1,]
write.table (data, file ="C:\\Users\\Administrator\\Desktop\\resultCsv.csv", sep =",", row.names =FALSE)
write.csv()



####*customize the PCA plot using the ggplot function but failed :(
pcaData <- plotPCA(vsd, intgroup=c("sample", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()







