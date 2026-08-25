####################
# Ruhland dataset
###################

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(pheatmap)
#BiocManager::install("apeglm")
library(apeglm)

df <- read.csv("../052_r_featureCounts_trimmed_Ruhland.counts.extraAttributes.txt", sep = "\t", skip = 1)
head(df)
names(df)



idx.num <- 14:19
sample_names <- c("EtOH_1", "EtOH_2", "EtOH_3", "TAM_1", "TAM_2", "TAM_3")
cbind(names(df[idx.num]), sample_names)
colnames(df)[idx.num] <- sample_names
names(df)

dim(df)
table(df$gene_biotype)

# Remove extra variance by subsetting for protein_coding
keep_gene <- which(df$gene_biotype == 'protein_coding') 
df2 <- df[keep_gene,]

dim(df2)

# Went from 78348 to 21818

# 1. Set your filtering rules
GroupSize <- 3
minReads <- 10

# 2. Check each row of your raw counts data frame
# (This counts how many samples in each row have at least 10 reads)
keep <- rowSums(df2[,idx.num] >= minReads) >= GroupSize

# 3. Filter the data frame to keep only those rows
filtered_df <- df2[keep, ]



#note: levels let's us define the reference levels
treatment <- factor( c(rep("EtOH",3), rep("TAM",3)), levels=c("EtOH", "TAM") )
colData <- data.frame(treatment, row.names = colnames(filtered_df)[idx.num])
colData

numM <- filtered_df[,idx.num]
rownames(numM) <- filtered_df$Geneid

dds <- DESeqDataSetFromMatrix(
  countData = numM, colData = colData, 
  design = ~ treatment)
dim(dds)

dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vsd, intgroup=c("treatment"))

#png(filename = "docs/assets/images/DESeq2_mouseMT/mouseMT_pca1.png", units = "in", height = 4, width = 5, res = 600)
pcaData + geom_label(aes(x=PC1,y=PC2,label=name))
#dev.off()

res <- results(dds)
summary(res)

## What to do with outliers?

outlier <- which(is.na(res$padj))

pheatmap::pheatmap(assay(vsd)[outlier,],
                   show_rownames = TRUE,
                   #labels_row = paste0(filtered_df$Geneid[outlier], " - ", filtered_df$gene_symbol[outlier]),
                   #filename = '../05_figures/05_heatCheck_outliers_scaled.png' ,
                   #width = 7, 
                   #height = 16,
                   scale = "row")  


# Remove outliers
filtered_df_noOutliers <- filtered_df[-outlier, ]

numM <- filtered_df_noOutliers[,idx.num]
rownames(numM) <- filtered_df_noOutliers$Geneid

dds <- DESeqDataSetFromMatrix(
  countData = numM, colData = colData, 
  design = ~ treatment)
dim(dds)

dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vsd, intgroup=c("treatment"))

#png(filename = "docs/assets/images/DESeq2_mouseMT/mouseMT_pca1.png", units = "in", height = 4, width = 5, res = 600)
pcaData + geom_label(aes(x=PC1,y=PC2,label=name))
#dev.off()

res <- results(dds)
summary(res)



head(coef(dds)) # the second column corresponds to the difference between the 2 conditions


FDRthreshold = 0.05
logFCthreshold = 1
# add a column of NAs

res.lfc <- lfcShrink(dds, coef=2, res=res)

res.lfc$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.01, set as "UP" 
res.lfc$diffexpressed[res.lfc$log2FoldChange > logFCthreshold & res.lfc$padj < FDRthreshold] <- "UP"
# if log2Foldchange < 1 and pvalue < 0.01, set as "DOWN"
res.lfc$diffexpressed[res.lfc$log2FoldChange < -logFCthreshold & res.lfc$padj < FDRthreshold] <- "DOWN"

plot <- ggplot( data = data.frame( res.lfc ) , aes( x=log2FoldChange , y = -log10(padj) , col =diffexpressed ) ) + 
  geom_point() + 
  geom_vline(xintercept=c(-logFCthreshold, logFCthreshold), col="red") +
  geom_hline(yintercept=-log10(FDRthreshold), col="red") +
  scale_color_manual(values=c("blue", "grey", "red"))
plot

table(res.lfc$diffexpressed)


#ggsave(plot, filename = "docs/assets/images/DESeq2_mouseMT/mouseMT_volcano.png", dpi = 600)



vsd.counts <- assay(vsd)

topVarGenes <- head(order(rowVars(vsd.counts), decreasing = TRUE), 20)
mat  <- vsd.counts[ topVarGenes, ] #scaled counts of the top genes

#png(filename = "docs/assets/images/DESeq2_mouseMT/mouseMT_heatmap.png", units = "in", height = 4, width = 6, res = 600)
pheatmap(mat,
         scale = "row")
#dev.off()



master <- cbind(filtered_df_noOutliers, res)
write.csv(master, file = "051_r_Ruhland.DESeq2.results.csv", row.names = FALSE)



library(clusterProfiler)
library(org.Mm.eg.db)

allGenes <- master$Geneid

genes_universe <- bitr(allGenes, fromType = "ENSEMBL",
                       toType = c("ENTREZID", "SYMBOL"),
                       OrgDb = "org.Mm.eg.db")
head(genes_universe)


sig <- master[master$padj < 0.05, ]
sigGenes <- sig$Geneid

genes_DE <- bitr(sigGenes, fromType = "ENSEMBL",
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = "org.Mm.eg.db")
head(genes_DE)


# GO "biological process (BP)" enrichment
ego_bp <- enrichGO(gene          = as.character(unique(genes_DE$ENTREZID)),
                   universe      = as.character(unique(genes_universe$ENTREZID)),
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
# couple of minutes to run

head(ego_bp)
dotplot(ego_bp)
# sample plot, but with adjusted p-value as x-axis
#dotplot(ego_bp, x = "p.adjust", showCategory = 20)

#BiocManager::install("ReactomePA")
library(ReactomePA)

# Reactome pathways enrichment
reactome.enrich <- enrichPathway(gene=as.character(unique(genes_DE$ENTREZID)),
                                 organism = "mouse",
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.9,
                                 readable=T,
                                 universe = genes_universe$ENTREZID)
# <1 minute to run


dotplot(reactome.enrich, x = "p.adjust")



