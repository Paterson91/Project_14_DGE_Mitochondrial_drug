#!/usr/bin/env Rscript
setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #14 Nanna_therapeutics")
args = commandArgs(trailingOnly = TRUE)

#Variables
#gene_counts_in = args[1]
#metadata_in = args[2]
#design_in = args[3]

#source("config_phasecomp.R")

#gene_counts_in = "Counts_stranded_gtf.gene_tidied"
#metadata_in = "metadata_omit.csv"
#design_in = "~breed + group"

####################################################################
########################## Input error check #######################
####################################################################

# test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("\n\nPlease provide a name for gene_counts_in as part of Config.R. \nThe config.R file should specify the following input files;\n\n1.\tGene Count table\n2.\tMetadata table\n3.\tDesign formula\n\n", call.=FALSE)
# } else {
#   source(as.character(args[1]))
# }

####################################################################
########################## Input formatting ########################
####################################################################

sampleFiles<- "Counts_stranded_gtf.gene_tidied"

countData <- as.matrix(read.table(sampleFiles, header=TRUE,
                                  row.names = 1, sep = "\t",
                                  as.is=TRUE,
                                  check.names=FALSE))


#countData = gsub(" ", "", countData)

#colClean <- function(x){ colnames(x) <- gsub("\\..*", "", colnames(x)); x } 
#countData <- colClean(countData)

cat("\n\nHead countData input: \n\n")

head(countData)

sampleNames <- as.data.frame(colnames(countData))
colnames(sampleNames)="SampleNames"

metadata = read.csv("metadata.csv", header = TRUE)

cat("\n\nMetadata input: \n\n")
metadata

order1_all = as.character(metadata$SampleNames)
iso_table_all = countData[, match(order1_all, colnames(countData))]

head(iso_table_all)


#colData = merge(x=sampleNames, y=metadata, by.x="SampleNames", by.y="Sample")

dir.create("DESeq2_output")
setwd("DESeq2_output")

design_in = "~ Group"
design_in_name = gsub("^[[:punct:]]", "", design_in)

# dir.create(design_in_name, showWarnings = FALSE)

####################################################################
########################## DGE - Pairwise ##########################
####################################################################
 
# setwd(paste0("./", as.character(design_in_name)))

library("DESeq2")
ddsHTSeq = DESeqDataSetFromMatrix(iso_table_all, metadata, design = as.formula(design_in), ignoreRank = FALSE)
#ddsHTSeq$phase2_cat=relevel(ddsHTSeq$phase2_cat, ref = reference_in)
levels(ddsHTSeq$Group)

 # Guts of DESeq2
dds <- DESeq(ddsHTSeq)
keep <- rowSums(counts(dds)) >= 1
ddskeep = dds[keep,]
# 
colData(dds)
# 
dir.create("DESeq2_Pairwise_OutputTables", showWarnings = FALSE)
# 
res1 = results(ddskeep,contrast=c("Group","NTX006","DMSO"), alpha = 0.05)
head(res1)
write.csv(res1[order(res1["padj"]),], file = "DESeq2_Pairwise_OutputTables/DESeq2_DMSO_vs_NTX006.csv")

res2 = results(ddskeep,contrast=c("Group","NTX031","DMSO"), alpha = 0.05)
head(res2)
write.csv(res2[order(res2["padj"]),], file = "DESeq2_Pairwise_OutputTables/DESeq2_DMSO_vs_NTX031.csv")

res3 = results(ddskeep,contrast=c("Group","NTX418","DMSO"), alpha = 0.05)
head(res3)
write.csv(res3[order(res3["padj"]),], file = "DESeq2_Pairwise_OutputTables/DESeq2_DMSO_vs_NTX418.csv")

res4 = results(ddskeep,contrast=c("Group","NTX421","DMSO"), alpha = 0.05)
head(res4)
write.csv(res4[order(res4["padj"]),], file = "DESeq2_Pairwise_OutputTables/DESeq2_DMSO_vs_NTX421.csv")

res5 = results(ddskeep,contrast=c("Group","NTX008","DMSO"), alpha = 0.05)
head(res5)
write.csv(res5[order(res5["padj"]),], file = "DESeq2_Pairwise_OutputTables/DESeq2_DMSO_vs_NTX008.csv")

res6 = results(ddskeep,contrast=c("Group","NTX011","DMSO"), alpha = 0.05)
head(res6)
write.csv(res6[order(res6["padj"]),], file = "DESeq2_Pairwise_OutputTables/DESeq2_DMSO_vs_NTX011.csv")


# send normalized counts to tab delimited file for GSEA, etc.
# 
outputPrefix = "All"
write.csv(as.data.frame(counts(ddskeep, normalized=T)), file = "All_normalized_counts_DESeq2.csv")



####################################################################
############################## Plot ################################
####################################################################

# # Plot expression of gene with smallest P-Adj value per sample over time - 0.05
# 
 library(ggplot2)

top = head(res1_lrt[order(res1_lrt$padj),], 1)
title = rownames(top)
fiss <- plotCounts(dds_lrt, gene=rownames(top), 
                   intgroup = "phase2_cat", returnData = T)
most_sig = ggplot(fiss,
            aes(x = phase2_cat, y = count, color = phase2_cat, group = 1)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  theme_bw() +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Output_Figures/Most_significant_LRT.jpg", plot = most_sig)
dev.off()


gene_to_test = "CD8A"
gene = plotCounts(dds,gene_to_test, returnData = T, normalized = T, intgroup = "phase2_cat")

plot = ggplot(gene,
          aes(x = phase2_cat, y = count, color = phase2_cat, group = 1)) + 
          geom_point() + stat_summary(fun=mean, geom="line")  + 
          theme_bw() +
          ggtitle(gene_to_test) +
          theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("Output_Figures/",gene_to_test,".jpg"), plot = plot)
dev.off()

####################################################################
############################## PCA #################################
####################################################################

dir.create("Output_Figures", showWarnings = FALSE)

# transform raw counts into normalized values. variance stabilization - variance stabilization is very good for heatmaps, etc.
vsd <- varianceStabilizingTransformation(dds, blind=T)

#PCA Plots
# ,shape=Phenotype
library(ggplot2)
PCA=plotPCA(vsd, intgroup = c("Group","Phenotype"), ntop = 100000, returnData = TRUE)
percentVar <- round(100 * attr(PCA, "percentVar"))
PCA2 = ggplot(data = PCA, aes(PC1, PC2, color = Group,shape=Phenotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + ggtitle("Principal Component Analysis \n") +
  theme(plot.title = element_text(hjust = 0.5))
PCA2

ggsave("Output_Figures/All_PCA_DESeq2.jpg", plot = PCA2)
dev.off()

###### Additional Analysis

####################################################################
######################## Volcano Plot ##############################
####################################################################

library(EnhancedVolcano)

#VOLCANO PLOT OF DIFFERENTIALLY ABUNDANT OTUS BETWEEN IC AND IN INFECTED PIGS
#LABELS TOP 20 MOST SIGNIFICANTLY DIFFERNTIALLY ABUNDANT OTUS BY OTU NAME; LABEL AT LOWEST KNOWN TAXONOMIC RANK (CAN GET AMBIGOUS AT LOW LEVELS)

res1[order(res1["padj"]),]

volcano1 = EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                legendLabels = c('NS', expression(Log[2]~FC),"adj p-value", expression(adj~p-value~and~log[2]~FC)),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                selectLab = rownames(res1)[res1$padj <= 0.0003], #label 20 most signifcant ASVs
                drawConnectors = TRUE, 
                title = "\nMost Significant Hits",
                subtitle = "",
                xlim = c(-5, 5))
volcano1

ggsave("Output_Figures/Volcano_LRT.jpg", plot = volcano1)
dev.off()

####################################################################
############################ Heatmap ###############################
####################################################################

ColourStrategy = metadata$Phenotype
library("RColorBrewer")
library(viridis)
library("gplots")
library(rafalib)

n = 1000

# 100 top expressed genes with heatmap.2
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:n]
cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy))]
heatmap.2(assay(vsd)[select,], col=viridis(256), Colv=TRUE,
          scale="row", key=T, keysize=1, symkey=T, labCol = metadata$Group, ColSideColors=cols,lhei = c(1.5,5),
          density.info="none", trace="none",margins=c(6,5),
          cexCol=1.1, labRow=F,
          main=paste(n," Top Expressed Genes Heatmap"))
dev.copy(png, width = 1000,height = 1000, paste("Output_Figures/",n,"-Heatmap_DESeq2.png", sep=""))
dev.off()

####################################################################
################# Hierarchical Cluster Generation ##################
####################################################################
# #
# # library(factoextra)
# #
# K_in = 5
# # clusters= hclust(dist(t(iso_table_all),method = "euclidean"), method = "complete")
# # plot(clusters, cex = 0.6, hang = -1)
# # rect.hclust(clusters, k = K_in, border = 2:5)
# # myplclust(clusters, labels=metadata$group, lab.col=as.fumeric(as.character(ColourStrategy)), cex=0.5)
# # dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-cluster_samples.png"))
# # #Optimise level to cluster
# # ClusterLevel = 5000
# # abline(h=ClusterLevel)
# # hclusters <- cutree(clusters, h=ClusterLevel)
# # table(true=metadata$group, cluster=hclusters)
# # hclusters <- cutree(clusters, k=3)
# # table(true=metadata$group, cluster=hclusters)
# # distance <- get_dist(t(iso_table_all))
# # fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
# #
# # select <- iso_table[order(rowMeans(iso_table),decreasing=T)[1:10],]
# # clusters= hclust(dist(select,method = "euclidean"), method = "complete")
# # plot(clusters, cex = 0.6, hang = -1)
# # rect.hclust(clusters, k = 6, border = 2:5)
# # dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-cluster_genes.png"))
# # dev.off()
# 
# #distance <- get_dist(iso_table)
# #fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
# 
####################################################################
################### K-Means Cluster Generation #####################
####################################################################
# library(factoextra)
# k = K_in
# iso_table2 = t(iso_table_all)
# iso_table3=iso_table2[,which(colSums(iso_table2)>1)]
# fviz_nbclust(iso_table3, kmeans, method = "gap_stat")
# #dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-K-means_optimisation.png"))
# dev.off()
# 
# kmeanscluster = kmeans(iso_table3, centers = k, nstart = 25)
# str(kmeanscluster)
# fviz_cluster(kmeanscluster, data = iso_table3)
# #dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-K-means_cluster.png"))
# dev.off()
# 
# k2 <- kmeans(iso_table3, centers = 2, nstart = 25)
# k3 <- kmeans(iso_table3, centers = 3, nstart = 25)
# k4 <- kmeans(iso_table3, centers = 4, nstart = 25)
# k5 <- kmeans(iso_table3, centers = 5, nstart = 25)
# 
# # plots to compare
# p1 <- fviz_cluster(k2, geom = "point", data = iso_table3) + ggtitle("k = 2")
# #dev.off()
# p2 <- fviz_cluster(k3, geom = "point",  data = iso_table3) + ggtitle("k = 3")
# #dev.off()
# p3 <- fviz_cluster(k4, geom = "point",  data = iso_table3) + ggtitle("k = 4")
# #dev.off()
# p4 <- fviz_cluster(k5, geom = "point",  data = iso_table3) + ggtitle("k = 5")
# #dev.off()
# p5 <- fviz_cluster(k5, geom = "point",  data = iso_table3) + ggtitle("k = 5")
# 
# library(gridExtra)
# grid.arrange(p1, p2, p3, p4, nrow = 2)




####################################################################
###################### Cluster Profiling ###########################
####################################################################
# 
# library(clusterProfiler)
# library(enrichplot)
# require("biomaRt")
# 
# #LOAD GENES AND ANNOTATE
# 
# n=11
# 
# mart <- useMart("ENSEMBL_MART_ENSEMBL")
# datasets = listDatasets(mart)
# datasets
#  
# mart <- useDataset("sscrofa_gene_ensembl", mart)
# ens = read.csv("Full_table.csv",header = T, row.names = 1)

# if (exists("annotLookup_table")) {
#   print("Lookup Table already created")
# } else {
#   annotLookup_table <- getBM(
#     mart=mart,
#     attributes=c("ensembl_gene_id", "external_gene_name","entrezgene_id"),
#     filter="ensembl_gene_id",
#     values=ens$gene,
#     uniqueRows=TRUE)
#   annotLookup_match = merge(x=ens, y=annotLookup_table, by.x="gene", by.y = "ensembl_gene_id", all.x = TRUE)
# }
# 
# 
# head(ens)
# head(annotLookup_table)
# head(annotLookup_match)
# 
# subset_pvalue = annotLookup_match[annotLookup_match$pvalue<0.05,]
# nrow(subset_pvalue)
# 
# GeneList_reads = subset_pvalue[,c(1,8:ncol(subset_pvalue))]
# GeneList_reads$external_gene_name=NULL
# GeneList_reads$gene_biotype=NULL
# GeneList_reads$gene=NULL
# 
# #ordering <- order(rowMeans(as.numeric(GeneList_reads)),decreasing=T)[1:n]
# #GeneList_reads_ordered = GeneList_reads[ordering,]
# #head(GeneList_reads_ordered)
# 
# 
# GeneList_fc = subset_pvalue[,c(1,3)]
# GeneList_fc$log2FoldChange = 2^GeneList_fc$log2FoldChange
# colnames(GeneList_fc)=c("Gene","FoldChange")
# GeneList_fc_order = GeneList_fc[order(-GeneList_fc$FoldChange),]
# head(GeneList_fc_order)
# nrow(GeneList_fc_order)
# 
# ### GO Classification
# 
# #Set levels
# #Set choice of; CC - Cellular Component, MF - Molecular Function, BP - Biological Process
# 
# ont_selection = "CC"
# level_selection = 7
# 
# ggo1 <- groupGO(gene     = as.character(GeneList_reads$entrezgene),
#                 keyType = "ENTREZID",
#                 OrgDb    = org.Mm.eg.db,
#                 ont      = ont_selection,
#                 level    = level_selection,
#                 readable = TRUE)
# barplot(ggo1, drop=TRUE, showCategory=12, border = c(10,10))
# dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_Go.png"))
# dev.off()
# 
# #Set levels
# #Set choice of; CC - Cellular Component, MF - Molecular Function, BP - Biological Process
# 
# ont_selection = "CC"
# level_selection = 7
# 
# de = GeneList_fc_order$Gene
# edo <- enrichDGN(de)
# 
# 
# 
# barplot(ggo1, drop=TRUE, showCategory=12, border = c(10,10))
# dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_Go.png"))
# dev.off()
# 
# 
# ### GO over-representation test
# 
# 
# ego <- enrichGO(gene          = as.character(annotLookup_match$ensembl_gene_id),
#                 OrgDb         = org.Mm.eg.db,
#                 keyType = "ENSEMBL",
#                 ont           = ont_selection,
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)
# head(ego)
# barplot(ego, showCategory=8)
# dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_OverRep.png"))
# dev.off()
# 
# ### GSEA test
# 
# geneList = as.numeric(GeneList_fc_order[,2])
# names(geneList) = as.character(GeneList_fc_order[,1])
# 
# ego3 <- gseGO(geneList     = geneList,
#               OrgDb        = org.Mm.eg.db,
#               ont          = "CC",
#               keyType ="ENSEMBL",
#               nPerm        = 1000,
#               minGSSize    = 100,
#               maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               verbose      = FALSE)
# 
# ridgeplot(ego3)
# 
# gseaplot(ego3, geneSetID = 2, title = ego3$Description[2])
# gseaplot2(ego3, geneSetID = 1:6)
# 
# 
# ### Pubmed Trend of enriched terms
# 
# terms <- ego3$Description[1:3]
# p <- pmcplot(terms, 2010:2017)
# p2 <- pmcplot(terms, 2010:2017, proportion=FALSE)
# plot_grid(p, p2, ncol=2)
# 
# pmcplot("Ifit1", 2000:2017)
# 
# # Visualisation
# 
# dotplot(ggo)
# emapplot(ego)
# goplot(ego)
# gseaplot(kk2, geneSetID = "hsa04145")
# ## categorySize can be scaled by 'pvalue' or 'geneNum'
# cnetplot(ego, categorySize="pvalue", foldChange=geneList)
# 
# 
# 
# 
# write.table(sessionInfo(), sep = "\t", file = "Session Info")


