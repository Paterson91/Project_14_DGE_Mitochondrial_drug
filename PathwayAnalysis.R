setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #14 Nanna_therapeutics/")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("biomaRt")
# BiocManager::install("UpsetR")

library(clusterProfiler)
library(enrichplot)
require("biomaRt")
library(org.Hs.eg.db)
library(UpSetR)
library(dplyr)

# Load up biomart for H_Sapiens
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

# Load up data tables
input_files = list.files(pattern="*.csv", path = "DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables/", full.names = T)
input_files
all_input = lapply(input_files, read.csv,header = T)
names(all_input) <- gsub(".csv","",
                         list.files(pattern="*.csv", path = "DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables/", full.names = F),
                         fixed = TRUE)

require("biomaRt")

#Further Annotation

filters = listFilters(mart)
filters

annots = lapply(all_input, function(x) {
  y = getBM(
    mart=mart,
    attributes=c("external_gene_name","hgnc_symbol","ensembl_gene_id","entrezgene_id","description"),
    filter="ensembl_gene_id",
    values=x$X,
    uniqueRows=TRUE)
})

#Merge datasets with annotations

merged = Map(merge, x=all_input, y=annots, by.x = "X", by.y = "ensembl_gene_id", all.x=TRUE)

dir.create("DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables_annoted/")
write.csv(merged$DESeq2_DMSO_vs_NTX006, file = "DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables_annoted/DESeq2_DMSO_vs_NTX006_annotated.csv")
write.csv(merged$DESeq2_DMSO_vs_NTX008, file = "DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables_annoted/DESeq2_DMSO_vs_NTX008_annotated.csv")
write.csv(merged$DESeq2_DMSO_vs_NTX011, file = "DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables_annoted/DESeq2_DMSO_vs_NTX011_annotated.csv")
write.csv(merged$DESeq2_DMSO_vs_NTX031, file = "DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables_annoted/DESeq2_DMSO_vs_NTX031_annotated.csv")
write.csv(merged$DESeq2_DMSO_vs_NTX418, file = "DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables_annoted/DESeq2_DMSO_vs_NTX418_annotated.csv")
write.csv(merged$DESeq2_DMSO_vs_NTX421, file = "DESeq2_output_gene_id/DESeq2_Pairwise_OutputTables_annoted/DESeq2_DMSO_vs_NTX421_annotated.csv")


# Simplify
simplified = lapply(merged, function(x) {
  x[,c("entrezgene_id","log2FoldChange","padj")]
})
new_col_names = c("Entrez","L2FC","padj")
simplified = lapply(simplified, setNames, nm = new_col_names)


#Convert to fold change
# simplified = lapply(simplified, function(x) {
#   x$FC =(2^x$L2FC);return(x)
# })

# #Select p-sig values
# 
# sig = lapply(simplified, function(x) {
#   x[which(x$padj<0.05),]
# })

#Omit NA 
sig_na_omit = lapply(simplified, function(x) na.omit(x))

#Select only fold change
genelist = lapply(sig_na_omit, function(x){
  x[,2]
})

#Name each gene to fold change
names(genelist$DESeq2_DMSO_vs_NTX006) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX006$Entrez)
names(genelist$DESeq2_DMSO_vs_NTX008) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX008$Entrez)
names(genelist$DESeq2_DMSO_vs_NTX011) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX011$Entrez)
names(genelist$DESeq2_DMSO_vs_NTX031) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX031$Entrez)
names(genelist$DESeq2_DMSO_vs_NTX418) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX418$Entrez)
names(genelist$DESeq2_DMSO_vs_NTX421) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX421$Entrez)

#Sort into decreasing FC
genelist_dec = lapply(genelist, function(x){
  sort(x, decreasing = TRUE)
})

genelist_dec_name = lapply(genelist_dec, function(x){
  names(x)
})

wp2gene <- read.gmt("wikipathways-20200610-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

# Enricher - BASED ON ARBITRARY CUT OFF OF >2FC


### Wiki pathways
genelist_dec_2 = lapply(genelist_dec, function(x){
  names(x)[abs(x) > 2]
})

all_enricher = lapply(genelist_dec_2, function(x){
          enricher(x,
           TERM2GENE = wpid2gene, 
           TERM2NAME = wpid2name,
           pvalueCutoff = 0.05)
})

all_enricher_named = lapply(all_enricher, function(x){
  setReadable(x, org.Hs.eg.db, keyType = "ENTREZID")
})


barplot(all_enricher_named$DESeq2_DMSO_vs_NTX008)


### KEGG Pathway Analysis
dir.create("GSEA_KEGG_output")

all_KEGG_enrich = lapply(genelist_dec_2, function(x){
  enrichKEGG(gene         = x,
             organism     = 'hsa',
             pvalueCutoff = 0.05)
})

head(all_KEGG_enrich)

all_KEGG_GSEA = lapply(genelist_dec, function(x){
  gseKEGG(geneList     = x,
          organism     = 'hsa',
          nPerm        = 10000,
          minGSSize    = 120,
          pvalueCutoff = 1,
          verbose      = FALSE)
})

head(all_KEGG_GSEA)

p1 <- dotplot(all_KEGG_GSEA$DESeq2_DMSO_vs_NTX006, showCategory=25) + ggtitle("NTX006")
p1
p2 <- dotplot(all_KEGG_GSEA$DESeq2_DMSO_vs_NTX008, showCategory=25) + ggtitle("NTX008")
p3 <- dotplot(all_KEGG_GSEA$DESeq2_DMSO_vs_NTX011, showCategory=25) + ggtitle("NTX011")
p4 <- dotplot(all_KEGG_GSEA$DESeq2_DMSO_vs_NTX031, showCategory=25) + ggtitle("NTX031")
p5 <- dotplot(all_KEGG_GSEA$DESeq2_DMSO_vs_NTX418, showCategory=25) + ggtitle("NTX418")
p6 <- dotplot(all_KEGG_GSEA$DESeq2_DMSO_vs_NTX421, showCategory=25) + ggtitle("NTX421")

plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, align = "hv", label_x = 0.5)
ggsave("GSEA_KEGG_output/GSEA_dotplot.png")


all_KEGG_named_GSEA = lapply(all_KEGG_GSEA, function(x){
  setReadable(x, org.Hs.eg.db, keyType = "ENTREZID")
})

all_KEGG_named_enrich = lapply(all_KEGG_enrich, function(x){
  setReadable(x, org.Hs.eg.db, keyType = "ENTREZID")
})

p1 <- cnetplot(all_KEGG_named_GSEA$DESeq2_DMSO_vs_NTX006,  foldChange=genelist_dec$DESeq2_DMSO_vs_NTX006, circular = F) + ggtitle("Network plot of GSEA terms - NTX006") +
  theme(plot.title = element_text(face = "bold"))
p1
ggsave("GSEA_KEGG_output/Network plot of GSEA terms - NTX006.png")

p2 <- cnetplot(all_KEGG_named_GSEA$DESeq2_DMSO_vs_NTX008, foldChange=genelist_dec$DESeq2_DMSO_vs_NTX008) + ggtitle("Network plot of GSEA terms - NTX008") +
  theme(plot.title = element_text(face = "bold"))
p2
ggsave("GSEA_KEGG_output/Network plot of GSEA terms - NTX008.png")


p3 <- cnetplot(all_KEGG_named_GSEA$DESeq2_DMSO_vs_NTX011, foldChange=genelist_dec$DESeq2_DMSO_vs_NTX011) + ggtitle("Network plot of GSEA terms - NTX011") +
  theme(plot.title = element_text(face = "bold"))
p3
ggsave("GSEA_KEGG_output/Network plot of GSEA terms - NTX011.png")


p4 <- cnetplot(all_KEGG_named_GSEA$DESeq2_DMSO_vs_NTX031, foldChange=genelist_dec$DESeq2_DMSO_vs_NTX031) + ggtitle("Network plot of GSEA terms - NTX031") +
  theme(plot.title = element_text(face = "bold"))
p4
ggsave("GSEA_KEGG_output/Network plot of GSEA terms - NTX031.png")


p5 <- cnetplot(all_KEGG_named_GSEA$DESeq2_DMSO_vs_NTX418, foldChange=genelist_dec$DESeq2_DMSO_vs_NTX418) + ggtitle("Network plot of GSEA terms - NTX418") +
  theme(plot.title = element_text(face = "bold"))
p5
ggsave("GSEA_KEGG_output/Network plot of GSEA terms - NTX418.png")


p6 <- cnetplot(all_KEGG_named_GSEA$DESeq2_DMSO_vs_NTX421, foldChange=genelist_dec$DESeq2_DMSO_vs_NTX421) + ggtitle("Network plot of GSEA terms - NTX421") +
  theme(plot.title = element_text(face = "bold"))
p6
ggsave("GSEA_KEGG_output/Network plot of GSEA terms - NTX421.png")

upsetplot(all_KEGG_named_GSEA$DESeq2_DMSO_vs_NTX011)


gseaplot2(all_KEGG_named_GSEA$DESeq2_DMSO_vs_NTX418, geneSetID = 1:3, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")


















### Disease analysis

all_disease = lapply(genelist_dec_2, function(x){
  enrichDO(gene          = x,
            ont           = "DO",
            pvalueCutoff  = 0.05,
            pAdjustMethod = "BH",
            universe      = names(genelist_dec$DESeq2_DMSO_vs_NTX421),
            minGSSize     = 5,
            maxGSSize     = 500,
            qvalueCutoff  = 0.05,
            readable      = FALSE)
})


lapply(genelist_dec)

mapply(function(x,y){
  enrichDO(gene          = x,
           ont           = "DO",
           pvalueCutoff  = 0.05,
           pAdjustMethod = "BH",
           universe      = names(genelist_dec$DESeq2_DMSO_vs_NTX421),
           minGSSize     = 5,
           maxGSSize     = 500,
           qvalueCutoff  = 0.05,
           readable      = FALSE)
})



x <- enrichDO(gene          = genelist_dec_2$DESeq2_DMSO_vs_NTX421,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(genelist_dec$DESeq2_DMSO_vs_NTX421),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

x <- enrichDO(gene          = genelist_dec_2$DESeq2_DMSO_vs_NTX008,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(genelist_dec$DESeq2_DMSO_vs_NTX008),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)






















# GSEA - No significance found

all_GSEA = lapply(genelist_dec, function(x){
      GSEA(x,
       TERM2GENE = wpid2gene, 
       TERM2NAME = wpid2name,
       verbose=FALSE)
})

# enrich

edo = enrichDGN(genelist_dec_2$DESeq2_DMSO_vs_NTX006)

library(enrichplot)
barplot(edo, showCategory=20)












###
#UpsetR plot

id_only = lapply(sig_na_omit, function(x) {
  x[,c("Entrez")]
})

id_only$DESeq2_DMSO_vs_NTX008

lapply(Reduce(intersect, list(id_only)), head)


Reduce(intersect, list(id_only$DESeq2_DMSO_vs_NTX421, id_only$DESeq2_DMSO_vs_NTX421, id_only$DESeq2_DMSO_vs_NTX421))



upset(fromList(id_only), order.by = "freq",
      empty.intersections = "off",
      nsets = 6,
      main="TEST")
      # queries = list(list(query = intersects, params = list("DESeq2_DMSO_vs_NTX008")),
      #                     list(query = intersects, params = list("DESeq2_DMSO_vs_NTX011","DESeq2_DMSO_vs_NTX008")),
      #                     list(query = intersects, params = list("DESeq2_DMSO_vs_NTX011"))
                   
###



X=fromList(id_only)



venn.diagram(X, filename ="1.tiff", height = 1000, width = 1000)











#Gene ontology

ont_selection = "BP"
level_selection = 7

ggo = lapply(sig_na_omit, function(x) {
  groupGO(gene     = as.character(x$Entrez),
          keyType = "ENTREZID",
          OrgDb    = org.Hs.eg.db,
          ont      = ont_selection,
          level    = level_selection,
          readable = TRUE)
})

barplot(ggo$DESeq2_DMSO_vs_NTX421, drop=TRUE, showCategory=12, border = c(10,10))






### GO over-representation test

# Simplify
simplified_ensembl = lapply(merged, function(x) {
  x[,c("X","log2FoldChange","padj")]
})
new_col_names = c("Ensembl","L2FC","padj")
simplified_ensembl = lapply(simplified_ensembl, setNames, nm = new_col_names)
#Convert to fold change
simplified_ensembl = lapply(simplified_ensembl, function(x) {
  x$FC =(2^x$L2FC);return(x)
})
# #Select p-sig values
#
# sig = lapply(simplified, function(x) {
#   x[which(x$padj<0.05),]
# })
#Omit NA
sig_na_omit_ensembl = lapply(simplified_ensembl, function(x) na.omit(x))
#Select only fold change
genelist_ensembl = lapply(sig_na_omit_ensembl, function(x){
  x[,4]
})
#Name each gene to fold change
names(genelist_ensembl$DESeq2_DMSO_vs_NTX006) = as.character(sig_na_omit_ensembl$DESeq2_DMSO_vs_NTX006$Ensembl)
names(genelist_ensembl$DESeq2_DMSO_vs_NTX008) = as.character(sig_na_omit_ensembl$DESeq2_DMSO_vs_NTX008$Ensembl)
names(genelist_ensembl$DESeq2_DMSO_vs_NTX011) = as.character(sig_na_omit_ensembl$DESeq2_DMSO_vs_NTX011$Ensembl)
names(genelist_ensembl$DESeq2_DMSO_vs_NTX031) = as.character(sig_na_omit_ensembl$DESeq2_DMSO_vs_NTX031$Ensembl)
names(genelist_ensembl$DESeq2_DMSO_vs_NTX418) = as.character(sig_na_omit_ensembl$DESeq2_DMSO_vs_NTX418$Ensembl)
names(genelist_ensembl$DESeq2_DMSO_vs_NTX421) = as.character(sig_na_omit_ensembl$DESeq2_DMSO_vs_NTX421$Ensembl)
#Sort into decreasing FC
genelist_dec_ensembl = lapply(genelist_ensembl, function(x){
  sort(x, decreasing = TRUE)
})


genelist_dec_ensembl_up_reg = lapply(genelist_dec_ensembl, function(x){
  x[which(x>1.5)]
})


dir.create("DESeq2_output_gene_id/ClusterProfileR/Upreg/", showWarnings = F)

ont_selection = "CC"
level_selection = 7
dir.create(paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/", ont_selection), showWarnings = F)
ego_cc <-  lapply(genelist_dec_ensembl_up_reg, function(x) {
  enrichGO(gene          = names(x),
                OrgDb         = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont           = ont_selection,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
})

barplot(ego_cc$DESeq2_DMSO_vs_NTX006, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX006, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX008, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX008, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX011, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX011, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX031, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX031, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX418, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX418, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX421, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX421, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()



ont_selection = "MF"
level_selection = 7
dir.create(paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/", ont_selection), showWarnings = T)
ego_mf <-  lapply(genelist_dec_ensembl_up_reg, function(x) {
  enrichGO(gene          = names(x),
           OrgDb         = org.Hs.eg.db,
           keyType = "ENSEMBL",
           ont           = ont_selection,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

barplot(ego_mf$DESeq2_DMSO_vs_NTX006, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX006, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX008, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX008, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX011, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX011, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX031, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX031, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX418, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX418, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX421, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX421, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()


ont_selection = "BP"
level_selection = 7
dir.create(paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/", ont_selection), showWarnings = F)
ego_bp <-  lapply(genelist_dec_ensembl_up_reg, function(x) {
  enrichGO(gene          = names(x),
           OrgDb         = org.Hs.eg.db,
           keyType = "ENSEMBL",
           ont           = ont_selection,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

barplot(ego_bp$DESeq2_DMSO_vs_NTX006, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX006, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX008, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX008, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX011, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX011, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX031, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX031, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX418, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX418, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX421, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX421, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Upreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()




genelist_dec_ensembl_down_reg = lapply(genelist_dec_ensembl, function(x){
  x[which(x<1.5)]
})

dir.create("DESeq2_output_gene_id/ClusterProfileR/Downreg/", showWarnings = F)

ont_selection = "CC"
level_selection = 7
dir.create(paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/", ont_selection), showWarnings = F)
ego_cc <-  lapply(genelist_dec_ensembl_down_reg, function(x) {
  enrichGO(gene          = names(x),
           OrgDb         = org.Hs.eg.db,
           keyType = "ENSEMBL",
           ont           = ont_selection,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

barplot(ego_cc$DESeq2_DMSO_vs_NTX006, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX006, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX008, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX008, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX011, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX011, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX031, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX031, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX418, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX418, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_cc$DESeq2_DMSO_vs_NTX421, main="TEST", showCategory=10)
write.csv(ego_cc$DESeq2_DMSO_vs_NTX421, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()



ont_selection = "MF"
level_selection = 7
dir.create(paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/", ont_selection), showWarnings = T)
ego_mf <-  lapply(genelist_dec_ensembl_down_reg, function(x) {
  enrichGO(gene          = names(x),
           OrgDb         = org.Hs.eg.db,
           keyType = "ENSEMBL",
           ont           = ont_selection,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

barplot(ego_mf$DESeq2_DMSO_vs_NTX006, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX006, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX008, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX008, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX011, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX011, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX031, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX031, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX418, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX418, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_mf$DESeq2_DMSO_vs_NTX421, main="TEST", showCategory=10)
write.csv(ego_mf$DESeq2_DMSO_vs_NTX421, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()


ont_selection = "BP"
level_selection = 7
dir.create(paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/", ont_selection), showWarnings = F)
ego_bp <-  lapply(genelist_dec_ensembl_down_reg, function(x) {
  enrichGO(gene          = names(x),
           OrgDb         = org.Hs.eg.db,
           keyType = "ENSEMBL",
           ont           = ont_selection,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

barplot(ego_bp$DESeq2_DMSO_vs_NTX006, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX006, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX006_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX008, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX008, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX008_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX011, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX011, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX011_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX031, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX031, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX031_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX418, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX418, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX418_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()

barplot(ego_bp$DESeq2_DMSO_vs_NTX421, main="TEST", showCategory=10)
write.csv(ego_bp$DESeq2_DMSO_vs_NTX421, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.csv"))
dev.copy(png, width = 1000,height = 1000, paste0("DESeq2_output_gene_id/ClusterProfileR/Downreg/",  ont_selection, "/DESeq2_DMSO_vs_NTX421_", ont_selection, "_lvl_", level_selection, "_ENRICH.png"))
dev.off()




# Compare Clusters

#Select p-sig values

sig = lapply(simplified, function(x) {
  x[which(x$padj<0.05),]
})
#Omit NA
sig_na_omit = lapply(sig, function(x) na.omit(x))
#Select only fold change
genelist_entrez = lapply(sig_na_omit, function(x){
  x[,4]
})
#Name each gene to fold change
names(genelist_entrez$DESeq2_DMSO_vs_NTX006) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX006$Entrez)
names(genelist_entrez$DESeq2_DMSO_vs_NTX008) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX008$Entrez)
names(genelist_entrez$DESeq2_DMSO_vs_NTX011) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX011$Entrez)
names(genelist_entrez$DESeq2_DMSO_vs_NTX031) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX031$Entrez)
names(genelist_entrez$DESeq2_DMSO_vs_NTX418) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX418$Entrez)
names(genelist_entrez$DESeq2_DMSO_vs_NTX421) = as.character(sig_na_omit$DESeq2_DMSO_vs_NTX421$Entrez)
#Sort into decreasing FC
genelist_dec_ensembl = lapply(genelist_entrez, function(x){
  sort(x, decreasing = TRUE)
})

genelist_dec_clust = lapply(genelist_dec_ensembl, function(x){
  names(x)
})

TEST = genelist_dec_clust$DESeq2_DMSO_vs_NTX006

ck <- compareCluster(geneCluster = TEST, 
                     fun = "enrichKEGG",
                     organism = "hsa", 
                     pvalueCutoff=1)

ck
head(as.data.frame(ck))
dotplot(ck, showCategory = 20,)







genelist_up_reg = lapply(genelist_dec, function(x){
  x[which(x>1.5)]
})

genelist_down_reg = lapply(genelist_dec, function(x){
  x[which(x<0.5)]
})


data(geneList)
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))