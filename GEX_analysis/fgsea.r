#################################
#Load Required Packages and Data#
#################################

library(Seurat)
library(SeuratObject)
library(dplyr)
library(fgsea)
library(reactome.db)
library(msigdbr)
library(tibble)

#load data 
all.merge <- readRDS("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/synovial-b-cells/harmony.rds")

##############################################################
#Create Ranked List for the DN2 Cluster vs All Other Clusters#
##############################################################

dn2_volc <- FindMarkers(all.merge, ident.1 = "DN2", test.use = "MAST", only.pos = FALSE)

write.csv(dn2_volc, "~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/fgsea/ranks/dn2_volc.csv", row.names = TRUE)

#dn2_volc <- read.csv("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/fgsea/ranks/dn2_volc.csv")

dn2_volc_ranks <- dn2_volc$avg_log2FC
dn2_volc_gene <- dn2_volc$X
dn2_volc_ranks <- data.frame(dn2_volc_gene, dn2_volc_ranks, row.names = NULL)
dn2_volc_ranks <- dn2_volc_ranks %>% arrange(desc(dn2_volc_ranks))
dn2_volc_ranks <- deframe(dn2_volc_ranks)

###########
#Run fgsea#
###########

#for GO biological pathways and the ranked list of genes for the DN2 cluster vs all other clusters
set.seed(42)
dn2_bp <- fgsea(pathways=gmtPathways("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/fgsea/gmt/c5.go.bp.v7.5.1.symbols.gmt"), dn2_volc_ranks, nperm=1000)

#plot barplot for top 10 enriched pathways
pdf("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/fgsea/figures/dn2_bp_top10.pdf", width = 15, height = 5)
top_n(dn2_bp, n=10, NES) %>%
  ggplot(., aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO BP pathways NES from GSEA") + 
  theme_minimal()
dev.off()
