

#################################
#Load Required Packages and Data#
#################################

library(scater)
library(Seurat)
library(patchwork)
library(slingshot)
library(uwot)
library(RColorBrewer)
library(grDevices)
library(scales)
library(BiocGenerics)
library(tradeSeq)
library(dplyr)
library(pheatmap)


#load data 
all.merge <- readRDS("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/synovial-b-cells/harmony.rds")

#find variable features 
all.merge <- FindVariableFeatures(all.merge, selection.method = "vst", nfeatures = 1000)

#remove Id VDJ genes from variable features list
igtcr.index <- grep(pattern = "IG[HKL][VDJ]|AC233755.1|IGH[GMDEA]|IGKC|IGLC|IGLL|TR[ABGD][CV]", x = VariableFeatures(all.merge), value = FALSE) 
VariableFeatures(all.merge) <- VariableFeatures(all.merge)[-igtcr.index]

#convert Seurat object to Single Cell Experiment
all.sce <- as.SingleCellExperiment(all.merge)

##############
#Run Sligshot#
##############

all.sce <- slingshot(all.sce, clusterLabels = all.sce$clusters, reducedDim = 'UMAP', stretch = 2, approx_points = 150, extend = 'n')


###############
#Plot Lineages#
###############

#adjust cluster colours to match Seurat
seurat_colours <- c("#ED68ED", "#ABA300","#00BFC4", "#FF61CC", "#0CB702", "#00A9FF", "#00C19A", "#00B8E7", "#E68613", "#F8766D", "#7CAE00", "#C77CFF")
clus <- all.sce@colData@listData[["clusters"]]
colvec <- seurat_colours[factor(clus)]

#plot minimum spanning tree on UMAP
pdf("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/slingshot/dot_lineages_umap.pdf", width = 10, height = 8)
plot(reducedDims(all.sce)$UMAP, col = colvec, pch=16, asp = 1, cex = 0.4)
lines(SlingshotDataSet(all.sce), lwd=2,  type = 'lineages', col='black')
dev.off()

#plot lineage 1 as smoothed curve over the UMAP
pdf("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/slingshot/lineage1_umap.pdf", width = 9, height = 6)
plot(reducedDims(all.sce)$UMAP, col = colvec, pch=16, asp = 1, cex = 0.2)
lines(SlingshotDataSet(all.sce), linInd = 1, lwd=2, col='black')
title("Lineage 1")
dev.off()

#plot lineage 2 as smoothed curve over the UMAP
pdf("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/slingshot/lineage2_umap.pdf", width = 9, height = 6)
plot(reducedDims(all.sce)$UMAP, col = colvec, pch=16, asp = 1, cex = 0.2)
lines(SlingshotDataSet(all.sce), linInd = 2, lwd=2, col='black')
title("Lineage 2")
dev.off()

#plot lineage 1 as smoothed curve over the UMAP
pdf("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/slingshot/lineage3_umap.pdf", width = 9, height = 6)
plot(reducedDims(all.sce)$UMAP, col = colvec, pch=16, asp = 1, cex = 0.2)
lines(SlingshotDataSet(all.sce), linInd = 3, lwd=2, col='black')
title("Lineage 3")
dev.off()

#plot lineage 1 with pseudotime colour
pdf("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/slingshot/lineage1_pseudotime_umap.pdf", width = 6.5, height = 6.5)
plot(reducedDims(all.sce)$UMAP, asp=1, col = 'grey75', pch = .5)
points(reducedDims(all.sce)$UMAP, col = hcl.colors(100)[cut(-slingPseudotime(all.sce)[,1], 100)] )
lines(SlingshotDataSet(all.sce), linInd = 1, lwd=2, col='black')
title("Lineage 1")
dev.off()


##############
#Run tradeseq#
##############

#generate list of genes to consider from top variable features
counts <- all.merge@assays[["RNA"]]@var.features

# fit negative binomial GAM - takes several hours
all.sce <- fitGAM(all.sce, genes = counts)

# test for dynamic expression
ATres <- associationTest(all.sce, lineages=TRUE)

#save SCE
saveRDS(all.sce, "all_sce_fitGAM.rds")

all.sce <- readRDS("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/slingshot/all_sce_fitGAM.rds")


#plot heatmap
#code addapted from
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7612943/pdf/EMS146188.pdf

slingpseudotime1 <- all.sce@colData@listData[["slingPseudotime_1"]]
slingpseudotime1clusters <- all.sce@colData@listData[["clusters"]]
pseudotime_1 <- data.frame(slingpseudotime1, slingpseudotime1clusters, stringsAsFactors=FALSE)
pseudotime_1 <- na.omit(pseudotime_1)
pseudotime_1 <- arrange(pseudotime_1, slingpseudotime1)

# reduce number of points along the curves (via interpolation).
# for each pseudotime point set the cluster identity to the majority cluster for that time point
timepoints <- seq(0, max(pseudotime_1$slingpseudotime1), length.out=65)
timestep <- timepoints[2] - timepoints[1]

Mode <- function(x) {
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

clusters_interp <- data.frame(cluster=factor(rep(NA, length(timepoints)), levels=unique(pseudotime_1$slingpseudotime1clusters)), row.names=seq(length(timepoints)))

for (i in seq(length(timepoints))) {
  indx <- between(pseudotime_1$slingpseudotime1, timepoints[i]-timestep/2, timepoints[i]+timestep/2)
  if (any(indx)) {
    clusters_interp$cluster[i] <- Mode(pseudotime_1$slingpseudotime1clusters[indx])
  }
}

clusters_interp <- na.omit(clusters_interp)

loess_curves <- select(pseudotime_1, -slingpseudotime1clusters) %>% distinct()
curves_interp <- apply(loess_curves[,-1], 2, function(x) approx(pseudotime_1$slingpseudotime1, x, xout=timepoints)$y) %>% as.matrix() %>% t()
rownames(clusters_interp) <- colnames((-yhatSmooth1))[150:101]

# annotation colors for clusters
annoCol<-list(cluster=c('Naive 1'="#00B8E7", 'Naive 2'="#E68613", 'Non-Switched Memory'="#F8766D", 'Switched Memory'="#7CAE00", 'DN2'= "#0CB702", 'ASC 3' = "#00BFC4", 'ASC 2'="#ABA300"))


#identify genes to plot
lineage1Genes <-  rownames(ATres)[
  which(p.adjust(ATres$pvalue_1, "fdr") <= 0.05)]
topgenes1 <- lineage1Genes[1:250]

#construct heatmap
yhatSmooth1 <- predictSmooth(all.sce, gene = topgenes1, nPoints = 50, tidy = FALSE)
yhatSmooth1 <- yhatSmooth1[,order(ncol(yhatSmooth1):1)]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth1[, 101:150]))),
                       cluster_cols = FALSE,
                       show_rownames = TRUE, 
                       show_colnames = FALSE,
                       cutree_rows = 5,
                       annotation=clusters_interp,
                       annotation_colors = annoCol)

#save heatmap
pdf("~/Documents/Edinburgh/Year 2/scRNA_seq/analysis/slingshot/heatmap_1_smooth_with_clusters.pdf", width=7, height=7)
heatSmooth
dev.off()











