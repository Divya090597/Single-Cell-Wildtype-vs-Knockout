
#########  SingleR

### Cd36_ko

library(SingleR)

ref <- celldex::MouseRNAseqData()

results1 <- SingleR(test = as.SingleCellExperiment(Cd36_ko), ref = ref, labels = ref$label.main)

colData(ref)

Cd36_ko$SingleR_labels <- results1$labels

view(results1)

Cd36_ko[[]]

DimPlot(Cd36_ko, reduction = 'umap', group.by = 'SingleR_labels',label = TRUE)

##### Cd36_wt

results2 <- SingleR(test = as.SingleCellExperiment(Cd36_wt), ref = ref, labels = ref$label.main)

colData(ref)

Cd36_wt$SingleR_labels <- results2$labels

view(results2)

Cd36_wt[[]]

DimPlot(Cd36_wt, reduction = 'umap', group.by = 'SingleR_labels',label = TRUE)


########## mergeddata

DimPlot(mergeddata, reduction = 'umap', group.by = 'SingleR_labels',label = TRUE)

############ integrated_data

DimPlot(integrated_data, reduction = 'umap', group.by = 'SingleR_labels', label = TRUE)

####### Annotation Diagnostics ########

### Based on Scores within cells
results1$scores
plotScoreHeatmap(results1)

### Based on deltas across the cells

plotDeltaDistribution(results1)
