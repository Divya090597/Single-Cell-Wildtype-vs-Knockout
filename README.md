# Single-Cell
Single_cell

**Step 1: Read and Prepare Data**
Read the scRNA-seq data files for both wild-type (WT) and knockout (KO) samples. This includes barcode, feature, and matrix files. The files are then renamed to match the format expected by the Read10X function.

**Step 2: Create Seurat Objects**
Create Seurat objects for the KO and WT data. Seurat objects are used to store both the expression data and associated metadata.
CreateSeuratObject: Initializes a Seurat object. Parameters:
counts: Raw gene expression data.
project: Name for the project.
min.cells: Minimum number of cells expressing a gene for it to be included.
min.features: Minimum number of genes expressed in a cell for it to be included.

**Step 3: Quality Control**
Calculate quality control metrics, visualize them using violin plots, and filter out low-quality cells.
PercentageFeatureSet: Computes the percentage of mitochondrial gene expression, which is a common QC metric.
VlnPlot: Creates violin plots to visualize the distribution of QC metrics.
FeatureScatter: Creates scatter plots to visualize relationships between QC metrics.

<img width="367" alt="Cd36_wt QC metrics" src="https://github.com/Divya090597/Single-Cell/assets/156469276/34f9e3de-f9a2-4522-a3e9-6bbdd4c151e5">

<img width="367" alt="Cd36_ko QC metrics" src="https://github.com/Divya090597/Single-Cell/assets/156469276/65874a81-afff-457f-a67a-f48d3d525659">

<img width="367" alt="FeatureScatter_wt" src="https://github.com/Divya090597/Single-Cell/assets/156469276/20abba79-29b9-493f-8e55-3014445827f5">
<img width="367" alt="FeatureScatter_ko" src="https://github.com/Divya090597/Single-Cell/assets/156469276/2673b8ea-77a8-4c8b-842c-2b295ce318e9">



**Step 4: Data Normalization**
Normalize the data to account for differences in sequencing depth across cells.
NormalizeData: Normalizes the expression data.
normalization.method = "LogNormalize": Normalizes gene expression values by log-transformation.
scale.factor = 10000: Scaling factor for normalization.

**Step 5: Identify Highly Variable Features**
Identify the most variable genes, which are likely to be biologically significant.
FindVariableFeatures: Identifies highly variable genes.
selection.method = "vst": Variance Stabilizing Transformation method.
nfeatures = 2000: Number of variable genes to identify.

<img width="367" alt="Variable features_wt" src="https://github.com/Divya090597/Single-Cell/assets/156469276/973b6f58-2d54-4dca-9739-7c78fa088e67">
<img width="367" alt="Variable features_ko" src="https://github.com/Divya090597/Single-Cell/assets/156469276/314de807-b9ff-411b-86a1-7b8ee58cfedc">


**Step 6: Scaling the Data**
Scale the data to ensure that each gene has a mean of zero and variance of one.
ScaleData: Centers and scales the data for each gene.
features = all.genes: Scales all genes in the dataset.

**Step 7: PCA for Dimensionality Reduction**
Perform Principal Component Analysis (PCA) to reduce dimensionality and identify major sources of variation.

<img width="367" alt="PCA_wt" src="https://github.com/Divya090597/Single-Cell/assets/156469276/b8ab520a-4746-4bba-a42e-76b3f9201a9b">
<img width="367" alt="PCA_wt Dimplot" src="https://github.com/Divya090597/Single-Cell/assets/156469276/7921e74d-8af3-4c27-a39d-72c8235e373f">
<img width="367" alt="PCA_ko" src="https://github.com/Divya090597/Single-Cell/assets/156469276/d4486174-d032-41d9-aaad-1ac408d83b3c">
<img width="367" alt="PCA_ko Dimplot" src="https://github.com/Divya090597/Single-Cell/assets/156469276/5a044b4d-25ca-4710-b24e-e75780cc4bb7">

**Step 8: Determine the Dimensionality**
Use ElbowPlot to determine the appropriate number of principal components (PCs) to retain for downstream analyses.
ElbowPlot: Helps in choosing the number of significant PCs by plotting the variance explained by each PC.

<img width="367" alt="Elbowplot_wt" src="https://github.com/Divya090597/Single-Cell/assets/156469276/29bc490b-31bd-46ec-9b58-82f7bb9ed48c">
<img width="367" alt="Elbowplot_ko" src="https://github.com/Divya090597/Single-Cell/assets/156469276/297321f8-36b1-4b75-9190-2c0f4da89974">


**Step 9: Clustering the Cells**
Cluster the cells based on PCA scores and visualize the clusters using UMAP.
FindNeighbors: Computes a k-nearest neighbors graph.
dims = 1:10: Uses the first 10 PCs.
FindClusters: Identifies clusters of cells.
resolution = 0.5: Determines the granularity of the clustering.
RunUMAP: Conducts Uniform Manifold Approximation and Projection (UMAP) for visualization.
DimPlot: Visualizes the clusters in a 2D space.
![CD36KO_UMAP_LBL](https://github.com/Divya090597/Single-Cell/assets/156469276/c8a039dd-12fb-4a07-b17e-51bb74dee2c8)
![CD36WT_UMAP_LBL](https://github.com/Divya090597/Single-Cell/assets/156469276/1b32c5ae-38b1-4a54-9b9e-57121d340c57)








