# Single-Cell

#--------------------------**Seurat Clustering**---------------------------------

**Step 1: Read and Prepare Data**
   
      ~ Read the scRNA-seq data files for both wild-type (WT) and knockout (KO) samples. 
   
      ~ This includes barcode, feature, and matrix files. 
   
      ~ The files are then renamed to match the format expected by the Read10X function.

**Step 2: Create Seurat Objects**

     Create Seurat objects for the KO and WT data. Seurat objects are used to store both the expression data and associated 
     metadata.
     CreateSeuratObject: Initializes a Seurat object. Parameters:
     counts: Raw gene expression data.
     project: Name for the project.
     min.cells: Minimum number of cells expressing a gene for it to be included.
     min.features: Minimum number of genes expressed in a cell for it to be included.

**Step 3: Quality Control**

     Calculate quality control metrics, visualize them using violin plots, and filter out low-quality cells.
     Filter out low quality cells and genes and normalize the data
     Calculate mitochondrial QC metrics
      QC metrics: "nfeature_RNA", "nCount_RNA", "percent.mt"
     Low quality cells or empty droplets often have very few genes
     Cell doublets or multiplets have high values of nfeature_RNA & nCount_RNA
     Low quality cells often have high percentage of mitochondrial genes

     The [[ operator can add columns to object metadata. This is a great place to stash QC stats
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

**Step 10: Finding differentially expressed features (cluster biomarkers)**

**Finding Cluster-Specific Markers**
1. Finding Markers for Cluster 2:

          we use the 'FindMarkers' function to find genes that are differentially expressed in cluster 2 compared to all other clusters.
          The 'head' function displays the top 5 markers.
2. Finding Markers that Distinguish Cluster 5 from Clusters 0 and 3:

          we identify markers that distinguish cluster 5 from clusters 0 and 3.
          The 'ident.1' parameter specifies the cluster of interest (5),
          while 'ident.2' lists the clusters for comparison (0 and 3).
3. Finding markers for every cluster compared to all remaining cells:

          The ,FindAllMarkers' function finds markers for all clusters.
          The 'only.pos' parameter ensures that only 'upregulated (positive)' markers are returned.
4. Filtering and Selecting Top Markers:

          Code groups markers by cluster and filters them to retain those with an average log2 fold-change greater than 1, indicating strong differential 
          expression.
6. Selecting the Top 10 Markers Per Cluster
   
<img width="367" alt="Vlnplot_features_ko" src="https://github.com/Divya090597/Single-Cell/assets/156469276/2e486d51-00de-4c68-b26b-de5077b00f96">
<img width="367" alt="Vln_features_wt" src="https://github.com/Divya090597/Single-Cell/assets/156469276/5f030bd7-8454-4444-bbf8-795dce5548d8">

# Step 10: Merging the Seurat Odjects
**STANDARD PREPROCESSING WORKFLOW**
<img width="462" alt="MD_Vln variablefeatureplot" src="https://github.com/Divya090597/Single-Cell/assets/156469276/5433d0bc-3c3d-4e6d-b6c8-b75504b42d5d">

![MD_PCA2](https://github.com/Divya090597/Single-Cell/assets/156469276/aefc99f0-8d08-433f-ac9b-ea7358e1aa7c)
![MD_PCA1](https://github.com/Divya090597/Single-Cell/assets/156469276/3e64be7c-3f0e-46b3-a252-315d07892622)
<img width="462" alt="Elbow_MD" src="https://github.com/Divya090597/Single-Cell/assets/156469276/c22303b2-aabd-4fcb-b994-3790914bd29d">

![,M_UMAP](https://github.com/Divya090597/Single-Cell/assets/156469276/d98bd49b-5505-46db-b725-27010a8e3264)

<img width="462" alt="tsne_MD" src="https://github.com/Divya090597/Single-Cell/assets/156469276/28621fd5-d096-4caf-9ac1-1159c2aafc8c">

# Step 11 : Integration

    Integrate the KO and WT datasets for combined analysis.
    SelectIntegrationFeatures: Identifies features to be used for integration.
    FindIntegrationAnchors: Finds anchor points between datasets.
    IntegrateData: Integrates multiple Seurat objects into a single object.

**Step 12:** 

    1.DefaultAssay: Sets the default assay.
    2.ScaleData: Scales the integrated data.
    3.RunPCA: Performs PCA on the integrated data.
    4.FindNeighbors: Computes a k-nearest neighbors graph on the integrated data.
    5.FindClusters: Identifies clusters on the integrated data.
    RunUMAP: Conducts UMAP on the integrated data.
    DimPlot: Visualizes clusters in UMAP space.
<img width="367" alt="Elbow_ID" src="https://github.com/Divya090597/Single-Cell/assets/156469276/74982062-16d0-46e7-821c-602f488b27f6">

![Int_UMAP](https://github.com/Divya090597/Single-Cell/assets/156469276/b76b0d6c-5e6f-4a7c-bae4-f72e7694503d)

**Step 12: Differential Expression on Integrated Data**

    Perform differential expression analysis to identify genes differentially expressed between conditions.
    FindMarkers: Identifies differentially expressed genes between specified conditions or clusters.
![Feature plot integrated data](https://github.com/Divya090597/Single-Cell/assets/156469276/9a94a9a3-0b77-4b0d-838e-430d37642b74)

**-----------------------------**Single R:**---------------------------------------**
 **Annotation Diagnostics**
 
These steps involve evaluating the quality of the SingleR annotations. SingleR assigns cell types to individual cells based on reference datasets.
1. Based on Scores within Cells
    Access the scores assigned by SingleR for each cell. Each cell has scores indicating the similarity to each reference cell type. Higher scores indicate a higher similarity.

![Score_Dgns](https://github.com/Divya090597/Single-Cell/assets/156469276/d2bc5b86-27f2-49f3-87bd-3cb04244f7b7)

2. Based on Deltas Across Cells
   SingleR uses deltas, which represent the difference between the highest score and the second-highest score for each cell. This metric helps in assessing the confidence of the cell type assignment.

![Deltas_Dgns](https://github.com/Divya090597/Single-Cell/assets/156469276/b4bf74f0-8051-4a3a-b1ed-0ae9a39c0e46)

# Pseudo-Bulk

 Create pseudo-bulk samples by aggregating expression profiles for downstream bulk RNA-seq style analyses.
 
 Save the counts data and metadata to files.
 
 These files can be used for further pseudo-bulk analyses outside of Seurat, such as differential expression analysis using 
 tools designed for bulk RNA-seq data.













