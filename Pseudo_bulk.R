
################### PSEUDO BULK ###########################################

### Aggregate Expression Data

agg = AggregateExpression(mergeddata, return.seurat = T, group.by = c('orig.ident', 'seurat_clusters','SingleR_labels'), normalization.method = "LogNormalize",scale.factor = 10000)


### View and Extract Data

View(agg@assays$RNA$data)
cnames = colnames(agg@assays$RNA$data)

###Parse Column Names and Create Metadata

m1 = do.call(rbind, strsplit(cnames, "_"))
print(head(m1))
meta_text = cbind(cnames, m1)
colnames(meta_text) = c("ID", "Sample", "Cluster", "CellType")
head(meta_text)

### Create Sample_CellType Column

Sample_celltype = paste(meta_text[, "Cluster"], meta_text[,"CellType"], sep="_")
meta_text = cbind(meta_text, Sample_celltype);
colnames(meta_text)
head(meta_text)

### Define Output Folder and File Paths
output_folder <- "~/DataAnalysis/R"
Project_Name <- "Bulk_data"
matrix_file = paste(output_folder, Project_Name, "Bulk_data.txt", sep="");

### Write Data to Files
write.table(2^ agg@assays$RNA$data, file = matrix_file, sep="\t", row.names = TRUE,col.names=NA)
head(agg@assays$RNA$data)
meta_file = paste(output_folder, Project_Name, "_meta.txt", sep="");
write.table(meta_text, file = meta_file, sep="\t", row.names=F)

### Write out the matrix ###
outputmatrix = matrix_file
outputmeta = meta_file

### Generation of the EASY app ###

output_dir <- file.path(output_folder,paste(Project_Name,"Mouse_EASY_App", sep = ""))

outputmatrix = paste(output_dir,"/", Project_Name,"Bulk_data.txt",sep = "")

outputmeta = paste(output_dir,"/", Project_Name,"_meta.txt", sep = "")


