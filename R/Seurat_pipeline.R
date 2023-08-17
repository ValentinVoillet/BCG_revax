library(Seurat)
library(here)
library(tidyverse)
library(data.table)



################################
################################
###--- Set up Seurat object
#- GEX
mat.GEX <- Matrix::readMM(here("NCBI", "GEX", "sparsematrix.mtx"))
mat.GEX <- as.data.frame(mat.GEX)
colnames(mat.GEX) <- as.data.frame(data.table::fread(here("NCBI", "GEX", "colnames.tsv"), header = FALSE))[, 1] %>% as.character()
rownames(mat.GEX) <- as.data.frame(data.table::fread(here("NCBI", "GEX", "rownames.tsv"), header = FALSE))[, 1] %>% as.character()

#- ADT
mat.ADT <- Matrix::readMM(here("NCBI", "ADT", "sparsematrix.mtx"))
mat.ADT <- as.data.frame(mat.ADT)
colnames(mat.ADT) <- as.data.frame(data.table::fread(here("NCBI", "ADT", "colnames.tsv"), header = FALSE))[, 1] %>% as.character()
rownames(mat.ADT) <- as.data.frame(data.table::fread(here("NCBI", "ADT", "rownames.tsv"), header = FALSE))[, 1] %>% as.character()

#- meta.data
meta.data <- read.table(file = here("NCBI", "meta-data.csv"), sep = ",")
rownames(meta.data) <- meta.data$rowname
# identical(rownames(meta.data), colnames(mat.GEX)) # TRUE
# identical(rownames(meta.data), colnames(mat.ADT)) # TRUE

#- Create the Seurat object.
merged.seurat <- CreateSeuratObject(counts = as(as.matrix(mat.GEX), "sparseMatrix"))
merged.seurat[["ADT"]] <- CreateAssayObject(counts = as(as.matrix(mat.ADT), "sparseMatrix"))
# identical(meta.data$rowname, rownames(merged.seurat@meta.data)) # TRUE
# identical(rownames(meta.data), rownames(merged.seurat@meta.data)) # TRUE
merged.seurat@meta.data <- meta.data
# saveRDS(object = merged.seurat, file = here("data", "001_merged-seurat.rds"))



################################
################################
###--- Data Integration (GEX) - AgSpec & AgSpecNot - T cells
#- Data (subset)
# merged.seurat <- readRDS(file = here("data", "001_merged-seurat.rds"))
merged.seurat <- subset(merged.seurat, sort_population != "CD3+")
merged.seurat@meta.data$batch <- sapply(X = merged.seurat@meta.data$orig.ident, FUN = function(x) str_split(string = x, pattern = "_")[[1]][1])
merged.seurat@meta.data$sample <- paste(merged.seurat@meta.data$batch, merged.seurat@meta.data$pubid, merged.seurat@meta.data$arm, merged.seurat@meta.data$stim, merged.seurat@meta.data$visit, sep = "_") # using as combination factor of batch, pubid, arm, stimulation & visit
DefaultAssay(object = merged.seurat) <- "RNA"
seurat.list <- SplitObject(merged.seurat, split.by = "sample")

#- Normalization
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <- SCTransform(seurat.list[[i]], verbose = TRUE)
}

#- Data Integration (with features, reference & rpca)
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 2000, verbose = TRUE)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features, verbose = TRUE)
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <- RunPCA(seurat.list[[i]], features = features, verbose = TRUE)
}
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", dims = 1:50, anchor.features = features, reduction = "rpca", reference = 1:2, verbose = TRUE) # One sample (pubid & Stimulation) is used as reference
integrated.seurat <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", dims = 1:50, verbose = TRUE)
DefaultAssay(integrated.seurat) <- "integrated"

#- Normalization
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "RNA", normalization.method = "LogNormalize", verbose = TRUE)
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "ADT", normalization.method = "CLR", verbose = TRUE)

#- PCA
# TCR-genes are removed
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDV"))]
integrated.seurat <- RunPCA(object = integrated.seurat, npcs = 100, verbose = TRUE)
chosen.nPCs <- 50

#- UMAP
integrated.seurat <- RunUMAP(object = integrated.seurat, dims = 1:chosen.nPCs, verbose = TRUE, seed.use = 1234, min.dist = .01)

#- Output
# saveRDS(object = integrated.seurat, file = here("data", "001_integrated-seurat_GEX.rds"))



################################
################################
###--- Data Integration (ADT) - AgSpec & AgSpecNot - T cells
#- Data (subset)
# merged.seurat <- readRDS(file = here("data", "001_merged-seurat_TEST.rds"))
merged.seurat <- subset(merged.seurat, sort_population != "CD3+")
merged.seurat@meta.data$batch <- sapply(X = merged.seurat@meta.data$orig.ident, FUN = function(x) str_split(string = x, pattern = "_")[[1]][1])
merged.seurat@meta.data$sample <- paste(merged.seurat@meta.data$batch, merged.seurat@meta.data$pubid, merged.seurat@meta.data$arm, merged.seurat@meta.data$stim, merged.seurat@meta.data$visit, sep = "_") # using as combination factor of batch, pubid, arm, stimulation & visit
DefaultAssay(object = merged.seurat) <- "ADT"
seurat.list <- SplitObject(merged.seurat, split.by = "sample")

#- Normalization
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <- NormalizeData(seurat.list[[i]], assay = "ADT", normalization.method = "CLR", verbose = TRUE)
  seurat.list[[i]] <- FindVariableFeatures(object = seurat.list[[i]], assay = "ADT")
}

#- Data Integration (with reference & rpca)
features <- SelectIntegrationFeatures(object.list = seurat.list)
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <- ScaleData(object = seurat.list[[i]], assay = "ADT", features = features, vars.to.regress = "nCount_ADT", verbose = TRUE)
  seurat.list[[i]] <- RunPCA(object = seurat.list[[i]], assay = "ADT", features = features, verbose = TRUE)
}
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:50, reduction = "rpca", reference = 1:2, verbose = TRUE) # One sample (pubid & Stimulation) is used as reference
integrated.seurat <- IntegrateData(anchorset = seurat.anchors, dims = 1:50, verbose = TRUE)
DefaultAssay(integrated.seurat) <- "integrated"

#- Normalization
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "RNA", normalization.method = "LogNormalize", verbose = TRUE)
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "ADT", normalization.method = "CLR", verbose = TRUE)

#- Scaling & PCA
integrated.seurat <- ScaleData(object = integrated.seurat, vars.to.regress = "nCount_ADT",verbose = TRUE)
integrated.seurat <- RunPCA(object = integrated.seurat, npcs = 100, verbose = TRUE)
chosen.nPCs <- 50

#- UMAP
integrated.seurat <- RunUMAP(object = integrated.seurat, dims = 1:chosen.nPCs, verbose = TRUE, seed.use = 1234, min.dist = .01)

#- Output
# saveRDS(object = integrated.seurat, file = here("data", "001_integrated-seurat_ADT.rds"))



################################
################################
###--- Data Integration - WNN - AgSpec & AgSpecNot - T cells
# integrated.seurat_GEX <- readRDS(file = here("data", "001_integrated-seurat_GEX.rds"))
# integrated.seurat_ADT <- readRDS(file = here("data", "001_integrated-seurat_ADT.rds"))

#- Add PCA results to one Seurat object
integrated.seurat_GEX[["pca_ADT"]] <- CreateDimReducObject(embeddings = Embeddings(object = integrated.seurat_ADT, reduction = "pca"), assay = "integrated")

#- Run WNN
integrated.seurat_WNN <- FindMultiModalNeighbors(object = integrated.seurat_GEX,
                                                 reduction.list = list("pca", "pca_ADT"),
                                                 dims.list = list(1:50, 1:50),
                                                 modality.weight.name = c("RNA.weight", "ADT.weight"),
                                                 prune.SNN = 1/20)

#- UMAP
integrated.seurat_WNN <- RunUMAP(object = integrated.seurat_WNN, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = 1234, min.dist = .01, verbose = TRUE)

#- Clustering
set.seed(1234)
resolutions <- c(seq(from = .1, to = 1, by = .1), seq(from = 1.5, to = 2, by = .5))
for(i in resolutions) {
  integrated.seurat_WNN <- FindClusters(object = integrated.seurat_WNN, graph.name = "wsnn", algorithm = 3, resolution = i, n.start = 10, verbose = TRUE)
}

#- Output
# saveRDS(object = integrated.seurat_WNN, file = here("data", "001_integrated-seurat_WNN.rds"))



################################
################################
###---  Data Integration (GEX) - AgSpec & AgSpecNot - CD4+ T
#- Data (subset)
# integrated.seurat <- readRDS(file = here("data", "001_integrated-seurat_GEX.rds"))
integrated.seurat <- subset(integrated.seurat, sort_population == "AgSpec")
integrated.seurat <- subset(integrated.seurat, T_phenotype == "CD4+")
DefaultAssay(integrated.seurat)

#- Data Integration
# Done - including all T cells (see code above)

#- Normalization
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "RNA", normalization.method = "LogNormalize", verbose = TRUE)
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "ADT", normalization.method = "CLR", verbose = TRUE)

#- PCA
# TCR-genes are removed
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDV"))]
integrated.seurat <- RunPCA(object = integrated.seurat, npcs = 100, verbose = TRUE)
chosen.nPCs <- 50

#- UMAP
integrated.seurat <- RunUMAP(object = integrated.seurat, dims = 1:chosen.nPCs, verbose = TRUE, seed.use = 1234, min.dist = .01)

#- Output
# saveRDS(object = integrated.seurat, file = here("data", "002_integrated-seurat_GEX.rds"))



################################
################################
###---  Data Integration (ADT) - AgSpec & AgSpecNot - CD4+ T
#- Data (subset)
# integrated.seurat <- readRDS(file = here("data", "001_integrated-seurat_ADT.rds"))
integrated.seurat <- subset(integrated.seurat, sort_population == "AgSpec")
integrated.seurat <- subset(integrated.seurat, T_phenotype == "CD4+")
DefaultAssay(integrated.seurat)

#- Data Integration
# Done - including all T cells (see code above)

#- Normalization
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "RNA", normalization.method = "LogNormalize", verbose = TRUE)
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "ADT", normalization.method = "CLR", verbose = TRUE)

#- Scaling & PCA
integrated.seurat <- ScaleData(object = integrated.seurat, vars.to.regress = "nCount_ADT",verbose = TRUE)
integrated.seurat <- RunPCA(object = integrated.seurat, npcs = 100, verbose = TRUE)
chosen.nPCs <- 25

#- UMAP
integrated.seurat <- RunUMAP(object = integrated.seurat, dims = 1:chosen.nPCs, verbose = TRUE, seed.use = 1234, min.dist = .01)

#- Output
# saveRDS(object = integrated.seurat, file = here("data", "002_integrated-seurat_ADT.rds"))



################################
################################
###--- Data Integration - WNN - AgSpec & AgSpecNot - CD4+ T
#- Data
# integrated.seurat_GEX <- readRDS(file = here("data", "002_integrated-seurat_GEX.rds"))
# integrated.seurat_ADT <- readRDS(file = here("data", "002_integrated-seurat_ADT.rds"))

#- Add PCA results to one Seurat object
integrated.seurat_GEX[["pca_ADT"]] <- CreateDimReducObject(embeddings = Embeddings(object = integrated.seurat_ADT, reduction = "pca"), assay = "integrated")

#- Run WNN
integrated.seurat_WNN <- FindMultiModalNeighbors(object = integrated.seurat_GEX,
                                                 reduction.list = list("pca", "pca_ADT"),
                                                 dims.list = list(1:50, 1:25),
                                                 modality.weight.name = c("RNA.weight", "ADT.weight"),
                                                 prune.SNN = 1/20)

#- UMAP
integrated.seurat_WNN <- RunUMAP(object = integrated.seurat_WNN, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = 1234, min.dist = .1, n.neighbors = 15, verbose = TRUE)

#- Clustering
set.seed(1234)
resolutions <- c(seq(from = .1, to = 1, by = .1), seq(from = 1.2, to = 2, by = .2))
for(i in resolutions) {
  integrated.seurat_WNN <- FindClusters(object = integrated.seurat_WNN, graph.name = "wsnn", algorithm = 3, resolution = i, n.start = 10, verbose = TRUE)
}

#- Output
saveRDS(object = integrated.seurat_WNN, file = here("data", "002_integrated-seurat_WNN.rds"))



################################
################################
###---  Data Integration (GEX) - AgSpec & AgSpecNot - CD8+ T
#- Data (subset)
# integrated.seurat <- readRDS(file = here("data", "001_integrated-seurat_GEX.rds"))
integrated.seurat <- subset(integrated.seurat, sort_population == "AgSpec")
integrated.seurat <- subset(integrated.seurat, T_phenotype == "CD8+")
DefaultAssay(integrated.seurat)

#- Data Integration
# Done - including all T cells (see code above)

#- Normalization
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "RNA", normalization.method = "LogNormalize", verbose = TRUE)
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "ADT", normalization.method = "CLR", verbose = TRUE)

#- PCA
# TCR-genes are removed
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDJ"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDD"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRAV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRBV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRGV"))]
VariableFeatures(integrated.seurat) <- VariableFeatures(integrated.seurat)[which(!str_detect(string = VariableFeatures(integrated.seurat), pattern = "^TRDV"))]
integrated.seurat <- RunPCA(object = integrated.seurat, npcs = 100, verbose = TRUE)
chosen.nPCs <- 50

#- UMAP
integrated.seurat <- RunUMAP(object = integrated.seurat, dims = 1:chosen.nPCs, verbose = TRUE, seed.use = 1234, min.dist = .01)

#- Output
# saveRDS(object = integrated.seurat, file = here("data", "003_integrated-seurat_GEX.rds"))



################################
################################
###---  Data Integration (ADT) - AgSpec & AgSpecNot - CD8+ T
#- Data (subset)
# integrated.seurat <- readRDS(file = here("data", "001_integrated-seurat_ADT.rds"))
integrated.seurat <- subset(integrated.seurat, sort_population == "AgSpec")
integrated.seurat <- subset(integrated.seurat, T_phenotype == "CD8+")
DefaultAssay(integrated.seurat)

#- Data Integration
# Done - including all T cells (see code above)

#- Normalization
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "RNA", normalization.method = "LogNormalize", verbose = TRUE)
integrated.seurat <- NormalizeData(object = integrated.seurat, assay = "ADT", normalization.method = "CLR", verbose = TRUE)

#- Scaling & PCA
integrated.seurat <- ScaleData(object = integrated.seurat, vars.to.regress = "nCount_ADT",verbose = TRUE)
integrated.seurat <- RunPCA(object = integrated.seurat, npcs = 100, verbose = TRUE)
chosen.nPCs <- 25

#- UMAP
integrated.seurat <- RunUMAP(object = integrated.seurat, dims = 1:chosen.nPCs, verbose = TRUE, seed.use = 1234, min.dist = .01)

#- Output
# saveRDS(object = integrated.seurat, file = here("data", "003_integrated-seurat_ADT.rds"))



################################
################################
###--- Data Integration - WNN - AgSpec & AgSpecNot - CD8+ T
#- Data
# integrated.seurat_GEX <- readRDS(file = here("data", "003_integrated-seurat_GEX.rds"))
# integrated.seurat_ADT <- readRDS(file = here("data", "003_integrated-seurat_ADT.rds"))

#- Add PCA results to one Seurat object
integrated.seurat_GEX[["pca_ADT"]] <- CreateDimReducObject(embeddings = Embeddings(object = integrated.seurat_ADT, reduction = "pca"), assay = "integrated")

#- Run WNN
integrated.seurat_WNN <- FindMultiModalNeighbors(object = integrated.seurat_GEX,
                                                 reduction.list = list("pca", "pca_ADT"),
                                                 dims.list = list(1:50, 1:25),
                                                 modality.weight.name = c("RNA.weight", "ADT.weight"),
                                                 prune.SNN = 1/20)

#- UMAP
integrated.seurat_WNN <- RunUMAP(object = integrated.seurat_WNN, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = 1234, min.dist = .1, n.neighbors = 15, verbose = TRUE)

#- Clustering
set.seed(1234)
resolutions <- c(seq(from = .1, to = 1, by = .1), seq(from = 1.2, to = 2, by = .2))
for(i in resolutions) {
  integrated.seurat_WNN <- FindClusters(object = integrated.seurat_WNN, graph.name = "wsnn", algorithm = 3, resolution = i, n.start = 10, verbose = TRUE)
}

#- Output
saveRDS(object = integrated.seurat_WNN, file = here("data", "003_integrated-seurat_WNN.rds"))

