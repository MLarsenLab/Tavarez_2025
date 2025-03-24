#packages were loaded using the pacman package, you will need to download the packages prior to using pacman, they are below
pacman::p_load(SingleCellExperiment, celda, Seurat, ggplot2, cowplot, scater, patchwork, 
               dplyr, sctransform, glmGamPoi, stringr, scCustomize, scDblFinder, singleCellTK, biomaRt,
               clustree, hdf5r, rhdf5, devtools, sctransform)

newFunc <- function(aSeurat, clusterResolution = 0.5, numFeatures = 500){
  if(0){
    aSeurat <- agg.merge_3day
  }
  
  aSeurat <- NormalizeData(aSeurat)
  aSeurat = FindVariableFeatures(aSeurat, selection.method = "vst", nfeatures = numFeatures)
  aSeurat = RunPCA(aSeurat, npcs = 35)
  aSeurat = FindNeighbors(aSeurat, dims = 1:35)
  aSeurat = FindClusters(aSeurat, resolution = clusterResolution)
  
  return(aSeurat)
}
Read_CellBender_h5_Mat <- function(file_name, use.names = TRUE, unique.features = TRUE) {
  # Check hdf5r installed
  if(0){
    file_name <- "Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_020_21_PDGFRaTdTom_ControlHarvest/outputDirectory_T41"
  }
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    cli_abort(message = c("Please install hdf5r to read HDF5 files",
                          "i" = "`install.packages('hdf5r')`")
    )
  }
  # Check file
  if (!file.exists(file_name)) {
    stop("File not found")
  }
  
  if (use.names) {
    feature_slot <- 'features/name'
  } else {
    feature_slot <- 'features/id'
  }
  
  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")
  
  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]
  
  
  sparse.mat <- Matrix::sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )
  
  if (unique.features) {
    features <- make.unique(names = features)
  }
  
  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
  
  infile$close_all()
  
  return(sparse.mat)
}
removeDoublets <- function(seuratObject, pctDbl = 10){
  #uses scDblFinder and DoubletFinder to find droplets that may contain mutiple cells
  #Returns a seruat object where cells that both scDblFinder and DoubletFinder agree on as doublets have been removed
  
  #scDblFinder needs a single cell experement as input and outputs a score and class (doublet or singlet)
  #Assume 20% doublets, this is an over estimation of number of doublets as we are looking for agreement 
  aSCE <- as.SingleCellExperiment(seuratObject)
  aSCE <- scDblFinder(aSCE, dbr = (pctDbl/100))
  
  #copy the cell asignemtns to new metadata in the seurat object
  seuratObject$scDblFinder.score <- aSCE$scDblFinder.score
  seuratObject$scDblFinder.class <- aSCE$scDblFinder.class
  
  #DoubletFinder requires data to be normalized before it can identify doublets
  #The object passed in may not be ready to be normalized or may need to be normalized in a differt way so
  #a copy is made and the copy is normalized for DoubletFinder and the cell assignments copied to the original
  dblSeurat <- seuratObject
  dblSeurat <- NormalizeData(dblSeurat)
  dblSeurat <- FindVariableFeatures(dblSeurat, selection.method = "vst", nfeatures = 2000)
  geneList <- rownames(dblSeurat)
  dblSeurat <- ScaleData(dblSeurat, feature = geneList)
  dblSeurat <- RunPCA(dblSeurat, features = VariableFeatures(object = dblSeurat))
  dblSeurat <- RunUMAP(dblSeurat, dims = 1:40)
  
  #DoubletFinder needs to know how many doublets to look for, using 20% of total number of cells
  dblCount <- as.integer(length(rownames(dblSeurat@meta.data))*(pctDbl/100))
  dblSeurat <- doubletFinder(dblSeurat, PCs = 1:40, pN = 0.25, pK = 0.1, nExp = dblCount, reuse.pANN = FALSE, sct = FALSE)
  #DoubletFinder returns the classification in a new metadata with the name based on the parameters used
  DF <- sprintf("DF.classifications_0.25_0.1_%i", dblCount)
  DFIndex <- which(colnames(dblSeurat@meta.data) == DF)
  
  #Copy cell assigment to the seurat object metadata
  seuratObject$DF.classifications <- dblSeurat@meta.data$DF
  
  #Find cells that are clasified as a doublet by both scDblFinder and doubletFinder
  doublets1 <- subset(seuratObject, scDblFinder.class == "doublet")
  trueDoublets <- subset(doublets1, DF.classifications == "Doublet")
  seuratObject$True_Doublet <- ifelse(colnames(seuratObject)%in% colnames(trueDoublets), "Doublet", "Singlet")
  
  #Print the percentage of doublets found
  print(sprintf("%i doublets found in %i total cells (%.2f%%)", ncol(trueDoublets), ncol(seuratObject), 
                (ncol(trueDoublets)/ncol(seuratObject)*100)))
  
  #Make a subset of only singlets and return the cleaned object
  seuratObject <- subset(seuratObject, True_Doublet == "Singlet")
  return(seuratObject)
}


#####pre-processing and merging, need to use the file locations where you download the files on GEO, for GEO number check Tavarez_2025
#agg.reps 1-3 are homeostatic, the other reps are labeled by the timepoints. The order which you load the data should not matter, just make
#sure that you put the right timepoints into the right names, you can change the names if you want. The files ending in h5 have undergone 
#cell bender cleanup but you can also use the raw files using the code found in the seurat vignette by the Satija Lab.
afileName <- ""
aMatrix <- Read_CellBender_h5_Mat(afileName)

agg.rep1 <- CreateSeuratObject(aMatrix, project = "Homeostatic", min.cells = 3, min.features = 200)

bfileName <- ""
bMatrix <- Read_CellBender_h5_Mat(bfileName)

agg.rep2 <- CreateSeuratObject(bMatrix, project = "Homeostatic", min.cells = 3, min.features = 200)

cfileName <- ""
cMatrix <- Read_CellBender_h5_Mat(cfileName)

agg.rep3 <- CreateSeuratObject(cMatrix, project = "Homeostatic", min.cells = 3, min.features = 200)

dfileName <- ""
dMatrix <- Read_CellBender_h5_Mat(dfileName)

agg.3day.rep1 <- CreateSeuratObject(dMatrix, project = "3day", min.cells = 3, min.features = 200)

efileName <- ""
eMatrix <- Read_CellBender_h5_Mat(efileName)

agg.3day.rep2 <- CreateSeuratObject(eMatrix, project = "3day", min.cells = 3, min.features = 200)

ffilename =""
fmatrix = Read_CellBender_h5_Mat(ffilename)

agg.7day.rep1 <- CreateSeuratObject(fmatrix, project = "3day", min.cells = 3, min.features = 200)

gfilename =""
gmatrix = Read_CellBender_h5_Mat(gfilename)

agg.7day.rep2 <- CreateSeuratObject(gmatrix, project = "7day", min.cells = 3, min.features = 200)

hfileName <- ""
hMatrix <- Read_CellBender_h5_Mat(hfileName)

agg.3day.rep3 <- CreateSeuratObject(hMatrix, project = "7day", min.cells = 3, min.features = 200)

rm(afileName,bfileName,cfileName,dfileName,ffilename,gfilename,hfileName,aMatrix,bMatrix,cMatrix,dMatrix,eMatrix,fmatrix,gmatrix,hMatrix)

afileName <- ""
aMatrix <- Read_CellBender_h5_Mat(afileName)

agg.1day.rep1 <- CreateSeuratObject(aMatrix, project = "1day", min.cells = 3, min.features = 200)

afileName <- ""
aMatrix <- Read_CellBender_h5_Mat(afileName)

agg.1day.rep2 <- CreateSeuratObject(aMatrix, project = "1day", min.cells = 3, min.features = 200)

afileName <- ""
aMatrix <- Read_CellBender_h5_Mat(afileName)

agg.1day.rep3 <- CreateSeuratObject(aMatrix, project = "1day", min.cells = 3, min.features = 200)

afileName <- ""
aMatrix <- Read_CellBender_h5_Mat(afileName)

agg.7day.rep3 <- CreateSeuratObject(aMatrix, project = "7day", min.cells = 3, min.features = 200)

agg.rep1[["percent.mt"]] <- PercentageFeatureSet(agg.rep1, pattern = "^mt-")
agg.rep2[["percent.mt"]] <- PercentageFeatureSet(agg.rep2, pattern = "^mt-")
agg.rep3[["percent.mt"]] <- PercentageFeatureSet(agg.rep3, pattern = "^mt-")
agg.3day.rep1[["percent.mt"]] <- PercentageFeatureSet(agg.3day.rep1, pattern = "^mt-")
agg.3day.rep2[["percent.mt"]] <- PercentageFeatureSet(agg.3day.rep2, pattern = "^mt-")
agg.3day.rep3[["percent.mt"]] <- PercentageFeatureSet(agg.3day.rep3, pattern = "^mt-")
agg.7day.rep1[["percent.mt"]] <- PercentageFeatureSet(agg.7day.rep1, pattern = "^mt-")
agg.7day.rep2[["percent.mt"]] <- PercentageFeatureSet(agg.7day.rep2, pattern = "^mt-")
agg.1day.rep1[["percent.mt"]] <- PercentageFeatureSet(agg.1day.rep1, pattern = "^mt-")
agg.1day.rep2[["percent.mt"]] <- PercentageFeatureSet(agg.1day.rep2, pattern = "^mt-")
agg.1day.rep3[["percent.mt"]] <- PercentageFeatureSet(agg.1day.rep3, pattern = "^mt-")
agg.7day.rep3[["percent.mt"]] <- PercentageFeatureSet(agg.7day.rep3, pattern = "^mt-")

agg.1day.rep1 <- subset(agg.1day.rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.1day.rep2 <- subset(agg.1day.rep2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.1day.rep3 <- subset(agg.1day.rep3, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.7day.rep3 <- subset(agg.7day.rep3, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.rep1 <- subset(agg.rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.rep2 <- subset(agg.rep2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.rep3 <- subset(agg.rep3, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.3day.rep1 <- subset(agg.3day.rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.3day.rep2 <- subset(agg.3day.rep2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.3day.rep3 <- subset(agg.3day.rep3, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.7day.rep1 <- subset(agg.7day.rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
agg.7day.rep2 <- subset(agg.7day.rep2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
#adding identifiable metadata
agg.1day.rep1@meta.data[,"replicate"] = "1-day1"
agg.1day.rep2@meta.data[,"replicate"] = "1-day2"
agg.1day.rep3@meta.data[,"replicate"] = "1-day3"
agg.7day.rep3@meta.data[,"replicate"] = "7-day3"
agg.1day.rep1@meta.data[,"dataset"] = "1-day"
agg.1day.rep2@meta.data[,"dataset"] = "1-day"
agg.1day.rep3@meta.data[,"dataset"] = "1-day"
agg.7day.rep3@meta.data[,"dataset"] = "7-day"
agg.rep1@meta.data[,"replicate"] = "Homeostatic1"
agg.rep2@meta.data[,"replicate"] = "Homeostatic2"
agg.rep3@meta.data[,"replicate"] = "Homeostatic3"
agg.3day.rep1@meta.data[,"replicate"] = "3-day1"
agg.3day.rep2@meta.data[,"replicate"] = "3-day2"
agg.3day.rep3@meta.data[,"replicate"] = "3-day3"
agg.7day.rep1@meta.data[,"replicate"] = "7-day1"
agg.7day.rep2@meta.data[,"replicate"] = "7-day2"
agg.rep1@meta.data[,"dataset"] = "Homeostatic"
agg.rep2@meta.data[,"dataset"] = "Homeostatic"
agg.rep3@meta.data[,"dataset"] = "Homeostatic"
agg.3day.rep1@meta.data[,"dataset"] = "3-day"
agg.3day.rep2@meta.data[,"dataset"] = "3-day"
agg.3day.rep3@meta.data[,"dataset"] = "3-day"
agg.7day.rep1@meta.data[,"dataset"] = "7-day"
agg.7day.rep2@meta.data[,"dataset"] = "7-day"
#doublet removal
agg.1day.rep1 = removeDoublets(agg.1day.rep1)
agg.1day.rep2 = removeDoublets(agg.1day.rep2)
agg.1day.rep3 = removeDoublets(agg.1day.rep3)
agg.rep1 = removeDoublets(agg.rep1)
agg.rep2 = removeDoublets(agg.rep2)
agg.rep3 = removeDoublets(agg.rep3)
agg.3day.rep1 = removeDoublets(agg.3day.rep1)
agg.3day.rep2 = removeDoublets(agg.3day.rep2)
agg.3day.rep3 = removeDoublets(agg.3day.rep3)
agg.7day.rep1 = removeDoublets(agg.7day.rep1)
agg.7day.rep2 = removeDoublets(agg.7day.rep2)
agg.7day.rep3 = removeDoublets(agg.1day.rep3)

#The datasets to be integrated need to be in a list
#Make sure you add any meta data to each data set if you need to identify cells by data set (see the lines just above doublet removal)
seuratList <- c(agg.rep1, agg.rep2, agg.rep3, agg.1day.rep1, agg.1day.rep2, agg.1day.rep3, agg.3day.rep1, agg.3day.rep2, agg.3day.rep3, agg.7day.rep1, agg.7day.rep2, agg.7day.rep3)
#SCTransform has the same purpose as NormalizeData(), FindVariableFeatures() and ScaleData() 
rm(agg.1day.rep1,agg.1day.rep2,agg.1day.rep3,agg.3day.rep1,agg.3day.rep2,agg.3day.rep3,agg.7day.rep1,agg.7day.rep2,agg.7day.rep3,
   agg.rep1,agg.rep2,agg.rep3)
gc()
seuratList <- lapply(X = seuratList, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})



features <- SelectIntegrationFeatures(object.list = seuratList, nfeatures = 3000)
seuratList <- PrepSCTIntegration(object.list = seuratList, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features, normalization.method = "SCT")
rm(seuratList)
gc()
agg.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(anchors)
gc()

agg.integrated = RunPCA(agg.integrated, npcs = 35)
agg.integrated = FindNeighbors(agg.integrated, dims = 1:35)
agg.integrated = FindClusters(agg.integrated, resolution = .5)
agg.integrated <- RunUMAP(agg.integrated, reduction = "pca", dims = 1:35)
DefaultAssay(agg.integrated) = "SCT"
agg.integrated = PrepSCTFindMarkers(agg.integrated)

#sometimes the metadata is lost for some cells, this ensures that all timepoints are labeled by 'dataset'

agg.integrated$dataset[agg.integrated$replicate == 'Homeostatic1'] <- "Homeostatic"
agg.integrated$dataset[agg.integrated$replicate == 'Homeostatic2'] <- "Homeostatic"
agg.integrated$dataset[agg.integrated$replicate == 'Homeostatic3'] <- "Homeostatic"

agg.integrated$dataset[agg.integrated$replicate == '1-day1'] <- "1-day"
agg.integrated$dataset[agg.integrated$replicate == '1-day2'] <- "1-day"
agg.integrated$dataset[agg.integrated$replicate == '1-day3'] <- "1-day"

agg.integrated$dataset[agg.integrated$replicate == '3-day1'] <- "3-day"
agg.integrated$dataset[agg.integrated$replicate == '3-day2'] <- "3-day"
agg.integrated$dataset[agg.integrated$replicate == '3-day3'] <- "3-day"

agg.integrated$dataset[agg.integrated$replicate == '7-day1'] <- "7-day"
agg.integrated$dataset[agg.integrated$replicate == '7-day2'] <- "7-day"
agg.integrated$dataset[agg.integrated$replicate == '7-day3'] <- "7-day"

#this is sometimes required depending on seurat version
agg.integrated = JoinLayers(agg.integrated)

#Integrated Cell Identification Process
#Use markers from figure 1C to help guide cluster identification, depending on your seurat version, some cluster numbers may change but the 
#gene expression should be similar enough that you can adjust the clusters below 

agg.integrated$Celltype <- "Missing"
table(agg.integrated$Celltype)

agg.integrated$Celltype[agg.integrated$seurat_clusters == '1'] <- "Fibroblast"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '24'] <- "Fibroblast"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '9'] <- "Fibroblast"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '23'] <- "Fibroblast"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '0'] <- "Endothelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '6'] <- "Endothelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '8'] <- "Endothelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '11'] <- "Endothelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '18'] <- "Endothelial"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '2'] <- "Macrophage/Monocyte"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '7'] <- "Macrophage/Monocyte"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '14'] <- "Macrophage/Monocyte"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '19'] <- "Macrophage/Monocyte"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '25'] <- "Macrophage/Monocyte"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '3'] <- "Immune-T"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '13'] <- "Immune-T"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '4'] <- "Immune-NK"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '5'] <- "Epithelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '10'] <- "Epithelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '15'] <- "Epithelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '21'] <- "Epithelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '26'] <- "Epithelial"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '29'] <- "Epithelial"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '12'] <- "Pericyte"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '16'] <- "Dividing"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '17'] <- "Mixed"
agg.integrated$Celltype[agg.integrated$seurat_clusters == '28'] <- "Mixed"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '20'] <- "Glia"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '22'] <- "Immune"

agg.integrated$Celltype[agg.integrated$seurat_clusters == '27'] <- "Nerve"

#Figure 1E, this code allows you to count the number of cells per timepoint that falls within the celltype categories

table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Fibroblast" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Endothelial" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Macrophage/Monocyte" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Immune-T" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Immune-NK" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Epithelial" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Pericyte" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Dividing" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Mixed" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Glia" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Immune" )
table(agg.integrated$dataset == "Homeostatic" & agg.integrated$Celltype == "Nerve" )

table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Fibroblast" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Endothelial" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Macrophage/Monocyte" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Immune-T" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Immune-NK" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Epithelial" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Pericyte" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Dividing" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Mixed" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Glia" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Immune" )
table(agg.integrated$dataset == "1-day" & agg.integrated$Celltype == "Nerve" )

table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Fibroblast" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Endothelial" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Macrophage/Monocyte" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Immune-T" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Immune-NK" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Epithelial" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Pericyte" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Dividing" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Mixed" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Glia" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Immune" )
table(agg.integrated$dataset == "3-day" & agg.integrated$Celltype == "Nerve" )

table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Fibroblast" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Endothelial" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Macrophage/Monocyte" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Immune-T" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Immune-NK" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Epithelial" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Pericyte" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Dividing" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Mixed" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Glia" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Immune" )
table(agg.integrated$dataset == "7-day" & agg.integrated$Celltype == "Nerve" )

#creates a new metadata "CT_D" which combines celltype and dataset allowing for timepoint analysis for each cell type
agg.integrated$CT_D <- sprintf("%s_%s", agg.integrated$Celltype, agg.integrated$dataset)

Idents(agg.integrated) = 'Celltype'

#creates figure 1b

print(DimPlot(agg.integrated, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.35, repel = TRUE, raster = FALSE)) + scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) + aes(UMAP1,UMAP2) + guides(x = axis, y = axis) + theme(axis.line = element_line
                                                                                            (arrow = arrow(type = 'closed',
                                                                                                           length = unit(10,'pt'))),
                                                                                            axis.title = element_text(hjust = 0))
#creates figure 1d
agg.integrated$dataset <- factor(x = agg.integrated$dataset, levels = c('Homeostatic','1-day','3-day','7-day'))
Idents(agg.integrated) = 'dataset'

png('Umap_integrated_dataset.png', units = 'in', width = 18, height = 7, res = 600)
DimPlot(agg.integrated,reduction = "umap",group.by = 'Celltype',label = FALSE,label.size = 6,pt.size = .5,repel = TRUE, split.by = 'dataset',raster = FALSE) + 
  labs(title = '') +aes(UMAP1,UMAP2)
dev.off()

#to generate fibroblast subsets, select fibroblasts and reintegrate all datasets, if you need more space, remove
#the integrated datset after subsetting with the 'rm(agg.integrated)' function and then use 'gc()'

agg.tdtom.fibroblasts = subset(agg.integrated, subset = Celltype == "Fibroblast" | tdTomato > 1.0)
Idents(agg.tdtom.fibroblasts) = 'replicate'
levels(agg.tdtom.fibroblasts)
seuratList = c()

for(aIdent in levels(agg.tdtom.fibroblasts)){
  aseurat = subset(agg.tdtom.fibroblasts, idents = aIdent)
  seuratList = c(seuratList, aseurat)
}
rm(agg.tdtom.fibroblasts, aseurat)
gc()
seuratList <- lapply(X = seuratList, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seuratList, nfeatures = 3000)
seuratList <- PrepSCTIntegration(object.list = seuratList, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features, normalization.method = "SCT")
rm(seuratList)
gc()
agg.tdtom.fibroblasts <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(anchors)
gc()
DefaultAssay(agg.tdtom.fibroblasts) = "integrated"
agg.tdtom.fibroblasts = RunPCA(agg.tdtom.fibroblasts, npcs = 35)
agg.tdtom.fibroblasts = FindNeighbors(agg.tdtom.fibroblasts, dims = 1:35)
agg.tdtom.fibroblasts = FindClusters(agg.tdtom.fibroblasts, resolution = (seq(0, 1, by = .05)))
clustree(agg.tdtom.fibroblasts, prefix = "integrated_snn_res.")
agg.tdtom.fibroblasts = FindClusters(agg.tdtom.fibroblasts, resolution = .65)
agg.tdtom.fibroblasts <- RunUMAP(agg.tdtom.fibroblasts, reduction = "pca", dims = 1:35)
DefaultAssay(agg.tdtom.fibroblasts) = "SCT"

#again, depending on packages, clusters may move around a bit. Clusters were chosen to eliminate contaminating tdtom+ cells which were obviously
#immune or endothelial cells, look for Ptprc and Pecam1
#subset was reintegrated 

agg.tdtom.fibroblasts2 = subset(agg.tdtom.fibroblasts2, idents = c('4','16','14','0','8','1','2','3','12','16','7'))
DefaultAssay(agg.tdtom.fibroblasts2) = "integrated"
agg.tdtom.fibroblasts2 = RunPCA(agg.tdtom.fibroblasts2,npcs = 35)
agg.tdtom.fibroblasts2 = FindNeighbors(agg.tdtom.fibroblasts2, dims = 1:35)
agg.tdtom.fibroblasts2 = FindClusters(agg.tdtom.fibroblasts2, resolution = (seq(0, 1, by = .05)))
clustree(agg.tdtom.fibroblasts2, prefix = "integrated_snn_res.",show_axis = TRUE)
agg.tdtom.fibroblasts2 = FindClusters(agg.tdtom.fibroblasts2, resolution = .65)
agg.tdtom.fibroblasts2 <- RunUMAP(agg.tdtom.fibroblasts2, reduction = "pca", dims = 1:35)
DefaultAssay(agg.tdtom.fibroblasts2) = "SCT"
agg.tdtom.fibroblasts2 = PrepSCTFindMarkers(agg.tdtom.fibroblasts2)

#Creates figure 2a
png('Umap_subset_.png', units = 'in', width = 8, height = 7, res = 600)
DimPlot(agg.tdtom.fibroblasts2,reduction = "umap",label = TRUE,label.size = 6,pt.size = 1.5,repel = TRUE) + 
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  aes(UMAP1,UMAP2) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow(type = 'closed', length = unit(10,'pt'))), axis.title = element_text(hjust = 0))+
  labs(title = "UMAP",subtitle = "Fibroblasts and TdTomato Subset")
dev.off()

#this creates a metadata, not necessary for figure generation but seurat does not like it when the active ident starts with a number

agg.tdtom.fibroblasts2$SC = 'missing'
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '0'] <- "Cluster 0"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '1'] <- "Cluster 1"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '2'] <- "Cluster 2"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '3'] <- "Cluster 3"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '4'] <- "Cluster 4"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '5'] <- "Cluster 5"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '6'] <- "Cluster 6"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '7'] <- "Cluster 7"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '8'] <- "Cluster 8"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '9'] <- "Cluster 9"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '10'] <- "Cluster 10"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '11'] <- "Cluster 11"
agg.tdtom.fibroblasts2$SC[agg.tdtom.fibroblasts2$seurat_clusters == '12'] <- "Cluster 12"

agg.tdtom.fibroblasts2$SC <- factor(x = agg.tdtom.fibroblasts2$SC, levels = c('Cluster 0','Cluster 3','Cluster 11','Cluster 9','Cluster 12',
                                                                              'Cluster 1','Cluster 2','Cluster 8','Cluster 10','Cluster 6',
                                                                              'Cluster 4','Cluster 5','Cluster 7'))
Idents(agg.tdtom.fibroblasts2) = 'SC'

#creates figure 2C

ID.averages <- AggregateExpression(agg.tdtom.fibroblasts2, return.seurat = TRUE)
Figure_genes = c('Adam12','Adamts2','Ctsh','Col1a1','Col5a1','Col12a1','Col14a1','Mmp23','Loxl2','Loxl3','Dpt','Cthrc1','Tgfbr2','Acta2','Postn','Pi16',
                 'Adamts12','Ctsl','Col4a1','Col15a1','Plod1','Serpine2','Bpgm','Cilp','Col8a1','Dbp','Fzd1','Il34',
                 'Mmp3','Piezo2','Spp1','Thbs4','Itgbl1','Phlda3','Fgf10','Adh1','Ltbp4','Smoc2')
DoHeatmap(ID.averages, features = Figure_genes,angle = 45,size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE) + guides(colour = FALSE)

#creates figure 2B
agg.tdtom.fibroblasts2$Groups = 'missing'
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '0'] <- "G1"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '1'] <- "G2"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '2'] <- "G2"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '3'] <- "G1"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '4'] <- "G3"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '5'] <- "G3"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '6'] <- "G3"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '7'] <- "G3"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '8'] <- "G2"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '9'] <- "G1"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '10'] <- "G2"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '11'] <- "G1"
agg.tdtom.fibroblasts2$Groups[agg.tdtom.fibroblasts2$seurat_clusters == '12'] <- "G1"
Idents(agg.tdtom.fibroblasts2) = 'Groups'

DimPlot(agg.tdtom.fibroblasts2,reduction = "umap",label = TRUE,label.size = 8,pt.size = 1.75,repel = TRUE) + 
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  aes(UMAP1,UMAP2) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow(type = 'closed', length = unit(10,'pt'))), axis.title = element_text(hjust = 0))+
  labs(title = "UMAP")

#adjusts the order of the groups
agg.tdtom.fibroblasts2$Groups <- factor(x = agg.tdtom.fibroblasts2$Groups, levels = c('G1','G2','G3'))
Idents(agg.tdtom.fibroblasts2) = 'Groups'

#figure 3 code
#this generates differential expressed genes by group
fibroblast_expression_group_1 = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G1", ident.2 = c("G2",'G3'), only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_2 = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G2", ident.2 = c("G1",'G3'), only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_3 = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G3", ident.2 = c("G2",'G1'), only.pos = TRUE, min.pct = .25, logfc.threshold = .25)

all.markers <- FindAllMarkers(object = agg.tdtom.fibroblasts2, only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
ID.averages <- AggregateExpression(agg.tdtom.fibroblasts2, return.seurat = TRUE)
png('heatmap_G1_.png', units = 'in', width = 8, height = 7, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_1)[1:50],angle = 45,size = 6,
          group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()
png('heatmap_G2_.png', units = 'in', width = 8, height = 7, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_2)[1:50],angle = 45,size = 6,
          group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()
png('heatmap_G3_.png', units = 'in', width = 8, height = 7, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_3)[1:50],angle = 45,size = 6,
          group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()
#this creates the CSVs used for GoTerm analysis, set your own path
write.csv(fibroblast_expression_group_1, file = "")
write.csv(fibroblast_expression_group_2, file = "")
write.csv(fibroblast_expression_group_3, file = "")

#figure 4a
agg.tdtom.fibroblasts2$dataset <- factor(x = agg.tdtom.fibroblasts2$dataset, levels = c('Homeostatic','1-day','3-day','7-day'))

png('umap_groups_dataset.png', units = 'in', width = 8, height = 7, res = 600)
DimPlot(agg.tdtom.fibroblasts2,reduction = "umap",group.by = 'Groups',label = FALSE,label.size = 8,pt.size = 1.75,repel = TRUE, split.by = 'dataset', ncol = 2) + 
  labs(title = '') +aes(UMAP1,UMAP2) +theme(text = element_text(size = 24))
dev.off()

#Creates metadata dividing groups by timepoint
agg.tdtom.fibroblasts2$G_D <- sprintf("%s_%s", agg.tdtom.fibroblasts2$Groups, agg.tdtom.fibroblasts2$dataset)
Idents(agg.tdtom.fibroblasts2) = 'G_D'
agg.tdtom.fibroblasts2$G_D <- factor(x = agg.tdtom.fibroblasts2$G_D, levels = c('G1_Homeostatic','G1_1-day',
                                                                                'G1_3-day',
                                                                                'G1_7-day','G2_Homeostatic',
                                                                                'G2_1-day','G2_3-day',
                                                                                'G2_7-day','G3_Homeostatic','G3_1-day',
                                                                                'G3_3-day', 'G3_7-day'
))
Idents(agg.tdtom.fibroblasts2) = 'G_D'

#Cell numbers for figure 4B
table(Idents(agg.tdtom.fibroblasts2))

#Figure 4C
fibroblast_expression_group_1_1day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G1_1-day", ident.2 = "G1_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_1_3day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G1_3-day", ident.2 = "G1_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_1_7day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G1_7-day", ident.2 = "G1_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_2_1day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G2_1-day", ident.2 = "G2_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_2_3day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G2_3-day", ident.2 = "G2_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_2_7day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G2_7-day", ident.2 = "G2_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_3_1day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G3_1-day", ident.2 = "G3_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_3_3day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G3_3-day", ident.2 = "G3_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)
fibroblast_expression_group_3_7day = FindMarkers(agg.tdtom.fibroblasts2, ident.1 = "G3_7-day", ident.2 = "G3_Homeostatic", only.pos = TRUE, min.pct = .25, logfc.threshold = .25)

write.csv(fibroblast_expression_group_1_1day, file = "")
write.csv(fibroblast_expression_group_1_3day, file = "")
write.csv(fibroblast_expression_group_1_7day, file = "")
write.csv(fibroblast_expression_group_2_1day, file = "")
write.csv(fibroblast_expression_group_2_3day, file = "")
write.csv(fibroblast_expression_group_2_7day, file = "")
write.csv(fibroblast_expression_group_3_1day, file = "")
write.csv(fibroblast_expression_group_3_3day, file = "")
write.csv(fibroblast_expression_group_3_7day, file = "")

agg.tmp = subset(agg.tdtom.fibroblasts2, idents = c('G1_Homeostatic','G1_1-day','G1_3-day','G1_7-day'))
ID.averages <- AggregateExpression(agg.tmp, return.seurat = TRUE)
png('heatmap_G1_1day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_1_1day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()
png('heatmap_G1_3day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_1_3day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()
png('heatmap_G1_7day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_1_7day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()

agg.tmp = subset(agg.tdtom.fibroblasts2, idents = c('G2_Homeostatic','G2_1-day','G2_3-day','G2_7-day'))
ID.averages <- AggregateExpression(agg.tmp, return.seurat = TRUE)

png('heatmap_G2_1day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_2_1day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()

png('heatmap_G2_3day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_2_3day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()

png('heatmap_G2_7day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_2_7day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()

agg.tmp = subset(agg.tdtom.fibroblasts2, idents = c('G3_Homeostatic','G3_1-day','G3_3-day','G3_7-day'))
ID.averages <- AggregateExpression(agg.tmp, return.seurat = TRUE)

png('heatmap_G3_1day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_3_1day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()

png('heatmap_G3_3day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_3_3day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()

png('heatmap_G3_7day.png', units = 'in', width = 8, height = 8, res = 600)
DoHeatmap(ID.averages, features = rownames(fibroblast_expression_group_3_7day)[1:50],angle = 45,
          size = 6,group.bar.height = .04,draw.lines = FALSE,label = TRUE)+guides(colour = FALSE)
dev.off()

#Figure 5 codes

geneList <- c()
geneList[["Inflammation"]] <- c("Abca1",'Acvr1b','Acvr2a','Adgre1','Adm','Adora2b','Adrm1','Ahr','Aplnr',
                                'Aqp9','Atp2a2','Atp2b1','Atp2c1','Axl','Bdkrb1','Best1','Bst2','Btg2','C3ar1',
                                'Calcrl','Ccl17','Ccl5','Ccl7','Ccr7','Ccrl2','Cd14','Cd40','Cd48','Cd69','Cd82',
                                'Cdkn1a','Chst2','Clec5a','Cmklr1','Csf1','Csf3','Csf3r','Cx3cl1','Cxcl10','Cxcl11',
                                'Cxcl5','Cxcl9','Cxcr6','Cybb','Dcbld2','Ebi3','Edn1','Eif2ak2','Emp3','Ereg','F3',
                                'Ffar2','Fpr1','Fzd5','Gabbr1','Gch1','Gna15','Gnai3','Gp1ba','Gpc3','Gpr132','Gpr183',
                                'Has2','Hbegf','Hif1a','Hpn','Hrh1','Icam1','Icam4','Icosl','Ifitm1','Ifnar1','Ifngr2',
                                'Il10','Il10ra','Il15','Il15ra','Il18','Il18r1','Il18rap','Il1a','Il1b','Il1r1','Il2rb',
                                'Il4ra','Il6','Il7r','Inhba','Irak2','Irf1','Irf7','Itga5','Itgb3','Itgb8','Kcna3','Kcnj2',
                                'Kif1b','Klf6','Lck','Lcp2','Ldlr','Lif','Lpar1','Lta','Ly6e','Lyn','Mefv','Met',
                                'Mmp14','Msr1','Mxd1','Myc','Nampt','Ndp','Nfkb1','Nfkbia','Nlrp3','Nmi','Nod2',
                                'Olr1','Osm','Osmr','P2rx4','P2rx7','P2ry2','Pcdh7','Pde4b','Pdpn','Pik3r5'
                                ,'Plaur','Psen1','Ptafr','Ptger2','Ptger4','Ptgir','Ptpre','Pvr','Raf1','Rasgrp1',
                                'Rela','Rgs1','Rgs16','Rhog','Ripk2','Rnf144b','Rtp4','Scarf1','Scn1b','Sele','Selenos'
                                ,'Sell','Sema4d','Serpine1','Sgms2','Slamf1','Slc11a2','Slc1a2','Slc28a2','Slc31a1','Slc31a2'
                                ,'Slc4a4','Slc7a1','Slc7a2','Sphk1','Sri','Stab1','Tacr1','Tacr3','Tapbp','Timp1','Tlr1','Tlr2'
                                ,'Tlr3','Tnfaip6','Tnfrsf1b','Tnfrsf9','Tnfsf10','Tnfsf15','Tnfsf9','Tpbg')
geneList[["Collagen"]] <- c('Col10a1','Col11a1','Col11a2','Col12a1','Col13a1','Col14a1','Col15a1','Col16a1','Col17a1',
                            'Col18a1','Col1a1','Col1a2','Col20a1',
                            'Col23a1','Col24a1','Col25a1','Col26a1','Col27a1','Col28a1','Col3a1',
                            'Col4a1','Col4a2','Col4a3','Col4a4','Col4a5','Col4a6','Col5a1','Col5a2','Col5a3','Col6a1',
                            'Col6a2','Col6a3','Col6a4','Col6a5','Col6a6','Col7a1','Col8a1','Col8a2','Col9a2')
geneList[["Fibril_collagen"]] <- c('Col1a1','Col1a2','Col2a1','Col3a1','Col5a1','Col5a2','Col11a1','Col11a2','Col23a1','Col27a1')
geneList[["Network_collagen"]] <- c('Col4a1','Col4a2','Col4a3','Col4a4','Col4a5','Col4a6','Col8a1','Col8a2','Col10a1')
geneList[["ADAM"]] <- c('Adam10','Adam11','Adam12','Adam15','Adam17','Adam18','Adam19',
                        'Adam1a','Adam1b','Adam2','Adam20','Adam21','Adam22','Adam23','Adam24','Adam25',
                        'Adam26a','Adam26b','Adam28','Adam29','Adam30','Adam32','Adam33','Adam34',
                        'Adam39','Adam4','Adam5','Adam6a','Adam6b','Adam7','Adam8','Adam9','Adamdec1',
                        'Adamts1','Adamts10','Adamts12','Adamts13','Adamts14','Adamts15','Adamts16','Adamts17'
                        ,'Adamts18','Adamts19','Adamts2','Adamts20','Adamts3','Adamts4','Adamts5','Adamts6','Adamts7'
                        ,'Adamts8','Adamts9','Adamtsl1','Adamtsl2','Adamtsl3','Adamtsl4','Adamtsl5')
geneList[["Cystatin"]] <- c('Cst10','Cst11','Cst12','Cst13','Cst3','Cst6','Cst7','Cst8','Cst9','Csta1',
                            'Cstb','Cstl1')
geneList[["Cathepsin"]] <- c('Cts3','Cts6','Cts7','Cts8','Ctsa','Ctsb','Ctsc','Ctsd','Ctse','Ctsf','Ctsg',
                             'Ctsh','Ctsj','Ctsk','Ctsl','Ctsll3','Ctsm','Ctso','Ctsq','Ctsr','Ctss','Ctsw','Ctsz')
geneList[["Mmps"]] <- c('Mmp10','Mmp11','Mmp12'
                        ,'Mmp13','Mmp14','Mmp15','Mmp16','Mmp17','Mmp19','Mmp1a','Mmp1b','Mmp2','Mmp20','Mmp21','Mmp23',
                        'Mmp25','Mmp27','Mmp28','Mmp3','Mmp7','Mmp8','Mmp9')
geneList[['Serpin']] = c('Serpina10','Serpina11',
                         'Serpina12','Serpina1a','Serpina1b','Serpina1c','Serpina1d','Serpina1e','Serpina1f','Serpina3a'
                         ,'Serpina3b','Serpina3c','Serpina3f','Serpina3g','Serpina3k','Serpina3m','Serpina3n','Serpina5'
                         ,'Serpina6','Serpina7','Serpina9','Serpinb10','Serpinb11','Serpinb12','Serpinb13','Serpinb1a'
                         ,'Serpinb1b','Serpinb1c','Serpinb2','Serpinb3a','Serpinb3b','Serpinb3c','Serpinb3d','Serpinb5'
                         ,'Serpinb6a','Serpinb6b','Serpinb6c','Serpinb6d','Serpinb7','Serpinb8','Serpinb9','Serpinb9b',
                         'Serpinb9c','Serpinb9d','Serpinb9e','Serpinb9f','Serpinb9g','Serpinc1','Serpind1','Serpine1',
                         'Serpine2','Serpine3','Serpinf1','Serpinf2','Serping1','Serpinh1','Serpini1','Serpini2')
geneList[['Regulators']] = c('A2m','Adam10','Adam11','Adam12','Adam15','Adam17','Adam18','Adam19',
                             'Adam1a','Adam1b','Adam2','Adam20','Adam21','Adam22','Adam23','Adam24','Adam25',
                             'Adam26a','Adam26b','Adam28','Adam29','Adam30','Adam32','Adam33','Adam34',
                             'Adam39','Adam4','Adam5','Adam6a','Adam6b','Adam7','Adam8','Adam9','Adamdec1',
                             'Adamts1','Adamts10','Adamts12','Adamts13','Adamts14','Adamts15','Adamts16','Adamts17'
                             ,'Adamts18','Adamts19','Adamts2','Adamts20','Adamts3','Adamts4','Adamts5','Adamts6','Adamts7'
                             ,'Adamts8','Adamts9','Adamtsl1','Adamtsl2','Adamtsl3','Adamtsl4','Adamtsl5','Agt','AI182371',
                             'Ambp','BC048546','BC051665','BC100530','BC117090','Bmp1','Cd109','Cela1','Cela2a',
                             'Cela3b','Cpn2','Cst10','Cst11','Cst12','Cst13','Cst3','Cst6','Cst7','Cst8','Cst9','Csta1',
                             'Cstb','Cstl1','Cts3','Cts6','Cts7','Cts8','Ctsa','Ctsb','Ctsc','Ctsd','Ctse','Ctsf','Ctsg',
                             'Ctsh','Ctsj','Ctsk','Ctsl','Ctsll3','Ctsm','Ctso','Ctsq','Ctsr','Ctss','Ctsw','Ctsz','Egln1',
                             'Egln2','Egln3','Elane','F10','F13a1','F13b','F2','F7','F9','Fam20a','Fam20b','Fam20c',
                             'Gm10334','Gm13011','Gm4758','Gm4787','Gm5483','Gm5771','Habp2','Hpse'
                             ,'Hpse2','Hrg','Htra1','Htra3','Htra4','Hyal1','Hyal2','Hyal3','Hyal4','Hyal5','Hyal6','Itih1'
                             ,'Itih2','Itih3','Itih4','Itih5','Kazald1','Kng1','Kng2','Lepre1','Leprel1','Leprel2',
                             'Lox','Loxl1','Loxl2','Loxl3','Loxl4','Masp1','Masp2','Mep1a','Mep1b','Mmp10','Mmp11','Mmp12'
                             ,'Mmp13','Mmp14','Mmp15','Mmp16','Mmp17','Mmp19','Mmp1a','Mmp1b','Mmp2','Mmp20','Mmp21','Mmp23',
                             'Mmp25','Mmp27','Mmp28','Mmp3','Mmp7','Mmp8','Mmp9','Mug2','Ngly1','Ogfod1','Ogfod2',
                             'P4ha1','P4ha2','P4ha3','P4htm','Pamr1','Pappa','Pappa2','Pcsk5','Pcsk6','Plat','Plau',
                             'Plod1','Plod2','Plod3','Prss12','Pzp','Serpina10','Serpina11',
                             'Serpina12','Serpina1a','Serpina1b','Serpina1c','Serpina1d','Serpina1e','Serpina1f','Serpina3a'
                             ,'Serpina3b','Serpina3c','Serpina3f','Serpina3g','Serpina3k','Serpina3m','Serpina3n','Serpina5'
                             ,'Serpina6','Serpina7','Serpina9','Serpinb10','Serpinb11','Serpinb12','Serpinb13','Serpinb1a'
                             ,'Serpinb1b','Serpinb1c','Serpinb2','Serpinb3a','Serpinb3b','Serpinb3c','Serpinb3d','Serpinb5'
                             ,'Serpinb6a','Serpinb6b','Serpinb6c','Serpinb6d','Serpinb7','Serpinb8','Serpinb9','Serpinb9b',
                             'Serpinb9c','Serpinb9d','Serpinb9e','Serpinb9f','Serpinb9g','Serpinc1','Serpind1','Serpine1',
                             'Serpine2','Serpine3','Serpinf1','Serpinf2','Serping1','Serpinh1','Serpini1','Serpini2','Slpi'
                             ,'Spam1','Spinkl','St14','Stfa1','Stfa2','Stfa2l1','Stfa3','Sulf1','Sulf2','Tgm1','Tgm2','Tgm3',
                             'Tgm4','Tgm5','Tgm6','Tgm7','Timp1','Timp2','Timp3','Timp4','Tll1','Tll2','Tmprss15','Tpbpa',
                             'Tpbpb','Try10','Try4','Try5')
geneList[['Secreted']] = c('Amh','Angpt1','Angpt2','Angpt4','Angptl1','Angptl2','Angptl3',
                           'Angptl4','Angptl6','Angptl7','Areg','Artn','Bdnf','Bmp10','Bmp15','Bmp2',
                           'Bmp3','Bmp4','Bmp5','Bmp6','Bmp7','Bmp8a','Bmp8b','Brinp2','Brinp3','Btc',
                           'Cbln1','Cbln2','Cbln3','Cbln4','Ccbe1','Ccl1','Ccl11','Ccl12','Ccl17','Ccl19',
                           'Ccl2','Ccl20','Ccl21a','Ccl21b','Ccl21c','Ccl22','Ccl24','Ccl25','Ccl26','Ccl27a'
                           ,'Ccl28','Ccl3','Ccl4','Ccl5','Ccl6','Ccl7','Ccl8','Ccl9','Cfc1','Chrd','Chrdl1'
                           ,'Chrdl2','Clcf1','Cntf','Crhbp','Crlf1','Crlf3','Crnn','Csf1','Csf2','Csf3','Ctf1'
                           ,'Ctf2','Cx3cl1','Cxcl1','Cxcl10','Cxcl11','Cxcl12','Cxcl13','Cxcl14','Cxcl15',
                           'Cxcl2','Cxcl3','Cxcl5','Cxcl9','Dhh','Ebi3','Eda','Egf','Egfl6','Egfl7','Egfl8',
                           'Epgn','Epo','Ereg','Fam132a','Fam132b','Fasl','Fgf1','Fgf10','Fgf11','Fgf12','Fgf13'
                           ,'Fgf14','Fgf15','Fgf16','Fgf17','Fgf18','Fgf2','Fgf20','Fgf21','Fgf22','Fgf23'
                           ,'Fgf3','Fgf4','Fgf5','Fgf6','Fgf7','Fgf8','Fgf9','Fgfbp1','Fgfbp3','Figf','Flg',
                           'Flg2','Flt3l','Frzb','Fst','Fstl1','Fstl3','Gdf1','Gdf10','Gdf11','Gdf15','Gdf2',
                           'Gdf3','Gdf5','Gdf6','Gdf7','Gdf9','Gdnf','Gh','Gm13271','Gm13272','Gm13275','Gm13276'
                           ,'Gm13277','Gm13278','Gm13279','Gm13283','Gm13285','Gm13287','Gm13288','Gm13289'
                           ,'Gm13290','Gm13306','Gm5849','Hbegf','Hcfc1','Hcfc2','Hgf','Hgfac','Hhip','Hrnr'
                           ,'Ifna1','Ifna11','Ifna12','Ifna13','Ifna14','Ifna15','Ifna16','Ifna2','Ifna4'
                           ,'Ifna5','Ifna6','Ifna7','Ifna9','Ifnab','Ifnb1','Ifne','Ifng','Ifnk','Ifnz','Igf1'
                           ,'Igf2','Ihh','Il10','Il11','Il12a','Il12b','Il13','Il15','Il16','Il17a','Il17b',
                           'Il17c','Il17d','Il17f','Il18','Il19','Il1a','Il1b','Il1f10','Il1f5','Il1f6','Il1f8'
                           ,'Il1f9','Il1rn','Il2','Il20','Il22','Il23a','Il24','Il25','Il3','Il34','Il4','Il5',
                           'Il6','Il7','Il9','Inha','Inhba','Inhbb','Inhbc','Inhbe','Ins1','Ins2','Insl3',
                           'Insl5','Insl6','Ism1','Ism2','Kitl','Lefty1','Lefty2','Lep','Lif','Lta','Ltb','Mdk'
                           ,'Megf10','Megf11','Megf6','Megf8','Megf9','Mst1','Mstn','Ngf','Nodal','Nrg1','Nrg2'
                           ,'Nrg3','Nrg4','Nrtn','Ntf3','Ntf5','Osm','Pdgfa','Pdgfb','Pdgfc','Pdgfd','Pf4','Pgf'
                           ,'Pik3ip1','Ppbp','Prl','Prl2a1','Prl2b1','Prl2c1','Prl2c2','Prl2c3','Prl2c4','Prl2c5'
                           ,'Prl3a1','Prl3b1','Prl3c1','Prl3d1','Prl3d2','Prl3d3','Prl4a1','Prl5a1','Prl6a1',
                           'Prl7a1','Prl7a2','Prl7b1','Prl7c1','Prl7d1','Prl8a1','Prl8a2','Prl8a6','Prl8a8',
                           'Prl8a9','Pspn','Ptn','Rptn','S100a1','S100a10','S100a11','S100a13','S100a14',
                           'S100a16','S100a2','S100a3','S100a4','S100a5','S100a6','S100a7a','S100a8','S100a9'
                           ,'S100b','S100g','S100z','Scube1','Scube2','Scube3','Sfrp1','Sfrp2','Sfrp4','Sfrp5',
                           'Shh','Tchh','Tchhl1','Tdgf1','Tgfa','Tgfb1','Tgfb2','Tgfb3','Thpo','Tnf','Tnfsf10'
                           ,'Tnfsf11','Tnfsf12','Tnfsf12tnfsf13','Tnfsf13','Tnfsf13b','Tnfsf14','Tnfsf15',
                           'Tnfsf18','Tnfsf4','Tnfsf8','Tnfsf9','Tpo','Vegfa','Vegfb','Vegfc','Vwc2','Vwc2l'
                           ,'Wfikkn1','Wfikkn2','Wif1','Wnt1','Wnt10a','Wnt10b','Wnt11','Wnt16','Wnt2','Wnt2b'
                           ,'Wnt3','Wnt3a','Wnt4','Wnt5a','Wnt5b','Wnt6','Wnt7a','Wnt7b','Wnt8a','Wnt8b','Wnt9a'
                           ,'Wnt9b','Xcl1')
geneList[['Glycoproteins']] = c('5430419D17Rik','Abi3bp','Adipoq','Aebp1','Agrn','Ambn','Amelx',
                                'AW551984','Bglap2','Bglap3','Bmper','Bsph1','Bsph2','Cdcp2','Cilp','Cilp2',
                                'Coch','Colq','Comp','Creld1','Creld2','Crim1','Crispld1','Crispld2','Ctgf','Cthrc1',
                                'Cyr61','Ddx26b','Dmbt1','Dmp1','Dpt','Dspp','Ecm1','Ecm2','Edil3','Efemp1','Efemp2'
                                ,'Egfem1','Egflam','Eln','Emid1','Emilin1','Emilin2','Emilin3','Fbln1','Fbln2','Fbln5'
                                ,'Fbln7','Fbn1','Fbn2','Fga','Fgb','Fgg','Fgl1','Fgl2','Fn1','Fndc1','Fndc7','Fndc8',
                                'Fras1','Gas6','Gldn','Hmcn1','Hmcn2','Ibsp','Igfals','Igfbp1','Igfbp2','Igfbp3','Igfbp4'
                                ,'Igfbp5','Igfbp6','Igfbp7','Igfbpl1','Igsf10','Kcp','Lama1','Lama2','Lama3','Lama4',
                                'Lama5','Lamb1','Lamb2','Lamb3','Lamc1','Lamc2','Lamc3','Lgi1','Lgi2','Lgi3','Lgi4',
                                'Lrg1','Ltbp1','Ltbp2','Ltbp3','Ltbp4','Matn1','Matn2','Matn3','Matn4','Mepe','Mfap1a'
                                ,'Mfap1b','Mfap2','Mfap3','Mfap4','Mfap5','Mfge8','Mgp','Mmrn1','Mmrn2','Ndnf','Nell1'
                                ,'Nell2','Nid1','Nid2','Nov','Npnt','Ntn1','Ntn3','Ntn4','Ntn5','Ntng1','Ntng2','Oit3',
                                'Otog','Otogl','Otol1','Papln','Pcolce','Pcolce2','Postn','Pxdn','Reln','Rspo1','Rspo2'
                                ,'Rspo3','Rspo4','Sbspon','Slamf6','Slit1','Slit2','Slit3','Smoc1','Smoc2','Sned1','Sparc'
                                ,'Sparcl1','Spon1','Spon2','Spp1','Srpx','Srpx2','Sspo','Svep1','Tecta','Tectb','Tgfbi',
                                'Thbs1','Thbs2','Thbs3','Thbs4','Thsd4','Tinag','Tinagl1','Tnc','Tnfaip6','Tnn','Tnr','Tnxb'
                                ,'Tsku','Tspear','Vit','Vtn','Vwa1','Vwa2','Vwa3a','Vwa3b','Vwa5a','Vwa5b1','Vwa5b2','Vwa7',
                                'Vwa9','Vwce','Vwde','Vwf','Wisp1','Wisp2','Wisp3','Zp1','Zp2','Zp3','Zp3r','Zpld1')
geneList[['Bmp']] = c('Bmp10','Bmp15','Bmp2',
                      'Bmp3','Bmp4','Bmp5','Bmp6','Bmp7','Bmp8a','Bmp8b')
geneList[['Chemokine']] = c('Ccl1','Ccl11','Ccl12','Ccl17','Ccl19',
                            'Ccl2','Ccl20','Ccl21a','Ccl21b','Ccl21c','Ccl22','Ccl24','Ccl25','Ccl26','Ccl27a'
                            ,'Ccl28','Ccl3','Ccl4','Ccl5','Ccl6','Ccl7','Ccl8','Ccl9','Cx3cl1','Cxcl1','Cxcl10','Cxcl11','Cxcl12','Cxcl13','Cxcl14','Cxcl15',
                            'Cxcl2','Cxcl3','Cxcl5','Cxcl9')
geneList[['FGFs']] = c('Fgf1','Fgf10','Fgf11','Fgf12','Fgf13'
                       ,'Fgf14','Fgf15','Fgf16','Fgf17','Fgf18','Fgf2','Fgf20','Fgf21','Fgf22','Fgf23'
                       ,'Fgf3','Fgf4','Fgf5','Fgf6','Fgf7','Fgf8','Fgf9')
geneList[['Interferons']] = c('Ifna1','Ifna11','Ifna12','Ifna13','Ifna14','Ifna15','Ifna16','Ifna2','Ifna4'
                              ,'Ifna5','Ifna6','Ifna7','Ifna9','Ifnab','Ifnb1','Ifne','Ifng','Ifnk','Ifnz')
geneList[['Interleukins']] = c('Il10','Il11','Il12a','Il12b','Il13','Il15','Il16','Il17a','Il17b',
                               'Il17c','Il17d','Il17f','Il18','Il19','Il1a','Il1b','Il1f10','Il1f5','Il1f6','Il1f8'
                               ,'Il1f9','Il1rn','Il2','Il20','Il22','Il23a','Il24','Il25','Il3','Il34','Il4','Il5',
                               'Il6','Il7','Il9')
geneList[['Tnf']] = c('Tnf','Tnfsf10'
                      ,'Tnfsf11','Tnfsf12','Tnfsf12tnfsf13','Tnfsf13','Tnfsf13b','Tnfsf14','Tnfsf15',
                      'Tnfsf18','Tnfsf4','Tnfsf8','Tnfsf9')
geneList[['Wnts']] = c('Wnt1','Wnt10a','Wnt10b','Wnt11','Wnt16','Wnt2','Wnt2b'
                       ,'Wnt3','Wnt3a','Wnt4','Wnt5a','Wnt5b','Wnt6','Wnt7a','Wnt7b','Wnt8a','Wnt8b','Wnt9a'
                       ,'Wnt9b')
geneList[['Proteoglycans']] = c('Acan','Aspn','Bcan','Bgn','Chad','Chadl','Dcn','Epyc','Esm1',
                                'Fmod','Hapln1','Hapln2','Hapln3','Hapln4','Hspg2','Impg1','Impg2','Kera','Lum',
                                'Ncan','Nepn','Nyx','Ogn','Omd','Optc','Podn','Podnl1','Prelp','Prg2','Prg3','Prg4',
                                'Spock1','Spock2','Spock3','Srgn','Vcan')

#loop below use the lists above and counts how many genes are expressed by a cell and adds a value of 1 per cell per gene expressed.
#this meta data can be used to generate figures

for (alistname in names(geneList)) {
  clean.genes = c()
  for (aIdent in geneList[[alistname]]) {
    if(length(grep(aIdent, rownames(agg.tdtom.fibroblasts2@assays$SCT)))){
      if(aIdent == rownames(agg.tdtom.fibroblasts2@assays$SCT)[grep(aIdent, rownames(agg.tdtom.fibroblasts2@assays$SCT))[[1]]])
        clean.genes = c(clean.genes,aIdent)
    }
    
  }
  metaDataName = sprintf('%s.count',alistname)
  agg.tdtom.fibroblasts2@meta.data[[metaDataName]] <- 0
  for(aIdent in clean.genes){
    agg.tdtom.fibroblasts2@meta.data[[metaDataName]] = ifelse((agg.tdtom.fibroblasts2@assays$SCT@counts[aIdent, ] >= 1),
                                                              (agg.tdtom.fibroblasts2@meta.data[[metaDataName]] + 1),
                                                              (agg.tdtom.fibroblasts2@meta.data[[metaDataName]]))
  }
}

#loop below use the lists above and counts how many transcripts are expressed by a cell for each gene and 
#adds a value of all transcripts per cell per gene expressed.
#this meta data can be used to generate figures

for (alistname in names(geneList)) {
  clean.genes = c()
  for (aIdent in geneList[[alistname]]) {
    if(length(grep(aIdent, rownames(agg.tdtom.fibroblasts2@assays$SCT)))){
      if(aIdent == rownames(agg.tdtom.fibroblasts2@assays$SCT)[grep(aIdent, rownames(agg.tdtom.fibroblasts2@assays$SCT))[[1]]])
        clean.genes = c(clean.genes,aIdent)
    }
    
  }
  metaDataName = sprintf('%s.total',alistname)
  agg.tdtom.fibroblasts2@meta.data[[metaDataName]] <- 0
  for(aIdent in clean.genes){
    agg.tdtom.fibroblasts2@meta.data[[metaDataName]] = 
      ifelse((agg.tdtom.fibroblasts2@assays$SCT@counts[aIdent, ] >= 1),
             (agg.tdtom.fibroblasts2@meta.data[[metaDataName]] + agg.tdtom.fibroblasts2@assays$SCT@counts[aIdent, ]),
             (agg.tdtom.fibroblasts2@meta.data[[metaDataName]]))
  }
}

#Figure 6A
VlnPlot(agg.tdtom.fibroblasts2, features = 
          c('Collagen.total','Inflammation.total','Regulators.total','Glycoproteins.total','Secreted.total','Proteoglycans.total'),
        stack = TRUE,pt.size = 0,log = FALSE, flip = TRUE)+theme(legend.position =  'none')+
  labs(title = "Gene Total By Cluster",subtitle = "Fibroblasts and TdTomato Subset")
#Figure 6B
DotPlot(agg.tdtom.fibroblasts2, features = 
          c('Proteoglycans.count','Secreted.count','Glycoproteins.count','Regulators.count','Inflammation.count','Collagen.count'),
        cols = c("lightblue", "red3"),
        scale.min = 0 ,scale.max = 100)+coord_flip()+
  theme(text = element_text(size = 15),axis.text = element_text(size = 15)) + theme_bw()+RotatedAxis() +
  labs(title = "Genes Expressed By Cluster",subtitle = "Fibroblasts and TdTomato Subset") 

#Figure 6 data takes longer to compute, I recommend trimming the geneList to the specific lists you are interested in. Otherwise
#it will take longer to compute and gice you a huge dataframe where 90% of it would be unused
Idents(agg.tdtom.fibroblasts2) = 'G_D'
#the dataframe columns must be generated before hand
#gene is the gene in said list
gene = c('Yes')
#group is the fibroblast group, because the ident is G_D, the group will be the group divided by timepoint, if you want to look at timepoint
#or the group seperately, you can add the values or you can change the ident in Idents code above
group = c('Yes')
#transcripts shows how many transcripts expressed in total by the group
transcripts = c('Yes')
#cell number is the number of cells per identity
Cell_Number = c('Yes')
#to normalize for cell number differences, I divide the transcripts by the number of cells
Per = c('Yes')
#category is the genelist category
Category = c('Yes')
df = data.frame(Group = group, Category = Category, Genes = gene, Transcripts = transcripts, Cell_Number = Cell_Number, Per = Per)
for (adataset in unique(agg.tdtom.fibroblasts2$G_D)) {
  agg.tmp = subset(agg.tdtom.fibroblasts2, idents = adataset)
  cells = ncol(agg.tmp)
  for (alistname in names(geneList)) {
    clean.genes = c()
    for (aIdent in geneList[[alistname]]) {
      if(length(grep(aIdent, rownames(agg.tmp@assays$SCT)))){
        if(aIdent == rownames(agg.tdtom.fibroblasts2@assays$SCT)[grep(aIdent, rownames(agg.tdtom.fibroblasts2@assays$SCT))[[1]]])
          clean.genes = c(clean.genes,aIdent)
      }}
    for (agene in clean.genes) {
      label = alistname
      #gene = c(gene, agene)
      #group = c(group, adataset)
      transcripts <- sum(FetchData(agg.tmp, vars = agene, layer = 'counts'))
      index <- nrow(df) + 1
      df[index,] <- c(adataset, label, agene, transcripts, cells, (transcripts/cells))
    }}}
write.csv(df, file = "")