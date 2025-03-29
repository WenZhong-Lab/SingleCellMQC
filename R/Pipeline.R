

#' @title Run Single-Cell Analysis Pipeline on Each Group
#'
#' @description
#' This function splits a Seurat object into multiple groups based on a specified metadata column (e.g., orig.ident),
#' and then applies the RunPipeline function to each group independently. This is useful for batch-specific analysis
#' or when you want to preprocess and cluster each group separately.
#'
#' @param object A Seurat object.
#' @param split.by A metadata column name used to split the Seurat object into separate groups (e.g., orig.ident).
#' Default is "orig.ident".
#' @param ... Additional arguments passed to the `RunPipeline` function. See `RunPipeline` documentation for details.
#'
#' @return A list of Seurat objects, where each element corresponds to the result of the `RunPipeline` function
#' applied to a specific group.
#' @export
RunPipelinePerGroup <- function(object, split.by = "orig.ident", ...){
  split_object <- splitObject(object, split.by = split.by)

  split_out <- lapply(split_object, function(x){
    RunPipeline(x, ...)
  })
  names(split_out) <- names(split_object)
  return(split_out)
}


#' @title Initial Clustering of Single-Cell Data
#'
#' @description
#' Default single-cell data preprocessing, including normalization, feature selection, scaling, dimensionality reduction, and clustering.
#'
#' @param object A Seurat object.
#' @param preprocess The method of preprocessing and clustering to use. Options include:
#' \itemize{
#' \item "rna.pca": Perform the analysis to PCA step on RNA assay.
#' \item "rna.umap": Perform the analysis to UMAP step on RNA assay.
#' \item "rna.umap.har": Perform the analysis to UMAP step on RNA assay with Harmony integration.
#' \item "adt.pca": Perform the analysis to PCA step on ADT assay.
#' \item "adt.umap": Perform the analysis to UMAP step on ADT assay.
#' \item "adt.umap.har": Perform the analysis to UMAP step on ADT assay with Harmony integration.
#' \item "wnn.umap": Perform the analysis to UMAP step using weighted nearest neighbor (WNN) method on RNA+ADT assays.
#' \item "wnn.umap.har": Perform the analysis to UMAP step using weighted nearest neighbor (WNN) method on RNA+ADT assays with Harmony integration.
#' }
#' The default setting: "rna.umap".
#' @param resolution The resolution for clustering. Default is 1.
#' @param har.group.by The column name in meta.data to group by for harmony integration. Default is "orig.ident".
#' @param nfeatures Number of features to select as top variable features. Default is 2000.
#' @param exclude_patterns Regular expression patterns to exclude from variable features. Default is:`^IGH[VDJ]*`, `^IGK[JV]*`, `^IGL[JLV]*`, `^IGL[ON].*`, `^TRA[JV]*`, `^TRB[VDJ]*`
#' @param dims The dimensions to use for analysis. Default is 1:20.
#' @param threads.pca Number of threads to use for PCA computation for BPCells. Default is NULL (automatic).
#' @param threads.knn_hnsw Number of threads to use for KNN/HNSW computation for BPCells. Default is 1.
#' @param BPtmpdir Temporary directory for BPCells processing. Default is "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/".
#' @param BP.cutoff The min cutoff for the number of cells to convert to BPCells processing data. Default is 200000.
#' @param BPdir Directory for BPCells data storage. Default is "./BPCellData/".
#' @param wnn.upper.cutoff The max cutoff for the number of cells to perform WNN analysis. Default is 200000.
#' @param wnn.knn.range KNN range for WNN processing. Default is 300.
#' @param wnn.prune.SNN Pruning parameter for WNN SNN graph. Default is 1/20.
#' @param har.seed Seed for harmony integration. Default is 1.
#' @param add.BPscale_data Add BPCells scale data. Default is FALSE.
#'
#' @return The Seurat object with preprocessing and clustering results included, such as PCA scores, UMAP coordinates, and cluster assignments.
#'
#' @details
#' For more details on these methods, please see the \code{\link{Seurat}} and BPCells documentation.
#' If the number of cells is greater than the \code{BP.cutoff}, the function will convert data to BPCells type.
#' If the number of cells is greater than the \code{wnn.upper.cutoff}, the function will perform WNN analysis.
#'
#' @examples
#' \dontrun{
#' ##do rna umap analysis
#' seuratObject <- RunPipeline(seuratObject, preprocess = "rna.umap", resolution = 1)
#'
#' ##do adt pca analysis
#' seuratObject <- RunPipeline(seuratObject, preprocess = "adt.pca")
#'
#' ##do wnn umap analysis
#' seuratObject <- RunPipeline(seuratObject, preprocess = "wnn.umap", resolution = 1)
#'
#' ##do rna harmony umap analysis
#' seuratObject <- RunPipeline(seuratObject, preprocess = "rna.umap.har", resolution = 1)
#'
#' }
#'
#' @export
#' @importFrom dplyr %>%
RunPipeline <- function(object, threads.pca = NULL, threads.knn_hnsw = 1, BPtmpdir = "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/",
                        BP.cutoff=200000,BPdir="./BPCellData/",add.BPscale_data=F,
                    preprocess = "rna.umap", resolution = 1, har.group.by = "orig.ident",
                    dims = 1:20, wnn.upper.cutoff = 200000, wnn.knn.range = 300, wnn.prune.SNN = 1/20,
                    har.seed = 1, nfeatures = 2000, exclude_patterns = c("^IGH[VDJ].*", "^IGK[JV].*",
                                                                         "^IGL[JLV].*", "^IGL[ON].*", "^TRA[JV].*", "^TRB[VDJ].*")) {
  if (!is(object, "Seurat")) {
    stop("Input object must be a Seurat object")
  }

  # check preprocess
  if (!preprocess %in% c("rna.pca", "rna.umap", "rna.umap.har", "adt.pca", "adt.umap", "adt.umap.har", "wnn.umap", "wnn.umap.har")) {
    stop("Invalid preprocess value: must be 'rna.pca', 'rna.umap', 'rna.umap.har', 'adt.pca', 'adt.umap', 'adt.umap.har', 'wnn.umap', or 'wnn.umap.har' ")
  }

  if(preprocess %in% c("rna.pca", "rna.umap", "rna.umap.har", "wnn.umap", "wnn.umap.har") & BP.cutoff < nrow(object)  ) {
    if (!("BPCells" %in% attr(class( SeuratObject::GetAssayData(object, assay ="RNA", slot="counts") ), "package")) ) {
      object <- ConvertToBPCells(object = object, BPdir = BPdir, assay = "RNA", slot = "counts")
    }
  }

  if(preprocess %in% c("adt.pca", "adt.umap", "adt.umap.har", "wnn.umap", "wnn.umap.har") & BP.cutoff < nrow(object)  ) {
    if (!("BPCells" %in% attr(class( SeuratObject::GetAssayData(object, assay ="ADT", slot="counts") ), "package")) ) {
      object <- ConvertToBPCells(object = object, BPdir = BPdir, assay = "ADT", slot = "counts")
    }
  }

  # Helper function for PCA processing
  processPCA <- function(object, assay, nfeatures, exclude_patterns = NULL, threads = NULL, tmpdir = NULL, reduction.name) {
    if ( "BPCells" %in% attr(class(Seurat::GetAssayData(object, assay = assay, slot = "counts")), "package") ) {
      object <- StepBPToPCA(object = object, assay = assay, nfeatures = nfeatures, exclude_patterns = exclude_patterns, threads = threads,
                            tmpdir = tmpdir, output_pca_name = reduction.name, add.scale_data=add.BPscale_data)
    } else {
      object <- .stepSeuratToPCA(object = object, assay = assay, nfeatures = nfeatures, exclude_patterns = exclude_patterns, reduction.name = reduction.name)
    }
    return(object)
  }

  # Helper function for UMAP processing
  processUMAP <- function(object, pca_name, assay, output_umap_name, dims, resolution, output_cluster_name, threads = NULL) {
    if ( "BPCells" %in% attr(class(Seurat::GetAssayData(object, assay = assay, slot = "counts")), "package") ) {
      object <- StepBPPCAToUMAP(object = object, pca_name = pca_name, assay = assay, output_umap_name = output_umap_name, dims = dims)
      object <- StepBPPCAToCluster(object = object, assay = assay, pca_name = pca_name, output_cluster_name = output_cluster_name, threads = threads, resolution = resolution, dims = dims)
    } else {
      object <- .stepSeuratToUMAP(object = object, pca_name = pca_name, assay = assay, output_umap_name = output_umap_name, graph.name = paste0(assay, "_"), dims = dims, resolution = resolution, output_cluster_name = output_cluster_name)
    }
    return(object)
  }

  # Helper function for Harmony processing
  processHarmony <- function(object, group.by.vars, assay.use, reduction.use, reduction.save) {
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Harmony"))
    object <- harmony::RunHarmony(object = object, group.by.vars = group.by.vars, assay.use = assay.use, reduction.use = reduction.use, verbose = TRUE, project.dim = FALSE, reduction.save = reduction.save)
    return(object)
  }

  # Helper function for WNN processing
  processWNN <- function(object, reduction.list, dims.list, modality.weight.name, knn.range, prune.SNN, reduction.name, cluster.name, resolution) {
    object <- Seurat::FindMultiModalNeighbors(object = object, reduction.list = reduction.list, dims.list = dims.list, modality.weight.name = modality.weight.name, knn.range = knn.range, prune.SNN = prune.SNN)
    object <- Seurat::RunUMAP(object = object, nn.name = "weighted.nn", reduction.name = reduction.name, reduction.key = "wnnUMAP_")
    cluster_name <- if (length(cluster.name) != length(resolution)) paste0(cluster.name, "_", resolution) else cluster.name
    object <- Seurat::FindClusters(object = object, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE, cluster.name = cluster_name)
    object$seurat_clusters <- object@meta.data[[cluster_name[length(cluster_name)]]]
    return(object)
  }

  # Main logic based on preprocess value
  switch(preprocess,
         "rna.pca" = {
           object <- processPCA(object = object, assay = "RNA", nfeatures = nfeatures, exclude_patterns = exclude_patterns, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "rna.pca")
         },
         "rna.umap" = {
           Seurat::DefaultAssay(object) <- "RNA"
           object <- processPCA(object = object, assay = "RNA", nfeatures = nfeatures, exclude_patterns = exclude_patterns, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "rna.pca")
           object <- processUMAP(object = object, pca_name = "rna.pca", assay = "RNA", output_umap_name = "rna.umap", dims = dims, resolution = resolution, output_cluster_name = "rna_cluster", threads = threads.knn_hnsw)
         },
         "rna.umap.har" = {
           Seurat::DefaultAssay(object) <- "RNA"
           object <- processPCA(object = object, assay = "RNA", nfeatures = nfeatures, exclude_patterns = exclude_patterns, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "rna.pca")
           object <- processHarmony(object = object, group.by.vars = har.group.by, assay.use = "RNA", reduction.use = "rna.pca", reduction.save = "rna.har")
           object <- processUMAP(object = object, pca_name = "rna.har", assay = "RNA", output_umap_name = "rna.umap.har", dims = dims, resolution = resolution, output_cluster_name = "rna_har_cluster", threads = threads.knn_hnsw)
         },
         "adt.pca" = {
           object <- processPCA(object = object, assay = "ADT", nfeatures = nfeatures, exclude_patterns = NULL, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "adt.pca")
         },
         "adt.umap" = {
           Seurat::DefaultAssay(object) <- "ADT"
           object <- processPCA(object = object, assay = "ADT", nfeatures = nfeatures, exclude_patterns = NULL, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "adt.pca")
           object <- processUMAP(object = object, pca_name = "adt.pca", assay = "ADT", output_umap_name = "adt.umap", dims = dims, resolution = resolution, output_cluster_name = "adt_cluster", threads = threads.knn_hnsw)
         },
         "adt.umap.har" = {
           Seurat::DefaultAssay(object) <- "ADT"
           object <- processPCA(object = object, assay = "ADT", nfeatures = nfeatures, exclude_patterns = NULL, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "adt.pca")
           object <- processHarmony(object = object, group.by.vars = har.group.by, assay.use = "ADT", reduction.use = "adt.pca", reduction.save = "adt.har")
           object <- processUMAP(object = object, pca_name = "adt.har", assay = "ADT", output_umap_name = "adt.umap.har", dims = dims, resolution = resolution, output_cluster_name = "adt_har_cluster", threads = threads.knn_hnsw)
         },
         "wnn.umap" = {
           if(wnn.upper.cutoff < nrow(object)){
             stop("wnn.upper.cutoff must be greater than the number of cells in the object, large datasets may require more memory! ")
           }
           object <- processPCA(object = object, assay = "RNA", nfeatures = nfeatures, exclude_patterns = exclude_patterns, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "rna.pca")
           object <- processPCA(object = object, assay = "ADT", nfeatures = nfeatures, exclude_patterns = NULL, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "adt.pca")
           object <- processWNN(object = object, reduction.list = list("rna.pca", "adt.pca"), dims.list = list(1:20, 1:10), modality.weight.name = c("RNA.weight", "ADT.weight"), knn.range = wnn.knn.range, prune.SNN = wnn.prune.SNN, reduction.name = "wnn.umap", cluster.name = "wnn_clusters", resolution = resolution)
         },
         "wnn.umap.har" = {
           if(wnn.upper.cutoff < nrow(object)){
             stop("wnn.upper.cutoff must be greater than the number of cells in the object, large datasets may require more memory! ")
           }
           object <- processPCA(object = object, assay = "RNA", nfeatures = nfeatures, exclude_patterns = exclude_patterns, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "rna.pca")
           object <- processPCA(object = object, assay = "ADT", nfeatures = nfeatures, exclude_patterns = NULL, threads = threads.pca, tmpdir = BPtmpdir, reduction.name = "adt.pca")
           object <- processHarmony(object = object, group.by.vars = har.group.by, assay.use = "RNA", reduction.use = "rna.pca", reduction.save = "rna.har")
           object <- processHarmony(object = object, group.by.vars = har.group.by, assay.use = "ADT", reduction.use = "adt.pca", reduction.save = "adt.har")
           object <- processWNN(object = object, reduction.list = list("rna.har", "adt.har"), dims.list = list(1:20, 1:10), modality.weight.name = c("RNA.weight", "ADT.weight"), knn.range = wnn.knn.range, prune.SNN = wnn.prune.SNN, reduction.name = "wnn.umap.har", cluster.name = "wnn_har_clusters", resolution = resolution)
         },
         stop("Invalid preprocess value: must be 'rna.pca', 'rna.umap', 'rna.umap.har', 'adt.pca', 'adt.umap', 'adt.umap.har', 'wnn.umap', or 'wnn.umap.har' ")
  )
  return(object)
}





.stepSeuratToUMAP<- function(object, pca_name="rna.pca",assay="RNA", output_umap_name="rna.umap", graph.name="RNA_",
                             dims=1:20,resolution=1,output_cluster_name="rna_cluster", ...){
  Seurat::DefaultAssay(object) <- assay
  object <- Seurat::RunUMAP(object, dims = dims, reduction = pca_name, reduction.name = output_umap_name, reduction.key = paste0(output_umap_name,'_') )

  graph.name <- c(paste0(graph.name, "_nn"), paste0(graph.name, "_snn"))
  object <- Seurat::FindNeighbors(object, dims = dims, reduction = pca_name, graph.name = graph.name )
  if(length(output_cluster_name)!= length(resolution) ){
    cluster_name <- paste0(output_cluster_name,"_", resolution)
  }else{
    cluster_name <- output_cluster_name
  }
  object <- Seurat::FindClusters(object, resolution = resolution, graph.name = graph.name[2], cluster.name = cluster_name )
  object$seurat_clusters <- object@meta.data[[cluster_name[length(cluster_name)] ]]
  return(object)
}

.stepSeuratToPCA <- function(object, assay = "RNA", nfeatures=2000, exclude_patterns=NULL, reduction.name="rna.pca") {
  if (assay == "RNA") {
    Seurat::DefaultAssay(object) <- "RNA"
    object <- Seurat::NormalizeData(object)
    object <- Seurat::FindVariableFeatures(object, nfeatures = nfeatures)
    if(!is.null(exclude_patterns)){
      HVFdata <- SeuratObject::HVFInfo(object = object, assay = assay)
      genes <- rownames(HVFdata)[order(HVFdata[, "variance.standardized"], decreasing = TRUE)]
      del_genes <- setdiff(grep(paste0(exclude_patterns, collapse = "|"), genes, value = TRUE), "TRAV1-2")
      genes <- as.character(setdiff(genes, del_genes))[1:nfeatures]
      Seurat::VariableFeatures(object) <- genes
    }

  } else if (assay == "ADT") {
    Seurat::DefaultAssay(object) <- "ADT"
    Seurat::VariableFeatures(object) <- rownames(object[["ADT"]])
    object <- Seurat::NormalizeData(object, normalization.method = 'CLR', margin = 2)
  } else {
    stop("Invalid assay value: must be 'RNA' or 'ADT'")
  }
  object <- Seurat::ScaleData(object)
  object <- Seurat::RunPCA(object, reduction.name = reduction.name)
  return(object)
}


#' @title Perform preprocessing and PCA using BPCells
#' @description Perform preprocessing and PCA on a Seurat or BP matrix using BPCells.
#'
#' @param object A Seurat object or BPCells matrix.
#' @param assay The Seurat object assay to use (default: "RNA").
#' @param nfeatures Number of features to select (default: 2000).
#' @param exclude_patterns Patterns to exclude from feature selection (default: NULL).
#' @param threads Number of threads to use (default: NULL).
#' @param tmpdir Temporary directory for intermediate files (default: "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/").
#' @param output_pca_name Name of the PCA output (default: "rna.pca").
#' @param add.norm_data Whether to add normalized data to the Seurat object (default: TRUE).
#' @param add.scale_data Whether to add scaled data to the Seurat object (default: FALSE).
#'
#' @return A Seurat object with PCA results or a PCA matrix.
#' @export
StepBPToPCA <- function(object, assay = "RNA", nfeatures=2000, exclude_patterns=NULL, threads=NULL, tmpdir ="./temp/SingleCellMQC_BPCellsStepBPToPCAScale/",
                        output_pca_name="rna.pca", add.norm_data=TRUE, add.scale_data=FALSE) {
  if("Seurat" %in% class(object)){
    mat <- Seurat::GetAssayData(object, assay = assay, slot = "counts")
  }else{
    mat <- object
  }
  if( !("BPCells" %in% attr(class(mat), "package"))){
    mat <- as(mat, "IterableMatrix")
  }
  if (assay == "RNA") {
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Normalization"))
    mat <- BPCells::multiply_cols(mat, 1/Matrix::colSums(mat))
    mat <- log1p(mat * 10000) # Log normalization
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Variance"))
    if(is.null(threads)){
      stats <- BPCells::matrix_stats(mat, row_stats="variance")
    }else{
      stats <- BPCells::matrix_stats(mat, row_stats="variance", threads = threads)
    }

    variable_genes <- order(stats$row_stats["variance",], decreasing=TRUE)
    genes <- colnames(stats$row_stats)[variable_genes]
    if(is.null(exclude_patterns)){
      variable_genes <- variable_genes[1:nfeatures]
    }else{
      del_genes <- setdiff(grep(paste0(exclude_patterns, collapse = "|"), genes, value = TRUE), "TRAV1-2")
      variable_genes <- variable_genes[genes %in% setdiff(genes, del_genes)][1:nfeatures]
    }

    mat_norm <- mat[variable_genes,]
  } else if (assay == "ADT") {
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Normalization"))

    log_sums <- BPCells::colSums(log1p(mat ), na.rm = TRUE)
    scale_factors <- exp(log_sums / nrow(mat))
    mat <- log1p(  BPCells::multiply_cols(mat ,1/ scale_factors) )
    mat_norm <- mat
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Variance"))

    if(is.null(threads)){
      stats <- BPCells::matrix_stats(mat_norm, row_stats="variance")
    }else{
      stats <- BPCells::matrix_stats(mat_norm, row_stats="variance", threads = threads)
    }
    variable_genes <-rownames(mat)
  } else {
    stop("Invalid assay value: must be 'RNA' or 'ADT'")
  }
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Save Temp to ", tmpdir))
  mat_norm <- mat_norm %>% BPCells::write_matrix_dir(tempfile("mat", tmpdir =tmpdir  ))
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Scale"))
  gene_means <- stats$row_stats["mean",variable_genes]
  gene_vars <- stats$row_stats["variance", variable_genes]
  mat_norm <- (mat_norm - gene_means) / gene_vars
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- PCA"))
  if(is.null(threads)){
    k <- min(nrow(mat_norm), 50)
    svd <- BPCells::svds(mat_norm, k=k)
  }else{
    k <- min(nrow(mat_norm), 50)
    svd <- BPCells::svds(mat_norm, k=k, threads=threads)
  }
  pca <- BPCells::multiply_cols(svd$v, svd$d)
  colnames(pca) <- paste0("PC_", 1:ncol(pca))
  rownames(pca) <- colnames(mat)

  if("Seurat" %in% class(object)){
    pca_dim <- Seurat::CreateDimReducObject(embeddings =pca,  key =  paste0(assay, "_PC_"),assay = assay)
    if(add.norm_data){
      object <- SeuratObject::SetAssayData(object = object,
                                               assay = assay,
                                               slot = "data",
                                               new.data = mat)
    }
    if(add.scale_data){
      object <- SeuratObject::SetAssayData(object = object,
                                               assay = assay,
                                               slot = "scale.data",
                                               new.data = mat_norm)
    }else{
      unlink(mat_norm@matrix@dir,recursive = TRUE)
    }
    object[[output_pca_name]] <- pca_dim
    return(object)
  }else{
    return(pca)
  }
}


#' @title Perform clustering on PCA results using BPCells
#' @description Perform clustering on PCA results using BPCells.
#'
#' @param object A Seurat object or PCA matrix.
#' @param assay The Seurat assay to use (default: "RNA").
#' @param pca_name Name of the PCA object (default: "rna.pca").
#' @param output_cluster_name Name of the cluster output (default: "rna_cluster").
#' @param threads Number of threads to use (default: 1).
#' @param resolution Clustering resolution (default: 1).
#' @param min_val Minimum value for SNN graph (default: 0.05).
#' @param dims Dimensions to use for clustering (default: 1:20).
#'
#' @return A Seurat object with clustering results or a cluster data.frame.
#' @export
StepBPPCAToCluster <- function(object, assay="RNA", pca_name="rna.pca", output_cluster_name="rna_cluster", threads = 1,
                               resolution=1 ,min_val=0.05, dims=1:20){
  if("Seurat" %in% class(object)){
    pca <-  Seurat::Embeddings(object, pca_name)
  }else{
    pca <- object
  }
  max_dims <- ncol(pca)
  if(max(dims) > max_dims) {
    dims <- 1:max_dims
    message(sprintf("Using all available PCs (1:%d) as input dimensions", max_dims))
  }
  pca <- pca[,dims]

  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- knn_hnsw & knn_to_snn_graph"))
  knn_result <-  BPCells::knn_hnsw(pca, ef=500, threads=threads) %>% # Find approximate nearest neighbors
    BPCells::knn_to_snn_graph( min_val = min_val) # Find shared nearest neighbors

  cluster_out <- lapply(resolution, function(x){
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- cluster_graph_leiden with resolution ", x))
    temp <- knn_result %>% # Convert to a SNN graph
      BPCells::cluster_graph_leiden(resolution = x) # Perform graph-based clustering
    return(temp)
  })
  cluster_out <- do.call(cbind, cluster_out)

  if(length(output_cluster_name)!= length(resolution) ){
    colnames(cluster_out) <- paste0(output_cluster_name,"_", resolution)
  }else{
    colnames(cluster_out) <- output_cluster_name
  }
  rownames(cluster_out) <- rownames(pca)

  if("Seurat" %in% class(object)){
    object <- Seurat::AddMetaData(object, cluster_out)
    object$seurat_clusters <- cluster_out[,dim(cluster_out)[2], drop=T]
    return(object)
  }else{
    return(data.frame(cluster_out))
  }
}


#' @title  Perform UMAP on PCA results using BPCells
#' @description Perform UMAP on PCA results using BPCells.
#'
#' @param object A Seurat object or PCA result.
#' @param pca_name Name of the PCA object (default: "rna.pca").
#' @param assay The Seurat assay to use (default: "RNA").
#' @param output_umap_name Name of the UMAP output (default: "rna.umap").
#' @param dims Dimensions to use for UMAP (default: 1:20).
#' @param ... Additional arguments passed to uwot::umap.
#'
#' @return A Seurat object with UMAP results or a UMAP matrix.
#' @export
StepBPPCAToUMAP <- function(object, pca_name="rna.pca",assay="RNA", output_umap_name="rna.umap",dims=1:20, ...){
  if("Seurat" %in% class(object)){
    pca <-  Seurat::Embeddings(object, pca_name)
  }else{
    pca <- object
  }
  max_dims <- ncol(pca)
  if(max(dims) > max_dims) {
    dims <- 1:max_dims
    message(sprintf("Using all available PCs (1:%d) as input dimensions", max_dims))
  }
  pca <- pca[,dims]
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- UMAP"))
  set.seed(12341512)
  umap <- uwot::umap(pca,...)
  colnames(umap) <- paste0(assay, c("_UMAP_1", "_UMAP_2"))
  if("Seurat" %in% class(object)){
    umap_dim <- Seurat::CreateDimReducObject(embeddings =umap,  key =  paste0(output_umap_name,"_"),assay = assay)
    object[[output_umap_name]] <- umap_dim
    return(object)
  }else{
    return(umap)
  }
}



