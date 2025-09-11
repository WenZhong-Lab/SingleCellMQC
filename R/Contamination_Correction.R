#' Detect the comtamination causing genes by scCDC package
#'
#' This function identifies potential contamination features in a given Seurat object. It utilizes the
#' `scCDC::ContaminationDetection` method to detect contamination, making it easier to identify unwanted or erroneous features that might affect
#' the downstream analysis.
#'
#' @param object A Seurat object or a list of Seurat objects. This is the primary input that contains single-cell data in the form
#' of a Seurat object. The function will analyze the specified assay within this object.
#' @param assay Character. The name of the assay to be analyzed. Default is "RNA". This assay contains
#' the expression data to be used for contamination detection.
#' @param group.by Character. The metadata column used for grouping cells. Default is "seurat_clusters".
#' This parameter is used to identify groups within the cells that may have different contamination profiles.
#' @param restriction_factor Numeric. A parameter to control the sensitivity of contamination detection.
#' Values should be between 0 and 1, where a higher value implies a stricter detection. Default is 0.1
#' @param ... Additional arguments passed to `scCDC::ContaminationDetection` for more customized analysis.
#' @param split.by A character string specifying the sample column in the Seurat object's metadata (e.g., "orig.ident"). Default is "orig.ident".
#'
#' @return A list of features or metrics that indicate potential contamination, as returned by
#' `scCDC::ContaminationDetection`.
#' @export
#'
#' @references
#' Wang, W., Cen, Y., Lu, Z., Xu, Y., Sun, T., Xiao, Y., Liu, W., Li, J. J., & Wang, C. (2024). scCDC: a computational method for gene-specific contamination detection and correction in single-cell and single-nucleus RNA-seq data. Genome biology, 25(1), 136. https://doi.org/10.1186/s13059-024-03284-wIF: 1
#'
FindContaminationFeature <- function(object, assay= "RNA", group.by="seurat_clusters", split.by="orig.ident", restriction_factor=0.1, ...){
  object_list <- splitObject(object, split.by = split.by)

  GCGname_list <- lapply(object_list, function(x){
    x <- SeuratObject::CreateSeuratObject(counts = x[[assay]] , assay = assay, meta.data = x@meta.data)

    if(class(x[[assay]]) %in% "Assay5"){
      x[[assay]] <- as(x[[assay]], "Assay")
    }

    SeuratObject::Idents(x) <- paste0("g_",x@meta.data[[group.by]])

    GCGs <- scCDC::ContaminationDetection(x, restriction_factor=restriction_factor, ...)
    return(rownames(GCGs))
  })
  names(GCGname_list) <- names(object_list)

  if(length(GCGname_list)==1){
    GCGname_list <- GCGname_list[[1]]
  }

  return(GCGname_list)
}


#' Perform Contamination Correction Using scCDC
#'
#' This function applies the scCDC contamination correction method to a Seurat object or a list of Seurat objects.
#' It supports multi-sample analysis and returns a Seurat object with a new "Corrected" assay.
#'
#' @param object A Seurat object or a list of Seurat objects.
#' @param assay Character. The name of the assay to be analyzed. Default is "RNA".
#' @param features Vector or a list of features (genes) to correct. If NULL, all features are used. Default is NULL.
#' @param split.by Character. The column name in the metadata used to split the object into multiple samples (e.g., "orig.ident"). Default is "orig.ident"
#' @param do.merge Logical. Whether to merge the corrected objects into a single Seurat object. Default is TRUE.
#' @param ... Additional arguments passed to `scCDC::ContaminationCorrection` for customized analysis.
#' @param group.by A character string specifying the metadata column used for grouping cells (default: `"seurat_clusters"`).
#' @param saveBP_dir Directory to save BPCells data. Defaults to `./scCDC_BP/`
#'
#' @return A Seurat object containing the corrected data in a new assay named "Corrected". If `do.merge` is TRUE and multiple samples are provided, a merged Seurat object is returned.
#' @export
RunCorrection_scCDC <- function(object, assay= "RNA", features=NULL, split.by="orig.ident", do.merge=T, group.by="seurat_clusters",
                                saveBP_dir="./scCDC_BP/",
                                ...){

  object_list <- splitObject(object, split.by = split.by)
  tmpdir= "./temp/SingleCellMQC_tempBPCellSplitSeurat/"

  object_name <- names(object_list)

  if(!requireNamespace("scCDC", quietly = TRUE)){
    stop("Please install the scCDC package first.")
  }

  if(packageVersion("scCDC") >= "1.4"){
    if(packageVersion("Seurat") < "5.0.0"){
      stop("scCDC package version >= 1.4 only support for Seurat package version >= 5.0.0, please install the Seurat package version >= 5.0.0 or downgrade the scCDC package version < 1.4.")
    }
    object_list <- lapply(object_name, function(x) {
      count_data = Seurat::GetAssayData(object_list[[x]], assay = assay, slot = "counts")
      index = inherits(count_data, "IterableMatrix")
      if (index) {
        count_data <- as(count_data, "dgCMatrix")

      }
      seuObject <- Seurat::CreateSeuratObject(
        counts = count_data,
        assay = assay,
        meta.data = object_list[[x]]@meta.data
      )
      SeuratObject::Idents(seuObject) <- seuObject@meta.data[[group.by]]
      cont_genes <- if (!is(features, "list")) features else features[[x]]
      seuObject <- scCDC::ContaminationCorrection(seuObject, cont_genes = cont_genes, ...)
      count <- SeuratObject::GetAssayData(seuObject, assay = "Corrected", slot = "counts")
      if (index) {
        count <- ConvertToBPCells(count, BPdir =paste0(saveBP_dir, "/", assay, "_split/", x) )
      }
      return(count)
    })

  }else if(packageVersion("scCDC") == "1.3"){

    object_list <- lapply(object_name, function(x) {
      count_data = Seurat::GetAssayData(object_list[[x]], assay = assay, slot = "counts")
      if (inherits(count_data, "IterableMatrix")) {
        count_data <- as(count_data, "dgCMatrix")
      }
      seuObject <- Seurat::CreateSeuratObject(
        counts = count_data,
        assay = assay,
        meta.data = object_list[[x]]@meta.data
      )
      if (inherits(seuObject[[assay]], "Assay5")) {
        seuObject[[assay]] <- as(seuObject[[assay]], "Assay")
      }
      cont_genes <- if (!is(features, "list")) features else features[[x]]
      seuObject <- scCDC::ContaminationCorrection(seuObject, cont_genes = cont_genes, ...)
      count <- SeuratObject::GetAssayData(seuObject, assay = "Corrected", slot = "counts")
    })
  }else{
    stop("Please install the scCDC package version >= 1.3.")
  }

  names(object_list) <- object_name

  if (dir.exists(tmpdir)) {
    unlink(tmpdir, recursive = TRUE)
  }

  if (length(object_list) == 1) {
    # If input is a single Seurat object, add the corrected assay to the original object
    object_list <- object_list[[1]]

    if (inherits(object, "Seurat")) {
      if(packageVersion("Seurat") < "5.0.0"){
        object[[paste0("scCDC_", assay)]] <- SeuratObject::CreateAssayObject(counts = object_list)
      }else{
        object[[paste0("scCDC_", assay)]] <- SeuratObject::CreateAssay5Object(counts = object_list)
      }
      return(object)
    }
  } else {
    if (do.merge) {
      # Merge the corrected objects
      object_list <- MergeMatrix(object_list, prefix=F)
      if (inherits(object_list, "IterableMatrix")) {
        object_list <- ConvertToBPCells(
          object_list,
          BPdir = paste0(saveBP_dir,"/", assay)
        )
      }

      if (dir.exists( paste0(saveBP_dir, "/", assay, "_split/") ) ) {
        unlink(paste0(saveBP_dir, "/", assay, "_split/"), recursive = TRUE)
      }

      # If the original input is a Seurat object, embed the corrected assay into it
      if (inherits(object, "Seurat")) {
        if(packageVersion("Seurat") < "5.0.0"){
          object[[paste0("scCDC_", assay)]] <- SeuratObject::CreateAssayObject(counts = object_list)
        }else{
          object[[paste0("scCDC_", assay)]] <- SeuratObject::CreateAssay5Object(counts = object_list)
        }
        return(object)
      } else {
        # If the original input is a list, return the merged object
        return(object_list)
      }
    } else {
      # If do.merge is FALSE, return the list of corrected objects
      name=names(object)
      object <- lapply(names(object), function(x){
        if(packageVersion("Seurat") < "5.0.0"){
          object[[x]][[paste0("scCDC_", assay)]] <- SeuratObject::CreateAssayObject(counts = object_list[[x]])
        }else{
          object[[x]][[paste0("scCDC_", assay)]] <- SeuratObject::CreateAssay5Object(counts = object_list[[x]])
        }
        return(object[[x]])
      })
      names(object) <- name
      return(object)
    }
  }

  return(object_list)
}


#' Perform the contamination correction by DecontX
#'
#' @param object A Seurat object.
#' @param split.by A character string specifying the metadata column in the Seurat object to split the data into batches.
#'                 Default is "orig.ident".
#' @param ... Additional arguments passed to `decontX::decontX` for more customized analysis.
#'
#' @return A Seurat object contains additional corrected assay.
#' @export

RunCorrection_DecontX <- function(object, split.by = "orig.ident", ...){
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  if(!requireNamespace("decontX", quietly = TRUE)){
    stop("Please install the decontX package first.")
  }

  counts <- Seurat::GetAssayData(object, assay = "RNA", slot = "count")
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  sce <- decontX::decontX(sce, batch=object@meta.data[[split.by]], ...)
  object[["DecontX_RNA"]] <- Seurat::CreateAssayObject(counts = decontX::decontXcounts(sce))
  return(object)
}
