# Background --------------------------------------------------------------
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



# RunCorrection_scCDC <- function(object, assay= "RNA", features=NULL, split.by="orig.ident", do.merge=T, group.by="seurat_clusters", ...){
#
#   object_list <- splitObject(object, split.by = split.by)
#   object_name <- names(object_list)
#
#   if(!requireNamespace("scCDC", quietly = TRUE)){
#     stop("Please install the scCDC package first.")
#   }
#
#   if(packageVersion("scCDC") >= "1.4"){
#     if(packageVersion("Seurat") < "5.0.0"){
#       stop("scCDC package version >= 1.4 only support for Seurat package version >= 5.0.0, please install the Seurat package version >= 5.0.0 or downgrade the scCDC package version < 1.4.")
#     }
#
#     object_list <- lapply(object_name, function(x) {
#       count_data = Seurat::GetAssayData(object_list[[x]], assay = assay, slot = "counts")
#       if (inherits(count_data, "IterableMatrix")) {
#         count_data <- as(count_data, "dgCMatrix")
#       }
#       seuObject <- Seurat::CreateSeuratObject(
#         counts = count_data,
#         assay = assay,
#         meta.data = object_list[[x]]@meta.data
#       )
#       SeuratObject::Idents(seuObject) <- seuObject@meta.data[[group.by]]
#       cont_genes <- if (!is(features, "list")) features else features[[x]]
#       seuObject <- scCDC::ContaminationCorrection(seuObject, cont_genes = cont_genes, ...)
#       Seurat::CreateSeuratObject(
#         counts = SeuratObject::GetAssayData(seuObject, assay = "Corrected", slot = "counts"),
#         assay = "Corrected"
#       )
#     })
#
#   }else if(packageVersion("scCDC") == "1.3"){
#
#     object_list <- lapply(object_name, function(x) {
#       count_data = Seurat::GetAssayData(object_list[[x]], assay = assay, slot = "counts")
#       if (inherits(count_data, "IterableMatrix")) {
#         count_data <- as(count_data, "dgCMatrix")
#       }
#       seuObject <- Seurat::CreateSeuratObject(
#         counts = count_data,
#         assay = assay,
#         meta.data = object_list[[x]]@meta.data
#       )
#       if (inherits(seuObject[[assay]], "Assay5")) {
#         seuObject[[assay]] <- as(seuObject[[assay]], "Assay")
#       }
#       cont_genes <- if (!is(features, "list")) features else features[[x]]
#       seuObject <- scCDC::ContaminationCorrection(seuObject, cont_genes = cont_genes, ...)
#       Seurat::CreateSeuratObject(
#         counts = SeuratObject::GetAssayData(seuObject, assay = "Corrected", slot = "counts"),
#         assay = "Corrected"
#       )
#     })
#   }else{
#     stop("Please install the scCDC package version >= 1.3.")
#   }
#
#   names(object_list) <- object_name
#
#
#   if (length(object_list) == 1) {
#     # If input is a single Seurat object, add the corrected assay to the original object
#     object_list <- object_list[[1]]
#
#     if (inherits(object, "Seurat")) {
#       object[[paste0("scCDC_", assay)]] <- object_list@assays$Corrected
#       return(object)
#     } else {
#       # If the original input is a list, return the merged object
#       return(object_list)
#     }
#   } else {
#     if (do.merge) {
#       # Merge the corrected objects
#       Corrected_data <- merge(object_list[[1]], y = object_list[2:length(object_list)])
#       if (inherits(object_list[[1]]@assays$Corrected, "Assay5")) {
#         Corrected_data <- SeuratObject::JoinLayers(Corrected_data, assay = "Corrected")
#       }
#
#       # If the original input is a Seurat object, embed the corrected assay into it
#       if (inherits(object, "Seurat")) {
#         object[[paste0("scCDC_", assay)]] <- Corrected_data@assays$Corrected
#         return(object)
#       } else {
#         # If the original input is a list, return the merged object
#         return(Corrected_data)
#       }
#     } else {
#       # If do.merge is FALSE, return the list of corrected objects
#       return(object_list)
#     }
#   }
#
#   return(object_list)
# }
#
# RunCorrection_scCDC <- function(object, assay= "RNA", features=NULL, split.by="orig.ident", do.merge=T, group.by="seurat_clusters", ...){
#
#   object_list <- splitObject(object, split.by = split.by)
#   object_name <- names(object_list)
#
#   if(!requireNamespace("scCDC", quietly = TRUE)){
#     stop("Please install the scCDC package first.")
#   }
#
#   if(packageVersion("scCDC") >= "1.4"){
#     if(packageVersion("Seurat") < "5.0.0"){
#       stop("scCDC package version >= 1.4 only support for Seurat package version >= 5.0.0, please install the Seurat package version >= 5.0.0 or downgrade the scCDC package version < 1.4.")
#     }
#
#
#     if(!is(features, "list")){
#       object_list <- lapply( object_name, function(x){
#         seuObject <- Seurat::CreateSeuratObject(counts = Seurat::GetAssayData(object_list[[x]], assay=assay, slot= "counts") , assay = assay, meta.data = object_list[[x]]@meta.data)
#         SeuratObject::Idents(seuObject) <- seuObject@meta.data[[group.by]]
#         seuObject <- scCDC::ContaminationCorrection(seuObject, cont_genes=features, ...)
#         seuObject <- Seurat::CreateSeuratObject(SeuratObject::GetAssayData(seuObject, assay="Corrected", slot="counts"), assay = "Corrected")
#       })
#     }else{
#       object_list <- lapply( object_name, function(x){
#         seuObject <- Seurat::CreateSeuratObject(counts = object_list[[x]][[assay]] , assay = assay, meta.data = object_list[[x]]@meta.data)
#         SeuratObject::Idents(seuObject) <-  seuObject@meta.data[[group.by]]
#         seuObject <- scCDC::ContaminationCorrection(seuObject, cont_genes=features[[x]], ...)
#         seuObject <- Seurat::CreateSeuratObject(SeuratObject::GetAssayData(seuObject, assay="Corrected", slot="counts"), assay = "Corrected")
#       })
#     }
#
#
#   }else if(packageVersion("scCDC") == "1.3"){
#     if(!is(features, "list")){
#       object_list <- lapply(object_name, function(x){
#         seuObject <- Seurat::CreateSeuratObject(counts = object_list[[x]][[assay]] , assay = assay, meta.data = object_list[[x]]@meta.data)
#         if(class(seuObject[[assay]]) %in% "Assay5"){
#           seuObject[[assay]] <- as(seuObject[[assay]], "Assay")
#         }
#         seuObject <- scCDC::ContaminationCorrection(seuObject, cont_genes=features, ...)
#         seuObject <- Seurat::CreateSeuratObject(SeuratObject::GetAssayData(seuObject, assay="Corrected", slot="counts"), assay = "Corrected")
#       })
#     }else{
#       object_list <- lapply(object_name, function(x){
#         seuObject <- Seurat::CreateSeuratObject(counts = object_list[[x]][[assay]] , assay = assay, meta.data = object_list[[x]]@meta.data)
#         if(class(seuObject[[assay]]) %in% "Assay5"){
#           seuObject[[assay]] <- as(seuObject[[assay]], "Assay")
#         }
#         seuObject <- scCDC::ContaminationCorrection(seuObject, cont_genes=features[[x]], ...)
#         seuObject <- Seurat::CreateSeuratObject(SeuratObject::GetAssayData(seuObject, assay="Corrected", slot="counts"), assay = "Corrected")
#       })
#     }
#
#   }else{
#     stop("Please install the scCDC package version >= 1.3.")
#   }
#
#   names(object_list) <- object_name
#
#
#   if (length(object_list) == 1) {
#     # If input is a single Seurat object, add the corrected assay to the original object
#     object_list <- object_list[[1]]
#
#     if (inherits(object, "Seurat")) {
#       object[[paste0("scCDC_", assay)]] <- object_list@assays$Corrected
#       return(object)
#     } else {
#       # If the original input is a list, return the merged object
#       return(object_list)
#     }
#   } else {
#     if (do.merge) {
#       # Merge the corrected objects
#       Corrected_data <- merge(object_list[[1]], y = object_list[2:length(object_list)])
#       if (inherits(object_list[[1]]@assays$Corrected, "Assay5")) {
#         Corrected_data <- SeuratObject::JoinLayers(Corrected_data, assay = "Corrected")
#       }
#
#       # If the original input is a Seurat object, embed the corrected assay into it
#       if (inherits(object, "Seurat")) {
#         object[[paste0("scCDC_", assay)]] <- Corrected_data@assays$Corrected
#         return(object)
#       } else {
#         # If the original input is a list, return the merged object
#         return(Corrected_data)
#       }
#     } else {
#       # If do.merge is FALSE, return the list of corrected objects
#       return(object_list)
#     }
#   }
#
#   return(object_list)
# }


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

#' Perform the contamination correction by DecontPro
#'
#' @param object A Seurat object.
#' @param group.by The column name in metadata that contains the clusters information.
#' @param ... Additional arguments passed to `decontX::decontPro` for more customized analysis.
#'
#' @return A Seurat object contains additional corrected assay.
#' @export

RunCorrection_DecontPro <- function(object, group.by="seurat_clusters",...){

  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  if(!requireNamespace("decontX", quietly = TRUE)){
    stop("Please install the decontX package first.")
  }

  counts <- Seurat::GetAssayData(object, assay = "ADT", slot = "count")
  clusters <- as.integer(object@meta.data[[group.by]])
  counts <- as.matrix(counts)
  out <- decontX::decontPro(counts, clusters, ...)
  object[["DecontPro_ADT"]] <- Seurat::CreateAssayObject(counts = out$decontaminated_counts )
  return(object)
}

#' A wrapper function for contamination correction
#'
#' This function performs contamination correction on a Seurat object using the specified method. The available methods include `DecontX`, `DecontPro`, and `scCDC`.
#'
#' @param object A Seurat object. The input data must be a valid Seurat object to proceed with correction.
#' @param method Character string specifying the correction method to use. Options include:
#'   \itemize{
#'     \item `"DecontX"`: Runs the `RunCorrection_DecontX()` function.
#'     \item `"DecontPro"`: Runs the `RunCorrection_DecontPro()` function.
#'     \item `"scCDC"`: Runs the `RunCorrection_scCDC()` function.
#'   }
#'   If not specified, the function will return an error indicating that a valid method must be chosen.
#' @param group.by Character string specifying the grouping variable (default is `"seurat_clusters"`). This parameter is mainly used when the method is `"DecontPro"`, allowing you to perform the correction in a grouped manner.
#' @param features Optional vector of feature names. Used specifically by the `scCDC` method to indicate which features to consider during correction.
#' @param ... Additional arguments passed to the respective correction function. These arguments will depend on the selected correction method.
#'
#' @return A Seurat object with corrected data, depending on the chosen method.
#'
#' @export
RunCorrection <- function(object, method=NULL, group.by="seurat_clusters", features=NULL, ... ){
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }
  out <- switch (method,
                 "DecontX" = RunCorrection_DecontX(object, ...),
                 "DecontPro" = RunCorrection_DecontPro(object, group.by=group.by, ...),
                 "scCDC" = RunCorrection_scCDC(object,features=features, ...),
                 stop("Invalid `method`, only `DecontX`, `DecontPro`, `scCDC`. "))
  return(out)

}



# batch feature -----------------------------------------------------------




#' @title Plot variance explained per feature for specified variables
#'
#' @description This function visualizes the percentage of variance explained per feature (e.g., gene) for specified variables in a Seurat object.
#' It supports creating density plots and bar plots to represent the distribution or highlight key features with high variance explained by the selected variables.
#'
#' @param object A Seurat object or a data frame. If a Seurat object is provided, the variance explained for each feature will be computed using the `RunVarExplained` function. If a data frame is provided, it should contain precomputed variance explained values.
#' @param assay Character string specifying which assay to use. Options include `"RNA"` (default) or other assays present in the Seurat object.
#' @param variables Character vector specifying the variables of interest from the metadata (e.g., `"condition"`, `"batch"`). These variables will be used to explain the variance for each feature.
#' @param plot.type Character string specifying the type of plot to generate. Options are `"density"` to visualize the distribution of variance explained or `"bar"` to visualize the top features. Default is `"density"`.
#' @param color.density Optional vector of colors for density plots. If `NULL`, default colors are generated.
#' @param color.bar Character string specifying the color of bars in the bar plot. Default is `"#56B4E9"`.
#' @param ntop Numeric value indicating the number of top features to highlight in the bar plot based on their variance explained. Default is `10`.
#' @param return.type Character string or vector indicating the type of output to return. Options are `"plot"`, `"interactive_table"`, or both. Default is `"plot"`.
#' @param csv.name Character. Name for interactive_table csv file.
#'
#' @return The function returns either a `ggplot` object, an interactive table, or both, depending on the `return.type` specified.
#'   \itemize{
#'     \item If `"plot"` is specified, a density plot or bar plot (or both) showing the variance explained per feature is returned.
#'     \item If `"interactive_table"` is specified, an interactive table summarizing the variance explained per feature is returned.
#'     \item If both are specified, a list containing both the plot and the interactive table is returned.
#'   }
#'
#' @export
PlotVEPerFeature <- function(object, assay="RNA", variables=NULL, plot.type="density", color.density=NULL, color.bar="#56B4E9",ntop=10, return.type="plot",
                             csv.name=paste0(assay, "_variance_explained")){

  if( length(setdiff(plot.type, c("density", "bar")))!=0 ){
    stop("Invalid `plot.type`, only: `density` or/and `bar` ")
  }

  ##

  if (is(object, "Seurat")) {
    matrix_r2 <- RunVarExplained(object, assay="RNA", variables=variables)
    rsquared_mat <- data.table::data.table(feature=rownames(matrix_r2), matrix_r2)
  }
  else {
    rsquared_mat <- data.table::data.table(feature=rownames(object), object)
    variables <- colnames(rsquared_mat)[-1]
  }
  out <- list()
  if("interactive_table"  %in% return.type){
    top_rows <- lapply(variables, function(col) {
      rsquared_mat[order(-get(col))][1:ntop, feature]
    })
    unique_rows <- unique(unlist(top_rows))
    dt <- rsquared_mat[feature %in% unique_rows]
    dt <- as.data.frame(dt)
    out$interactive_table <- .re_table(dt, csv.name =csv.name, elementId = paste0(csv.name, "_table"), maxWidth = NULL,
                                       subtitle = paste0("Variance Explained per Feature (top ",ntop, ")") )
  }

  rsquared_mat <- data.table::melt(rsquared_mat, id.vars="feature")
  if(plot.type=="density"){
    if(is.null(color.density)){
      color.density <- get_colors(length( unique(rsquared_mat$variable)))
    }
  }

  if("plot" %in% return.type){
    plot_out <- switch (plot.type,
                        "density"= ggplot2::ggplot(rsquared_mat, ggplot2::aes(x = value, colour = variable)) +
                          ggplot2::geom_line(stat = "density",alpha = 0.8, linewidth = 1) +
                          ggplot2::geom_vline(xintercept = 1, linetype = 2, linewidth=0.8) +
                          ggplot2::scale_x_log10(breaks = 10^(-3:2), labels = c(0.001,0.01, 0.1, 1, 10, 100)) +
                          ggplot2::labs(subtitle = paste0("Explanatory Variables in ",assay), x="% variance explained", y="Density", color = "Variable") +
                          ggplot2::coord_cartesian(xlim = c(10^(-3),
                                                   100))+
                          ggplot2::scale_color_manual(values =color.density )+
                          ggplot2::theme_bw(base_size = 15)+
                          ggplot2::theme(#text = element_text(face = "bold"),
                            panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),

                            axis.text = ggplot2::element_text(color = "black"),
                            axis.text.x = ggplot2::element_text(angle = 45,hjust = 1,vjust = 1)##
                          ),
                        "bar"= {
                          lapply(variables, function(x){
                            bar_data <- rsquared_mat[rsquared_mat$variable %in% x, ]
                            bar_data <- bar_data[order(-(bar_data$value) ),][1:ntop,]
                            bar_data$feature <- factor(bar_data$feature, levels = bar_data$feature[length(bar_data$feature):1])
                            ggplot2::ggplot(data = bar_data) +
                              ggplot2::geom_col(ggplot2::aes(x=feature, y = .data[["value"]]),alpha=0.7,fill=color.bar)+
                              ggplot2::coord_flip()+
                              ggplot2::theme_bw(base_size = 14)+
                              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                                    panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),

                                    axis.text = ggplot2::element_text(color = "black"),
                                    axis.title.y = ggplot2::element_blank()
                              )+
                              # ggplot2::theme(panel.grid.major=ggplot2::element_blank(),panel.grid.minor=ggplot2::element_blank())+
                              ggplot2::labs(subtitle = x, y = "% variance explained")
                          })
                        },
                        stop("Invalid plot.type, only `heatmap` or `bar`  ")
    )
    if(plot.type=="bar"){
      names(plot_out) <- variables
    }
    out$plot <- plot_out
  }

  if(length(out)==1){
    out <- out[[1]]
  }

  return(out)
}



