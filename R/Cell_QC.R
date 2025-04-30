
.replace_params <- function(known_params, filtered_args) {
  matched_args <- intersect(names(filtered_args), names(known_params))
  for (arg in matched_args) {
    known_params[[arg]] <- filtered_args[[arg]]
  }
  return(known_params)
}

# Low Quality Cell Detection ----------------------------------------------

#' Low-Quality Cells Detection
#'
#' This function handles multiple quality assessment methods to identify low-quality cells in single-cell data.
#'
#' @param object A Seurat object.
#' @param methods A character vector specifying the quality assessment methods to apply (one or more). Available method options are:
#' \itemize{
#'   \item "MAD": Outliers based on the median absolute deviation. Suitable for RNA and ADT.
#'   \item "ddqc": Only RNA. An unsupervised adaptive QC framework to perform flexible and data-driven QC at the level of cell types while retaining critical biological insights and improved power for downstream analysis. (Subramanian, Ayshwarya et al., Genome biology, 2022, PMID: 36575523)
#'   \item "miQC": Only RNA. An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data. (Hippen, Ariel A et al., PLoS computational biology, 2021, PMID: 34428202)
#'   \item "fixed": Fixed cutoffs for RNA and ADT metrics.
#'   }
#'   The default is all methods.
#' @param add.Seurat Logical; if TRUE, adds the low-quality cell detection results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the low-quality cell detection results.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param ... Additional arguments.
#'
#' @references
#' 1.Hippen, Ariel A et al. “miQC: An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data.” PLoS computational biology vol. 17,8 e1009290. 24 Aug. 2021, doi:10.1371/journal.pcbi.1009290
#'
#' 2.Subramanian, Ayshwarya et al. “Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics.” Genome biology vol. 23,1 267. 27 Dec. 2022, doi:10.1186/s13059-022-02820-w
#'
#' @return If `add.Seurat` is TRUE, a Seurat object with appended metadata; otherwise, a data frame with the low-quality cell detection results.
#' @details The 'Fail' value represents low-quality cells.
#'
#' @examples
#' \dontrun{
#' # Assuming 'seuratObj' is a Seurat object prepared for low-quality cell detection:
#' seuratObj <- RunLQ(object = seuratObj, methods = c("MAD", "ddqc", "miQC"), add.Seurat = TRUE)
#' }
#'
#' @seealso \code{\link{RunLQ_MAD}}, \code{\link{RunLQ_miQC}}, \code{\link{RunLQ_ddqc}}
#' @export
RunLQ <- function(object, methods=c("ddqc", "miQC", "MAD"), add.Seurat = T,split.by=NULL, ...) {
  if (!is.character(methods)) {
    stop("methods should be a character vector.")
  }

  valid_methods <- c("ddqc", "miQC", "MAD", "fixed")
  if (!all(methods %in% valid_methods)) {
    stop("Unknown method(s) found. Valid methods are: ", paste(valid_methods, collapse = ", "))
  }

  results <- lapply(methods, function(method) {
    func_name <- paste0("RunLQ_", method)
    func <- match.fun(func_name)
    if(method == "fixed"){
      out <- func(object, add.Seurat=F, ...)
    }else{
      out <- func(object ,add.Seurat=F, split.by=split.by, ...)
    }
    out <- data.table(cell=rownames(out), out)
    return(out)
  }
  )

  if (length(methods) > 1) {
    out <- Reduce(function(x, y){
      merge_data <- merge(x, y, by = "cell", all = TRUE)
    }, results)
  } else {
    out <- results[[1]]
  }
  out <- data.frame(out)
  rownames(out) <- out[, 1]
  out <- out[, -1, drop=F]

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = out)
    return(object)
  }
  return(out)
}


#' Low quality cells detection based on fixed cutoff
#'
#' Identifies low-quality cells based on the fixed cutoff on a Seurat object.
#'
#' @param object A Seurat object to be filtered.
#' @param min.nFeature_RNA (Optional) Minimum number of RNA features per cell.
#' @param min.nCount_RNA (Optional) Minimum RNA count per cell.
#' @param percent.mt (Optional) Maximum percentage of mitochondrial genes per cell.
#' @param percent.rb (Optional) Maximum percentage of ribosomal genes per cell.
#' @param percent.hb (Optional) Maximum percentage of hemoglobin genes per cell.
#' @param min.nFeature_ADT (Optional) Minimum number of ADT features per cell.
#' @param min.nCount_ADT (Optional) Minimum ADT count per cell.
#' @param add.Seurat (Optional) Logical. If TRUE, the final filtering result is added
#'                   to the Seurat object as a new metadata column named `db_fixed`.
#'                   Default is TRUE.
#' @param ... Additional arguments.
#'
#' @return If `add.Seurat` is FALSE, returns a vector of "Pass" or "Fail" for each cell
#'         based on the filtering criteria. If `add.Seurat` is TRUE, returns the modified
#'         Seurat object with the final filtering result added as a new metadata column
#'         named `db_fixed`.
#'
#' @export
RunLQ_fixed <- function(object, min.nFeature_RNA=NULL,  min.nCount_RNA=NULL,
                        percent.mt=NULL, percent.rb=NULL, percent.hb=NULL,
                        min.nFeature_ADT=NULL,  min.nCount_ADT=NULL,
                        add.Seurat = TRUE,...
){
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  # Initialize the output data frame
  out <- data.frame(row.names = colnames(object))

  # Filter based on nFeature_RNA
  if(!is.null(min.nFeature_RNA)){
    out$min.nFeature_RNA <- ifelse(object$nFeature_RNA >= min.nFeature_RNA, "Pass", "Fail")
  }

  # Filter based on nCount_RNA
  if(!is.null(min.nCount_RNA)){
    out$min.nCount_RNA <- ifelse(object$nCount_RNA >= min.nCount_RNA, "Pass", "Fail")
  }


  # Filter based on percent.mt
  if(!is.null(percent.mt)){
    out$percent.mt <- ifelse(object$percent.mt <= percent.mt, "Pass", "Fail")
  }

  # Filter based on percent.rb
  if(!is.null(percent.rb)){
    out$percent.rb <- ifelse(object$percent.rb <= percent.rb, "Pass", "Fail")
  }

  # Filter based on percent.hb
  if(!is.null(percent.hb)){
    out$percent.hb <- ifelse(object$percent.hb <= percent.hb, "Pass", "Fail")
  }

  # Filter based on nFeature_ADT
  if(!is.null(min.nFeature_ADT)){
    out$min.nFeature_ADT <- ifelse(object$nFeature_ADT >= min.nFeature_ADT, "Pass", "Fail")
  }

  # Filter based on nCount_ADT
  if(!is.null(min.nCount_ADT)){
    out$min.nCount_ADT <- ifelse(object$nCount_ADT >= min.nCount_ADT, "Pass", "Fail")
  }

  # Combine all results into a single "Pass" or "Fail" column
  out$lq_fixed <- ifelse(rowSums(out == "Fail") == 0, "Pass", "Fail")

  # If add.Seurat is TRUE, add the lq_fixed result to the Seurat object

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = out[,"lq_fixed", drop = FALSE])
    return(object)
  } else {
    return(out[,"lq_fixed", drop = FALSE])
  }
}





RunOutlier <-
  function(object,
           method = "mad",#sd
           nmads  = 3,
           nsds = 3 ,
           constant = 1.4826,
           log = F,
           group.by = NULL,
           share.median = F,
           share.mean=F,
           share.mad = F,
           share.sd = F,
           metrics = NULL,
           type = c("both", "lower", "higher"),
           ...) {

    if (!method %in% c("mad", "sd")) {
      stop("Invalid method. Please choose either 'mad' or 'sd'.")
    }
    median_dif=mean_dif=NULL
    if (is.null(metrics)) {
      object <- data.frame(metrics = object)
      metrics <- "metrics"
    }
    if (log) {
      object[[metrics]] <- log(object[[metrics]])
    }

    if(is.null(group.by)){
      data <- split(object,as.character( rep("metrics", length(object[[metrics]]))))
    }else{
      data <- split(object,as.character( object[[group.by]]))
    }

    out <- lapply(data, function(x){
      if("mad" %in% method){
        if( !is.null(group.by) &  share.median){
          stat_value <- median(object[[metrics]],na.rm =T)
        }else{
          stat_value <- median(x[[metrics]],na.rm =T)
        }
        if(!is.null(group.by) & share.mad){
          stat_dif <- data.table(object[,c(metrics, group.by),drop=F])[, median_dif:=  get(metrics)- median(get(metrics)),by=list(get(group.by)) ]$median_dif
        }else{
          stat_dif <- x[[metrics]]-stat_value
        }
        cal_value <- median(abs( na.omit(stat_dif) )) * constant

        nds = nmads
      }else{
        if(!is.null(group.by) & share.mean){
          stat_value <- mean(object[[metrics]],na.rm=T)
        }else{
          stat_value <- mean(x[[metrics]],na.rm =T)
        }
        if(!is.null(group.by) & share.sd){
          stat_dif <- data.table(object[,c(metrics, group.by),drop=F])[, mean_dif:=  (get(metrics)- mean(get(metrics)))^2 ,by=list(get(group.by)) ]$mean_diff
        }else{
          stat_dif <- na.omit(x[[metrics]]-stat_value)
        }
        cal_value <- sqrt(sum(stat_dif)/ (length(stat_dif)-1))

        nds = nsds
      }
      cutoff <- nds * cal_value
      cut_low <- stat_value - cutoff
      cut_high <-  stat_value + cutoff


      if(type[1] == "lower"){
        filtered_x <- x[[metrics]] < cut_low
        if(log){
          cut_low <- exp(cut_low)
        }
        cut_high <- Inf
      }else if(type[1] == "higher"){
        filtered_x <- x[[metrics]] > cut_high
        cut_low <- Inf
        if (log) {
          cut_high <- exp(cut_high)
        }
      }else{
        filtered_x <- abs(x[[metrics]] - stat_value) < cutoff
        if (log) {
          cut_low <- exp(cut_low)
          cut_high <- exp(cut_high)
        }
      }

      value <- c("lower"=cut_low, "higher"=cut_high)
      attr(filtered_x, "thresholds") <- value
      names(filtered_x) <- rownames(x)
      return(filtered_x)
    })

    thresholds <- do.call(rbind, lapply(out, function(x){
      attr(x, "thresholds")
    }))
    names(out) <- NULL
    out  <- do.call(c, out)
    attr(out, "thresholds") <- thresholds
    return(out)
  }



.get_named_column <- function(df, col_name) {
  result <- df[[col_name]]
  names(result) <- rownames(df)
  return(result)
}

#' Low quality cells detection with MAD
#'
#' Identifies low-quality cells based on the median absolute deviation (MAD) method on a Seurat object.
#'
#' @param object Seurat object or data frame containing metrics for quality control. RNA metrics should contain "nCount_RNA", "nFeature_RNA" and "percent.mt". If you have ADT data, it should contain "nCount_ADT", "nFeature_ADT" and "percent.isotype". These metrics can all be generated via the \code{\link{CalculateMetrics}} function.
#' @param nmads Numeric. The number of MADs to use for identifying outliers. Default is 3.
#' @param split.by Optional; Column name by which to group the data before calculating statistics.
#' @param share.medians Logical; if TRUE, forces all groups to have the same median. Default is FALSE.
#' @param share.mads Logical; if TRUE, forces all groups to have the same MAD. Default is FALSE.
#' @param add.Seurat Logical; if TRUE, adds the MAD results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the MAD results.
#' @param ... Additional arguments to be passed to \code{\link[scater]{isOutlier}} function from the `scater` package.
#'
#' @return If `add.Seurat` is TRUE, a Seurat object with appended metadata; otherwise, a data frame with outlier results.
#'
#' @export
RunLQ_MAD <- function(object, nmads  = 3, split.by=NULL, add.Seurat=T, share.medians=F, share.mads=F,  ...) {
  if("Seurat" %in% class(object)){
    metadata <- object@meta.data
  }else{
    metadata <- object
  }

  known_params <- names(formals(scater::isOutlier))
  filtered_args <- list(...)
  filtered_args <- setdiff(filtered_args[names(filtered_args) %in% known_params],c("log", "batch", "type", "nmads") )

  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "--------Low-quality method: MAD --------"))
  cat("--------RNA --------")
  if( length(intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"), colnames(metadata)) )!=3 ){
    stop("metadata colnames not obtain `nCount_RNA`, `nFeature_RNA`, `percent.mt`")
  }

  nCount_RNA <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="nCount_RNA"), "batch"=metadata[[split.by]], "log"=T, "nmads"=nmads), filtered_args))
  nFeature_RNA <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="nFeature_RNA"), "batch"=metadata[[split.by]], "log"=T, "nmads"=nmads), filtered_args))
  mt <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="percent.mt"), "batch"=metadata[[split.by]], "log"=F, "nmads"=nmads), filtered_args))

  mad_result <- data.frame(nCount_RNA, nFeature_RNA, mt)
  colnames(mad_result) <- c( paste0(c("lq_RNA_umi_","lq_RNA_gene_","lq_RNA_mt_"), nmads, c("mad")) )

  if( length(intersect(c("nFeature_ADT", "nCount_ADT", "percent.isotype"), colnames(metadata))==2) ){
    cat("--------ADT --------")
    nCount_ADT <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="nCount_ADT"), "batch"=metadata[[split.by]], "log"=T, "nmads"=nmads), filtered_args))
    nFeature_ADT <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="nFeature_ADT"), "batch"=metadata[[split.by]], "log"=T, "nmads"=nmads), filtered_args))
    isotype <- do.call(scater::isOutlier, c(list("metric" = .get_named_column(metadata, col_name="percent.isotype"), "batch"=metadata[[split.by]], "log"=F, "nmads"=nmads), filtered_args))
    mad_result <- cbind(mad_result, nCount_ADT, nFeature_ADT, isotype)
    colnames(mad_result) <- c( paste0(c("lq_RNA_umi_","lq_RNA_gene_","lq_RNA_mt_", "lq_ADT_umi_","lq_ADT_pro_", "lq_ADT_isotype_"), nmads, c("mad")) )
  }

  mad_result <- as.data.frame(lapply(mad_result, function(x) {
    ifelse(x == TRUE, "Fail", "Pass")
  }))
  rownames(mad_result) <- names(nCount_RNA)
  mad_result <- mad_result[match(rownames(metadata), names(nCount_RNA) ), ]

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = mad_result)
    return(object)
  }

  return(mad_result)
}


#' Low quality cells detection with miQC
#'
#' Performs quality control using the miQC method on a Seurat object, an adaptive probabilistic framework for quality control of single-cell RNA-sequencing data.
#'
#' @param object A Seurat object.
#' @param model_type (character) What type of model to generate. A linear mixture model ("linear") based on mitochondrial percentage and library complexity is recommended. B-spline ("spline") and two-degree polynomial ("polynomial") models are also supported. For a simpler model, a one-dimensional gaussian mixture model ("one_dimensional") based on mitochondrial percentage only is available. Default = "linear".
#' @param name.nFeature_RNA The name of the column in metadata representing the number of detected features per cell.
#' @param name.percent.mt The name of the column in metadata representing MT pct per cell.
#' @param add.Seurat Logical; if TRUE, adds the miQC results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the miQC results.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param random.seed Random seed for reproducibility.
#' @param ... Additional arguments to be passed to \code{\link[miQC]{filterCells}} function from the `miQC` package.
#'
#' @return A Seurat object with appended miQC metadata if `add.Seurat` is TRUE; otherwise, a data frame with the miQC results.
#' @details The 'Fail' value represents low-quality cells.
#' @examples
#' \dontrun{
#' # Assuming `seuratObject` is a Seurat object
#' seuratObject <- RunLQ_miQC(seuratObject, add.Seurat = TRUE)
#' }
#' @seealso \code{\link{RunLQ}}
#' @export
#' @references
#' Hippen, Ariel A et al. “miQC: An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data.” PLoS computational biology vol. 17,8 e1009290. 24 Aug. 2021, doi:10.1371/journal.pcbi.1009290
#'
RunLQ_miQC <- function(object, model_type = "linear",  split.by= "orig.ident",
                       name.nFeature_RNA = "nFeature_RNA", name.percent.mt = "percent.mt", add.Seurat = TRUE,
                       random.seed=1 , ...) {
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  if(!requireNamespace("miQC", quietly = TRUE)){
    stop("Error: Please install the miQC package first.")
  }

  known_params <- names(formals(miQC::filterCells))
  filtered_args <- list(...)
  filtered_args <- filtered_args[names(filtered_args) %in% known_params]


  metadata <- object@meta.data[, c(name.nFeature_RNA, name.percent.mt), drop = FALSE]
  colnames(metadata)[1:2] <- c("detected", "subsets_mito_percent")
  object <- Seurat::AddMetaData(object, metadata = metadata)

  if(!is.null(split.by)){
    if(  length(unique(object@meta.data[[split.by]]))>1  ){
      message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Split Seurat Object"))
      split_object <- Seurat::SplitObject(object, split.by = split.by)
    }else{
      split_object <- list(Sample = object)
      names(split_object) <- as.character(unique(object@meta.data[[split.by]]))
    }
  }else{
    split_object <- list(Sample = object)
  }

  out <- lapply( names(split_object), function(x){
    set.seed(random.seed)
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Low-quality detection (miQC): ","Sample ",x))
    sce <- suppressWarnings(Seurat::as.SingleCellExperiment(split_object[[x]], assay = "RNA"))
    model <- miQC::mixtureModel(sce,model_type=model_type)

    out <- do.call(miQC::filterCells, c(list("sce" = sce, "model"=model, "verbose"=F), filtered_args))
    index <- colnames(out)

    lq_miQC <- data.frame( lq_miQC = rep("Fail", dim(split_object[[x]])[2] ) )
    rownames(lq_miQC) <- rownames(split_object[[x]]@meta.data)
    lq_miQC$lq_miQC[ match(index,rownames(lq_miQC))  ] <- "Pass"
    cat(" ", sum( lq_miQC$lq_miQC=="Fail" ), "fail cells \n", "", sum( lq_miQC$lq_miQC=="Pass" ), "pass cells \n" )
    return(lq_miQC)
  }
  )
  names(out) <- NULL
  out <- do.call(rbind, out)

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = out)
    return(object)
  }
  return(out)
}



#' Low quality cells detection with ddqc
#'
#' Performs quality control using the ddqc method on a Seurat object, an unsupervised adaptive QC framework.
#'
#' @param object A Seurat object.
#' @param add.Seurat Logical; if TRUE, adds the ddqc results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the ddqc results.
#' @param ... Additional arguments to be passed to ddqcR::ddqc.metrics function.
#' @param split.by Column name by which to group the data before calculating statistics.
#'
#' @return A Seurat object with appended ddqc metadata if `add.Seurat` is TRUE; otherwise, a data frame with the ddqc results.
#' @details The 'Fail' value represents low-quality cells.
#'
#' @examples
#' \dontrun{
#' # Assuming `seuratObject` is a Seurat object
#' seuratObject <- RunLQ_ddqc(seuratObject, add.Seurat = TRUE)
#' }
#' @seealso \code{\link{RunLQ}}
#' @export
#' @references
#' Subramanian, Ayshwarya et al. “Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics.” Genome biology vol. 23,1 267. 27 Dec. 2022, doi:10.1186/s13059-022-02820-w
#'

RunLQ_ddqc <- function(object,add.Seurat=TRUE,split.by=NULL,  ...) {
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  if(!requireNamespace("ddqcR", quietly = TRUE)){
    stop("Error: Please install the ddqcR package first.")
  }

  if(!is.null(split.by)){
    if(  length(unique(object@meta.data[[split.by]]))>1  ){
      message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Split Seurat Object"))
      split_object <- Seurat::SplitObject(object, split.by = split.by)
    }else{
      split_object <- list(Sample = object)
      names(split_object) <- as.character(unique(object@meta.data[[split.by]]))
    }
  }else{
    split_object <- list(Sample = object)
  }

  known_params <- formals(ddqcR::ddqc.metrics)
  known_params$data <-NULL
  filtered_args <- list(...)
  updated_params <- .replace_params(known_params, filtered_args)

 out <- lapply( names(split_object) , function(x){
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Low-quality detection (ddqc): ","Sample ",x))
    object_sample <- split_object[[x]]
    result <- data.frame(lq_ddqc=rep(FALSE, length(rownames(object_sample@meta.data)) ))
    rownames(result) <- rownames(object_sample@meta.data)
    object_sample <- ddqcR::initialQC(object_sample)
    # out <- !(do.call(ddqcR::ddqc.metrics, c(list("data" = object_sample), filtered_args))[, "passed.qc", drop = F])
    out <- !(ddqcR::ddqc.metrics(data = object_sample, res = known_params$res, threshold = known_params$threshold,
                                 do.counts = known_params$do.counts, do.genes = known_params$do.genes, do.mito = known_params$do.mito,
                                 do.ribo = known_params$do.ribo,n.genes.lower.bound = known_params$n.genes.lower.bound,
                                 percent.mito.upper.bound = known_params$percent.mito.upper.bound, random.state = known_params$random.state
                                 )[, "passed.qc", drop = F])
    out <- as.data.frame(out)
    colnames(out) <- "lq_ddqc"

    result[match(rownames(out),rownames(result)), 1] <- out$lq_ddqc

    result$lq_ddqc <- ifelse(result$lq_ddqc == TRUE , "Fail", "Pass")

    cat(" ", sum( result$lq_ddqc=="Fail" ), "fail cells \n", "", sum( result$lq_ddqc=="Pass" ), "pass cells \n" )

    return(result)

  } )

  names(out) <- NULL
  out <- do.call(rbind, out)

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "LowQuality",metadata = out)
    return(object)
  }
  return(out)
}

# Doublet Detection ---------------------------------------------------------------------

#' @title Doublet Detection in Single-Cell Data
#' @description A wrapper function that applies various doublet detection methods to a single-cell dataset and optionally adds results to Seurat object.
#'
#' @param object A Seurat object.
#' @param methods A character vector specifying which doublet detection methods to apply. Supported methods include: "hybrid", "bcds", "cxds", "scDblFinder", "DoubletFinder", "VDJ", and "ADT".
#' @param add.Seurat  Logical; if TRUE, adds the doublet detection results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the doublet detection results.
#' @param split.by A character string specifying the column in the Seurat object's metadata to use for stratifying the dataset before applying doublet detection (e.g., "orig.ident").
#' @param db.rate The expected doublet rate, if known; otherwise, it will be estimated 0.08*nCell/10000. Not applicable to VDJ and ADT methods.
#' @param feature.RNA1 Character vector specifying RNA features for ADT doublet detection method.
#' @param feature.RNA2 Character vector specifying additional RNA features for ADT doublet detection method.
#' @param feature.ADT1 Character vector specifying ADT features for ADT doublet detection method.
#' @param feature.ADT2 Character vector specifying additional ADT features for ADT doublet detection method.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008, (so 0.016 among 1600 cells).
#' @param ... Additional arguments passed to individual doublet detection methods. \code{\link{RunDbt_hybrid}}; \code{\link{RunDbt_bcds}}; \code{\link{RunDbt_cxds}}; \code{\link{RunDbt_scDblFinder}}; \code{\link{RunDbt_DoubletFinder}}; \code{\link{RunDbt_VDJ}}; \code{\link{RunDbt_ADT}};
#'
#' @return If `add.Seurat` is TRUE, returns the modified Seurat object with doublet detection results added. Otherwise, returns a list or dataframe with the doublet scores and classifications.
#' @seealso \code{\link{RunDbt_hybrid}}; \code{\link{RunDbt_bcds}}; \code{\link{RunDbt_cxds}}; \code{\link{RunDbt_scDblFinder}}; \code{\link{RunDbt_DoubletFinder}}; \code{\link{RunDbt_VDJ}}; \code{\link{RunDbt_ADT}};
#' @export
RunDbt <- function( object, methods=c("hybrid", "bcds", "cxds", "scDblFinder"),
                    add.Seurat=T, split.by ="orig.ident", db.rate=NULL, split.db.rate.1000=0.008, feature.RNA1 = NULL, feature.RNA2 = NULL,
                    feature.ADT1 =NULL, feature.ADT2 =NULL, ...){
  if (!is.character(methods)) {
    stop("methods should be a character vector.")
  }

  valid_methods <- c("hybrid", "bcds", "cxds", "scDblFinder","DoubletFinder", "VDJ", "ADT")
  if (!all(methods %in% valid_methods)) {
    stop("Unknown method(s) found. Valid methods are: ", paste(valid_methods, collapse = ", "))
  }
  object_list <- splitObject(object = object, split.by = split.by)
  results <- lapply(methods, function(method){
    func_name <- paste0("RunDbt_", method)
    func <- match.fun(func_name)
    out <- func(
      object_list,
      split.by = split.by,
      db.rate = db.rate,
      add.Seurat = F,
      ...
    )
    out <- data.table::data.table(cell = rownames(out), out)
    return(out)
  })

  if (length(methods) > 1) {
    out <- Reduce(function(x, y){
      merge_data <- data.table::merge.data.table(x, y, by = "cell", all = TRUE)
    }, results)
  } else {
    out <- results[[1]]
  }

  out <- data.frame(out)
  rownames(out) <- out[, 1]
  out <- out[, -1, drop=F]

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "Doublet",metadata = out)
    return(object)
  }
  return(out)
}


#' Doublet Detection with hybrid
#'
#' Performs doublet detection using a hybrid method in the `scds` package on a Seurat object.
#'
#' @param object A Seurat object.
#' @param db.rate The expected doublet rate, if known; otherwise, it will be estimated 0.08*nCell/10000.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param add.Seurat Logical; if TRUE, adds the hybrid results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the hybrid results.
#' @param ... Additional arguments to be passed to \code{\link[scds]{cxds_bcds_hybrid}} function from the `scds` package.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008, (so 0.016 among 1600 cells).
#' @param BPtmpdir Temporary directory for BPCells matrix processing. Default is "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/".
#' @return If `add.Seurat` is TRUE, a Seurat object with appended metadata; otherwise, a data frame with the hybrid results.
#' @details The results contained hybrid method scores for each cell and classification as 'Pass' or 'Fail'.
#' @export
#' @seealso \code{\link{RunDbt}}
#' @references
#' Bais, Abha S, and Dennis Kostka. “scds: computational annotation of doublets in single-cell RNA sequencing data.” Bioinformatics (Oxford, England) vol. 36,4 (2020): 1150-1158. doi:10.1093/bioinformatics/btz698
RunDbt_hybrid <- function(object, db.rate=NULL, split.db.rate.1000=0.008,
                          split.by="orig.ident",add.Seurat=T, BPtmpdir= "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/", ...){
  if(!requireNamespace("scds", quietly = TRUE)){
    stop("Error: Please install the scds package first.")
  }

  split_object <- splitObject(object, split.by = split.by, assay="RNA", tmpdir= BPtmpdir)

  # split_object <- Seurat::SplitObject(object = object, split.by = split.by)
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Doublet dection (hybrid)"))

  known_params <- names(formals(scds::cxds_bcds_hybrid))
  filtered_args <- list(...)
  filtered_args <- filtered_args[names(filtered_args) %in% known_params]

  gc()

  ## Run parallel processing
  p <- progressr::progressor(along = 1:length(split_object))
  db_hybrid <- smart_lapply(split_object, function(object_temp){

    x <- ConvertToSCE(object_temp, assay = "RNA")
    if( is(x@assays@data@listData[["counts"]], "IterableMatrix" ) ){
      x@assays@data@listData[["counts"]] <- as(object = x@assays@data@listData[["counts"]], Class = "dgCMatrix")
    }

    if(is.null(db.rate)){
      db.rate <- split.db.rate.1000*x@colData@nrows/1000
    }
    x <- do.call(scds::cxds_bcds_hybrid, c(list("sce"=x), filtered_args))
    names(x[["hybrid_score"]]) <- x@colData@rownames
    ntop <- db.rate * length(x[["hybrid_score"]])
    db_qc <- names(sort(x[["hybrid_score"]], decreasing = TRUE)[1:ntop])

    db_hybrid <- data.frame(db_hybrid_score =x[["hybrid_score"]])
    db_hybrid$db_hybrid <- "Pass"
    db_hybrid[db_qc, 2] <- "Fail"

    message( paste0(format(Sys.time(), "%H:%M:%S"), "-- ", getSampleName(object_temp, sample.by = split.by), ": ",
                    sum( db_hybrid$db_hybrid=="Fail" ), " doublets, ", sum( db_hybrid$db_hybrid=="Pass" ), " singlets" ))
    p()
    return(db_hybrid)
  },future.seed = 1)

  names(db_hybrid) <- NULL
  db_hybrid <- do.call(rbind, db_hybrid)

  if (dir.exists(BPtmpdir)) {
    unlink(BPtmpdir, recursive = TRUE)
  }

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "Doublet",metadata = db_hybrid)
    return(object)
  }
  return(db_hybrid)
}


#' @title Doublet Detection with bcds
#'
#' @description Performs doublet detection using a bcds method in the `scds` package on a Seurat object.
#'
#' @param object A Seurat object.
#' @param db.rate The expected doublet rate, if known; otherwise, it will be estimated 0.08*nCell/10000.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param add.Seurat Logical; if TRUE, adds the bcds results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the bcds results.
#' @param ... Additional arguments to be passed to \code{\link[scds]{bcds}} function from the `scds` package.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008, (so 0.016 among 1600 cells).
#' @param BPtmpdir Temporary directory for BPCells matrix processing. Default is "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/".
#'
#' @return If `add.Seurat` is TRUE, a Seurat object with appended metadata; otherwise, a data frame with the bcds results.
#' @details The results contained bcds method scores for each cell and classification as 'Pass' or 'Fail'.
#' @export
#' @seealso \code{\link{RunDbt}}
#' @references
#' Bais, Abha S, and Dennis Kostka. “scds: computational annotation of doublets in single-cell RNA sequencing data.” Bioinformatics (Oxford, England) vol. 36,4 (2020): 1150-1158. doi:10.1093/bioinformatics/btz698

RunDbt_bcds <- function(object, db.rate=NULL, split.db.rate.1000=0.008, split.by="orig.ident", add.Seurat=T,
                        BPtmpdir= "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/", ...){
  if(!requireNamespace("scds", quietly = TRUE)){
    stop("Error: Please install the scds package first.")
  }

  split_object <- splitObject(object, split.by = split.by, assay="RNA", tmpdir= BPtmpdir)
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Doublet dection (bcds)"))

  known_params <- names(formals(scds::bcds))
  filtered_args <- list(...)
  filtered_args <- filtered_args[names(filtered_args) %in% known_params]


  ## Run parallel processing

  p <- progressr::progressor(along = 1:length(split_object))
  db_bcds <- smart_lapply(split_object, function(object_temp){

    x <- ConvertToSCE(object_temp, assay = "RNA")
    if( is(x@assays@data@listData[["counts"]], "IterableMatrix" ) ){
      x@assays@data@listData[["counts"]] <- as(object = x@assays@data@listData[["counts"]], Class = "dgCMatrix")
    }

    if(is.null(db.rate)){
      db.rate <- split.db.rate.1000*x@colData@nrows/1000
    }
    x <- do.call(scds::bcds, c(list("sce"=x), filtered_args))
    names(x[["bcds_score"]]) <- x@colData@rownames
    ntop <- db.rate * length(x[["bcds_score"]])
    db_qc <- names(sort(x[["bcds_score"]], decreasing = TRUE)[1:ntop])

    db_bcds <- data.frame(db_bcds_score =x[["bcds_score"]])

    db_bcds$db_bcds <- "Pass"
    db_bcds[db_qc, 2] <- "Fail"

    message( paste0(format(Sys.time(), "%H:%M:%S"), "-- ", getSampleName(object_temp, sample.by = split.by), ": ", sum( db_bcds$db_bcds=="Fail" ), " doublets, ", sum( db_bcds$db_bcds=="Pass" ), " singlets" ))
    p()
    return(db_bcds)
  }, future.seed=1)
  names(db_bcds) <- NULL
  db_bcds <- do.call(rbind, db_bcds)

  if (dir.exists(BPtmpdir)) {
    unlink(BPtmpdir, recursive = TRUE)
  }

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "Doublet",metadata = db_bcds)
    return(object)
  }

  return(db_bcds)
}


#' @title Doublet Detection with cxds
#' @description Performs doublet detection using a cxds method in the `scds` package on a Seurat object.
#'
#' @param object A Seurat object.
#' @param db.rate The expected doublet rate, if known; otherwise, it will be estimated 0.08*nCell/10000.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param add.Seurat Logical; if TRUE, adds the cxds results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the cxds results.
#' @param ... Additional arguments to be passed to \code{\link[scds]{cxds}} function from the `scds` package.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008, (so 0.016 among 1600 cells).
#' @param BPtmpdir Temporary directory for BPCells matrix processing. Default is "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/".
#' @return If `add.Seurat` is TRUE, a Seurat object with appended metadata; otherwise, a data frame with the cxds results.
#' @details The results contained cxds method scores for each cell and classification as 'Pass' or 'Fail'.
#' @export
#' @seealso \code{\link{RunDbt}}
#' @references
#' Bais, Abha S, and Dennis Kostka. “scds: computational annotation of doublets in single-cell RNA sequencing data.” Bioinformatics (Oxford, England) vol. 36,4 (2020): 1150-1158. doi:10.1093/bioinformatics/btz698

RunDbt_cxds <- function(object, db.rate=NULL,  split.db.rate.1000=0.008, split.by="orig.ident",add.Seurat=T, BPtmpdir= "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/", ...){
  if(!requireNamespace("scds", quietly = TRUE)){
    stop("Error: Please install the scds package first.")
  }

  split_object <- splitObject(object, split.by = split.by, assay="RNA", tmpdir= BPtmpdir)
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Doublet dection (cxds)"))
  known_params <- names(formals(scds::cxds))
  filtered_args <- list(...)
  filtered_args <- filtered_args[names(filtered_args) %in% known_params]

  ## Run parallel processing

  p <- progressr::progressor(along = 1:length(split_object))
  db_cxds <- smart_lapply(split_object, function(object_temp){
    x <- ConvertToSCE(object_temp, assay = "RNA")
    if( is(x@assays@data@listData[["counts"]], "IterableMatrix" ) ){
      x@assays@data@listData[["counts"]] <- as(object = x@assays@data@listData[["counts"]], Class = "dgCMatrix")
    }

    if(is.null(db.rate)){
      db.rate <- split.db.rate.1000*x@colData@nrows/1000
    }
    x <- do.call(scds::cxds, c(list("sce"=x), filtered_args))
    names(x[["cxds_score"]]) <- x@colData@rownames
    ntop <- db.rate * length(x[["cxds_score"]])
    db_qc <- names(sort(x[["cxds_score"]], decreasing = TRUE)[1:ntop])

    db_cxds <- data.frame(db_cxds_score =x[["cxds_score"]])

    db_cxds$db_cxds <- "Pass"
    db_cxds[db_qc, 2] <- "Fail"

    message( paste0(format(Sys.time(), "%H:%M:%S"), "-- ", getSampleName(object_temp, sample.by = split.by), ": ", sum( db_cxds$db_cxds=="Fail" ), " doublets, ", sum( db_cxds$db_cxds=="Pass" ), " singlets" ))
    p()
    return(db_cxds)
  },future.seed=1)
  names(db_cxds) <- NULL
  db_cxds <- do.call(rbind, db_cxds)

  if (dir.exists(BPtmpdir)) {
    unlink(BPtmpdir, recursive = TRUE)
  }

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "Doublet",metadata = db_cxds)
    return(object)
  }
  return(db_cxds)
}


#' Doublet Detection with scDblFinder
#'
#' Performs doublet detection using a scDblFinder method in the `scDblFinder` package on a Seurat object.
#'
#' @param object A Seurat object. The input must be a valid Seurat object containing single-cell RNA-seq data.
#' @param do.topscore Logical; if TRUE, uses a top scoring method to identify doublets based on scDblFinder scores. Default is FALSE.
#' @param db.rate Numeric; the expected doublet rate for the top score method. Default is calculated as 0.08 * number of cells / 10000.
#' @param split.by A character string specifying the column in the Seurat object's metadata to use for stratifying the dataset before applying doublet detection (e.g., "orig.ident"). Default is "orig.ident".
#' @param seed Random seed for reproducibility. Default is 1.
#' @param add.Seurat Logical; if TRUE, adds the doublet detection results as metadata to the Seurat object and stores it in the `SingleCelMQC` slot in `misc`. Otherwise, returns a data.frame of the doublet detection results. Default is FALSE.
#' @param ... Additional arguments to be passed to \code{\link[scDblFinder]{scDblFinder}} function from the `scDblFinder` package.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008, (so 0.016 among 1600 cells).
#' @param BPtmpdir Temporary directory for BPCells matrix processing. Default is "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/".
#' @return If `add.Seurat` is TRUE, returns the modified Seurat object with doublet detection results added. Otherwise, returns a list or data.frame with the doublet scores and classifications.
#' @details The results contained scDblFinder method scores for each cell and classification as 'Pass' or 'Fail'. 'Fail' represents doublets.
#' @export
#' @seealso \code{\link{RunDbt}}
#' @references
#' Germain, Pierre-Luc et al. “Doublet identification in single-cell sequencing data using scDblFinder.” F1000Research vol. 10 979. 28 Sep. 2021, doi:10.12688/f1000research.73600.2
#'
RunDbt_scDblFinder <- function(object, do.topscore=F, db.rate=NULL, split.db.rate.1000=0.008, split.by="orig.ident", seed=1,add.Seurat=T,BPtmpdir= "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/", ...){
  if(!requireNamespace("scDblFinder", quietly = TRUE)){
    stop("Error: Please install the scDblFinder package first.")
  }

  split_object <- splitObject(object, split.by = split.by, assay="RNA", tmpdir= BPtmpdir)
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Doublet dection (scDblFinder)"))
  known_params <- union(setdiff(names(formals(scDblFinder::scDblFinder)), c("sce", "...", "verbose", "dbr") ), names(formals(scDblFinder::getArtificialDoublets)))
  filtered_args <- list(...)
  filtered_args <- filtered_args[names(filtered_args) %in% known_params]

  ## Run parallel processing
  p <- progressr::progressor(along = 1:length(split_object))
  db_scDblFinder <- smart_lapply(split_object, function(object_temp){

    x <- ConvertToSCE(object_temp, assay = "RNA")
    if( is(x@assays@data@listData[["counts"]], "IterableMatrix" ) ){
      x@assays@data@listData[["counts"]] <- as(object = x@assays@data@listData[["counts"]], Class = "dgCMatrix")
    }

    ####
    if( !do.topscore ){
      if(is.null(db.rate) ){
        dbr <- split.db.rate.1000*ncol(x)/1000
      }
      set.seed(seed)
      suppressWarnings(x <- do.call(scDblFinder::scDblFinder, c(list("sce"=x, "dbr"=db.rate, "verbose"=F  ), filtered_args)))
      db_scDblFinder <- data.frame(db_scDblFinder_score=as.numeric(x[["scDblFinder.score"]]),
                                   db_scDblFinder=as.character(x$scDblFinder.class))
      rownames(db_scDblFinder) <- rownames(x@colData)
      db_scDblFinder$db_scDblFinder <- ifelse(db_scDblFinder$db_scDblFinder=="singlet", "Pass", "Fail")


    }else{
      if(is.null(db.rate) ){
        db.rate <- split.db.rate.1000*ncol(x)/1000
      }
      set.seed(seed)
      x <- do.call(scDblFinder::scDblFinder, c(list("sce"=x, "dbr"=db.rate, "verbose"=F  ), filtered_args))
      db_scDblFinder <- data.frame(db_scDblFinder_score=as.numeric(x[["scDblFinder.score"]]),
                                   db_scDblFinder="Pass")
      rownames(db_scDblFinder) <- rownames(x@colData)
      names(x[["scDblFinder.score"]]) <- x@colData@rownames
      ntop <- db.rate * length(x[["scDblFinder.score"]])
      db_qc <- names(sort(x[["scDblFinder.score"]], decreasing = TRUE)[1:ntop])
      db_scDblFinder[db_qc, 2] <- "Fail"
    }

    message( paste0(format(Sys.time(), "%H:%M:%S"), "-- ", getSampleName(object_temp, sample.by = split.by), ": ",
                    sum( db_scDblFinder$db_scDblFinder=="Fail" ), " doublets ", "", sum( db_scDblFinder$db_scDblFinder=="Pass" ), " singlets" ))
    p()
    return(db_scDblFinder)
  }, future.seed=1)
  names(db_scDblFinder) <- NULL
  db_scDblFinder <- do.call(rbind, db_scDblFinder)
  if (dir.exists(BPtmpdir)) {
    unlink(BPtmpdir, recursive = TRUE)
  }

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "Doublet",metadata = db_scDblFinder)
    return(object)
  }
  return(db_scDblFinder)
}





#' VDJ Doublet Detection in single cell TCR and/or BCR data.
#'
#' Performs doublet detection based on VDJ sequencing information.
#' @param ... Additional arguments
#' @inheritParams RunDbt
#' @return If `add.Seurat` is TRUE, returns the modified Seurat object with doublet detection results added. Otherwise, returns a list or data.frame with the doublet classifications.
#' @details
#' The classification as 'Pass' or 'Fail'. 'Fail' represents doublets.
#' 'multichain' and 'ambiguous' in the `chain_pair` column of the \code{\link{CalculateMetricsPerCell}} function are identified as doublets.
#'
#' @examples
#' \dontrun{
#' # Assuming 'seuratObj' is a Seurat object
#' results_VDJ <- RunDbt_VDJ(object = seurat_obj, add.Seurat = TRUE)
#' }
#' @seealso \code{\link{RunDbt}}
#' @export
#'
RunDbt_VDJ <- function(object,add.Seurat=T,...){
  if("Seurat" %in% class(object)){
    metadata <- object@meta.data
  }else{
    metadata <- object
  }

  if(!("chain_pair" %in% colnames(metadata)) ){
    stop("Error: Please run `CalculateMetrics` function first !!")
  }

  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Doublet dection (VDJchain) "))

  db_VDJ <- data.frame(db_VDJ=metadata$chain_pair)
  rownames(db_VDJ) <- rownames(metadata)
  db_VDJ$db_VDJ <- ifelse(metadata$chain_pair %in% c("multichain", "ambiguous"), "Fail",
                          ifelse(!is.na(metadata$chain_pair),"Pass",NA))
  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "Doublet",metadata = db_VDJ)
    return(object)
  }
  return(db_VDJ)
}


#' Doublet Detection with DoubletFinder
#'
#' This function uses the DoubletFinder package to detect doublets within a Seurat object.
#'
#' @param object A Seurat object. The input must be a valid Seurat object containing single-cell RNA-seq data.
#' @param split.by Character. The name of the metadata column used to split the Seurat object. Default is "orig.ident".
#' @param add.Seurat Logical. If TRUE, the function will add the doublet information back to the input Seurat object
#'        as new metadata fields. If FALSE, the function will return a data frame containing the doublet detection results. Default is TRUE.
#' @param db.rate Numeric. The doublet rate score used to estimate the number of expected doublets.
#'        Default is calculated as 0.08 * ncol(object)/10000. If NULL, this default formula will be used.
#' @param ... Additional arguments.
#' @param split.db.rate.1000 Numeric. The expected doublet rate per 1000 cells. Default is 0.008, (so 0.016 among 1600 cells).
#' @param PCs Numeric. The number of principal components. Default is 10.
#' @param BPtmpdir Temporary directory for BPCells matrix processing. Default is "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/".
#' @return If add.Seurat is TRUE, returns the input Seurat object with additional metadata fields for doublet detection.
#'         Otherwise, returns a data frame with doublet scores and classifications for each cell.
#'
#' @export
RunDbt_DoubletFinder <- function(object,split.by="orig.ident",add.Seurat=TRUE, db.rate=NULL,
                                 split.db.rate.1000=0.008, PCs=10,
                                 BPtmpdir= "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/",...){
  if(!requireNamespace("DoubletFinder", quietly = TRUE)){
    stop("Error: Please install the DoubletFinder package first.")
  }

  if ( !(packageVersion("DoubletFinder") > "2.0.5")) {
    stop("Please update the `DoubletFinder` package to version >= 2.0.6")
  }

  split_object <- splitObject(object, split.by = split.by, assay="RNA", tmpdir= BPtmpdir)
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Doublet dection (DoubletFinder)"))


  ## Run parallel processing

  p <- progressr::progressor(along = 1:length(split_object))
  dbt_DoubletFinder <- smart_lapply(split_object, function(object_temp){
    x <- ConvertToSeurat(object_temp)
    ####
    SeuratObject::DefaultAssay(x) <- "RNA"

    if( is(SeuratObject::GetAssayData(x, assay = "RNA", slot = "counts"), "IterableMatrix" ) ){
      x <- SeuratObject::SetAssayData(x, assay = "RNA", slot = "counts",
                        new.data = as(object = SeuratObject::GetAssayData(x, assay = "RNA", slot = "counts"), Class = "dgCMatrix") )
    }

    currentSample <- quiet(RunPipeline(x, preprocess = "rna.pca"))
    # pK Identification (no ground-truth)
    sweep.res.list <- quiet(DoubletFinder::paramSweep(currentSample , PCs = 1:PCs, sct = FALSE))
    sweep.stats <- quiet(DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE))
    # sweep.stats <- .summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- quiet(DoubletFinder::find.pK(sweep.stats))
    # bcmvn <- .find.pK(sweep.stats)
    # best PK
    mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    # Run doublet finder
    if(is.null(db.rate)){
      db.rate <- split.db.rate.1000  * ncol(x)/1000
    }
    nExp_poi <- db.rate * ncol(x)
    out_doubletFinder <- quiet(
          DoubletFinder::doubletFinder(
            currentSample,
            PCs = 1:PCs,
            pN = 0.25,
            pK = mpK,
            nExp = nExp_poi,
            reuse.pANN = NULL,
            sct = FALSE
          ))

    index <- c(grep("pANN_", colnames(out_doubletFinder@meta.data)) ,
               grep("DF.classifications", colnames(out_doubletFinder@meta.data))  )
    dbt_DoubletFinder <- out_doubletFinder@meta.data[,index, drop=F ]
    colnames(dbt_DoubletFinder) <- c("db_DoubletFinder_score", "db_DoubletFinder")
    dbt_DoubletFinder$db_DoubletFinder <- ifelse(dbt_DoubletFinder$db_DoubletFinder == "Singlet", "Pass", "Fail")
    message(paste0(format(Sys.time(), "%H:%M:%S"), "-- ", getSampleName(object_temp, sample.by = split.by), ": ", sum( dbt_DoubletFinder$db_DoubletFinder=="Fail" ), " doublets, ", sum( dbt_DoubletFinder$db_DoubletFinder=="Pass" ), " singlets" ))

    p()
    return(dbt_DoubletFinder)
  }, future.seed = 1)

  names(dbt_DoubletFinder) <- NULL
  dbt_DoubletFinder <- do.call(rbind, dbt_DoubletFinder)

  if (dir.exists(BPtmpdir)) {
    unlink(BPtmpdir, recursive = TRUE)
  }
  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "Doublet",metadata = dbt_DoubletFinder)
    return(object)
  }
  return(dbt_DoubletFinder)
}



#' Doublet Detection (ADT) with mutually exclusive CITE-seq markers in CITE-seq data.
#'
#' This function based on mutually exclusive CITE-seq markers for doublet detection in CITE-seq data.
#'
#' @param object A Seurat object containing assays: ADT or RNA+ADT.
#' @param add.Seurat  Logical; if TRUE, adds the doublet detection results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the doublet detection results.
#' @param split.by A character string specifying the column in the Seurat object's metadata to use for stratifying the dataset before applying doublet detection (e.g., "orig.ident").
#' @param preprocess Character; name of the preprocessing method to use. Default is "adt.umap".
#' @param return.cutoff Logical; if TRUE, returns the cutoff values used for doublet detection. Default is FALSE.
#' @param feature1 A character vector of the names of the mutually exclusive CITE-seq markers to use for doublet detection.
#' @param feature2 A character vector of the names of the mutually exclusive CITE-seq markers to use for doublet detection.
#' @param method.predict Character; the method used to predict doublets. Default is "NB". Options are "NB", "ROC", "Logit", and "LDA".
#' @param resolution Numeric; the resolution parameter used for clustering. Default is 1.5.
#' @param dims Numeric; the dimensions to use for clustering. Default is 1:20.
#' @param do.correct Logical; if TRUE, performs correction for batch effects before doublet detection. Default is FALSE.
#' @param ... Additional arguments.
#' @param BPtmpdir Temporary directory for BPCells matrix processing. Default is "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/".
#' @return  If add.Seurat is TRUE, returns the modified Seurat object with doublet detection results added. Otherwise, returns a data frame with the doublet classifications.
#' @export
#'
RunDbt_ADT <- function(object, split.by="orig.ident",feature1 = NULL, feature2 = NULL, method.predict=c("NB","ROC", "Logit", "LDA"),
                       resolution=1.5, dims=1:20, do.correct=F, preprocess="adt.umap",
                       return.cutoff=F, add.Seurat=T, BPtmpdir= "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/",...){

  if(!requireNamespace("COSG", quietly = TRUE)){
    stop("Error: Please install the COSG package first.")
  }

  if(! "ADT" %in% SeuratObject::Assays(object)){
    stop("Error: No ADT assay in the object!!")
  }

  if(is.null(feature1) | is.null(feature2)){
    stop("Error: Please provide feature1 and feature2!!")
  }

  split_object <- splitObject(object, split.by = split.by, tmpdir= BPtmpdir)
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Doublet dection (ADT)"))

  ## Run parallel processing
  p <- progressr::progressor(along = 1:length(split_object))
  out <- smart_lapply(split_object, function(object_temp){

    suppressWarnings(suppressMessages(seu <- RunPipeline(object_temp, preprocess = preprocess, resolution = resolution, dims = dims)))
    if ( "BPCells" %in% attr(class(Seurat::GetAssayData(seu, assay = "ADT", slot = "counts")), "package") ) {
      expADT <- Seurat::GetAssayData(seu, assay = "ADT", slot = "counts")
      # log_sums <- BPCells::colSums(log1p(expADT ), na.rm = TRUE)
      # scale_factors <- exp(log_sums / nrow(expADT))
      # expADT <- log1p(  BPCells::multiply_cols(expADT ,1/ scale_factors) )
      seu <- SeuratObject::SetAssayData(object = seu,
                                           assay = "ADT",
                                           slot = "counts",
                                           new.data = as(expADT, "dgCMatrix") )
      Seurat::DefaultAssay(seu) <- "ADT"
      seu <- Seurat::NormalizeData(seu, normalization.method = 'CLR', margin = 2, verbose = F)

    }

      Seurat::DefaultAssay(seu) <- "ADT"
      Seurat::VariableFeatures(seu) <- rownames(seu[["ADT"]])


    c_num <- length(unique(seu$seurat_clusters))
    # If correction is requested, perform the correction step
    if (do.correct) {
      message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- scCDC correction for features"))
      seu <- RunCorrection_scCDC(seu, assay = "ADT", features = unique(c(feature1, feature2)) )
      Seurat::DefaultAssay(seu) <- "scCDC_ADT"
      Seurat::VariableFeatures(seu) <- rownames(seu[["scCDC_ADT"]])
      assay.pos <- "scCDC_ADT"
    } else {
      assay.pos <- "ADT"
    }

    seu <- Seurat::NormalizeData(seu, normalization.method = 'CLR', margin = 2, verbose = F)

    # Calculate the cutoff for each feature1
    cutoff1 <- lapply(feature1, function(x) {
      out <- getPosNegCutoff(seu, method = method.predict, assay.pos = assay.pos, assay.predict = "ADT",
                             feature = x, slot.pos = "data", slot.predict = "data", cluster.by = "seurat_clusters")
      cat(">>>>> Cutoff for ", x, " is ", out, "\n")
      return(out)
    })
    cutoff1 <- do.call(c, cutoff1)

    # Calculate the cutoff for each feature2
    cutoff2 <- lapply(feature2, function(x) {
      out <- getPosNegCutoff(seu, method = method.predict, assay.pos = assay.pos, assay.predict = "ADT",
                             feature = x, slot.pos = "data", slot.predict = "data", cluster.by = "seurat_clusters")
      cat(">>>>> Cutoff for ", x, " is ", out, "\n")

      return(out)
    })
    cutoff2 <- do.call(c, cutoff2)

    p()

    if(return.cutoff){
      cutout <- data.frame(Sample= getSampleName(object_temp, sample.by = split.by),cluster_num = c_num,resolution = resolution, feature1=feature1, cutoff1=cutoff1 ,feature2=feature2, cutoff2=cutoff2 )
      return(cutout)
    }

    ADT_data <- Seurat::GetAssayData(seu, assay="ADT", slot = "data")
    name_index <- colnames(ADT_data)
    name_list1 <- mapply(function(x, y){
      if(is.infinite(y)){
        return(NULL)
      }
      name_index[ ADT_data[x,]>y]
    }, x=feature1, y=cutoff1, SIMPLIFY = F)

    name_list2 <- mapply(function(x, y){
      if(is.infinite(y)){
        return(NULL)
      }
      name_index[ ADT_data[x,]>y]
    }, x=feature2, y=cutoff2, SIMPLIFY = F)

    name_list <- mapply(function(x, y){
      intersect(x,y)
    }, x=name_list1, y=name_list2, SIMPLIFY = F)
    name_list <- do.call(c, name_list)
    db_ADT <- data.frame(db_ADT= colnames(seu))
    db_ADT$db_ADT <- ifelse(db_ADT$db_ADT %in% name_list, "Fail", "Pass")
    rownames(db_ADT) <- colnames(seu)

    message(paste0(format(Sys.time(), "%H:%M:%S"), "-- ", getSampleName(object_temp, sample.by = split.by), ": ", sum( db_ADT$db_ADT=="Fail" ), " doublets, ", sum( db_ADT$db_ADT=="Pass" ), " singlets" ))
    return(db_ADT)
  }, future.seed = 1)

  ##
  ##

  names(out) <- NULL
  db_ADT <- do.call(rbind, out)

  if (dir.exists(BPtmpdir)) {
    unlink(BPtmpdir, recursive = TRUE)
  }

  if(return.cutoff){
    return(db_ADT)
  }

  if(add.Seurat){
    object  <- AddSingleCellMQCData(object, listname = "Doublet",metadata = db_ADT)
    return(object)
  }
  return(db_ADT)
}
getPosNegCutoff <- function(object, method = c("NB","ROC", "Logit", "LDA"), assay.pos = "ADT", assay.predict = "ADT", feature = NULL, slot.pos = "data", slot.predict = "data", cluster.by = "seurat_clusters"){
  # Check if the object is a Seurat object
  if (!("Seurat" %in% class(object))) {
    stop("object must be a Seurat object")
  }

  # Set the identity of the object to the specified cluster
  Seurat::Idents(object) <- object@meta.data[[cluster.by]]

  # Get the positive and negative data
  out_data <- getPosNeg(object, assay = assay.pos, feature = feature, slot = slot.pos)

  # Get the expression data for the prediction assay
  expdata <- Seurat::GetAssayData(object, assay = assay.predict, slot = slot.predict)

  # Match the feature to the expression data
  out_data$feature <- expdata[match(feature, rownames(expdata)), ]

  # Determine the cutoff based on the specified method
  method <- match.arg(method)
  cutoff <- switch(method,
                   "ROC" = .roc_cutoff(out_data),
                   "Logit" = .logit_cutoff(out_data),
                   "NB" = .naiveBayes_cutoff(out_data),
                   "LDA" = .lda_cutoff(out_data)
  )

  return(cutoff)
}
# Helper functions for each method
.roc_cutoff <- function(object) {
  roc <- cutoff::roc(object$feature, object$class)
  cutoff <- roc$cutoff
  return(cutoff)
}

.logit_cutoff <- function(object) {
  object$class <- factor(object$class)
  model <- stats::glm(formula = class ~ feature, family = binomial, data = object)
  pro <- stats::predict(model, object, type = "response")
  pos <- object$feature[pro >= 0.5]
  neg <- object$feature[pro < 0.5]
  cutoff <- mean(c(min(pos), max(neg)))
  return(cutoff)
}

.naiveBayes_cutoff <- function(object) {
  object$class <- factor(object$class)
  model <- e1071::naiveBayes(class ~ feature, data = object)
  pro <- stats::predict(model, object)
  pos <- object$feature[pro == "Pos"]
  neg <- object$feature[pro == "Neg"]
  cutoff <- mean(c(min(pos), max(neg)))
  return(cutoff)
}

.lda_cutoff <- function(object) {
  object$class <- factor(object$class)
  model <- MASS::lda(class ~ feature, data = object)
  pro <- stats::predict(model, object)$class
  pos <- object$feature[pro == "Pos"]
  neg <- object$feature[pro == "Neg"]
  cutoff <- mean(c(min(pos), max(neg)))
  return(cutoff)
}

getPosNeg<- function(object, feature =NULL, assay="ADT",slot="data"){
  genexcell <- Seurat::GetAssayData(object = object, assay=assay,  slot = slot)
  index <- match(feature, rownames(genexcell))
  genexcell <- genexcell[c(index,index[1]),,drop=F]
  rownames(genexcell) <- c(feature, "tempfeature")

  colnames(genexcell) <- Seurat::Idents(object)
  count_positive_by_id <- as.character(unique(Seurat::Idents(object)))

  #pct
  pct <- lapply(count_positive_by_id, function(x){
    #print(x)
    data <- genexcell[1,colnames(genexcell)==x]
    sum(data>0)/length(data)
  })
  pct <- do.call(cbind,pct)
  colnames(pct) <- count_positive_by_id

  #mean
  aver_mean <- lapply(count_positive_by_id, function(x){
    data <- genexcell[1,colnames(genexcell)==x,drop=F]
    Matrix::mean(data)
  })
  aver_mean <- do.call(cbind,aver_mean)
  colnames(aver_mean) <- count_positive_by_id

  ##penalty mean
  pen_mean <- aver_mean * pct
  gene_name = rownames(pen_mean)

  pen_sort <- sort(pen_mean[,], decreasing = T)
  name <- names(pen_sort)
  pct_0 <- pct[,name]
  i = 1
  pre_score=(-1)
  score = 0
  while((i < length(name)) && (pre_score < score)){
    if(pct_0[i] < 0.2){
      break
    }
    pre_score <- score
    cluster <- name[1:i]
    score <- suppressWarnings(.cosg_test(genexcell, cluster)$scores)
    i = i+1
  }
  if(pre_score > score){
    i = i-2
    score = pre_score
  }else{
    i = i-1
  }
  index_name <- colnames(object)[object@meta.data$seurat_clusters %in% name[1:i]]
  out_data <- data.frame(feature= genexcell[1,], cellname=colnames(object))
  out_data$class <- ifelse(out_data$cellname %in% index_name, "Pos", "Neg")
  return(out_data)
}
.cosg_test <- function(data, cluster){
  temp <- colnames(data)
  temp[temp %in% cluster] = "aim"
  names(temp) <- seq_along(temp)
  colnames(data) <- seq_along(temp)
  data <- Seurat::CreateSeuratObject(data, min.cells = 0, min.features = 0)
  data <- Seurat::AddMetaData(data, metadata = temp, col.name = "classtype")
  Seurat::Idents(data) <- data$classtype
  suppressWarnings(marker_cosg <- COSG::cosg(
    data,slot="counts",
    mu=1,
    n_genes_user=1,
    remove_lowly_expressed = F))
  marker_cosg1 <- do.call(data.frame, lapply(marker_cosg, function(x){
    x[["aim"]]
  }))
  out<-na.omit(marker_cosg1)
  return(out)
}






# Re-clustering --------------------------------------------------------------------

#' LQ cells re-clustering
#'
#' This function performs re-clustering on low quality cells.
#' Cells marked as 'Fail' in any of the LQ columns are extracted and preprocessed.
#'
#' @param object A Seurat object. The input must be a valid Seurat object.
#' @param method_columns Character. The name of the metadata column used for specifying the method.
#' @param cluster_resolution Numeric. The resolution used for clustering the low-quality cells. Default is 0.5.
#' @param preprocess Character. The name of the preprocessing method to use. Default is "rna.umap".
#'
#' @return A Seurat object containing only the low-quality cells that were re-clustered, with updated cluster assignments.
#'
#' @export
RunLQReClustering <- function(object, method_columns=NULL, cluster_resolution=0.5,preprocess="rna.umap"){
  if( !("Seurat" %in% is(object)) ){
    stop("Error: Input must be Seurat object.")
  }

  if(is.null(method_columns)){
    stop("Error: method_columns is required.")
  }

  metadata <- object@meta.data
  metaname <- colnames(metadata)
  # grep_col <- grep("^lq_(?!.*score$)", metaname, perl = TRUE, value = TRUE)
  meta_table <- metadata[, match( method_columns, metaname), drop=F]


  lq_union <- apply(meta_table, 1, function(row) {
    if(sum(row == "Fail", na.rm = T) >= 1) {
      return("Fail")
    } else {
      return("Pass")
    }
  })

  data_low <- object[, lq_union  %in% 'Fail']
  message(paste("--------Preprocess--------"))
  data_low <- RunPipeline(data_low, preprocess=preprocess, resolution=cluster_resolution)
  data_low$seurat_clusters <- paste0("Fail_",data_low$seurat_clusters)
  return(data_low)
}


#' Find markers by COSG Method
#'
#' This function identifies cluster-specific marker genes using the COSG package. It also calculates log fold changes,
#' the percentage of cells expressing each gene, and other related statistics for each cluster.
#'
#' @param object A Seurat object. The input must be a valid Seurat object.
#' @param cluster.by Character. The metadata column used for clustering the cells. Default is "seurat_clusters".
#' @param top_gene Numeric. The number of top marker genes to identify for each cluster. Default is 50.
#' @param assay Character. The name of the assay to use. Default is "RNA".
#'
#' @return A data frame containing the marker genes, cluster information, COSG scores, log fold changes, mean expressions, and percentages.
#'
#' @export
FindAllMarkerCOSG <- function(object, cluster.by= "seurat_clusters", top_gene = 50,assay="RNA"){

  if( !("Seurat" %in% is(object)) ){
    stop("Error: Input must be Seurat object.")
  }

  if(!requireNamespace("COSG", quietly = TRUE)){
    stop("Error: Please install the COSG package first.")
  }


  Seurat::Idents(object) <- cluster.by
  marker_cosg <- COSG::cosg(object, mu=1, n_genes_user = top_gene, assay=assay)
  unique_gene <- unique(do.call(c,as.vector(marker_cosg$names)))
  unique_data <- SeuratObject::GetAssayData(object, assay=assay,slot = "data")
  clutser_label <- object@meta.data[[cluster.by]]
  unique_label <- unique(clutser_label)
  out <- lapply(unique_label, function(x){
    data <- unique_data[marker_cosg$names[,x], ]
    pct_1 = Matrix::rowSums(data[, clutser_label %in% x]>0) / length(data[1, clutser_label %in% x])
    pct_2 = Matrix::rowSums(data[, !(clutser_label %in% x)]>0) / length(data[1, !(clutser_label %in% x)])
    mean_1 = log(Matrix::rowMeans(expm1(data[, clutser_label %in% x])) +1 )
    mean_2 = log(Matrix::rowMeans(expm1(data[,!(clutser_label %in% x)]))  +1 )
    logFC = mean_1-mean_2
    data.frame(gene=rownames(data), cluster=x , cosg= round(marker_cosg$scores[,x],3) ,
               avg_logFC = round(logFC,3) ,
               meanlog_1 = round(mean_1,4),
               meanlog_2 = round(mean_2,4),
               pct_1 = round(pct_1,4),
               pct_2 = round(pct_2,4)   )
  })
  out <- do.call(rbind, out)
  return(out)
}

#' Plot KEGG enrichment by clusterProfiler package.
#'
#' This function plots KEGG enrichment results by clusterProfiler package for marker genes identified by COSG.
#'
#' @param cosg_marker A data frame containing marker genes, with columns for gene symbols, log fold change, and cluster information.
#'
#' @return A plot representing a dot plot of KEGG pathway enrichment results for each cluster.
#'
#' @export
PlotMarkerEnrichKegg <- function(cosg_marker){
  if(!requireNamespace("clusterProfiler", quietly = TRUE)){
    stop("Error: Please install the clusterProfiler package first.")
  }

  cosg_marker <- cosg_marker[cosg_marker$avg_logFC>0,]
  gid <-clusterProfiler::bitr(unique(cosg_marker$gene), 'SYMBOL',"ENTREZID", OrgDb= 'org.Hs.eg.db')
  cosg_marker$ENTREZID <- gid$ENTREZID[match( cosg_marker$gene, gid$SYMBOL)]
  cosg_marker <- na.omit(cosg_marker)
  gene_list <- split(cosg_marker$ENTREZID, cosg_marker$cluster)
  com <- clusterProfiler::compareCluster(gene = gene_list, fun = "enrichKEGG", pvalueCutoff =0.05, qvalueCutoff =0.9)
  p1 <- clusterProfiler::dotplot(com, label_format=60,color="p.adjust", by="Count") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))+
    ggplot2::ggtitle("KEGG")
  return(p1)
}


# DropletUtils ------------------------------------------------------------

#The NA value represents an empty droplet with count < 100, which can be customized with the lower value
#FDR>0.05 is considered to be the captured droplet
#' Run empty droplet identification
#'
#' This function identifies empty droplets from single-cell RNA-seq data using the DropletUtils package.
#'
#'
#' @param dir_name Character. The path to the directory containing 10X count data.
#' @param lower Numeric. The lower UMI count threshold to consider a droplet as non-empty. Default is 100.
#' @param FDR Numeric. The false discovery rate threshold for identifying retained droplets. Default is 0.01.
#'
#' @return A SingleCellExperiment object containing only the retained droplets that passed the thresholding criteria.
#' @details
#' The NA value represents an empty droplet with count < 100, which can be customized with the `lower` parameter.
#' FDR > 0.05 is considered to be the captured droplet.
#'
#' @export
RunEmptyDroplet<-function(dir_name=NULL,lower=100,FDR=0.01){

  if(!requireNamespace("DropletUtils", quietly = TRUE)){
    stop("Error: Please install the DropletUtils package first.")
  }

  sce <- DropletUtils::read10xCounts(dir_name)
  bcrank1 <- DropletUtils::barcodeRanks(SingleCellExperiment::counts(sce))

  # Only showing unique points for plotting speed.
  uniq1 <- !duplicated(bcrank1$rank)
  plot(bcrank1$rank[uniq1], bcrank1$total[uniq1], log="xy",
       xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
  graphics::abline(h=S4Vectors::metadata(bcrank1)$inflection, col="darkgreen", lty=2)
  graphics::abline(h=S4Vectors::metadata(bcrank1)$knee, col="dodgerblue", lty=2)
  graphics::legend("bottomleft", legend=c("Inflection", "Knee"),
         col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

  set.seed(100)
  e.out <- DropletUtils::emptyDrops(SingleCellExperiment::counts(sce),lower = lower)
  summary_table<-summary(e.out$FDR <= FDR)


  retained <- sce[,which(e.out$FDR <= FDR)]
  bcrank2 <- DropletUtils::barcodeRanks(SingleCellExperiment::counts(retained))
  uniq2 <- !duplicated(bcrank2$rank)
  plot(bcrank2$rank[uniq2], bcrank2$total[uniq2], log="xy",
       xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
  print(summary_table)
  return(retained)
}


# FilterCells --------------------------------------------------------------

#' Filter Cells in a Seurat Object Based on Multiple Criteria
#'
#' This function filters cells in a Seurat object based on multiple criteria columns.
#' Each column in the Seurat object's metadata should contain "Pass" or "Fail" values,
#' and optionally NA values. The function can filter cells based on two logical conditions:
#' 1. "and": Filter cells where all specified columns are "Fail".
#' 2. "or": Filter cells where at least one of the specified columns is "Fail".
#'
#' @param object A Seurat object to be filtered.
#' @param filter_columns A character vector of column names in the Seurat object's metadata
#'                      that contain "Pass" or "Fail" values.
#' @param filter_logic A string specifying the filtering logic. Must be either "and" or "or".  Default is "or".
#'                     - `and`: Filter cells where all specified columns are "Fail".
#'                     - `or`: Filter cells where at least one of the specified columns is "Fail".
#' @param split.by (Optional) A column name in the Seurat object's metadata to split the data by.
#'                The function will report the number of cells before and after filtering for each group.
#'                Default is "orig.ident".
#' @param return.table Logical. If TRUE, returns a data.frame summarizing the number of cells
#'                     before and after filtering, including cases with and without split.by.
#'                     If FALSE, returns the filtered Seurat object. Default is FALSE.
#'
#' @return If `return.table` is TRUE, returns a data.frame summarizing the number of cells
#'         before and after filtering. If `return.table` is FALSE, returns the filtered Seurat object.
#'
#' @export
FilterCells <- function(object, filter_columns = NULL, filter_logic = "or", split.by = "orig.ident",return.table = FALSE ){
  if(!("Seurat" %in% class(object)) ){
    stop("Error: Seurat object must be as input!!")
  }

  if(!all(filter_columns %in% colnames(object@meta.data))){
    stop("Error: The filter_columns must be in the Seurat object's metadata!!")
  }

  if(!filter_logic %in% c("and", "or")){
    stop("Error: The filter_logic must be 'and' or 'or'!!")
  }

  if(!split.by %in% colnames(object@meta.data)){
    stop("Error: The split.by column must be in the Seurat object's metadata!!")
  }


  filter_data <- object@meta.data[, filter_columns, drop = FALSE]
  filter_data[is.na(filter_data)] <- "Pass"

  if(filter_logic == "and"){
    # Filter cells where all specified columns are "Fail"
    out <- !(rowSums(filter_data == "Fail") == length(filter_columns))
  } else {
    # Filter cells where at least one of the specified columns is "Fail"
    out <- !(rowSums(filter_data == "Fail") > 0)
  }

  summary_table <- data.frame(
    Sample = character(),
    Before_Filtering = integer(),
    After_Filtering = integer(),
    stringsAsFactors = FALSE
  )
  # Print filtering stats for the entire object
  cat("Total cells before filtering:", ncol(object), "\n")
  cat("Total cells after filtering:", sum(out), "\n")

  summary_table <- rbind(summary_table, data.frame(
    Sample = "Total",
    Before_Filtering = ncol(object),
    After_Filtering = sum(out),
    stringsAsFactors = FALSE
  ))

  if(!is.null(split.by)){
    split_groups <- unique(object@meta.data[[split.by]])
    for(group in split_groups){
      group_cells <- object@meta.data[[split.by]] == group
      cat(">>>> Sample:", group, "\n")
      cat("  Cells before filtering:", sum(group_cells), "\n")
      cat("  Cells after filtering:", sum(group_cells & out), "\n")
      summary_table <- rbind(summary_table, data.frame(
        Sample = group,
        Before_Filtering = sum(group_cells),
        After_Filtering = sum(group_cells & out),
        stringsAsFactors = FALSE
      ))
    }
  }

  if(return.table){
    return(summary_table)
  } else {
    # Return the filtered Seurat object
    object <- object[, out]
    return(object)
  }
}

