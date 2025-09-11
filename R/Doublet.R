
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
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008.
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

    if("VDJ" %in% valid_methods){
      out = func(object, add.Seurat=F)
    }else{
      out <- func(
        object_list,
        split.by = split.by,
        db.rate = db.rate,
        add.Seurat = F,
        ...
      )
    }
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
#' @param object A single-cell object (e.g., Seurat object).
#' @param db.rate The expected doublet rate, if known; otherwise, it will be estimated 0.08*nCell/10000.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param add.Seurat Logical; if TRUE, adds the hybrid results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the hybrid results.
#' @param ... Additional arguments to be passed to \code{\link[scds]{cxds_bcds_hybrid}} function from the `scds` package.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008.
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

    rm(x)
    gc()
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
#' @param object A single-cell object (e.g., Seurat object).
#' @param db.rate The expected doublet rate, if known; otherwise, it will be estimated 0.08*nCell/10000.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param add.Seurat Logical; if TRUE, adds the bcds results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the bcds results.
#' @param ... Additional arguments to be passed to \code{\link[scds]{bcds}} function from the `scds` package.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008.
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
    rm(x)
    gc()
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
#' @param object A single-cell object (e.g., Seurat object).
#' @param db.rate The expected doublet rate, if known; otherwise, it will be estimated 0.08*nCell/10000.
#' @param split.by Column name by which to group the data before calculating statistics.
#' @param add.Seurat Logical; if TRUE, adds the cxds results as metadata to the Seurat object and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a data.frame of the cxds results.
#' @param ... Additional arguments to be passed to \code{\link[scds]{cxds}} function from the `scds` package.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008.
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
    rm(x)
    gc()
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
#' @param object A single-cell object (e.g., Seurat object).
#' @param do.topscore Logical; if TRUE, uses a top scoring method to identify doublets based on scDblFinder scores.
#' @param db.rate Numeric; the expected doublet rate for the top score method. Default is calculated as 0.08 * number of cells / 10000.
#' @param split.by A character string specifying the column in the Seurat object's metadata to use for stratifying the dataset before applying doublet detection (e.g., "orig.ident"). Default is "orig.ident".
#' @param seed Random seed for reproducibility. Default is 1.
#' @param add.Seurat Logical; if TRUE, adds the doublet detection results as metadata to the Seurat object and stores it in the `SingleCelMQC` slot in `misc`. Otherwise, returns a data.frame of the doublet detection results. Default is FALSE.
#' @param ... Additional arguments to be passed to \code{\link[scDblFinder]{scDblFinder}} function from the `scDblFinder` package.
#' @param split.db.rate.1000 The expected doublet rate per 1000 cells. Default is 0.008.
#' @param BPtmpdir Temporary directory for BPCells matrix processing. Default is "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/".
#' @return If `add.Seurat` is TRUE, returns the modified Seurat object with doublet detection results added. Otherwise, returns a list or data.frame with the doublet scores and classifications.
#' @details The results contained scDblFinder method scores for each cell and classification as 'Pass' or 'Fail'. 'Fail' represents doublets.
#' @export
#' @seealso \code{\link{RunDbt}}
#' @references
#' Germain, Pierre-Luc et al. “Doublet identification in single-cell sequencing data using scDblFinder.” F1000Research vol. 10 979. 28 Sep. 2021, doi:10.12688/f1000research.73600.2
#'
RunDbt_scDblFinder <- function(object, do.topscore=T, db.rate=NULL, split.db.rate.1000=0.008, split.by="orig.ident", seed=1,add.Seurat=T,BPtmpdir= "./temp/SingleCellMQC_BPCellsStepBPToPCAScale/", ...){
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
    rm(x)
    gc()
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
#' @param split.db.rate.1000 Numeric. The expected doublet rate per 1000 cells. Default is 0.008.
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

    rm(currentSample)
    gc()

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



.replace_params <- function(known_params, filtered_args) {
  matched_args <- intersect(names(filtered_args), names(known_params))
  for (arg in matched_args) {
    known_params[[arg]] <- filtered_args[[arg]]
  }
  return(known_params)
}

