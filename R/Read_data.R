

#' Read V(D)J Sequencing Data From 10X.
#'
#' Reads V(D)J sequencing data for TCR or BCR provided by 10X and combines them into a list.
#'
#' @param dir10x  Directory containing filtered_contig_annotations.csv. e.g. "/data/SingleCellMQC/CellRanger/TP1/vdj_t/filtered_contig_annotations.csv"
#' @param sample  Optional; sample identifiers. If not specified, defaults to sequence along dir10x.
#' @return        A list of contig tables, each list associated with a sample.
#' @export
ReadVDJ <- function(dir10x, sample=seq_along(dir10x)){
  if(is.null(sample)){
    if (is.null(names(dir10x))) {
      sample = paste0("Sample_",  seq_along(dir10x))
    }else{
      sample = names(dir10x)
    }
  }
  names(dir10x) <- .standardizeNames(sample)
  sample=names(dir10x)

  contig_list <- mapply(function(x, y){
    out <- data.table::fread(y, data.table =F)
    out$id <- paste0(x,"_",out$barcode)
    return(out)
  }, x= sample, y=dir10x, SIMPLIFY = F)
  contig_list <- structure(contig_list, class = c(class(contig_list), "VDJ"))
  return(contig_list)
}



#' Load multi-omics data from 10x
#'
#' Read gene/protein expression (GEX), T-cell receptor (TCR), and B-cell receptor (BCR) data provided by 10X genomics.
#'
#' @param dir_GEX  Directory for Gene/Protein Expression data provided by 10X. Directory containing the matrix.mtx, features.tsv, and barcodes.tsv files provided by 10X. e.g. "/data/SingleCellMQC/CellRanger/TP1/sample_filtered_feature_bc_matrix/".
#' A vector or named vector can be given in order to load several data directories.
#' @param dir_TCR  Directory for TCR data provided by 10X. Directory containing filtered_contig_annotations.csv file provided by 10X. A vector or named vector can be given in order to load several data directories.
#' @param dir_BCR  Directory for BCR data provided by 10X. Directory containing filtered_contig_annotations.csv file provided by 10X. A vector or named vector can be given in order to load several data directories.
#' @param dir_BPCells Directory to save BPCells data. Used if is not NULL. Defaults : NULL.
#' @param gene.column Integer; the column of the features.tsv file that contains gene names. Defaults to `2`.
#' @param sample Optional; sample identifiers. If not provided, default names will be generated.
#' @param return.type Character; specifies the format of the returned GEX data. Options are `"Seurat"` (returns a Seurat object) or `"matrix"` (returns a dgCMatrix or IterableMatrix ).
#'   Defaults to `"Seurat"`.
#'
#' @return A list containing the following components:
#'   - `GEX`: Gene/Protein Expression data. If `return.type` is `"Seurat"`, this is a merged Seurat object; if `"matrix"`, it is a merged matrix.
#'   - `TCR`: T-cell receptor data, returned as a list of data frames.
#'   - `BCR`: B-cell receptor data, returned as a list of data frames.
#'
#'   If only one component is loaded (e.g., only GEX data), the function returns that component directly instead of a list.
#' @export
Read10XData <- function(dir_GEX=NULL, dir_TCR=NULL, dir_BCR=NULL, sample=NULL,
                         dir_BPCells=NULL,  gene.column=2, return.type="Seurat"){
  out_list <- list()
  saveBPCells = !is.null(dir_BPCells)
  #exp
  if(length(dir_GEX) != 0){
    if ("Seurat" %in% return.type & saveBPCells) {
      if (packageVersion("Seurat") < "5.0.0") {
        stop("BPCells in Seurat is only supported for versions >= 5.0.0.\n",
             "You can set `return.type` to 'matrix' to output as a matrix instead of a Seurat object.")
      }
    }
    if(!is.null(sample)){
      names(dir_GEX) <- sample
    }else{
      if( is.null( names(dir_GEX) )){
        names(dir_GEX) <- paste0("Sample",seq_along( dir_GEX) )
      }
    }
    names(dir_GEX) <- .standardizeNames(names(dir_GEX))
    sample =names(dir_GEX)

    # Parallelize the processing of each matrix
    message("Reading GEX data...")
    p <- progressr::progressor(along = seq_along(names(dir_GEX)))
    matrix_list <- smart_lapply(names(dir_GEX), function(x) {
      suppressMessages(GEX <- Seurat::Read10X(dir_GEX[x], gene.column=gene.column))
      if ("dgCMatrix" %in% class(GEX)) {
        GEX <- list("Gene Expression" = GEX)
      }
      if(saveBPCells){
        # check BPCells
        if(!requireNamespace("BPCells", quietly = TRUE)){
          stop("BPCells is not installed. Please install BPCells to use BPCells.")
        }
        if( !is.null(GEX[["Gene Expression"]]) ){
          GEX[["Gene Expression"]] <- ConvertToBPCells(GEX[["Gene Expression"]], BPdir = paste0(dir_BPCells, "/", x, "_RNA"))
        }
        if( !is.null(GEX[["Antibody Capture"]]) ){
          GEX[["Antibody Capture"]] <- ConvertToBPCells(GEX[["Antibody Capture"]], BPdir = paste0(dir_BPCells, "/", x, "_ADT"))
        }
      }
      p()
      return(GEX)
    })

    # check type
    types <- lapply(matrix_list, function(x){
      names(x)
    })
    types <- unique(do.call(c, types ))

    out <- list()

    # merge matrix
    if("Gene Expression" %in% types){
      out[["Gene Expression"]] <- Filter(function(x) length(x) > 0, lapply(matrix_list, function(x) x[["Gene Expression"]]))
      out[["Gene Expression"]] <- MergeMatrix(out[["Gene Expression"]])
    }

    if("Antibody Capture" %in% types){
      out[["Antibody Capture"]] <- Filter(function(x) length(x) > 0, lapply(matrix_list, function(x) x[["Antibody Capture"]]))
      out[["Antibody Capture"]] <- MergeMatrix(out[["Antibody Capture"]])
    }
    rm(matrix_list)

    # BPCells save
    if("Gene Expression" %in% types & saveBPCells){
      out[["Gene Expression"]] <- ConvertToBPCells(
        out[["Gene Expression"]],
        BPdir = paste0(dir_BPCells,"/", "RNA")
      )
      unlink(paste0(dir_BPCells, "/", names(dir_GEX), "_RNA") , recursive = TRUE)
    }
    if("Antibody Capture" %in% types & saveBPCells){
      out[["Antibody Capture"]] <- ConvertToBPCells(
        out[["Antibody Capture"]],
        BPdir = paste0(dir_BPCells,"/", "ADT")
      )
      unlink(paste0(dir_BPCells, "/", names(dir_GEX), "_ADT") , recursive = TRUE)
    }

    if("matrix" %in% return.type){
      out_list$GEX <- out
    }else if("Seurat" %in% return.type){
      out_list$GEX <- toSeuratObject(out)
    }
  }

  #tcr
  if(length(dir_TCR) != 0){
    out_list$TCR <- ReadVDJ(dir10x=dir_TCR, sample=sample)
  }

  #bcr
  if(length(dir_BCR) != 0){
    out_list$BCR <- ReadVDJ(dir_BCR, sample=sample)
  }
  if(length(out_list)==1){
    out_list <- out_list[[1]]
  }
  return(out_list)
}


#' Load multi-omics h5 data from 10x
#'
#' @param dir_GEX Directory for Gene/Protein Expression data provided by 10X. Directory containing h5 files provided by 10X. e.g. "/data/SingleCellMQC/CellRanger/TP1/sample_filtered_feature_bc_matrix.h5".
#' A vector or named vector can be given in order to load several data directories.
#' @param dir_TCR Directory for TCR data provided by 10X. Directory containing filtered_contig_annotations.csv file provided by 10X. A vector or named vector can be given in order to load several data directories.
#' @param dir_BCR Directory for BCR data provided by 10X. Directory containing filtered_contig_annotations.csv file provided by 10X. A vector or named vector can be given in order to load several data directories.
#' @param sample Optional; sample identifiers.
#' @param dir_BPCells Directory to save BPCells data. Used if is not NULL. Defaults : NULL.
#' @param gene.column Integer; the column of the features.tsv file that contains gene names. Defaults to `2`.
#' @param sample Optional; sample identifiers. If not provided, default names will be generated.
#' @param return.type Character; specifies the format of the returned GEX data. Options are `"Seurat"` (returns a Seurat object) or `"matrix"` (returns a dgCMatrix or IterableMatrix ).
#'   Defaults to `"Seurat"`.
#'
#' @return A list containing the following components:
#'   - `GEX`: Gene/Protein Expression data. If `return.type` is `"Seurat"`, this is a merged Seurat object; if `"matrix"`, it is a merged matrix.
#'   - `TCR`: T-cell receptor data, returned as a list of data frames.
#'   - `BCR`: B-cell receptor data, returned as a list of data frames.
#'
#'   If only one component is loaded (e.g., only GEX data), the function returns that component directly instead of a list.
#' @export
#'
Read10XH5Data <- function(dir_GEX = NULL, dir_TCR = NULL, dir_BCR = NULL, sample = NULL,
                           dir_BPCells = NULL, gene.column = 2,
                          return.type = "Seurat") {
  out_list <- list()
  saveBPCells = !is.null(dir_BPCells)

  # GEX data processing
  if (length(dir_GEX) != 0) {
    # Check Seurat version for BPCells
    if ("Seurat" %in% return.type & saveBPCells) {
      if (packageVersion("Seurat") < "5.0.0") {
        stop("BPCells in Seurat is only supported for versions >= 5.0.0.\n",
             "You can set `return.type` to 'matrix' to output as a matrix instead of a Seurat object.")
      }
    }

    # Assign sample names
    if (!is.null(sample)) {
      names(dir_GEX) <- sample
    } else if (is.null(names(dir_GEX))) {
      names(dir_GEX) <- paste0("Sample", seq_along(dir_GEX))
    }

    names(dir_GEX) <- .standardizeNames(names(dir_GEX))
    sample = names(dir_GEX)

    # Parallelize the processing of each matrix
    message("Reading GEX data...")
    p <- progressr::progressor(along = seq_along(names(dir_GEX)))
    matrix_list <- smart_lapply(names(dir_GEX), function(x) {
      suppressMessages(GEX <- Seurat::Read10X_h5(dir_GEX[x]))
      if ("dgCMatrix" %in% class(GEX)) {
        colnames(GEX) <- paste0(x, "_", colnames(GEX))
        GEX <- list("Gene Expression" = GEX)
      }

      # Save as BPCells if enabled
      if (saveBPCells) {
        if (!requireNamespace("BPCells", quietly = TRUE)) {
          stop("BPCells is not installed. Please install BPCells to use BPCells.")
        }
        if (!is.null(GEX[["Gene Expression"]])) {
          colnames(GEX[["Gene Expression"]]) <- paste0(x, "_", colnames(GEX[["Gene Expression"]]))
          GEX[["Gene Expression"]] <- ConvertToBPCells(
            GEX[["Gene Expression"]],
            BPdir = paste0(dir_BPCells, "/", x, "_RNA")
          )
        }
        if (!is.null(GEX[["Antibody Capture"]])) {
          colnames( GEX[["Antibody Capture"]]) <- paste0(x, "_", colnames( GEX[["Antibody Capture"]]))

          GEX[["Antibody Capture"]] <- ConvertToBPCells(
            GEX[["Antibody Capture"]],
            BPdir = paste0(dir_BPCells, "/", x, "_ADT")
          )
        }
      }
      p()
      return(GEX)
    })

    # Merge matrices
    types <- unique(do.call(c, lapply(matrix_list, names)))
    out <- list()

    if ("Gene Expression" %in% types) {
      out[["Gene Expression"]] <- Filter(function(x) length(x) > 0,
                                         lapply(matrix_list, function(x) x[["Gene Expression"]]))
      out[["Gene Expression"]] <- MergeMatrix(out[["Gene Expression"]],
                                              add_sample_ids = names(out[["Gene Expression"]]))
    }
    if ("Antibody Capture" %in% types) {
      out[["Antibody Capture"]] <- Filter(function(x) length(x) > 0,
                                          lapply(matrix_list, function(x) x[["Antibody Capture"]]))
      out[["Antibody Capture"]] <- MergeMatrix(out[["Antibody Capture"]],
                                               add_sample_ids = names(out[["Antibody Capture"]]))
    }
    rm(matrix_list)

    # Save merged BPCells
    if ("Gene Expression" %in% types & saveBPCells) {
      out[["Gene Expression"]] <- ConvertToBPCells(
        out[["Gene Expression"]],
        BPdir = paste0(dir_BPCells, "/RNA")
      )
      unlink(paste0(dir_BPCells, "/", names(dir_GEX), "_RNA"), recursive = TRUE)
    }
    if ("Antibody Capture" %in% types & saveBPCells) {
      out[["Antibody Capture"]] <- ConvertToBPCells(
        out[["Antibody Capture"]],
        BPdir = paste0(dir_BPCells, "/ADT")
      )
      unlink(paste0(dir_BPCells, "/", names(dir_GEX), "_ADT"), recursive = TRUE)
    }

    # Return type handling
    if ("matrix" %in% return.type) {
      out_list$GEX <- out
    } else if ("Seurat" %in% return.type) {
      out_list$GEX <- toSeuratObject(out)
    }
  }

  # TCR data processing
  if (length(dir_TCR) != 0) {
    out_list$TCR <- ReadVDJ(dir10x = dir_TCR, sample = sample)
  }

  # BCR data processing
  if (length(dir_BCR) != 0) {
    out_list$BCR <- ReadVDJ(dir_BCR, sample = sample)
  }

  # Simplify output if only one element
  if (length(out_list) == 1) {
    out_list <- out_list[[1]]
  }
  return(out_list)
}



#' Read 10X QC Metrics
#'
#' This function reads 10X Genomics QC metrics files and processes them based on the data type.
#' It supports multiple file types, including multi-library metrics, VDJ metrics, and count metrics.
#'
#' @param file_path A named vector or list of file paths to the QC metrics files.
#' @param sample_name An optional vector of sample names corresponding to the files in `file_path`.
#' If `file_path` is unnamed and `sample_name` is not provided, sample names will be generated automatically as "sample_1", "sample_2", etc.
#'
#' @return A data frame combining the processed QC metrics for all samples.
#' @export
#'
Read10XMetrics <- function(file_path, sample_name=NULL){
  if( is.null(names(file_path)) && is.null(sample_name) ){
    sample_name <- paste0( "sample_", seq_along(file_path) )
  }else{
    if( is.null(sample_name)){
      sample_name <- names(file_path)
    }
  }
  names(file_path) <- .standardizeNames(sample_name)
  sample_name = names(file_path)

  seq_list <- lapply( sample_name, function(x){
    seq <- data.table::fread( file_path[x])
    if( sum(grepl("Library Type", colnames(seq))) > 0 ){
      seq_table <- .metricsMulti(seq, x)
    }else if( sum(grepl("V-J", colnames(seq))) > 0 ){
      seq_table <- .metricsVDJ(seq, x)
    }else{
      seq_table <- .metricsCount(seq, x)
    }
  })

  out <- data.table::rbindlist(seq_list, fill = TRUE)
  out <- as.data.frame(out)
  return(out)
}
.metricsMulti <- function(object, sample_name=NULL){
  dt <- object[, list(`Metric Name`, `Metric Value`, `Library Type`)]
  dt <- unique(dt)
  dt$sample <- sample_name
  dt[grepl("%", `Metric Value`), `Metric Name` := paste0(`Metric Name`, " (%)")]
  dt[, `Metric Value` := gsub(",", "", `Metric Value`)]
  dt[, `Metric Value` := gsub("%", "", `Metric Value`)]
  dt[, `Metric Value` := as.numeric(`Metric Value`)]
  dt[, new_name := fifelse(`Library Type` == "Antibody Capture", "Antibody: ",
                           fifelse(`Library Type` == "Gene Expression", "RNA: ",
                                   fifelse(`Library Type` == "VDJ B", "VDJ B: ",
                                           fifelse(`Library Type` == "VDJ T", "VDJ T: ", NA_character_))))]
  dt$new_name <- paste0(dt$new_name, dt$`Metric Name`)
  dt <- data.table::dcast(dt, sample ~ new_name, value.var = "Metric Value")
  colnames(dt)[1] <- "sample"
  return(dt)
}
.metricsCount <- function(object, sample_name=NULL) {
  dt <- object
  data.table::setnames(dt, old = colnames(dt)[grepl("%", dt[1,] )],
                       new = paste0(colnames(dt)[grepl("%", dt[1,])], " (%)"))
  dt[, (names(dt)) := lapply(.SD, function(x) {
    x <- gsub(",", "", x)   #
    x <- gsub("%", "", x)   #
    as.numeric(x)           #
  })]
  dt <- data.table::data.table(sample=sample_name, dt)
  data.table::setnames(dt, old = colnames(dt)[ !grepl("Antibody:", colnames(dt)) ],
                       new = paste0("RNA: ",  colnames(dt)[ !grepl("Antibody:", colnames(dt)) ] ) )
  colnames(dt)[1] <- "sample"

  return(dt)
}
.metricsVDJ <- function(object, sample_name=NULL) {
  dt <- object
  data.table::setnames(dt, old = colnames(dt)[grepl("%", dt[1,] )],
                       new = paste0(colnames(dt)[grepl("%", dt[1,])], " (%)"))
  dt[, (names(dt)) := lapply(.SD, function(x) {
    x <- gsub(",", "", x)   #
    x <- gsub("%", "", x)   #
    as.numeric(x)           #
  })]
  dt <- data.table::data.table(sample=sample_name, dt)
  if( sum(grepl("TRA|TRB|TRG|TRD", colnames(dt)))!=0 ){
    data.table::setnames(dt, old = colnames(dt),
                         new = paste0("VDJ T: ",  colnames(dt)) )
  }else{
    data.table::setnames(dt, old = colnames(dt),
                         new = paste0("VDJ B: ",  colnames(dt)) )
  }
  colnames(dt)[1] <- "sample"
  return(dt)
}


















