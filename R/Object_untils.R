

toSeuratObject <- function(object, ...){
  if (!("Seurat" %in% class(object))) {
    if ("Antibody Capture" %in% names(object)) {
      rna <- Seurat::CreateSeuratObject(object$`Gene Expression`, min.cells = 0, ...)
      if( "Assay5" %in% class(rna@assays$RNA) ){
        adt <- SeuratObject::CreateAssay5Object(object$`Antibody Capture`, min.cells = 0, ...)
      }else{
        adt <- SeuratObject::CreateAssayObject(object$`Antibody Capture`, min.cells = 0, ...)
      }
      rna[["ADT"]] <- adt
      Seurat::DefaultAssay(rna) <- 'RNA'
    }else if( "dgCMatrix" %in% class(object) ){
      rna <- Seurat::CreateSeuratObject(object, min.cells = 0, ...)
    }else{
      rna <- Seurat::CreateSeuratObject(object$`Gene Expression`, min.cells = 0, ...)
    }
  }else{
    rna <- object
  }
  rm(object)
  return(rna)
}




#' @rdname ConvertToBPCells
#' @method ConvertToBPCells Seurat
#' @export
ConvertToBPCells.Seurat <- function(object, BPdir="./BPCellData/", assay = "RNA", slot="counts", ...) {
  BPdir <- paste0(BPdir, "/", assay)
  data <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  BPCells::write_matrix_dir(data, dir = BPdir, overwrite = TRUE)
  data <- BPCells::open_matrix_dir(BPdir)
  object <- SeuratObject::SetAssayData(object = object,
                                       assay = assay,
                                       slot = slot,
                                       new.data = data)
  return(object)
}

#' @rdname ConvertToBPCells
#' @method ConvertToBPCells dgCMatrix
#' @export
ConvertToBPCells.dgCMatrix <- function(object, BPdir="./BPCellData/", ...) {
  BPCells::write_matrix_dir(object, dir = BPdir, overwrite = TRUE)
  data <- BPCells::open_matrix_dir(BPdir)
  return(data)
}


#' @rdname ConvertToBPCells
#' @method ConvertToBPCells IterableMatrix
#' @export
ConvertToBPCells.IterableMatrix <- function(object, BPdir="./BPCellData/", ...) {
  BPCells::write_matrix_dir(object, dir = BPdir, overwrite = TRUE)
  data <- BPCells::open_matrix_dir(BPdir)
  return(data)
}


#' @rdname ConvertToSCE
#' @method ConvertToSCE Seurat
#' @export
ConvertToSCE.Seurat <- function(object, assay="RNA", slot="counts", ...) {
  x <- suppressWarnings(Seurat::as.SingleCellExperiment(object, assay = assay, slot = slot))
  return(x)
}

#' @rdname ConvertToSCE
#' @method ConvertToSCE SingleCellExperiment
#' @export
ConvertToSCE.SingleCellExperiment <- function(object, ...) {
  return(object)
}

#' @rdname ConvertToSCE
#' @method ConvertToSCE default
#' @export
ConvertToSCE.default <- function(object, ...) {
  x <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = object))
  return(x)
}


#' @method ConvertToSeurat list
ConvertToSeurat.list <- function(object, ...) {
  if (!("Gene Expression" %in% names(object))) {
    stop("Input list must contain a 'Gene Expression' component.")
  }

  # Create Seurat object from "Gene Expression"
  rna <- Seurat::CreateSeuratObject(object$`Gene Expression`, min.cells = 1, ...)

  # Add "Antibody Capture" as an additional assay if present
  if ("Antibody Capture" %in% names(object)) {
    if ("Assay5" %in% class(rna@assays$RNA)) {
      adt <- SeuratObject::CreateAssay5Object(object$`Antibody Capture`, min.cells = 1, ...)
    } else {
      adt <- SeuratObject::CreateAssayObject(object$`Antibody Capture`, min.cells = 1, ...)
    }
    rna[["ADT"]] <- adt
    Seurat::DefaultAssay(rna) <- "RNA"
  }

  return(rna)
}

#' @method ConvertToSeurat default
ConvertToSeurat.default <- function(object, ...) {
  Seurat::CreateSeuratObject(object, min.cells = 1, ...)
}

#' @method ConvertToSeurat Seurat
ConvertToSeurat.Seurat <- function(object, ...) {
  return(object)
}


#' @method getSampleName Seurat
getSampleName.Seurat <- function(object, sample.by="orig.ident", ...) {
  return(unique(object@meta.data[[sample.by]]) )
}

#' @method splitObject Seurat
splitObject.Seurat <- function(object, split.by = "orig.ident", assay=NULL,tmpdir= "./temp/SingleCellMQC_tempBPCellSplitSeurat/",  ...){
  if(!is.null(split.by)){
    if(  length(unique(object@meta.data[[split.by]]))>1  ){
      message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-------- Split Seurat Object"))
      split_object <- .splitObjectSeurat(object, split.by = split.by, assay=assay, tmpdir=tmpdir)
    }else{
      split_object <- list(Sample = object)
      names(split_object) <- as.character(unique(object@meta.data[[split.by]]))
    }
  }else{
    split_object <- list(Sample = object)
  }
  return(split_object)
}

#' @method splitObject Seurat
splitObject.list <- function(object, split.by = NULL, ...){
  return(object)
}

.splitObjectSeurat <- function(object, assay=NULL, split.by="orig.ident", tmpdir= "./temp/SingleCellMQC_tempBPCellSplitSeurat/"){
  if(is.null(assay)){
    assay_names <- SeuratObject::Assays(object)
  }else{
    assay_names = assay
  }

  metalist <- split(object@meta.data, as.character(object@meta.data[[split.by]]) )
  exp_list <- lapply(assay_names, function(assay_name) {
    counts <- SeuratObject::GetAssayData(object, assay = assay_name, slot = "counts")
    return(counts)
  })
  names(exp_list) <- assay_names

  if (is(exp_list[[1]], "IterableMatrix")) {
    # Check if the folder exists, and create it if it doesn't
    if (!dir.exists(tmpdir)) {
      dir.create(tmpdir, recursive = TRUE)
      message("BPCells Temp folder created: ", tmpdir)
    } else {
      message("BPCells Temp folder already exists: ", tmpdir)
    }
  }
  rm(object)
  ## Run parallel processing
  # progressr::with_progress({
  p <- progressr::progressor(along = 1:(length(metalist)+1) )
  object_list <- smart_lapply(metalist, function(x) {
    cell_name <- rownames(x)
    # first
    subexp <- exp_list[[ assay_names[1] ]][, cell_name, drop = FALSE]
    if (is(subexp, "IterableMatrix")) {
      subexp <- suppressWarnings(suppressMessages(
        BPCells::write_matrix_dir(subexp, tempfile("mat", tmpdir = tmpdir))
      ))
    }
    subObject <- suppressWarnings(suppressMessages(
      SeuratObject::CreateSeuratObject(subexp, meta.data = x, assay = assay_names[1])
    ))


    if (length(assay_names) > 1) {
      is_assay5 <- "Assay5" %in% class(subObject[[assay_names[1]]])

      exp_matrices <- lapply(seq_along(assay_names[-1]), function(i) {
        tempexp <- exp_list[[assay_names[-1][i]]][, cell_name, drop = FALSE]
        if (is(tempexp, "IterableMatrix")) {
          tempexp = BPCells::write_matrix_dir(tempexp, tempfile("mat", tmpdir = tmpdir))
        }
        return(tempexp)
      })
      for (i in seq_along(assay_names[-1])) {
        subObject[[assay_names[-1][i]]] <- if (is_assay5) {
          SeuratObject::CreateAssay5Object(exp_matrices[[i]])
        } else {
          SeuratObject::CreateAssayObject(exp_matrices[[i]])
        }
      }

      Seurat::DefaultAssay(subObject) <- assay_names[1]
    }
    p()
    return(subObject)
  }
  ,
  future.packages = c("BPCells"), future.seed = NULL
  )
  # }
  # ,handlers=  progressr::handlers(progressr::handler_progress(
  #   format = ":spin :bar :percent :elapsed :eta",
  #   clear = FALSE
  # ))
  # )
  rm(exp_list)
  return(object_list)
}



#' Merge Matrices
#'
#' This function merges a list of matrices (either `dgCMatrix` or `IterableMatrix`) into a single matrix.
#'
#' @param matrix_list A list of matrices to merge.
#' @param add_sample_ids Optional vector of sample IDs to append to column names.
#' @param prefix Logical, whether to add sample IDs as a prefix (default: TRUE).
#' @param cell_id_delimiter Character, delimiter to separate sample IDs and column names (default: "_").
#' @param ... Additional arguments passed to methods.
#'
#' @return A merged matrix of the same class as the input matrices.
#' @export
MergeMatrix <- function(matrix_list, add_sample_ids = NULL, prefix = TRUE, cell_id_delimiter = "_", ...) {
  # Check if input is a list
  if (!is.list(matrix_list)) {
    stop("Input must be a list of matrices.")
  }

  # Check if the list is empty
  if (length(matrix_list) == 0) {
    stop("Input list is empty.")
  }

  # Determine the class of the first matrix in the list without creating a new object
  if (inherits(matrix_list[[1]], "dgCMatrix")) {
    MergeMatrix_dgCMatrix(matrix_list, add_sample_ids, prefix, cell_id_delimiter, ...)
  } else if (inherits(matrix_list[[1]], "IterableMatrix")) {
    MergeMatrix_IterableMatrix(matrix_list, add_sample_ids, prefix, cell_id_delimiter, ...)
  } else {
    stop("Unsupported matrix class. Supported classes: dgCMatrix, IterableMatrix.")
  }
}


MergeMatrix_dgCMatrix <- function(
    matrix_list,
    add_sample_ids = NULL,
    prefix = TRUE,
    cell_id_delimiter = "_"
) {
  message("Merge dgCMatrix")
  # Check if input is a list and all elements are dgCMatrix
  if (!is.list(matrix_list) || !all(sapply(matrix_list, inherits, "dgCMatrix"))) {
    stop("Input must be a list of sparse matrices in dgCMatrix format.")
  }

  # If matrix_list has no names, generate default names (Sample1, Sample2, etc.)
  if (is.null(names(matrix_list))) {
    names(matrix_list) <- paste0("Sample", seq_along(matrix_list))
  }

  # Check for duplicate column names
  duplicated_barcodes <- any(duplicated(unlist(lapply(matrix_list, colnames))))

  # If add_sample_ids is NULL, handle duplicates only if they exist
  if (is.null(add_sample_ids)) {
    if (isTRUE(duplicated_barcodes)) {
      warning("Duplicate column names detected in input matrices. Unique column names will be generated dynamically using `matrix_list` names.")
      # Dynamically generate unique column names
      name_matrix = names(matrix_list)
      matrix_list <- lapply(seq_along(matrix_list), function(i) {
        mat <- matrix_list[[i]]
        if (prefix) {
          colnames(mat) <- paste0(name_matrix[i], cell_id_delimiter, colnames(mat))
        } else {
          colnames(mat) <- paste0(colnames(mat), cell_id_delimiter, name_matrix[i])
        }
        return(mat)
      })
    }
  } else {
    # If add_sample_ids is provided, check length and use it to generate unique column names
    if (length(add_sample_ids) != length(matrix_list)) {
      stop("The length of `add_sample_ids` must match the length of `matrix_list`.")
    }
    names(matrix_list) <- add_sample_ids
    # Dynamically generate unique column names
    name_matrix = names(matrix_list)
    matrix_list <- lapply(seq_along(matrix_list), function(i) {
      mat <- matrix_list[[i]]
      if (prefix) {
        colnames(mat) <- paste0(name_matrix[i], cell_id_delimiter, colnames(mat))
      } else {
        colnames(mat) <- paste0(colnames(mat), cell_id_delimiter, name_matrix[i])
      }
      return(mat)
    })
  }

  # Get the union of all row names
  all_rows <- unique(unlist(lapply(matrix_list, rownames), use.names = F))

  # Get all column names
  all_colnames <- unlist(lapply(matrix_list, colnames), use.names = F)

  # Check if dynamically generated column names are still duplicated
  if (any(duplicated(all_colnames))) {
    stop("Dynamically generated column names still contain duplicates. Ensure `matrix_list` names are unique.")
  }

  # Parallelize the processing of each matrix
  p <- progressr::progressor(along = seq_along(matrix_list))
  matrix_data <- smart_lapply(seq_along(matrix_list), function(k) {
    mat <- matrix_list[[k]]
    mat_i <- mat@i + 1  # Convert to R-style row indices (1-based)
    mat_p <- mat@p
    mat_x <- mat@x
    # Adjust column indices
    mat_j <- rep(seq_along(mat_p[-1]), diff(mat_p)) + ifelse(k == 1, 0, sum(sapply(matrix_list[1:(k-1)], ncol)))
    # Match row names to the union
    mat_rows <- match(rownames(mat)[mat_i], all_rows)
    p()
    # Return data
    list(i = mat_rows, j = mat_j, x = mat_x)
  }, future.seed = TRUE)

  # Combine all data
  all_i <- unlist(lapply(matrix_data, function(x) x$i))
  all_j <- unlist(lapply(matrix_data, function(x) x$j))
  all_x <- unlist(lapply(matrix_data, function(x) x$x))

  # Calculate total number of columns
  total_cols <- sum(sapply(matrix_list, ncol))

  # Create the new sparse matrix
  new_mat <- Matrix::sparseMatrix(i = all_i, j = all_j, x = all_x,
                                  dims = c(length(all_rows), total_cols),
                                  dimnames = list(all_rows, all_colnames),
                                  repr = "C")
  return(new_mat)
}


MergeMatrix_IterableMatrix <- function(
    matrix_list,
    add_sample_ids = NULL,
    prefix = TRUE,
    cell_id_delimiter = "_",
    ...
) {
  on.exit(expr = SeuratObject::CheckGC())  # Ensure garbage collection runs on exit
  message("Merge IterableMatrix")

  # Check if the input is a list
  if (!is.list(matrix_list)) {
    matrix_list <- list(matrix_list)
  }

  # Handle column name duplication
  # If matrix_list has no names, generate default names (Sample1, Sample2, etc.)
  if (is.null(names(matrix_list))) {
    names(matrix_list) <- paste0("Sample", seq_along(matrix_list))
  }

  # Check for duplicate column names
  duplicated_barcodes <- any(duplicated(unlist(lapply(matrix_list, colnames))))

  # If add_sample_ids is NULL, handle duplicates only if they exist
  if (is.null(add_sample_ids)) {
    if (isTRUE(duplicated_barcodes)) {
      warning("Duplicate column names detected in input matrices. Unique column names will be generated dynamically using `matrix_list` names.")
      # Dynamically generate unique column names
      name_matrix <- names(matrix_list)
      matrix_list <- lapply(seq_along(matrix_list), function(i) {
        mat <- matrix_list[[i]]
        if (prefix) {
          colnames(mat) <- paste0(name_matrix[i], cell_id_delimiter, colnames(mat))
        } else {
          colnames(mat) <- paste0(colnames(mat), cell_id_delimiter, name_matrix[i])
        }
        return(mat)
      })
    }
  } else {
    # If add_sample_ids is provided, check its length and use it to generate unique column names
    if (length(add_sample_ids) != length(matrix_list)) {
      stop("The length of `add_sample_ids` must match the length of `matrix_list`.")
    }
    names(matrix_list) <- add_sample_ids
    # Dynamically generate unique column names
    name_matrix <- names(matrix_list)
    matrix_list <- lapply(seq_along(matrix_list), function(i) {
      mat <- matrix_list[[i]]
      if (prefix) {
        colnames(mat) <- paste0(name_matrix[i], cell_id_delimiter, colnames(mat))
      } else {
        colnames(mat) <- paste0(colnames(mat), cell_id_delimiter, name_matrix[i])
      }
      return(mat)
    })
  }

  # Check if dynamically generated column names are still duplicated
  all_colnames <- unlist(lapply(matrix_list, colnames), use.names = FALSE)
  if (any(duplicated(all_colnames))) {
    stop("Dynamically generated column names still contain duplicates. Ensure `matrix_list` names are unique.")
  }

  # Get the union of all row names
  all_rows <- unique(unlist(lapply(matrix_list, rownames)))

  # Process each matrix to fill missing rows
  for (i in seq_along(along.with = matrix_list)) {
    missing_row <- setdiff(x = all_rows, y = rownames(x = matrix_list[[i]]))

    if (length(x = missing_row) > 0) {
      zero_i <- SeuratObject::SparseEmptyMatrix(
        nrow = length(x = missing_row),
        ncol = ncol(x = matrix_list[[i]]),
        colnames = colnames(x = matrix_list[[i]]),
        rownames = missing_row
      )
      zero_i <- as(zero_i, Class = 'IterableMatrix')
      matrix_list[[i]] <- rbind(matrix_list[[i]], zero_i)[all_rows, ]
    }
  }

  # Combine all matrices
  m <- Reduce(f = cbind, x = matrix_list)

  return(m)
}



quiet <- function(
    x,
    print_cat = TRUE,
    message = TRUE,
    warning = TRUE
) {
  stopifnot(
    is.logical(print_cat) && length(print_cat) == 1,
    is.logical(message) && length(message) == 1,
    is.logical(warning) && length(warning) == 1
  )

  if (print_cat) {
    tmp_output <- tempfile()
    sink(tmp_output, type = "output")
    on.exit({
      sink(type = "output")
      unlink(tmp_output)
    }, add = TRUE)
  }

  if (warning && message) {
    invisible(suppressMessages(suppressWarnings(force(x))))
  } else if (warning && !message) {
    invisible(suppressWarnings(force(x)))
  } else if (!warning && message) {
    invisible(suppressMessages(force(x)))
  } else {
    invisible(force(x))
  }
}


smart_lapply <- function(X, FUN, ..., .use_future = NULL, future.seed = TRUE) {
  if (is.null(.use_future)) {
    .use_future <- !inherits(future::plan(), "sequential")
  }
  if (.use_future) {
    future.apply::future_lapply(X, FUN, ..., future.seed = future.seed)
  } else {
    lapply(X, FUN)
  }
}
