

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


.standardizeNames <- function(x) {
  if (is.null(x)) {
    stop("Input object has no `names` attribute.")
  }

  old_names <- x
  new_names <- gsub("_", "-", old_names)

  if (!identical(old_names, new_names)) {
    warning(
      "Warning: names have underscores ('_'), replacing with dashes ('-').",
      call. = FALSE  # Suppresses the function call trace in the warning
    )
    x <- new_names
  }

  return(x)
}



.findDirValue <- function(obj) {
  if (isS4(obj)) {   # 这里直接用isS4，无需加methods::
    slots <- methods::slotNames(obj)
    if ("dir" %in% slots) {
      return(methods::slot(obj, "dir"))
    } else {
      for (s in slots) {
        result <- .findDirValue(methods::slot(obj, s))
        if (!is.null(result)) return(result)
      }
    }
  }
  return(NULL)
}



.bar_style <- function(width = 1, fill = "#e6e6e6", height = "100%",
                       align = c("left", "right"), color = NULL) {
  align <- match.arg(align)
  if (align == "left") {
    position <- paste0(width * 100, "%")
    image <- sprintf("linear-gradient(90deg, %1$s %2$s, transparent %2$s)", fill, position)
  } else {
    position <- paste0(100 - width * 100, "%")
    image <- sprintf("linear-gradient(90deg, transparent %1$s, %2$s %1$s)", position, fill)
  }
  list(
    # paddingLeft = "0.5%",
    backgroundImage = image,
    backgroundSize = paste("92%", height),
    backgroundRepeat = "no-repeat",
    backgroundPosition = "left",
    color = color
  )
}


.add_custom_text <- function(table, text, font_size = 16) {
  htmlwidgets::prependContent(
    table,
    htmltools::tags$p(
      text,
      style = paste0(
        "color:#000;",
        "background:#FFFFFF;",
        "text-align:left;",
        "font-size:", font_size, "px;",
        "font-style:normal;",
        "font-weight:bold;",
        "text-decoration:none;",
        "text-transform:none;",
        "letter-spacing:normal;",
        "word-spacing:normal;",
        "text-shadow:none;",
        "margin-top:0.5px;",
        "margin-right:0px;",
        "margin-bottom:0.5px;",
        "margin-left:0px;"
      )
    )
  )
}

.re_table <- function(object, csv.name="table", elementId="table", right.sparkline=F,
                      down.sparkline=F, label=1, first_name=colnames(object)[1],
                      other_sticky_column=NULL,
                      maxWidth=85,
                      number_type="auto",##auto, custom, initial
                      number_fmr="paste0(value,'%')",
                      subtitle=NULL
){
  count <- as.data.frame(object)
  colnames(count)[1] <- first_name
  filter_code <- "filterable = TRUE,
  filterMethod = htmlwidgets::JS(
    'function(rows, columnId, filterValue) {
                  return rows.filter(function(row) {
                    return row.values[columnId] >= filterValue
                  })
                }'
  ),"

  filter_input_code <- "
  filterable=TRUE,
    filterInput = function(values, name) {
      htmltools::tags$select(
        onchange = paste0('Reactable.setFilter(\\'', elementId, '\\', \\'', name, '\\', event.target.value || undefined)'),
        htmltools::tags$option(value = '', 'All'),
        lapply(unique(values), htmltools::tags$option),
        'aria-label' = sprintf('Filter %s', name),
        style = 'width: 100%; height: 28px;'
      )
    }
"
  cell_col <- switch(number_type,
                     "initial" = "",
                     "custom" = paste0("cell = function(value) { return(", number_fmr, ")}," ),
                     "auto" =""
  )




  count_table <- lapply(colnames(count), function(column) {
    if (is.numeric(count[[column]])) {
      if(number_type=="auto"){
        if(grepl("%", column)){
          cell_col <- "cell = function(value) {
            return(paste0( round(value,2),'%'))
          },"
        }else{
          cell_col <- "cell = function(value) {
              return(round(value,3))
            },"
        }
      }

      eval(parse(text = paste0("reactable::colDef(", cell_col,
                               filter_code,
                               "style = function(value) {",
                               "  if (!is.na(value)) {.bar_style(width = value / max(count$`", column, "` *1, na.rm =T), fill = '#A6CEE3', color = 'black', height = '80%')}",
                               "}, ",
                               "align = 'left')"
      )))
    } else {
      eval(parse(text = "reactable::colDef(
          cell = function(value) {
            return(value)
          },
          align = 'left'
        )"))
    }
  })
  names(count_table) <- colnames(count)

  count_table[[first_name]] <- eval(parse(text = paste0("reactable::colDef( sticky = 'left', style=list(background='#f7f7f7'),",filter_input_code ,")" )))


  if(!is.null(other_sticky_column)){
    other_sticky_list <- lapply(other_sticky_column, function(x){
      eval(parse(text = paste0("reactable::colDef( sticky = 'left', style=list(background='#f7f7f7'))" )))
    })
    names(other_sticky_list) <- colnames(count)[other_sticky_column]
    count_table <- c(count_table, other_sticky_list )
  }

  if(right.sparkline){
    count_table$Sparkline <-  eval(parse(text = paste0("reactable::colDef(sticky = 'right',style=list(background='#f7f7f7'), cell = function(values) {sparkline::sparkline(values, type = 'line')},",")" )))
    count$Sparkline <- apply(count[,-c(label),drop=F],1,function(x){
      c(as.numeric(x))
    },simplify = F)
  }
  if(down.sparkline){
    footer_function <- "function(values) {
      if (!is.numeric(values)) return()
      sparkline::sparkline(values, type = 'line')
    }"
  }else{
    footer_function<-NULL
  }

  defaultColDef_list <- eval(parse(text = paste0(
    "reactable::colDef(",",footer=", footer_function,
    ",align = 'left',
      maxWidth=",maxWidth,")")))


  out <- htmltools::browsable(
    htmltools::tagList(
      htmltools::tags$button("Download as CSV", onclick = paste0("Reactable.downloadDataCSV('",elementId, "', '", csv.name, "')")) ,
      reactable::reactable(
        count,
        columns = count_table,
        filterable = TRUE,
        theme = reactablefmtr::cosmo(header_font_size =14, font_size =14,  cell_padding =4),
        showPageSizeOptions = TRUE,
        elementId = elementId,
        defaultPageSize =10 ,
        defaultColDef = defaultColDef_list,
        showPagination = TRUE
      ) %>% .add_custom_text(subtitle, font_size = 16)
    )
  )
  return(out)
}


.re_table_fmtr <- function(object, csv.name="table", elementId="table", max="single", right.sparkline=F,
                           down.sparkline=F, label=1, first_name=colnames(object)[1],
                           other_sticky_column=NULL,
                           maxWidth=85,
                           subtitle=NULL,
                           bar_type= "data_bars(
                             count,
                             text_position = 'inside-base',
                             number_fmt = scales::percent,
                             fill_color 	='#A6CEE3',
                             animation='none'
                           )"
){
  count <- as.data.frame(object)
  colnames(count)[1] <- first_name
  filter_input_code <- "
  filterable=TRUE,
    filterInput = function(values, name) {
      htmltools::tags$select(
        onchange = paste0('Reactable.setFilter(\\'', elementId, '\\', \\'', name, '\\', event.target.value || undefined)'),
        htmltools::tags$option(value = '', 'All'),
        lapply(unique(values), htmltools::tags$option),
        'aria-label' = sprintf('Filter %s', name),
        style = 'width: 100%; height: 28px;'
      )
    }
"
  count_table <-list()


  count_table[[first_name]] <- eval(parse(text = paste0("reactable::colDef( sticky = 'left', style=list(background='#f7f7f7'),",filter_input_code ,")" )))
  if(!is.null(other_sticky_column)){
    other_sticky_list <- lapply(other_sticky_column, function(x){
      eval(parse(text = paste0("reactable::colDef( sticky = 'left', style=list(background='#f7f7f7'))" )))
    })
    names(other_sticky_list) <- colnames(count)[other_sticky_column]
    count_table <- c(count_table, other_sticky_list )
  }

  if(right.sparkline){
    count_table$Sparkline <-  eval(parse(text = paste0("reactable::colDef(sticky = 'right',style=list(background='#f7f7f7'), cell = function(values) {sparkline::sparkline(values, type = 'line')},",")" )))
    count$Sparkline <- apply(count[,-c(label),drop=F],1,function(x){
      c(as.numeric(x))
    },simplify = F)
  }
  if(down.sparkline){
    footer_function <- "function(values) {
      if (!is.numeric(values)) return()
      sparkline::sparkline(values, type = 'line')
    }"
  }else{
    footer_function<-NULL
  }


  defaultColDef_list <- eval(parse(text = paste0(
    "reactable::colDef(
    cell = ", bar_type,",footer=", footer_function,
    ",align = 'left',
      maxWidth=",maxWidth,")")))

  out <- htmltools::browsable(
    htmltools::tagList(
      htmltools::tags$button("Download as CSV", onclick = paste0("Reactable.downloadDataCSV('",elementId, "', '", csv.name, "')")) ,
      reactable::reactable(
        count,
        columns = count_table,
        filterable = TRUE,
        theme = reactablefmtr::cosmo(header_font_size =14, font_size =14, cell_padding =4),
        showPageSizeOptions = TRUE,
        elementId = elementId,
        defaultPageSize =10 ,
        defaultColDef = defaultColDef_list,
        showPagination = TRUE
      ) %>% .add_custom_text(subtitle, font_size = 16)
    )
  )
  return(out)
}


get_colors <- function(n) {
  base_colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FDBF6F", "#FF7F00",
                   "#FB9A99", "#E31A1C", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
  if (n <= length(base_colors)) {
    return(base_colors[1:n])
  } else {
    color_palette <- colorRampPalette(base_colors)
    return(color_palette(n))
  }
}



