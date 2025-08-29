

.lmR2 <- function(x, y) {
  if (ncol(x) != length(y)) {
    stop("The number of columns in x must match the length of y.")
  }

  mean_x <- BPCells::rowMeans(x)
  mean_y <- mean(y)

  # beta1
  ma_value <- BPCells::multiply_cols(BPCells::add_rows(x, - mean_x),  y - mean_y)
  numerator <- BPCells::rowSums( ma_value  )
  denominator <- BPCells::rowSums((x-mean_x)^2)
  beta1 <- numerator / denominator

  # beta0
  beta0 <- mean_y - beta1 * mean_x

  # pred value
  y_pred <- beta0 + x*beta1

  # R^2
  ss_total <- sum((y - mean_y)^2)
  ss_residual <- BPCells::rowSums( (BPCells::add_cols(y_pred, -y)) ^2  )
  R2 <- 1 - ss_residual / ss_total
  names(R2) <- rownames(x)
  return(R2)
}

.aovR2 <- function (x, y){
  if (ncol(x) != length(y)) {
    stop("The number of columns in x must match the length of y.")
  }

  #SST
  SST <- BPCells::rowVars(x)*( dim(x)[2]-1 )
  ##SSE
  group_indices <- split(colnames(x), y)
  result_list <- lapply(names(group_indices), function(group_name) {
    cols <- group_indices[[group_name]]
    if (length(cols) == 1) {
      return(0)
    }
    group_vars <- BPCells::rowVars(x[, cols, drop = FALSE]) * (length(cols) -1)
    return(group_vars)
  })
  SSE <- Matrix::rowSums( do.call(cbind, result_list) )
  R2 <- 1 - SSE/SST
  names(R2) <- rownames(x)
  return(R2)
}

#' @rdname RunVarExplained
#' @export
RunVarExplained.Seurat <- function(object, assay="RNA", variables=NULL, ...){

  if(is.null(variables)){
    stop(" `variables` is missing")
  }

  if(assay=="RNA"){
    SeuratObject::DefaultAssay(object) <-"RNA"
    object <- Seurat::NormalizeData(object,assay="RNA", verbose =F)
    exp <- Seurat::GetAssayData(object, assay = "RNA", slot = "data")
  }else if(assay=="ADT"){
    SeuratObject::DefaultAssay(object) <- 'ADT'
    if ( "BPCells" %in% attr(class(Seurat::GetAssayData(object, assay = "ADT", slot = "counts")), "package") ) {
      expADT <- Seurat::GetAssayData(object, assay = "ADT", slot = "counts")
      log_sums <- BPCells::colSums(log1p(expADT ), na.rm = TRUE)
      scale_factors <- exp(log_sums / nrow(expADT))
      exp <- log1p(  BPCells::multiply_cols(expADT ,1/ scale_factors) )
    }else{
      Seurat::DefaultAssay(object) <- "ADT"
      Seurat::VariableFeatures(object) <- rownames(object[["ADT"]])
      object <- Seurat::NormalizeData(object, normalization.method = 'CLR', margin = 2, verbose =F)
      exp <- Seurat::GetAssayData(object, assay = "ADT", slot = "data")
    }
  }else{
    stop("Invalid `assay`, only `RNA` or `ADT`.")
  }

  if(is(exp, "dgCMatrix") ){
    exp <- as(exp, "IterableMatrix")
  }else if(is(exp, "matrix") ){
    exp <- as(exp, "dgCMatrix")
    exp <- as(exp, "IterableMatrix")
  }

  rsquared_mat <- .getVarExplained(object=exp, metadata=object@meta.data, variables = variables)
  return(rsquared_mat)
}

#' @rdname RunVarExplained
#' @export
RunVarExplained.IterableMatrix <- function(object, metadata, variables=NULL, ...) {
  out <- .getVarExplained(object=object, metadata=metadata, variables=variables, ...)
  return(out)
}

#' @rdname RunVarExplained
#' @export
RunVarExplained.default <- function(object, metadata, variables=NULL, ...) {
  if(is(object, "dgCMatrix") ){
    object <- as(object, "IterableMatrix")
  }else if(is(object, "matrix") |  is(object, "data.frame")){
    object <- as(object, "dgCMatrix")
    object <- as(object, "IterableMatrix")
  }
  out <- .getVarExplained(object=object, metadata=metadata, variables=variables, ...)
  return(out)
}


.getVarExplained <- function(object, metadata, variables = NULL) {
  if (is.null(variables)) {
    variables <- colnames(metadata)
  }
  R2 <- lapply(variables, function(x) {
    if (is.character(metadata[[x]]) | is.factor(metadata[[x]])) {
      metadata[[x]] <- as.character(metadata[[x]])
      if (length(unique(metadata[[x]])) <= 1) {
        warning(sprintf("Variable '%s' has fewer than 2 unique levels. Skipping.", x))
        return(NULL)
      }
      if (length(unique(metadata[[x]])) == nrow(metadata)) {
        warning(sprintf(
          "Variable '%s' has a unique value for every observation (R2=100%% is trivial and uninformative). Skipping.",
          x
        ))
        return( NULL )
      }

      out <- .aovR2(object, metadata[[x]])

    } else if (is.numeric(metadata[[x]])) {
      if (nrow(metadata) <= 2) {
        warning(sprintf(
         " Variable '%s' has <= 2 value for lm. Skipping.",
          x
        ))
        return( NULL )
      }

      out <- .lmR2(object, metadata[[x]])

    } else {
      stop("Metadata must be numeric or character.")
    }

    return(out)
  })
  names(R2) <- variables
  R2 <- base::Filter(Negate(is.null), R2)
  R2 <- do.call(cbind, R2)
  rownames(R2) <- rownames(object)
  return(R2 * 100)
}






#' Aggregate Counts by Group Across Cells
#'
#' This function aggregates count data from a given matrix or data frame by grouping cells based on a specified metadata column. The output is a matrix where rows represent features and columns represent groups, with the values being the summed counts of the corresponding group.
#'
#' @param object A matrix or data frame containing the count data, where rows represent features (e.g., genes or proteins) and columns represent cells.
#' @param metadata A data frame containing metadata information for the cells. The number of rows should match the number of columns in `object`.
#' @param group.by A character string specifying the column in `metadata` to group cells by.
#'
#' @return A matrix where rows are features and columns are groups, with values representing the summed counts for each group.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # Aggregate counts by group
#' counts <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
#' metadata <- data.frame(group = rep(c("A", "B"), each = 5))
#' grouped_counts <- AggregateGroupAcrossCells(counts, metadata, group.by = "group")
#' }
#'
AggregateGroupAcrossCells <- function(object, metadata, group.by){
  split_object <- split(1:dim(object)[2], as.character(metadata[[group.by]]) )
  count_sum <- lapply(split_object, function(x){
    result <- Matrix::rowSums(object[, x, drop=F] , na.rm = T)
    return(result)
  })
  count_sum <- base::Filter(Negate(is.null), count_sum)
  count_sum <- do.call(cbind, count_sum)
  colnames(count_sum) <- names(split_object)
  return(count_sum)
}

.getVstExp <- function(count_sum, assay) {

  if (ncol(count_sum) < 2) {
    stop("Error: data must have at least 2 columns.")
  }

  # check each row
  if (all(rowSums(count_sum == 0) > 0)) {
    warning("Every gene has at least one zero. Adding pseudocount (+1) to avoid log(0) issues.")
    count_sum <- count_sum + 1  # +1
  }

  if (assay == "ADT") {
    if (nrow(count_sum) > 100) {
      vsd <- DESeq2::vst(count_sum, blind = TRUE, nsub = 100)
    } else {
      vsd <- DESeq2::vst(count_sum, blind = TRUE, nsub = nrow(count_sum) - 1)
    }
  } else {
    if (nrow(count_sum) > 1000) {
      vsd <- DESeq2::vst(count_sum, blind = TRUE, nsub = 1000)
    } else {
      vsd <- DESeq2::vst(count_sum, blind = TRUE)
    }
  }

  return(vsd)
}



.getPcaRes <- function(object, assay ="RNA", sample.by = "orig.ident",  ntop=2000, do.aggregate=T){
  if( !("Seurat" %in% is(object)) ){
    stop("Error: Input must be Seurat object.")
  }
  metadata <- object@meta.data
  exp <- Seurat::GetAssayData(object, slot = "counts", assay = assay)
  if(do.aggregate){
    count_sum <- AggregateGroupAcrossCells(exp, metadata, group.by = sample.by)
    vsd <- .getVstExp(count_sum, assay)
  }else{
    vsd <- as.matrix(exp)
  }
  var_result <- MatrixGenerics::rowVars(vsd)
  top_genes <- names(sort(var_result, decreasing = T))[1:ntop]
  vsd <- vsd[rownames(vsd) %in% top_genes,]
  t_vsd <- t(vsd)
  pca_res <- stats::prcomp(t_vsd)
  return(pca_res)
}
.pseudobulkMetadata <- function(object, sample.by, variable.by){
  index <- 0
  if( "Seurat" %in% is(object) ){
    sample_count <- object@misc[["SingleCellMQC"]][["perQCMetrics"]][["perSample"]][["count"]]
    object <- object@meta.data

    if(sum(variable.by %in% colnames(sample_count))!=0){
      count_table <- sample_count[, colnames(sample_count) %in% variable.by, drop=F]
      index <- 1
      variable.by <- setdiff(variable.by, colnames(count_table) )
    }

  }

  object[[sample.by]] <- as.character(object[[sample.by]])
  unique_sample_name <- unique(as.character(object[[sample.by]]))
  out_list <- lapply(variable.by, function(x){
    if( !is(object[[x]], "numeric") ){
      out  <- unique( object[, c(sample.by, x)] ,drop=F)[,2]
      out <- data.frame(out)
      colnames(out) <- x
    }else{
      out  <- data.table(object)[, lapply(.SD, function(y) mean( as.numeric(y), na.rm = TRUE) ),
                                 by = sample.by, .SDcols = x][[2]]
      out <- data.frame( out )
      if( nrow(object)==nrow(out) ){
        colnames(out) <- x
      }else{
        colnames(out) <- paste0("mean_",x)
      }

    }
    return(out)
  })
  out <- do.call(cbind, out_list)
  out <- data.frame(Sample=unique_sample_name, out, check.names = F)

  if(index==1){
    count_table<-count_table[match(out$Sample , count_table[[1]]), ]
    out <- data.frame(out, count_table, check.names = F)
  }
  return(out)
}
.pseudobulkExplanatory <- function(object, metadata){
  out <- RunVarExplained(t(object$x), variables= colnames(metadata)[-1], metadata= metadata)
  out <- data.frame(PCs= rownames(out), out, check.names = F)
  return(out)
}

.samplePCAResult <- function(object, assay ="RNA", sample.by = "orig.ident", variable.by=NULL, group.by= "orig.ident",
                             ntop=2000, color=NULL,size.point=3, size.text=3,Legend.lab=NULL, maxPCs=5,do.aggregate=T,
                             plot.type="pca", return.type="plot"){
  if( !("Seurat" %in% is(object)) ){
    stop("Error: Input must be Seurat object.")
  }

  ##
  if( length(setdiff(plot.type, c("pca", "svd", "explanatory")))!=0 ){
    stop("Invalid `plot.type`, only: `pca`, `svd` or/and `explanatory` ")
  }


  if( length(setdiff(return.type, c("plot", "interactive_table", "table")))!=0 ){
    stop("Invalid `return.type`, only: `plot`, `interactive_table`, or/and `table` ")
  }

  metadata <- object@meta.data
  pca_res <- suppressMessages(.getPcaRes(object=object, sample.by=sample.by, ntop = ntop,assay=assay, do.aggregate=do.aggregate))
  out <- list()

  if("pca" %in% plot.type){
    samples_groups <- .pseudobulkMetadata(object, sample.by = sample.by, variable.by =  group.by)
    samples_groups <- samples_groups[match(rownames(pca_res$x), samples_groups$Sample) , ]
    pca_out <- data.frame(pca_res$x, samples_groups[,2,drop=F])

    if("table" %in% return.type){
      pca_out_table <- data.frame(Sample=rownames(pca_res$x), pca_out)
      rownames(pca_out_table) <- NULL
      out$pca$table <- pca_out_table
    }

    if("interactive_table" %in% return.type){
      pca_out_table <- data.frame(Sample=rownames(pca_res$x), pca_out)
      rownames(pca_out_table) <- NULL
      out$pca$interactive_table <- .re_table_fmtr(pca_out_table, csv.name =paste0(assay, "_PCA"), elementId = paste0(assay,"_PCA"), right.sparkline = F, down.sparkline = F,
                                                  bar_type= "reactablefmtr::data_bars(
                             count,
                             text_position = 'inside-base',
                             number_fmt = scales::number,
                             fill_color 	='#A6CEE3',
                             animation='none'
                           )")
    }

    if("plot" %in% return.type){
      if(!is.null(Legend.lab)){
        colnames(pca_out)[length(colnames(pca_out))] <- Legend.lab
        group.by<-Legend.lab
      }

      if(is.null(color)){
        color <- get_colors( length(unique(pca_out[[group.by]])))
      }
      percentage <- round(pca_res$sdev/sum(pca_res$sdev) * 100,  2)
      percentage <- paste(colnames(pca_out), "(", paste(as.character(percentage), "%", ")", sep = ""))
      out$pca$plot <- ggplot2::ggplot(pca_out, mapping = ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data[[group.by]], shape = NULL,label = row.names(pca_out)))+
        ggplot2::geom_point(size = size.point, show.legend = TRUE,alpha = 0.9) +
        ggrepel::geom_text_repel(show.legend = F)+
        ggplot2::xlab(percentage[1]) +
        ggplot2::ylab(percentage[2]) +
        ggforce::geom_mark_ellipse(ggplot2::aes(fill = .data[[group.by]],label=NULL), color = grDevices::rgb(0, 0, 0, 0), alpha = 0.15, show.legend = F) +
        ggplot2::scale_color_manual(values = color ) +
        ggplot2::scale_fill_manual(values = color ) +
        ggplot2::theme_bw(base_size = 14)+
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8, fill = NA),
          axis.text = ggplot2::element_text(color = "black"))+
        ggplot2::ggtitle(assay)

    }
    if(length(out$pca)==1){
      out$pca <- out$pca[[1]]
    }
  }

  if("svd" %in% plot.type){
    samples_groups <- .pseudobulkMetadata(object, sample.by = sample.by, variable.by = variable.by )
    samples_groups <- samples_groups[match(rownames(pca_res$x), samples_groups$Sample) , ]
    variable_by <- colnames(samples_groups)[-1]
    sig_data <- .svdSig(pca_res, metadata = samples_groups, variables = variable_by, maxPCs = maxPCs)
    sig_data <- data.table(Variable=rownames(sig_data), sig_data, check.names = F)
    if("table" %in% return.type){
      out$svd$table <- sig_data
    }
    if("interactive_table" %in% return.type){
      out$svd$interactive_table <- .re_table(sig_data, csv.name = paste0(assay, "_SVD"), elementId = paste0(assay, "_SVD"), right.sparkline = F, down.sparkline = T, subtitle = "SVD significance pvalue")
    }

    if("plot" %in% return.type){
      sig_data <- melt.data.table(sig_data,id.vars="Variable")
      sig_data$Pvalue <- dplyr::case_when(
        sig_data$value >= 0.05 ~ "P >= 0.05",
        sig_data$value >= 0.01 ~ "P < 0.05",
        sig_data$value >= 0.001 ~ "P < 0.01",
        TRUE ~ "P < 0.001"
      )
      colnames(sig_data)[2] <- c("PCs")

      p_color <- colorRampPalette(c("white", "#E31A1C"))
      p_color <- p_color(4)
      out$svd$plot <- ggplot2::ggplot(sig_data)+
        ggplot2::geom_tile(
          ggplot2::aes(y=Variable,x=PCs,fill=Pvalue),na.rm=F,
          color="black",size=1)+
        ggplot2::theme_classic(base_size = 15)+
        shadowtext::geom_shadowtext(ggplot2::aes(y = Variable, x = PCs, label = round(sig_data[["value"]],3 )),  bg.colour = "white",color="black" )  +
        ggplot2::theme(axis.line=ggplot2::element_blank(),
              axis.ticks=ggplot2::element_blank(),
              axis.text.x=ggplot2::element_text(angle=60,hjust=0,vjust=0.5),
              legend.position="right")+
        ggplot2::guides(fill=ggplot2::guide_legend(title.position=NULL,byrow=TRUE))+
        ggplot2::labs(x=NULL,y=NULL)+
        ggplot2::scale_x_discrete(position="top")+
        ggplot2::scale_fill_manual(values=c("P >= 0.05"="white","P < 0.05"=p_color[2], "P < 0.01"= p_color[3],  "P < 0.001"= p_color[4]) )+
        ggplot2::coord_equal()
    }
    if(length(out$svd)==1){
      out$svd <- out$svd[[1]]
    }

  }

  if("explanatory" %in% plot.type){
    samples_groups <- .pseudobulkMetadata(object, sample.by = sample.by, variable.by = variable.by )
    samples_groups <- samples_groups[match(rownames(pca_res$x), samples_groups$Sample) , ]
    if(is.null(color)){
      color <- get_colors( dim(samples_groups)[1] )
    }
    #
    # return(list(object=pca_res, metadata=samples_groups))
    #

    plot_data <- .pseudobulkExplanatory(pca_res, samples_groups)
    plot_data <- plot_data[1:maxPCs, apply(plot_data, 2, function(x){sum(!is.na(x))}) > 0, drop=F]
    if(ncol(plot_data)==1){
      out$explanatory=NULL
    }else{

    rownames(plot_data) <- NULL
    if("table" %in% return.type){
      out$explanatory$table <- plot_data
    }
    if("interactive_table" %in% return.type){
      out$explanatory$interactive_table <- .re_table(plot_data, csv.name = paste0(assay,"_Explanatory"), elementId = paste0(assay,"_Explanatory"), right.sparkline = T, down.sparkline = F, subtitle = "Explanatory PCs for each variable",
                                                     number_type = "custom", number_fmr = "paste0( round(value,2),'%')")
    }

    if("plot" %in% return.type){
      plot_data <- melt.data.table(data.table(plot_data), id.vars = "PCs")
      # plot_data <- na.omit(plot_data)
      out$explanatory$plot <- ggplot2::ggplot(data=plot_data,ggplot2::aes(x=PCs,y=value, group = variable, color=variable ))+
        ggplot2::geom_point(size=2)+
        ggplot2::geom_line(size=1 )+
        ggplot2::theme_classic(base_size=15)+
        ggplot2::theme(panel.grid=ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank()
        )+
        ggplot2::scale_color_manual(values=color)+
        ggplot2::labs(y="% variance explained")
    }

    if(length(out$explanatory)==1){
      out$explanatory <- out$explanatory[[1]]
    }
    }
  }

  if(length(out$explanatory)==1){
    out$explanatory <- out$explanatory[[1]]
  }

  if(length(out)==1){
    out <- out[[1]]
  }
  return(out)
}


#' Visualize PCA Results on Pseudobulk Data
#'
#' This function generates PCA plots, SVD significance analysis, or explanatory variance plots using pseudobulk data. It allows for splitting data by a specified variable, aggregating counts, and visualizing the results.
#'
#' @param object A Seurat object containing the single-cell data to be analyzed.
#' @param assay A character string specifying the assay to use for PCA (default: `"RNA"`).
#' @param sample.by A character string specifying the metadata column used for defining samples (default: `"orig.ident"`).
#' @param svd.variable.by A character vector of metadata columns to use for SVD analysis (default: `NULL`).
#' @param pca.group.by A character string specifying the metadata column used for grouping cells in the PCA plot (default: `"orig.ident"`).
#' @param celltype.by A character string specifying the metadata column used to split the data into subsets (default: `NULL`).
#' @param ntop An integer specifying the number of top variable features to use for PCA (default: `2000`).
#' @param color A vector of colors to use for visualizations (default: `NULL`).
#' @param pca.size.point A numeric value specifying the point size for PCA plots (default: `3`).
#' @param pca.size.text A numeric value specifying the text size for PCA plots (default: `3`).
#' @param pca.Legend.lab A character string specifying the legend label for PCA plots (default: `NULL`).
#' @param svd.maxPCs An integer specifying the maximum number of PCs to include in the SVD plot (default: `3`).
#' @param plot.type A character vector specifying the types of plots to generate. Options include `"pca"`, `"svd"`, and `"explanatory"` (default: `"pca"`).
#' @param return.type A character vector specifying the output format. Options include `"plot"`, `"interactive_table"`, and `"table"` (default: `"plot"`).
#' @param do.aggregate Logical, whether to aggregate counts across groups before PCA (default: `TRUE`).
#'
#' @return A list containing plots, tables, or interactive tables based on the specified `plot.type` and `return.type`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate PCA plots
#' object <- CreateSeuratObject(counts = matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10))
#' pca_results <- PlotSamplePCA(object, sample.by = "orig.ident", plot.type = c("pca", "svd"))
#' }
#'
PlotSamplePCA <- function(object, assay = "RNA", sample.by = "orig.ident", svd.variable.by=NULL,
                          pca.group.by= "orig.ident", celltype.by = NULL,do.aggregate=T,
                          ntop=2000, color=NULL, pca.size.point=3, pca.size.text=3, pca.Legend.lab=NULL, svd.maxPCs=3,
                          plot.type="pca", return.type="plot"){
  if( !("Seurat" %in% is(object)) ){
    stop("Error: Input must be Seurat object.")
  }

  ##
  split_object <- list(all=object)
  tmpdir= "./temp/SingleCellMQC_tempBPCellSplitSeurat/"

  if(!is.null(celltype.by)){
    split_object <- c(split_object, splitObject(object, split.by = celltype.by))
  }
  out_list <- list()


  if(is.null(celltype.by)){
    out_list <- .samplePCAResult(split_object[[1]], assay = assay, sample.by = sample.by, variable.by = svd.variable.by,  maxPCs = svd.maxPCs,color=color,
                                 plot.type = plot.type, return.type = return.type, do.aggregate=do.aggregate, group.by = pca.group.by)
  }else{

    if("svd" %in% plot.type){
      out <- lapply(names(split_object), function(x){
        message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "svd ",  x))
        result <- tryCatch(
          {
            .samplePCAResult(
              split_object[[x]],
              assay = assay,
              sample.by = sample.by,
              variable.by = svd.variable.by,
              maxPCs = svd.maxPCs,
              plot.type = "svd",
              return.type = "table",
              do.aggregate = do.aggregate
            )
          },
          error = function(e) {
            warning("Warning in .samplePCAResult for element ", x, ": ", e$message)
            return(NULL)  # return NULL
          }
        )
      })
      names(out) <- names(split_object)
      out <- Filter(Negate(is.null), out)

      #
      plot_table <- lapply( names(out), function(x){
        data.frame(Group=x, data.table::melt.data.table(out[[x]], id.vars = "Variable"), check.names = F)
      })
      plot_table <- do.call(rbind, plot_table)
      plot_table <- split(plot_table, plot_table$variable)

      if("plot" %in% return.type){
        plot_list <- lapply(plot_table, function(sig_data){
          sig_data$Pvalue <- dplyr::case_when(
            sig_data$value >= 0.05 ~ "P >= 0.05",
            sig_data$value >= 0.01 ~ "P < 0.05",
            sig_data$value >= 0.001 ~ "P < 0.01",
            TRUE ~ "P < 0.001"
          )
          p_color <- colorRampPalette(c("white", "#E31A1C"))
          p_color <- p_color(4)

          p1 <- ggplot2::ggplot(sig_data)+
            ggplot2::geom_tile(
              ggplot2::aes(y=Group,x=Variable,fill=Pvalue),na.rm=F,
              color="black",size=1)+
            ggplot2::theme_classic(base_size = 15)+
            shadowtext::geom_shadowtext(ggplot2::aes(y = Group, x = Variable, label = round(.data[["value"]],3 )),  bg.colour = "white",color="black" )  +
            ggplot2::theme(axis.line=ggplot2::element_blank(),
                  axis.ticks=ggplot2::element_blank(),
                  axis.text.x=ggplot2::element_text(angle=60,hjust=0,vjust=0.5),
                  legend.position="right")+
            ggplot2::guides(fill=ggplot2::guide_legend(title.position=NULL,byrow=TRUE))+
            ggplot2::labs(x=NULL,y=NULL)+
            ggplot2::scale_x_discrete(position="top")+
            ggplot2::scale_fill_manual(values=c("P >= 0.05"="white","P < 0.05"=p_color[2], "P < 0.01"= p_color[3],  "P < 0.001"= p_color[4]) )
          return(p1)
        })
        names(plot_list) <- names(plot_table)
        out_list$svd$plot <- plot_list
      }

      if("interactive_table" %in% return.type){
        out_list$svd$interactive_table <- lapply( names(plot_table), function(x){
          in_table <- data.table::dcast(data.table(plot_table[[x]]), Group~Variable, value.var = "value")
          re_table <- .re_table(in_table, csv.name =paste0(assay,"_SVD_",x) ,elementId = paste0(assay, "_SVDtable_",x), down.sparkline = T, right.sparkline = T, maxWidth = NULL, subtitle = paste0("SVD significance pvalue (",x,")") )
          return(re_table)
        })
        names(out_list$svd$interactive_table) <- names(plot_table)
      }

      if(length(out_list$svd)==1){
        out_list$svd <- out_list$svd[[1]]
      }
    }

    if("explanatory" %in% plot.type){
      out <- lapply(names(split_object), function(x){
        message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "explanatory ",  x))
        result <- tryCatch(
          {
        .samplePCAResult(split_object[[x]], assay = assay, sample.by = sample.by, variable.by = svd.variable.by,  maxPCs = svd.maxPCs,
                         plot.type = "explanatory", return.type = "table",do.aggregate=do.aggregate)}
            ,error = function(e) {
              warning("Warning in .samplePCAResult for element ", x, ": ", e$message)
              return(NULL)  # return NULL
            }
        )
      })
      names(out) <- names(split_object)
      out <- Filter(Negate(is.null), out)

      if(length(out)==0){
        return(out_list)
      }

      plot_table <- lapply( names(out), function(x){
        data.frame(Group=x, data.table::melt.data.table( data.table(out[[x]]), id.vars = "PCs"), check.names = F)
      })
      plot_table <- do.call(rbind, plot_table)
      plot_table <- split(plot_table, plot_table$PCs)
      if("plot" %in% return.type){
        plot_list <- lapply(plot_table, function(plot_data){
          scientific_theme1 <- function(base_size = 12, base_family = "") {
            ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
              ggplot2::theme(
                text = ggplot2::element_text(color = "black"),
                # axis.title = element_text(size = rel(1.2), face = "bold"),
                axis.title.y = ggplot2::element_blank(),
                axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
                # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                axis.text.y = ggplot2::element_text(hjust = 1),
                axis.line = ggplot2::element_line(color = "black"),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.border = ggplot2::element_blank(),
                strip.background = ggplot2::element_rect(fill = "white", color = "black"),
                strip.text = ggplot2::element_text(face = "bold")
                # ,legend.position = "none"
              )}
          if(is.null(color)){
            color <- get_colors(length(unique(plot_data[["variable"]])))
          }
          p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(y = Group, x = value, fill = variable)) +
            ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.9)) +
            ggplot2::facet_wrap(~ variable, scales = "free_x", ncol =5 ) +
            ggplot2::scale_fill_manual(values = color ) +
            scientific_theme1() +
            ggplot2::labs(
              x = "% variance explained",
              subtitle = unique(plot_data$PCs)
            )

          return(p1)
        })
        names(plot_list) <- names(plot_table)
        out_list$explanatory$plot <- plot_list
      }

      if("interactive_table" %in% return.type){
        out_list$explanatory$interactive_table <- lapply( names(plot_table), function(x){
          in_table <- data.table::dcast(data.table(plot_table[[x]]), variable~Group, value.var = "value")
          re_table <- .re_table(in_table, csv.name =paste0(assay, "_explanatory",x) ,elementId = paste0(assay, "_explanatorytable_",x), down.sparkline = T, right.sparkline = T, maxWidth = NULL, subtitle = paste0("Explanatory PCs for each variable (",x,")"),
                                number_type = "custom", number_fmr = "paste0( round(value,2),'%')")
          return(re_table)
        })
        names(out_list$explanatory$interactive_table) <- names(plot_table)
      }

      if(length(out_list$explanatory)==1){
        out_list$explanatory <- out_list$explanatory[[1]]
      }
    }

  }
  if(length(out_list)==1){
    out_list <- out_list[[1]]
  }

  if (dir.exists(tmpdir)) {
    unlink(tmpdir, recursive = TRUE)
  }

  return(out_list)
}


#' @title Visualize Covariate Impact on Single-Cell Data
#' @description This function evaluates the impact of specified covariates on pseudobulk level. It supports PCA-based analyses (`"pseudobulk"`).
#'
#' @param object A Seurat object containing the single-cell data to be analyzed.
#' @param assay A character string specifying the assay to use for PCA (default: `"RNA"`).
#' @param variables A character vector specifying the metadata columns to use for variance explained analysis (default: `NULL`).
#' @param pseudobulk.sample.by A character string specifying the metadata column used for defining samples in the pseudobulk analysis (default: `"orig.ident"`).
#' @param return.type A character string specifying the output format. Options include `"plot"` and `"interactive_table"` (default: `"plot"`).
#' @param pseudobulk.celltype.by A character string specifying the metadata column used to split the data into subsets in the pseudobulk analysis (default: `NULL`).
#' @param color A vector of colors to use for visualizations (default: `NULL`).
#' @param ... Additional arguments to be passed to the plotting functions.
#' @return A plot or interactive table based on the specified `return.type`.
#' @export
#'
#'
#' @seealso  \code{\link{RunVarExplained}}, \code{\link{PlotVEPerFeature}}, \code{\link{PlotSamplePCA}}
#'
PlotCovariateImpact <- function(object, assay="RNA", variables=NULL, pseudobulk.sample.by="orig.ident",return.type="plot",
                                pseudobulk.celltype.by = NULL,
                                color=NULL ,...){
  if( !("Seurat" %in% is(object)) ){
    stop("Error: Input must be Seurat object.")
  }

  out <- PlotSamplePCA(object, assay = assay, svd.variable.by = variables, sample.by =pseudobulk.sample.by, celltype.by = pseudobulk.celltype.by,plot.type = c("pca", "svd", "explanatory"),
                       color = color, return.type = return.type,... )
  return(out)

}


#' @title Run Batch Test
#' @description This function performs a batch test on a Seurat object to evaluate the presence of batch effects. It calculates the Kendall W statistic, pairwise Kendall correlations, and the number of negative correlations between samples.
#' The function also checks for severe batch effects based on the results.
#'
#' @param object A Seurat object containing the single-cell data to be analyzed.
#' @param sample.by A character string specifying the metadata column used for defining samples (default: `"orig.ident"`).
#' @param cluster.by A character string specifying the metadata column used for grouping cells (default: `"seurat_clusters"`).
#' @param k An integer specifying the number of clusters to use for batch testing (default: `3`).
#' @param seed An integer specifying the seed for random sampling (default: `1`).
#' @param n An integer specifying the number of samples to use for batch testing (default: `100`).
#' @return A list containing the correlation matrix, Kendall W results, and the qbinom results.
#' @export
RunBatchTest <- function(object, sample.by="orig.ident", cluster.by="seurat_clusters", k=3, seed=1, n=100 ){
  test_data <- data.table(object@meta.data)[, .(count = .N), by = .(sample = get(sample.by), cluster = get(cluster.by))]
  test_data <- dcast(test_data, cluster ~ sample, value.var = "count", fill = 0)
  test_data <- as.data.frame(test_data)
  rownames(test_data) <- test_data[,1]
  test_data <- test_data[,-1]
  cor_data <- stats::cor(test_data,method ="kendall")
  dist_matrix <- stats::as.dist(1-cor_data)
  hc <- stats::hclust(dist_matrix)
  clusters <- stats::cutree(hc, k = k)

  k_list <- lapply(1:k, function(x){
    test_data[,colnames(test_data) %in% names(which(clusters==x)), drop=F]
  })

  sampled_list <- lapply(1:k, function(x){
    set.seed(seed)
    sampled_chars <- sample(colnames(k_list[[x]]), size = n, replace = TRUE)
    test_data[, match(sampled_chars, colnames(test_data)), drop=F]
  })

  kendallW_result <- lapply(1:n, function(x){
    merge_data <- lapply(sampled_list, function(y){
      y[, x,drop=F]
    })
    merge_data <- do.call(cbind, merge_data)
    merge_data <- merge_data[rowSums(merge_data) > 0, ]
    kendallW<-irr::kendall(merge_data)

    cor_data <- apply(utils::combn(dim(merge_data)[2], 2), 2, function(x){
      cor_temp <- merge_data[, x]
      cor_temp <- cor_temp[rowSums(cor_temp) > 0, ]
      stats::cor( as.numeric(cor_temp[,1]),  as.numeric(cor_temp[,2]), method ="kendall")
    })
    c(kendallW$value, kendallW$p.value, sum(cor_data<0) )
  })
  kendallW_result <- do.call(rbind, kendallW_result)
  colnames(kendallW_result) <- c("W", "pvalue", "InterNegNum")
  kendallW_result <- as.data.frame(kendallW_result)
  kendallW_result$Warning <- unlist(apply(kendallW_result, 1, function(x){
    (x[2] > 0.1) & (x[3] >= 1)
  }))

  if(sum(kendallW_result$Warning) >= stats::qbinom(0.05,n, 0.5)){
    message("The dataset may have severe clustering differences, which may be batch effect !!!")

  }
  return(list(cor_matrix=cor_data, kendallW_result=kendallW_result, qbinom_result=data.frame(sum=sum(kendallW_result$Warning), ntimes=n, qbinom_min= stats::qbinom(0.05,n, 0.5) ) ))
}






