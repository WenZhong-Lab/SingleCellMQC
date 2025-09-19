# ##############################################################################
# #
# #  VariancePartition (cell metadata)
# #
# ##############################################################################

#' Run Variance Partitioning on Cell Metadata
#'
#' @description
#' This function performs variance partition analysis on specified cell metadata
#' features using linear mixed models. It can analyze the entire dataset or
#' split it based on a grouping variable.
#'
#' @param object An object containing cell metadata.
#' @param features A character vector specifying the names of the metadata
#'   features to be analyzed. Common examples include "percent.mt", "percent.rb",
#'   "percent.hb", "percent.dissociation".
#' @param formula A `formula` object specifying the fixed and random effects
#'   to be included in the linear mixed model. This formula should only contain
#'   the explanatory variables (e.g., `~ (1|Batch)`). The
#'   response variable will be dynamically added within the function for each
#'   `feature`.
#' @param split.by (Optional) A character string specifying the name of a
#'   column in the metadata by which to split the data for separate variance
#'   partition analyses. If `NULL` (default), the analysis is performed on
#'   the entire object.
#' @param do.log1p A logical value indicating whether to apply a `log1p`
#'   transformation (i.e., `log(1 + x)`) to the specified `features` before
#'   performing the variance partition. Defaults to `TRUE`.
#'
#' @return A data frame summarizing the variance explained by each
#'   component in the `formula` for each `feature` and (if `split.by` is used)
#'   for each group. The data frame includes columns for 'Group', 'Feature',
#'   and columns for each variance component.
#'
#' @details
#' This function leverages `lmerTest::lmer` to fit linear mixed models and
#' `variancePartition::calcVarPart` to calculate the fraction of variance
#' explained by each covariate. If `split.by` is provided, the function
#' iterates through each unique value in the specified column, performs the
#' variance partition for that subset of data, and combines the results.
#'
#'
#'
#' @importFrom stats as.formula update
#' @importFrom lmerTest lmer
#' @importFrom variancePartition calcVarPart
#' @importFrom gtools smartbind
#' @export
#'
RunVarPartMetadata <- function(object,
                               features,
                               formula,
                               split.by =NULL,
                               do.log1p=T){
  object <- getMetaData(object)
  if(is.null(split.by)){
    out <- calMetadataVarPart(object, formula= formula, group_name = "No split.by",features=features, do.log1p=do.log1p )
  }else{
    split_object <- split(object, object[[split.by]])
    index_name <- names(split_object)
    out <- lapply(index_name, function(x){
      calMetadataVarPart(split_object[[x]], formula= formula, group_name = x,features=features, do.log1p=do.log1p )
    })
    out <- do.call(gtools::smartbind, out)
  }
  rownames(out) <- NULL
  return(out)
}

calMetadataVarPart <- function(object,
                               formula,
                               group_name,
                               cluster_name,
                               do.log1p=T,
                               features = c("percent.mt", "percent.rb", "percent.hb", "percent.dissociation")){

  metadata = getMetaData(object)
  in_name <- intersect(colnames(metadata), features )
  data = metadata[, c(in_name), drop=F ]
  if(do.log1p){
    data = log1p(data)
  }
  data <- data.frame(data, metadata[, cluster_name, drop=F])


  var_list <- lapply(in_name, function(x){
    formula <- stats::update(formula, as.formula( paste0(x, " ~ .")) )
    print(formula)
    fm <- lmerTest::lmer(formula, data = data)
    varPart <- variancePartition::calcVarPart(fm)
  })
  var_list <- do.call(rbind, var_list)
  rownames(var_list) <- in_name

  # varPart <- variancePartition::fitExtractVarPartModel(data, formula, metadata)
  varPart <- data.frame(Group = group_name, Feature=rownames(var_list), data.frame(var_list))
  print(11111)

  return(varPart)
}



# ##############################################################################
# #
# #  VariancePartition (Pseudobulk)
# #
# ##############################################################################

#' Run Variance Partitioning on Pseudobulk Data
#'
#' @description
#' This function performs variance partitioning analysis on pseudobulk
#' data. It can analyze either a single pseudobulk matrix or a list of matrices
#' (e.g., grouped by cell type), using `dreamlet` and `variancePartition` packages.
#'
#' @param object A list object containing pseudobulk data gets from `RunPseudobulkData` function.
#' @param variables (Optional) A character vector of column names from
#'   `pseudobulk_metadata` to be used as random effects in the variance
#'   partition model. If provided, the formula will be constructed using `(1|variable)`.
#'   Cannot be used with `formula`.
#' @param formula (Optional) A `formula` object specifying the fixed and random
#'   effects for the variance partition model (e.g., `~ (1|donor_id)`).
#'   If provided, `variables` is ignored. If neither `variables` nor `formula`
#'   are provided, an error will be thrown.
#'
#' @return A data frame summarizing the variance explained by each
#'   component in the `formula` for each gene. The data frame includes columns
#'   for 'Group' (indicating if multiple groups were processed), 'Feature' (gene name),
#'   and columns for each variance component.
#'
#' @details
#' This function internally uses `dreamlet::processOneAssay` to prepare the data
#' and `variancePartition::fitExtractVarPartModel` to perform the actual
#' variance partitioning.
#' If `pseudobulk_mat` is a list, the function iterates through each element,
#' performing variance partitioning independently for each group/cluster.
#'
#' @seealso
#' \code{\link[dreamlet]{processOneAssay}}, \code{\link[variancePartition]{fitExtractVarPartModel}}
#'
#' @importFrom stats as.formula
#' @importFrom dreamlet removeConstantTerms dropRedundantTerms
#' @importFrom variancePartition fitExtractVarPartModel
#' @importFrom gtools smartbind
#' @export
RunVarPartPseudobulk <- function(object, variables=NULL, formula=NULL){
  # Internal function to run variance partitioning for one assay
  run_varpart <- function(exp, metadat, ncells_vec, group_name) {
    metadat <- metadat[match(colnames(exp), rownames(metadat)), , drop=FALSE]
    ncells_vec <- ncells_vec[["ncells"]][match(colnames(exp), ncells_vec[,1])]

    if(!is.null(formula)){
      # if formula is not null
      formula <- stats::as.formula(formula)
    } else if(!is.null(variables)){
      # Build formula string from variables
      formula <- stats::as.formula(paste("~", paste0("(1|", variables, ")", collapse = " + ")))
    } else {
      stop("Either 'variables' or 'formula' must be provided!")
    }

    message(paste0(format(Sys.time(), "%H:%M:%S"), "-------- Processing expression data: ", group_name))
    print(formula)
    process_data <- dreamlet:::processOneAssay(exp, formula = formula, data = metadat, n.cells = ncells_vec)
    message(paste0(format(Sys.time(), "%H:%M:%S"), "-------- Variance partition: ", group_name))
    varPart_metadata <- metadat[match(colnames(process_data), rownames(metadat)), , drop=FALSE]
    formula <- dreamlet::removeConstantTerms(formula, varPart_metadata)
    formula <- dreamlet::dropRedundantTerms(formula, varPart_metadata)
    varPart <- variancePartition::fitExtractVarPartModel(process_data, formula, varPart_metadata)
    varPart <- data.frame(Group = group_name, Feature=rownames(varPart), data.frame(varPart))
    return(varPart)
  }

  pseudobulk_mat <- object[["pseudobulk_mat"]]
  metadat <- object[["pseudobulk_metadata"]]
  ncells_info <- object[["ncells"]]

  if(is.list(pseudobulk_mat)){
    # If multiple assays, process each one
    out <- lapply(names(pseudobulk_mat), function(x) {
      exp <- pseudobulk_mat[[x]]

      # Subset ncells for the current assay
      ncells_sub <- ncells_info[ncells_info[,2] %in% x, , drop=FALSE]

      run_varpart(exp, metadat, ncells_sub, group_name = x)
    })
    # Merge results across assays
    out <- do.call(gtools::smartbind, out)
  } else {
    # If only single assay
    ncells_sub <- ncells_info
    out <- run_varpart(pseudobulk_mat, metadat, ncells_sub, group_name = "No cluster")
  }
  rownames(out) <- NULL
  return(out)

}


#' Run Variance Partitioning on Pseudobulk PCA embeddings
#'
#' @description
#' This function performs variance partitioning analysis on the Principal Component
#' (PC) embeddings derived from pseudobulk expression data. It allows identifying
#' the factors contributing to the variance in the PCs.
#'
#' @param object A list object containing pseudobulk data gets from `RunPseudobulkData` function.
#' @param variables (Optional) A character vector of column names from
#'   `pseudobulk_metadata` to be used as random effects in the variance
#'   partition model. If provided, the formula will be constructed using `(1|variable)`.
#'   Cannot be used with `formula`.
#' @param nPCs An integer specifying the number of top principal components
#'   to analyze. Defaults to 2. The actual number of PCs analyzed will be
#'   the minimum of `nPCs` and the available number of PCs.
#' @param formula (Optional) A `formula` object specifying the fixed and random
#'   effects for the variance partition model (e.g., `~  (1|donor_id)`).
#'   If provided, `variables` is ignored. If neither `variables` nor `formula`
#'   are provided, an error will be thrown.
#'
#' @return A data frame summarizing the variance explained by each
#'   component in the `formula` for each PC. The data frame includes columns
#'   for 'Group' (indicating if multiple groups were processed), 'Feature' (PC name, e.g., 'PC_1'),
#'   and columns for each variance component.
#'
#' @details
#' This function first performs PCA on the pseudobulk expression data using
#' an internal `runPseudobulkPCA` function (assumed to return a Seurat-like
#' object). It then extracts the specified number of top PCs and treats each
#' PC's projection as a "feature" for a linear mixed model.
#' `lmerTest::lmer` is used to fit these models, and `variancePartition::calcVarPart`
#' calculates the variance explained by each covariate. If `pseudobulk_mat` is a
#' list, PCA and variance partitioning are performed independently for each
#' group/cluster.
#'
#'
#' @importFrom stats as.formula update
#' @importFrom lmerTest lmer
#' @importFrom variancePartition calcVarPart
#' @importFrom dreamlet removeConstantTerms
#' @importFrom Seurat Embeddings
#' @importFrom gtools smartbind
#' @export
RunVarPartPseudobulkPCA <- function(object, variables=NULL, nPCs=2, formula=NULL){
  mat <- object$pseudobulk_mat
  meta_data <- object$pseudobulk_metadata

  if(!is.null(formula)){
    # if formula is not null
    formula <- stats::as.formula( formula )
  } else if(!is.null(variables)){
    # Build formula string from variables
    formula <- stats::as.formula(paste("~", paste0("(1|", variables, ")", collapse = " + ")))
  } else {
    stop("Eitherd 'variables' or 'formula' must be provided!")
  }

  if(is.list(mat)){
    group_name <- names(mat)
    plot_list <- lapply(group_name, function(x){
      message(paste0(format(Sys.time(), "%H:%M:%S"), "-------- runPseudobulkPCA: ", x))
      seu <- runPseudobulkPCA(mat[[x]], metadata = meta_data)
      ## subset pca embeddings & metadata
      pca_res <- Seurat::Embeddings(seu, "pca")
      nPCs <- min(ncol(pca_res),nPCs)
      pca_res_meta <- meta_data[match(rownames(pca_res), rownames(meta_data)), , drop=F]
      formula <- dreamlet::removeConstantTerms(formula, pca_res_meta)
      formula <- stats::update(formula, as.formula("feature ~ ."))
      #formula <- dreamlet::dropRedundantTerms(formula, varPart_metadata)
      ## pca varpart
      # lmer
      out <- lapply(1:nPCs, function(y){
        res_temp <- data.frame(feature=pca_res[, y], pca_res_meta)
        fm <- lmerTest::lmer(formula, data = res_temp)
        varPart <- variancePartition::calcVarPart(fm)
      })
      out <- do.call(rbind, out)
      out <- data.frame(Group=x, Feature=colnames(pca_res)[1:nPCs], out)
      return(out)
    })
    names(plot_list) <- NULL
    plot_list <- do.call(gtools::smartbind, plot_list)
    return(plot_list)
  }else{
    message(paste0(format(Sys.time(), "%H:%M:%S"), "-------- runPseudobulkPCA "))
    seu <- runPseudobulkPCA(mat, metadata = meta_data)
    ## subset pca embeddings & metadata
    pca_res <- Seurat::Embeddings(seu, "pca")
    nPCs <- min(ncol(pca_res),nPCs)
    pca_res_meta <- meta_data[match(rownames(pca_res), rownames(meta_data)), , drop=F]
    formula <- dreamlet::removeConstantTerms(formula, pca_res_meta)
    formula <- stats::update(formula, as.formula("feature ~ ."))
    ## pca varpart
    # lmer
    out <- lapply(1:nPCs, function(y){
      res_temp <- data.frame(feature=pca_res[, y], pca_res_meta)
      fm <- lmerTest::lmer(formula, data = res_temp)
      varPart <- variancePartition::calcVarPart(fm)
    })
    out <- do.call(rbind, out)
    out <- data.frame(Group="No cluster", Feature=colnames(pca_res)[1:nPCs], out)
    return(out)
  }
}

#' Plot Variance Partition Results as Violin Plots
#'
#' @description
#' Generates violin plots to visualize the proportion of variance explained by
#' different factors across features (genes, PCs, or metadata features).
#'
#' @param object A data frame, typically the output of `RunVarPartMetadata`,
#'   `RunVarPartPseudobulk`, or `RunVarPartPseudobulkPCA`. It should contain
#'   at least 'Feature' and 'Group' columns (if applicable), and columns
#'   representing the variance components (e.g., 'donor_id', 'batch').
#' @param color (Optional) A character vector of colors to be used for filling
#'   the violin plots. If `NULL`, default colors are assigned.
#' @param do.split A logical value. If `TRUE`, it expects the `object` to
#'   already contain a 'Group' column and will facet the plots by this 'Group'
#'   column. If `FALSE` and multiple groups are present, it will generate a
#'   separate plot for each group.
#'
#' @return A `ggplot` object (or a list of `ggplot` objects if `do.split` is `FALSE`
#'   and multiple groups are present) representing the violin plots of variance explained.
#'
#'
#'
#' @importFrom data.table data.table melt
#' @importFrom ggplot2 ggplot geom_violin labs theme element_text element_blank facet_wrap
#' @importFrom ggplot2 coord_flip scale_y_reverse scale_fill_manual
#' @export
PlotVarPartVln <- function(object, color=NULL, do.split=F ){
  if(!do.split){
    group_len <- length(unique(object$Group))
    if(group_len!=1){
      object <- split(object, object$Group)
      index_name = names(object)
      out <- lapply(index_name, function(x){
        plotVarPartVln(object[[x]], color = color )
      })
      names(out) <- index_name
    }else{
      out <-  plotVarPartVln(object, color = color )
    }
  }else{
    out <-  plotVarPartVln(object, color = color, do.split=T )

  }
  return(out)

}

plotVarPartVln <- function(object, color=NULL, do.split=F ){
  if(!do.split){
    rsquared_mat <- data.table::data.table(object[, -1, drop=F])
    rsquared_long <- data.table::melt(rsquared_mat,id.vars="Feature")
  }else{
    rsquared_mat <- data.table::data.table(object)
    rsquared_long <- data.table::melt(rsquared_mat,id.vars=c("Feature", "Group") )
  }

  rsquared_long$value <- rsquared_long$value *100
  if(is.null(color)){
    color <- get_colors(length(unique(rsquared_long$variable))-1)
    color <- c(color, "#D3D3D3")
  }

  if(!do.split){
  varplot <- plotVln(object=rsquared_long, x="variable", y="value", log.y = F,combine = F,  group.by="variable", color = color)[[1]] +
    ggplot2::theme(legend.position = "none")+
    ggplot2::labs(y = "Variance explained (%)")
  }else{
    varplot <- plotVln(object=rsquared_long, x="variable", y="value", log.y = F,combine = F,  group.by="variable", color = color)[[1]] +
      ggplot2::theme(legend.position = "none")+
      ggplot2::labs(y = "Variance explained (%)")  +
      ggplot2::facet_wrap(~ Group,  ncol  = 3)

  }


  return(varplot)
}


#' Plot Variance Partition Results as Stacked Bar
#'
#' @description
#' Generates stacked bar charts to visualize the proportion of variance explained
#' for the top features, ordered by specified covariates.
#'
#' @param object A data frame, typically the output of `RunVarPartMetadata`,
#'   `RunVarPartPseudobulk`, or `RunVarPartPseudobulkPCA`. It should contain
#'   at least 'Feature' and 'Group' columns (if applicable), and columns
#'   representing the variance components (e.g., 'donor_id', 'batch').
#' @param sort.by (Optional) A character vector of column names (variance components)
#'   by which to sort the features. Features will be sorted descending by the
#'   sum of variance explained by these components. If `NULL`, features are
#'   plotted in their original order.
#' @param ntop An integer specifying the number of top features to display.
#'   Defaults to 10.
#' @param split.by A character string specifying a column in `object` to split
#'   the plotting by. Typically 'Group' when results from multiple groups are combined.
#'   Defaults to "Group".
#' @param color (Optional) A character vector of colors to be used for filling
#'   the stacked bars. If `NULL`, default colors are assigned.
#'
#' @return A `ggplot` object (or a list of `ggplot` objects if `split.by` leads to
#'   multiple facets/plots) representing the stacked bar charts of variance explained.
#'
#'
#' @importFrom data.table data.table melt
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text element_blank facet_wrap
#' @importFrom ggplot2 coord_flip scale_y_reverse scale_fill_manual theme_classic
#' @export
PlotVarPartStackBar <- function(object,
                                sort.by = NULL,
                                ntop = 10,
                                split.by = "Group",
                                color = NULL){

  group_len <- length(unique(object[[split.by]]))
  if(group_len!=1){
    index <- setdiff(c("Group", "Feature"),  split.by)
    object <- split(object, object[[split.by]])
    index_name = names(object)
    out <- lapply(index_name, function(x){
      object[[x]]$Feature <-object[[x]][[index]]
      plotVarPartStackBar(object[[x]], sort.by=sort.by, ntop=ntop, color = color )
    })
    names(out) <- index_name
  }else{
    out <-  plotVarPartStackBar(object, sort.by=sort.by, ntop=ntop, color = color )
  }
  return(out)
}


plotVarPartStackBar <- function(object,
                                sort.by = NULL,
                                ntop = 10,
                                color = NULL) {
  # Make sure object is a data.frame
  df <- object[, -1, drop=F]

  # Check column existence if sorting required
  if (!is.null(sort.by)) {
    if (!all(sort.by %in% colnames(df))) {
      stop("Some sort.by columns do not exist in the input object.")
    }
    ord <- do.call(order, c(df[sort.by], decreasing = TRUE))
    top_features <- df$Feature[ord][1:min(ntop, nrow(df))]

  } else {
    top_features <- utils::head(df$Feature, ntop)
  }

  # Filter for top_features
  df_top <- df[df$Feature %in% top_features, ]

  # Set Feature as factor to control facet order
  df_top$Feature <- factor(df_top$Feature, levels = rev(top_features))

  # Melt to long for ggplot
  df_top <-data.table::data.table(df_top)
  df_long <- data.table::melt(df_top, id.vars = "Feature")
  df_long$value <- df_long$value * 100
  df_long <- na.omit(df_long)
  if(is.null(color)){
    color <- get_colors(length(unique(df_long$variable))-1)
    color <- c(color, "#D3D3D3")
  }

  df_long$variable <- factor(df_long$variable, levels = unique(df_long$variable))


  plt <- ggplot2::ggplot(df_long, ggplot2::aes(x = Feature, y = value, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.9) +
    ggplot2::labs(y = "Variance explained (%)", y = NULL) +
    ggplot2::scale_y_reverse(breaks = seq(0, 100, by = 20), label = seq(100, 0, by = -20)) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      axis.title.y = ggplot2::element_blank()
    )+
    ggplot2::coord_flip()+
    ggplot2::scale_fill_manual(values = color, name = NULL)

  return(plt)
}



# ##############################################################################
# #
# #  VariancePartition (Cell single variable)
# #
# ##############################################################################


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






