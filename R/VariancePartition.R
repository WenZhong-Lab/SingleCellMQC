


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


PlotVarPartVln <- function(object, color=NULL ){
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
  return(out)
}

plotVarPartVln <- function(object, color=NULL ){
  rsquared_mat <- data.table::data.table(object[, -1, drop=F])
  rsquared_long <- data.table::melt(rsquared_mat,id.vars="Feature")
  rsquared_long$value <- rsquared_long$value *100
  if(is.null(color)){
    color <- get_colors(length(unique(rsquared_long$variable))-1)
    color <- c(color, "#D3D3D3")
  }
  varplot <- plotVln(object=rsquared_long, x="variable", y="value", log.y = F,combine = F,group.by="variable", color = color)[[1]] +
    ggplot2::theme(legend.position = "none")+
    ggplot2::labs(y = "Variance explained (%)")
  return(varplot)
}



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
    top_features <- head(df$Feature, ntop)
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
