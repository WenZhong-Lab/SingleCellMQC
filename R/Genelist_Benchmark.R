


# add_custom_features will add hvg features =  nfeatures + add_custom_features
# exclude_custom_features & exclude_patterns will contain always hvg features = nfeatures
stepRNAToPCA <- function(object,
                       assay = "RNA",
                       slot="counts",
                       layer=NULL,
                       nfeatures=2000,
                       exclude_custom_features=NULL,
                       add_custom_features = NULL,
                       exclude_patterns=NULL,
                       threads=NULL,
                       tmpdir ="./temp_SingleCellMQC/pca/",
                       output_pca_name="rna.pca",
                       vars.to.regress=NULL,
                       only.returnPCA =TRUE,
                       add.norm_data=FALSE,
                       add.scale_data=FALSE) {

  mat <- getMatrix(object,assay=assay, slot=slot, layer=layer)
  metadata <- getMetaData(object)
  feature_name <- rownames(mat)

  if( !("BPCells" %in% attr(class(mat), "package"))){
    mat <- as(mat, "IterableMatrix")
  }

  message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- Log Normalization"))
  # Log normalization
  mat <- BPCells::multiply_cols(mat, 1/Matrix::colSums(mat))
  mat <- log1p(mat * 10000)

  message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- VST & Select features"))
  vst_data <- Seurat::VST(mat)
  var_index <- order(vst_data$variance.standardized, decreasing=TRUE)
  var_feature <- feature_name[var_index]
  del_genes <- character(0)
  if (!is.null(exclude_patterns)) {
    del_genes <- union(del_genes, grep(paste0(exclude_patterns, collapse = "|"), feature_name, value = TRUE))
  }
  if (!is.null(exclude_custom_features)) {
    del_genes <- union(del_genes, exclude_custom_features)
  }
  mask <- !(var_feature %in% del_genes)
  nfeatures <- min(length(var_feature), nfeatures)
  var_index_top <- var_index[mask][1:nfeatures]

  mask_add <- var_feature %in% add_custom_features

  var_index_top <- union(var_index_top, var_index[mask_add])
  var_index_top <- na.omit(var_index_top)

  message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- nFeatures: ", length(var_index_top)))

  # scale
  message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- Scale"))
  # variance
  mat_norm <- mat[var_index_top,]
  mat_norm1 <- mat_norm

  # save temp
  message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- Save Tmp to ", tmpdir))
  mat_norm <- mat_norm %>% BPCells::write_matrix_dir(tempfile("mat", tmpdir =tmpdir  ))

  if(!is.null(vars.to.regress)){
    message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- vars.to.regress ", vars.to.regress))
    mat_norm <- BPCells::regress_out(mat_norm, metadata[, vars.to.regress,drop=F])
  }

  if(is.null(threads)){
    stats <- BPCells::matrix_stats(mat_norm, row_stats="variance")
  }else{
    stats <- BPCells::matrix_stats(mat_norm, row_stats="variance", threads = threads)
  }
  features.mean <- stats$row_stats['mean',]
  features.sd <- sqrt(stats$row_stats['variance',])
  features.sd[features.sd == 0] <- 0.01
  mat_norm <- BPCells::min_by_row(mat = mat_norm, vals = 10 * features.sd + features.mean)
  mat_norm <- (mat_norm - features.mean) / features.sd

  message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- PCA"))
  if(is.null(threads)){
    k <- min(nrow(mat_norm), 50)
    set.seed(1)
    svd <- BPCells::svds(mat_norm, k=k)
  }else{
    k <- min(nrow(mat_norm), 50)
    set.seed(1)
    svd <- BPCells::svds(mat_norm, k=k, threads=threads)
  }
  pca <- BPCells::multiply_cols(svd$v, svd$d)
  colnames(pca) <- paste0("PC_", 1:ncol(pca))
  rownames(pca) <- colnames(mat)


  ## unlink tmp
  message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- Unlink Tmp to ",.findDirValue(mat_norm)))
  unlink(.findDirValue(mat_norm),recursive = TRUE)

  if(only.returnPCA){
    return(pca)
  }


  ## save data
  if("Seurat" %in% class(object)){
    pca_dim <- Seurat::CreateDimReducObject(embeddings =pca,  key =  paste0(assay, "PC_"),assay = assay)
    if(add.norm_data){
      object <-setAssayData(object = object,assay = assay,slot = "data",new.data = mat)
    }
    if(add.scale_data){
      mat_norm1 <- BPCells::min_by_row(mat = mat_norm1, vals = 10 * features.sd + features.mean)
      mat_norm1 <- (mat_norm1 - features.mean) / features.sd
      object <- setAssayData(object = object,assay = assay,slot = "scale.data",new.data = mat_norm1)
    }
    object[[output_pca_name]] <- pca_dim
    return(object)
  }else{
    return(pca)
  }

}



stepPCAToCluster <- function(object, assay="RNA", pca_name="rna.pca", output_cluster_name="rna_cluster", threads = 1,
                               resolution=1 ,min_val=0.05, dims=1:20){
  if("Seurat" %in% class(object)){
    pca <-  Seurat::Embeddings(object, pca_name)
  }else{
    pca <- object
  }
  max_dims <- ncol(pca)
  if(max(dims) > max_dims) {
    dims <- 1:max_dims
    message(sprintf("Using all available PCs (1:%d) as input dimensions", max_dims))
  }
  pca <- pca[,dims]

  message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- knn_hnsw & knn_to_snn_graph"))
  knn_result <-  BPCells::knn_hnsw(pca, ef=500, threads=threads) %>% # Find approximate nearest neighbors
    BPCells::knn_to_snn_graph( min_val = min_val) # Find shared nearest neighbors

  cluster_out <- lapply(resolution, function(x){
    message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- cluster_graph_leiden with resolution ", x))
    temp <- knn_result %>% # Convert to a SNN graph
      BPCells::cluster_graph_leiden(resolution = x) # Perform graph-based clustering
    temp <- as.character(temp)
    return(temp)
  })
  cluster_out <- do.call(cbind, cluster_out)

  if(length(output_cluster_name)!= length(resolution) ){
    colnames(cluster_out) <- paste0(output_cluster_name,"_", resolution)
  }else{
    colnames(cluster_out) <- output_cluster_name
  }
  rownames(cluster_out) <- rownames(pca)

  if("Seurat" %in% class(object)){
    object <- Seurat::AddMetaData(object, cluster_out)
    object$seurat_clusters <- cluster_out[,dim(cluster_out)[2], drop=T]
    return(object)
  }else{
    return(data.frame(cluster_out))
  }
}




#' Plot Variance Explained by Noise Gene Cluster
#'
#' @description
#' This function visualizes the proportion of variance explained by different
#' factors. It generates
#' a stacked bar chart (or dodged bar chart in this specific re-implementation).
#'
#' @param metrics An object containing variance fraction metrics, typically
#'   the `Variance_fractions` component from the output of `CalNoiseGeneClusterMetrics` function.
#'
#' @return A `ggplot` object representing the bar plot of variance explained.
#'
#' @importFrom data.table data.table melt
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual labs theme_classic
#'   coord_flip theme element_text element_blank
#' @importFrom RColorBrewer brewer.pal
#' @export
PlotNoiseGeneClusterMetricsVar <- function(metrics ){
  metrics = metrics$Variance_fractions
    df <- data.table(metrics)
    df <- melt(df[,-1], id.vars=c("Feature"))
    df <- na.omit(df)
    df <- df[!(df$variable %in% "Residuals"), ]
    df$value <-df$value*100
    p1 <- ggplot(df, aes(x = Feature, y = value, fill = variable)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) + # 关键改变: position = "dodge"
      scale_fill_manual(values = get_colors(6)) + # 应用自定义颜色
      labs(
        x = "",
        y = "Variance explained (%)", # y轴标签可以根据你的数据含义修改，这里是"Value"
        fill = ""
      ) +
      theme_classic(base_size = 13)+
      ggplot2::coord_flip()+
      theme(
       # axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), # x轴文字旋转45度，右对齐，颜色黑色
        axis.text = element_text(color = "black"), # y轴文字颜色黑色
        axis.title.x = element_text(color = "black"), # x轴标题颜色黑色
        legend.text = element_text(color = "black"), # 图例文字颜色黑色
        plot.title = element_text(color = "black"), # 如果有主标题，也设置为黑色
        legend.position = "right", # 图例位置
      )
    return(p1)
}


#' Plot Gini Coefficient Differences for Noise Genes
#'
#' @description
#' This function visualizes the differences in Gini coefficients before and
#' after removing specific noise gene sets. It generates bar plots for the
#' top features (genes) with the largest Gini coefficient differences, and
#' accompanying dot plots showing the expression patterns of these genes across
#' clusters.
#'
#' @param object A Seurat object containing single-cell expression data
#'   and metadata.
#' @param metrics An object containing Gini coefficient differences, typically
#'   the `Diff_Gini` component from the output of `CalNoiseGeneClusterMetrics`.
#' @param cluster_res A data frame or matrix containing the clustering results from `RunNoiseGeneCluster` function.
#' @param ntop An integer specifying the number of top features (genes) to display
#'   in the bar plot based on Gini difference. Defaults to 10.
#' @param color.bar A character string specifying the color for the Gini
#'   difference bar plot. Defaults to "#1F78B4".
#'
#' @return A named list of lists. Each top-level element corresponds to a
#'   noise gene set (e.g., "mt_add"). Each inner list contains:
#'   \itemize{
#'     \item `diff_gini`: A `ggplot` object of the Gini difference bar plot.
#'     \item `add_dotplot`: A `ggplot` object of the DotPlot for features in
#'       clusters *after* adding back the specific noise gene set.
#'     \item `del_dotplot`: A `ggplot` object of the DotPlot for features in
#'       clusters *after* removing *all* noise genes (from `GetNoiseGeneList` function).
#'   }
#'
#' @importFrom data.table data.table melt
#' @importFrom ggplot2 ggplot aes geom_bar theme_classic coord_flip labs theme element_text
#' @importFrom stringr str_extract
#' @export
PlotNoiseGeneClusterMetricsGini <- function(object, metrics, cluster_res,  ntop=10, color.bar="#1F78B4"){
  metrics = metrics$Diff_Gini
  metrics_list <- split(metrics, metrics$variable)
  p_list <- lapply(names(metrics_list), function(x){
    df <- data.frame(metrics_list[[x]])
    df <- df[order(df$value, decreasing = T)[1:ntop], ]
    df$Feature <- factor(df$Feature, levels = rev(df$Feature))
    bar_p <- ggplot2::ggplot(df, ggplot2::aes(x=Feature, y=value)) +
      ggplot2::geom_bar(stat="identity", fill=color.bar) +
      ggplot2::theme_classic(base_size = 13)+
      ggplot2::coord_flip()+
      ggplot2::labs( x="", y="Gini(add) - Gini(del)")+
      ggplot2::theme(
        axis.text = ggplot2::element_text(color = "black")
      )

    ## dotplot
    object <- Seurat::AddMetaData(object, metadata = cluster_res)
    middle <- stringr::str_extract(x, "(?<=add_)[^_]+")
    add <- grep(middle,  colnames(cluster_res), value = T)
    del <- grep("del_all",  colnames(cluster_res), value = T)
    suppressMessages({add_dotplot <- Seurat::DotPlot(object, features = c( rev(as.character(df$Feature))), group.by = add)+ggplot2::coord_flip()+
      ggplot2::scale_size_area(max_size = 6, limits = c(0,100), breaks = c(0,25,50,75,100))  #ggplot2::scale_size_area()
})
    suppressMessages({del_dotplot <- Seurat::DotPlot(object, features =c(rev(as.character(df$Feature))), group.by = del)+ggplot2::coord_flip()+
      ggplot2::scale_size_area(max_size = 6, limits = c(0,100), breaks = c(0,25,50,75,100))  #ggplot2::scale_size_area()
    })
    return(list(diff_gini = bar_p, add_dotplot = add_dotplot, del_dotplot=del_dotplot))
  })
  names(p_list) <-names(metrics_list)
  return(p_list)
}

#' Run Noise Gene Clustering Benchmark
#'
#' @description
#' This function performs a benchmark to evaluate the impact of removing different
#' sets of "noise" genes (e.g., mitochondrial, ribosomal, dissociation-induced)
#' on the clustering of single-cell RNA-seq data. It compares clustering results
#' after removing all noise genes versus removing all noise genes but then adding
#' back specific subsets.
#'
#' @param object A single-cell object (e.g., Seurat object) from which
#'   expression data and metadata will be extracted.
#' @param gene_list (Optional) A named list of character vectors, where each
#'   vector contains gene names considered as a specific type of "noise" (e.g.,
#'   `list(mt = c("MT-G1", "MT-G2"))`). If `NULL`, `GetNoiseGeneList` is used
#'   to define common noise gene sets.
#' @param nfeatures An integer specifying the number of highly variable features
#'   to use for PCA. Defaults to 2000.
#' @param har.batch.by A character string or vector specifying the metadata
#'   columns to use for batch correction with Harmony.
#' @param resolution A numeric vector of resolution values for Leiden clustering.
#'   Defaults to `c(1, 2)`.
#' @param tmpdir A character string specifying a temporary directory for large
#'   intermediate files (particularly from `BPCells`). Defaults to
#'   "./temp_SingleCellMQC/BenchmarkDelGene/".
#' @inheritParams common_params
#'
#' @return A named list of data frames. Each element in the list corresponds
#'   to a clustering `resolution` (e.g., `resolution=1`, `resolution=2`).
#'   Each data frame contains clustering assignments for each cell, with columns
#'   like 'del_all' (clustering after removing all noise genes) and 'add_X'
#'   (clustering after removing all noise genes and then adding back noise set X).
#'
#' @details
#' The function first defines or takes specific noise gene sets. It then performs
#' the following steps:
#' 1. **"Del All" Baseline**: Runs PCA and Harmony (for batch correction) after
#'    excluding *all* genes specified in the combined `gene_list`. Clustering
#'    is then performed.
#' 2. **"Add Back" Scenarios**: For each individual noise gene set in `gene_list`,
#'    it repeats the PCA and Harmony process, but this time, after excluding
#'    all noise genes, it *adds back* the genes from the current specific noise set.
#'    Clustering is then performed.
#' All PCA, Harmony, and clustering steps are performed with fixed random seeds
#' for reproducibility across comparisons.
#'
#' @seealso
#' \code{\link{GetNoiseGeneList}}.
#'
#' @export
RunNoiseGeneCluster <-function(object,
                               assay="RNA",
                               slot="counts",
                               layer=NULL,
                                  gene_list =NULL,
                                  nfeatures=2000,
                                  har.batch.by,
                                  resolution=c(1, 2),
                                  tmpdir="./temp_SingleCellMQC/BenchmarkDelGene/"){

  metadata = getMetaData(object)
  ## gene list
  if(is.null(gene_list)){
    gene_list <- GetNoiseGeneList(object)
  }
  ## filter e
  gene_list <- Filter(function(x) length(x) > 0, gene_list)

  del_all_gene <- do.call(c, gene_list)
  if( length(del_all_gene)==0){
    stop("The list of genes to delete is empty. Please provide a valid gene_list or ensure GetNoiseGeneList(object) returns non-empty.")
  }

  # for in genelist to clustering (before & after)
  ### del all genes
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ">>>>>>>>>>>> del all "))
  set.seed(1)
  cluster_del_all <- stepRNAToPCA(object,
                                  assay = assay,
                                  slot= slot,
                                  layer=layer,
                                  nfeatures=nfeatures,
                                  exclude_custom_features=del_all_gene,
                                  tmpdir=tmpdir,
                                  only.returnPCA=T)
  set.seed(1)
  cluster_del_all <- harmony::RunHarmony(cluster_del_all,
    meta_data = metadata,
    vars_use = har.batch.by,
    return_object = FALSE
  )
  set.seed(1)
  cluster_del_all <- stepPCAToCluster(cluster_del_all, output_cluster_name= "del_all", resolution = resolution)


  ## split by
  cluster_dellist <- lapply( names(gene_list), function(x){
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ">>>>>>>>>>>> add ", x))
    set.seed(1)
    cluster_del <- stepRNAToPCA(object,
                                assay = assay,
                                slot= slot,
                                layer=layer,
                                nfeatures=nfeatures,
                                exclude_custom_features=del_all_gene,
                                add_custom_features=gene_list[[x]],
                                tmpdir=tmpdir)
    set.seed(1)
    cluster_del <- harmony::RunHarmony(cluster_del,
        meta_data = metadata,
        vars_use = har.batch.by,
        return_object = FALSE
      )
    set.seed(1)
    cluster_del <- stepPCAToCluster(cluster_del, output_cluster_name= paste0("add_", x) ,resolution=resolution)
    return(cluster_del)
  })
  cluster_del_all <- cbind(cluster_del_all, do.call(cbind, cluster_dellist))

  groups <- split(1:ncol(cluster_del_all), rep(resolution, ncol(cluster_del_all)/2))
  out <- lapply(groups, function(cols) cluster_del_all[, cols, drop=FALSE])
  names(out) <- resolution
  return(out)
}


#' Get Predefined Noise Gene Lists
#'
#' @description
#' This function identifies various categories of "noise" genes (Mitochondrial,
#' Ribosomal, Hemoglobin, Dissociation/Stress, TCR/BCR, and Sex-chromosome
#' linked genes).
#'
#' @details
#' For each category, the function identifies genes using two
#' **complementary** methods:
#' 1. **Grepping**: It searches for genes matching a regular expression (e.g., `mt.pattern`).
#' 2. **Filtering**: It looks for specific gene symbols provided in a list (e.g., `mt.genes`).
#'
#' The results from both methods are merged (union) and deduplicated (unique).
#' Finally, any genes specified in `exclude.genes` (e.g., HSPG2, TRAV1-2) are
#' removed from the final output. The function is case-insensitive to support
#' both Human and Mouse nomenclature.
#'
#' @param object A single-cell object.
#' @param mt.pattern Regex for mitochondrial genes (e.g., `^MT-`). Combined with `mt.genes`.
#' @param mt.genes List of specific mitochondrial symbols. Combined with `mt.pattern`.
#' @param rb.pattern Regex for ribosomal genes (e.g., `^RP[SL]`). Combined with `rb.genes`.
#' @param rb.genes List of specific ribosomal symbols. Combined with `rb.pattern`.
#' @param hb.pattern Regex for hemoglobin genes. Combined with `hb.genes`.
#' @param hb.genes List of specific hemoglobin symbols (e.g. "HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBZ"). Combined with `hb.pattern`.
#' @param dis.pattern Regex for stress/dissociation genes, Default is "^FOS|^JUN|^DNAJ|^HSP". Combined with `dis.genes`.
#' @param dis.genes List of specific stress symbols (e.g., FOS, JUN). Combined with `dis.pattern`.
#' @param tcr.pattern Regex for T-cell receptor genes, Default is alpha/beta `^TRA[JV].*|^TRB[VDJ].*`. Combined with `tcr.genes`.
#' @param tcr.genes List of specific TCR symbols. Combined with `tcr.pattern`.
#' @param bcr.pattern Regex for B-cell receptor genes, Default is `^IGH[VJ].*|^IGK[JV].*|^IGL[JV].*`. Combined with `bcr.genes`.
#' @param bcr.genes List of specific BCR symbols. Combined with `bcr.pattern`.
#' @param sex.genes List of sex-chromosome-linked genes (e.g., XIST, UTY).
#' @param exclude.genes A blacklist of genes to be removed from all results. Default: "HSPG2" and "TRAV1-2".
#'
#' @return A named list of gene vectors for each noise category.
#' @export
GetNoiseGeneList <- function(
    object,
    # 1. Mitochondria (MT)
    mt.pattern = "^MT-",
    mt.genes   = NULL,
    # 2. Ribosome (RB)
    rb.pattern = "^RP[SL]",
    rb.genes   = NULL,
    # 3. Hemoglobin (HB)
    hb.pattern = NULL,
    hb.genes   = c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBZ"),
    # 4. Dissociation / Stress
    dis.pattern = "^FOS|^JUN|^DNAJ|^HSP",
    dis.genes   = c("DUSP1", "IER2", "IER3", "BTG1", "BTG2", "ATF3", "EGR1", "CEBPB",
                    "CEBPB-AS1", "CEBPD", "CXCL1", "FOS","FOSB", "FOSL1", "FOSL1P1", "FOSL2",
                    "ID3","JUN", "JUNB", "JUND", "MT1A", "MT1B", "MT1E", "MT1F", "MT1G", "MT1H",
                    "MT1L", "MT1M", "MT1X", "MT2A", "NFKBIA", "NR4A1", "PPP1R15A", "SOCS3",
                    "UBC", "ZFP36"),

    # 5. TCR
    tcr.pattern = "^TRA[JV].*|^TRB[VDJ].*",
    tcr.genes   = NULL,
    # 6. BCR
    bcr.pattern = "^IGH[VJ].*|^IGK[JV].*|^IGL[JV].*",
    bcr.genes   = NULL,
    # 7. Sex Chromosome
    sex.genes   = c("ARSD", "CXorf15", "DDX3X", "HDHD1A", "KDM5C", "PNPLA4", "RIBC1",
                    "RPS4X", "KDM6A", "ZFX", "XIST", "EIF2S3X", "ZFY", "USP9Y", "UTY",
                    "PRKY", "CYorf15A", "CYorf15B", "RPS4Y1", "NCRNA00185", "KDM5D",
                    "EIF1AY", "DDX3Y"),
    # 8. Global Exclusion
    exclude.genes = c("HSPG2", "TRAV1-2")
) {

  # A. Extract gene names from the provided object
  gene_name <- tryCatch({
    if (is.matrix(object) || is(object, "dgCMatrix")) {
      rownames(object)
    } else if (exists("getMatrix")) {
      rownames(getMatrix(object))
    } else {
      rownames(object)
    }
  }, error = function(e) stop("Failed to extract gene names from the input object."))

  if (is.null(gene_name)) stop("Object rownames are empty.")

  # B. Helper function: Union of Grep & List, then Unique, then Exclude
  .collect <- function(all_genes, pa, ge, ex) {
    # 1. Matching by Pattern
    res_pa <- character(0)
    if (!is.null(pa) && pa != "") {
      res_pa <- grep(pa, all_genes, value = TRUE, ignore.case = TRUE)
    }

    # 2. Matching by Gene List (Case-insensitive)
    res_ge <- character(0)
    if (!is.null(ge)) {
      res_ge <- all_genes[tolower(all_genes) %in% tolower(ge)]
    }

    # 3. Union and Unique
    combined <- unique(c(res_pa, res_ge))

    # 4. Exclusion (Case-insensitive)
    if (length(combined) > 0 && !is.null(ex)) {
      combined <- combined[!(tolower(combined) %in% tolower(ex))]
    }

    return(combined)
  }

  # C. Assemble the named list for all categories
  gene_list <- list(
    "mt"              = .collect(gene_name, mt.pattern, mt.genes, exclude.genes),
    "rb"              = .collect(gene_name, rb.pattern, rb.genes, exclude.genes),
    "hb"              = .collect(gene_name, hb.pattern, hb.genes, exclude.genes),
    "dissociation"    = .collect(gene_name, dis.pattern, dis.genes, exclude.genes),
    "TCRab"           = .collect(gene_name, tcr.pattern, tcr.genes, exclude.genes),
    "BCR"             = .collect(gene_name, bcr.pattern, bcr.genes, exclude.genes),
    "Sex_chromosome"  = .collect(gene_name, NULL, sex.genes, exclude.genes)
  )

  return(gene_list)
}


#KIR23_gene <- grep("^KIR2|^KIR3", gene_name, value=TRUE)
# cell_cycle_s <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
#                       "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
#                       "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
#                       "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
#                       "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
#                       "BRIP1", "E2F8")
# cell_cycle_g2m <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
#                "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
#                "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
#                "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
#                "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
#                "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")



#' Calculate Noise Gene Cluster Metrics
#'
#' @description
#' This function calculates various metrics to evaluate the impact of noise
#' gene removal strategies on clustering and metadata variable associations.
#' It computes variance fractions explained by `del_all` and `add_X` clusterings,
#' and differences in Gini coefficients for relevant genes.
#'
#' @param object A single-cell object (e.g., Seurat object) from which
#'   expression data and metadata will be extracted for Gini calculations.
#' @param cluster_info A data frame or a named list of data frames containing
#'   clustering results, typically the output of `RunNoiseGeneCluster`.
#'   Each data frame should have columns 'del_all' and 'add_X' for different
#'   noise gene removal scenarios.
#' @inheritParams common_params
#'
#' @return If `cluster_info` is a single data frame, returns a list with:
#'   \itemize{
#'     \item `Variance_fractions`: A data frame summarizing the variance
#'       explained by clusterings for various metadata features (e.g., %MT, %RB, inferred cell cycle, dissociation genes).
#'     \item `Diff_Gini`: A data frame summarizing the differences in Gini
#'       coefficients for specific marker genes across 'add_X' vs 'del_all' scenarios.
#'   }
#'   If `cluster_info` is a list (e.g., for different resolutions), returns
#'   a named list of such lists, where names correspond to the resolutions.
#'
#' @details
#' The function uses two main types of metrics:
#' 1. **Variance Fractions**: For each metadata variable (e.g., percent.mt, percent.rb),
#'    it fits a mixed-effects model using `lmerTest::lmer` and
#' `variancePartition::calcVarPart` to calculate the fraction of variance
#' explained by the noise-aware clusterings
#'    (`del_all`, `add_X`). `calMetadataVarPart` (an internal helper) is used for this.
#' 2. **Gini Coefficient Differences**: For specific marker genes (e.g., from TCR, BCR, sex chromosomes),
#'    it calculates the Gini coefficient
#'    for gene expression proportion (pct) across clusters for both `del_all` and `add_X` scenarios.
#'    The differences (add_X - del_all) are then reported, indicating if adding back
#'    certain noise genes worsens the localized expression proportion of these markers.
#'
#' @importFrom stats as.formula
#' @importFrom gtools smartbind
#' @importFrom data.table data.table melt
#' @export
CalNoiseGeneClusterMetrics <- function(object,
                                       assay="RNA",
                                       slot="data",
                                       layer=NULL,
                                       cluster_info){
  if(is.list(cluster_info)){
    out_list <- lapply(cluster_info, function(x){
      calNoiseGeneMetrics(object, del_clusters=x, assay=assay,slot=slot, layer=layer)
    })
    names(out_list) <- names(cluster_info)
    return(out_list)
  }else{
    return(calNoiseGeneMetrics(object, del_clusters=cluster_info,assay=assay,slot=slot, layer=layer))
  }
}





calNoiseGeneMetrics <- function(object,
                                assay="RNA",
                                slot="data",
                                layer=NULL,
                                del_clusters){
  metadata = getMetaData(object)
  del_clusters <- del_clusters[match(rownames(metadata), rownames(del_clusters) ),]
  metadata <- cbind(metadata, del_clusters )
  cluster_name <- colnames(del_clusters)

  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ">>>>>>>>>>>> Variance Partitioning "))
  ## Variance Partitioning variance fractions
  del_all <- grep("del_all", cluster_name, value = T)
  del_all_out <- calMetadataVarPart(metadata, formula= stats::as.formula(paste0( "~ (1 |", del_all, ")")),
                            group_name = "del_all",
                            cluster_name=del_all,
                            features=c("percent.mt", "percent.rb", "percent.hb", "percent.dissociation"),
                            do.log1p=T )


  #
  varPartList <- lapply(c("mt", "rb", "hb", "dissociation"), function(x){
    in_name <- intersect(colnames(metadata), c( paste0("percent." , x) ))
    if(length(in_name)==0){
      return(NULL)
    }
    print(in_name)
    index_name <- grep(x, cluster_name, value = T)
    form = stats::as.formula( paste0("~  (1 | ", index_name, ")"))
    out <- calMetadataVarPart(metadata, formula= form,
                                      group_name ="add",
                              cluster_name=index_name,
                                      features=c(in_name, "nCount_RNA"),
                                      do.log1p=T )
    out <- out[out$Feature != "nCount_RNA", ,drop=F]
    return(out)
  })
  varPartList <- do.call(gtools::smartbind, c(list(del_all_out) , varPartList) )
  varPartList <- varPartList[, c(setdiff(names(varPartList), "Residuals"), "Residuals")]

  ## gini
  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ">>>>>>>>>>>> gini "))
  del_index <- lapply(c("mt", "rb", "hb", "dissociation"), function(x){grep(x, cluster_name, value = T)})
  cluster_name_gini <- setdiff(cluster_name, del_index)
  gini_list <- lapply(cluster_name_gini, function(x){
    message( paste0(format(Sys.time(), "%H:%M:%S"), "-------- runGini ", x))
    runGini(object,
            assay=assay,
            slot=slot,
            layer=layer,
            metadata = del_clusters,
            split.by = x)
  })
  names(gini_list) <- cluster_name_gini
  pct.gini <- do.call(cbind, gini_list )
  colnames(pct.gini) <- cluster_name_gini
  rownames(pct.gini) <- rownames(object)
  diff_pct.gini <- pct.gini[, -1] - pct.gini[, 1]
  noise_list <- GetNoiseGeneList(object)

  for (x in c("TCRab", "BCR", "Sex_chromosome")) {
    index_cols <- grep(x, colnames(diff_pct.gini), value = TRUE)
    if (length(index_cols) == 0) {
      next
    } else {
      diff_pct.gini[[index_cols]][!(rownames(diff_pct.gini) %in% noise_list[[x]]) ] <-NA # (C)
    }
  }
  diff_pct.gini <- data.table::data.table(Feature=rownames(diff_pct.gini), diff_pct.gini)
  diff_pct.gini <- data.table::melt(diff_pct.gini, id.vars="Feature")
  diff_pct.gini <- na.omit(diff_pct.gini)
  return(list(Variance_fractions=varPartList,
              Diff_Gini = diff_pct.gini ))
}



runGini <- function(object, assay="RNA", slot= "data", layer=NULL, metadata, split.by){
  mat <- getMatrix(object, assay=assay, slot=slot, layer=layer)
  noise_list <- GetNoiseGeneList(object)
  if( !("BPCells" %in% attr(class(mat), "package"))){
    mat <- as(mat, "IterableMatrix")
  }
  row_indices <- split(1:ncol(mat), metadata[[split.by]])
  split_matrices <- lapply(row_indices, function(idx) mat[, idx, drop = FALSE])

  out_list <- lapply( names(split_matrices), function(y){
    stats <- BPCells::matrix_stats(split_matrices[[y]], row_stats="nonzero")
    stats[["row_stats"]]["nonzero", ] <- stats[["row_stats"]]["nonzero", ]/dim(split_matrices[[y]])[2]
    #stats$pct <- stats$nonzero / dim(y)[2]
    return(stats)
  })
  names(out_list) <- names(split_matrices)
  pct_value <- do.call(cbind, lapply(out_list, function(x){x[["row_stats"]]["nonzero", ]}) )
  pct_metrics <- calGini(pct_value)
  return(pct_metrics)
}

# expr: genes x cells (matrix or data.frame)
calGini <- function(expr) {
  # Gini
  gini_coef <- function(x) {
    n <- length(x)
    x <- as.numeric(x)
    if (all(x == 0)) return(0) # 预防全零
    x_sorted <- sort(x)
    index <- seq_along(x_sorted)
    G <- (2 * sum(index * x_sorted) / (n * sum(x_sorted))) - ((n + 1) / n)
    return(G)
  }
  ginis <- apply(expr, 1, gini_coef)
  data.frame(Gini = ginis, row.names = rownames(expr))
}


