

# pipeline function -------------------------------------------------------

# add_custom_features will add hvg features =  nfeatures + add_custom_features
# exclude_custom_features & exclude_patterns will contain always hvg features = nfeatures
stepRNAToPCA <- function(object,
                       assay = "RNA",
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

  mat <- getMatrix(object)
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
    svd <- BPCells::svds(mat_norm, k=k)
  }else{
    k <- min(nrow(mat_norm), 50)
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
      object <- SeuratObject::SetAssayData(object = object,
                                           assay = assay,
                                           slot = "data",
                                           new.data = mat)
    }
    if(add.scale_data){
      mat_norm1 <- BPCells::min_by_row(mat = mat_norm1, vals = 10 * features.sd + features.mean)
      mat_norm1 <- (mat_norm1 - features.mean) / features.sd
      object <- SeuratObject::SetAssayData(object = object,
                                           assay = assay,
                                           slot = "scale.data",
                                           new.data = mat_norm1)
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





# benchmarking del feature ------------------------------------------------

PlotNoiseGeneClusterMetrics <- function(metrics, type="Variance_fractions"){
  if(type == "Variance_fractions" ){
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
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), # x轴文字旋转45度，右对齐，颜色黑色
        axis.text.y = element_text(color = "black"), # y轴文字颜色黑色
        axis.title.x = element_text(color = "black"), # x轴标题颜色黑色
        legend.text = element_text(color = "black"), # 图例文字颜色黑色
        plot.title = element_text(color = "black"), # 如果有主标题，也设置为黑色
        legend.position = "right", # 图例位置
      )
    return(p1)
  }


  if(type=="Diff_Gini"){

  }

}





RunNoiseGeneCluster <-function(object,
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
  cluster_del_all <- stepRNAToPCA(object, nfeatures=nfeatures, exclude_custom_features=del_all_gene, tmpdir=tmpdir, only.returnPCA=T) %>%
    harmony::RunHarmony(
      meta_data = metadata,
      vars_use = har.batch.by,
      return_object = FALSE
    ) %>%
    stepPCAToCluster(output_cluster_name= "del_all", resolution = resolution)


  # in_name <- intersect(colnames(metadata), c("percent.mt", "percent.rb", "percent.hb", "percent.dissociation") )
  # if(length(in_name)!=0){
  #   regress_cluster <- lapply(in_name, function(x){
  #     message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ">>>>>>>>>>>> del all regress ", x))
  #     out <- stepRNAToPCA(object, nfeatures=nfeatures, exclude_custom_features=del_all_gene, tmpdir=tmpdir, only.returnPCA=T, vars.to.regress = x) %>%
  #       harmony::RunHarmony(
  #         meta_data = metadata,
  #         vars_use = har.batch.by,
  #         return_object = FALSE
  #       ) %>%
  #       stepPCAToCluster(output_cluster_name= paste0("del_all_regress_", x), resolution = resolution)
  #     return(out)
  #   })
  #   cluster_del_all <- cbind(cluster_del_all, do.call(cbind, regress_cluster))
  # }

  ## split by
  cluster_dellist <- lapply( names(gene_list), function(x){
    message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ">>>>>>>>>>>> add ", x))
    cluster_del <- stepRNAToPCA(object, nfeatures=nfeatures, exclude_custom_features=del_all_gene, add_custom_features=gene_list[[x]], tmpdir=tmpdir) %>%
      harmony::RunHarmony(
        meta_data = metadata,
        vars_use = har.batch.by,
        return_object = FALSE
      ) %>%
      stepPCAToCluster(output_cluster_name= paste0("add_", x) ,resolution=resolution)
    return(cluster_del)
  })
  cluster_del_all <- cbind(cluster_del_all, do.call(cbind, cluster_dellist))

  groups <- split(1:ncol(cluster_del_all), rep(resolution, ncol(cluster_del_all)/2))
  out <- lapply(groups, function(cols) cluster_del_all[, cols, drop=FALSE])
  names(out) <- resolution
  return(out)
}



GetNoiseGeneList <- function(object){
  mat <- getMatrix(object)
  gene_name <- rownames(mat)
  mt_gene <- grep("^MT-", gene_name, value = T)
  rb_gene <- grep("^RP[SL]", gene_name, value = T)
  hb_gene <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBZ")

  dis_pa <- c("^FOS","^JUN","^DNAJ","^HSP")
  dis_gene <- c(setdiff(grep(paste0(dis_pa, collapse = "|"), gene_name, value = TRUE), "HSPG2"), "DUSP1", "IER2", "IER3", "BTG1", "BTG2", "ATF3","EGR1","ATF3","BTG2","CEBPB","CEBPB-AS1","CEBPD","CXCL1","EGR1","FOS","FOSB","FOSL1","FOSL1P1","FOSL2","ID3","IER2","JUN",
                "JUNB","JUND","MT1A","MT1B","MT1E","MT1F","MT1G","MT1H","MT1L","MT1M","MT1X","MT2A","NFKBIA","NR4A1","PPP1R15A","SOCS3",
                "UBC","ZFP36")

  TCRab_gene <- grep("^TRA[JV].*|^TRB[VDJ].*", gene_name, value = T)
  BCR_gene <- grep(paste0(c("^IGH[VDJ].*",
                              "^IGK[JV].*",
                              "^IGL[JLV].*",
                              "^IGL[ON].*"), collapse = "|"), gene_name, value = T)
  sex_gene <- c(
    # X chromosome
    "ARSD", "CXorf15", "DDX3X", "HDHD1A", "KDM5C", "PNPLA4",
    "RIBC1", "RPS4X", "KDM6A", "ZFX", "XIST",
    "EIF2S3X",
    # Y chromosome
    "ZFY", "USP9Y", "UTY", "PRKY", "CYorf15A", "CYorf15B",
    "RPS4Y1", "NCRNA00185", "KDM5D", "EIF1AY", "DDX3Y"
  )
  gene_list <- list(
    "mt" = intersect(mt_gene, gene_name),
    "rb" = intersect(rb_gene,gene_name),
    "hb" = intersect(hb_gene,gene_name),
    "dissociation" = intersect(dis_gene,gene_name),
    "TCRab" = intersect(TCRab_gene,gene_name),
    "BCR" = intersect(BCR_gene,gene_name),
    "Sex_chromosome" = intersect(sex_gene,gene_name)
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



# test --------------------------------------------------------------------
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
                               do.log1p=T,
                               features = c("percent.mt", "percent.rb", "percent.hb", "percent.dissociation")){
  metadata = getMetaData(object)
  in_name <- intersect(colnames(metadata), features )
  data = t(metadata[, c(in_name), drop=F ])
  if(do.log1p){
    data = log1p(data)
  }
  varPart <- variancePartition::fitExtractVarPartModel(data, formula, metadata)
  varPart <- data.frame(Group = group_name, Feature=rownames(varPart), data.frame(varPart))
  return(varPart)
}


test4 <- function(object){

}



CalNoiseGeneClusterMetrics <- function(object, cluster_info){
  out_list <- lapply(cluster_info, function(x){
    calNoiseGeneMetrics(object, x)
  })
  names(out_list) <- names(cluster_info)
  return(out_list)
}





calNoiseGeneMetrics <- function(object, del_clusters){
  metadata = getMetaData(object)
  del_clusters <- del_clusters[match(rownames(metadata), rownames(del_clusters) ),]
  metadata <- cbind(metadata, del_clusters )
  cluster_name <- colnames(del_clusters)

  message( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ">>>>>>>>>>>> fitExtractVarPartModel "))
  ## fitExtractVarPartModel variance fractions
  del_all <- grep("del_all", cluster_name, value = T)
  del_all_out <- calMetadataVarPart(metadata, formula= stats::as.formula(paste0("~ (1 |", del_all, ")")),
                            group_name = "del_all",
                            features=c("percent.mt", "percent.rb", "percent.hb", "percent.dissociation"),
                            do.log1p=T )

  #
  varPartList <- lapply(c("mt", "rb", "hb", "dissociation"), function(x){
    in_name <- intersect(colnames(metadata), c( paste0("percent." , x) ))
    if(length(in_name)==0){
      return(NULL)
    }
    index_name <- grep(x, cluster_name, value = T)
    form = stats::as.formula( paste0("~  (1 | ", index_name, ")"))
    out <- calMetadataVarPart(metadata, formula= form,
                                      group_name ="add",
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
    runGini(object, metadata = del_clusters, split.by = x)
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



runGini <- function(object, metadata, split.by){
  mat <- getMatrix(object, slot="data")
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


