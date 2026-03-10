

#' @title Calculate Multi-omics Metrics
#' @description
#' Calculates different metrics, including sample, cell, and feature metrics.
#' Utilizes `CalculateMetricsPerCell`, `CalculateMetricsPerSample`, and `CalculateMetricsPerFeature` to perform its tasks based on the `types` argument.
#'
#' @param object A Seurat object or a list containing GEX (Seurat object), TCR (list; optional), and BCR (list; optional) components.
#' @param sample.by Character string specifying the grouping variable for sample metrics. Defaults to "orig.ident".
#' @param types A character vector specifying the types of metrics to run: "Sample", "Cell", and/or "Feature". Defaults to c("Sample", "Cell").
#' @param exclude.genes A blacklist of genes to be excluded from the matched gene sets (e.g., MT, RB). Default: "HSPG2".
#' @inheritParams GetNoiseGeneList
#'
#' @return A Seurat object. All results are stored in `meta.data` and `SingleCelMQC` slot in `misc`.
#' @seealso \code{\link{CalculateMetricsPerCell}}, \code{\link{CalculateMetricsPerSample}}, \code{\link{CalculateMetricsPerFeature}}
#' @export
#'
CalculateMetrics <- function(object,
                             types=c("Sample", "Cell"),
                             sample.by="orig.ident",
                             mt.pattern = "^MT-",
                             mt.genes   = NULL,
                             rb.pattern = "^RP[SL]",
                             rb.genes   = NULL,
                             hb.pattern = NULL,
                             hb.genes   = c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBZ"),
                             dis.pattern = "^FOS|^JUN|^DNAJ|^HSP",
                             dis.genes   = c("DUSP1", "IER2", "IER3", "BTG1", "BTG2", "ATF3", "EGR1", "CEBPB",
                                             "CEBPB-AS1", "CEBPD", "CXCL1", "FOS","FOSB", "FOSL1", "FOSL1P1", "FOSL2",
                                             "ID3","JUN", "JUNB", "JUND", "MT1A", "MT1B", "MT1E", "MT1F", "MT1G", "MT1H",
                                             "MT1L", "MT1M", "MT1X", "MT2A", "NFKBIA", "NR4A1", "PPP1R15A", "SOCS3",
                                             "UBC", "ZFP36"),
                             exclude.genes = c("HSPG2")
                             ){
  if ("Cell" %in% types) {
    object <- CalculateMetricsPerCell(object,
                                      mt.pattern    = mt.pattern,
                                      mt.genes      = mt.genes,
                                      rb.pattern    = rb.pattern,
                                      rb.genes      = rb.genes,
                                      hb.pattern    = hb.pattern,
                                      hb.genes      = hb.genes,
                                      dis.pattern   = dis.pattern,
                                      dis.genes     = dis.genes,
                                      exclude.genes =exclude.genes)
  }
  if ("Sample" %in% types) {
    object <- CalculateMetricsPerSample(object, add.Seurat = TRUE, sample.by=sample.by )
  }
  if ("Feature" %in% types) {
    object <- CalculateMetricsPerFeature(object, add.Seurat = TRUE)
  }
  return(object)
}




#' Calculate Cell QC Metrics
#'
#' @description
#' A function for Calculating cell-specific metrics based on the data type.
#'
#' @param object The object is obtained by the `Read10XData` function, or a list containing GEX (single-cell object), TCR (list; optional), and BCR (list; optional) components, or a single-cell object.
#' @param exclude.genes A blacklist of genes to be excluded from the matched gene sets (e.g., MT, RB). Default: "HSPG2".
#'
#' @inheritParams GetNoiseGeneList
#'
#' @return A single-cell object. All results are stored in \code{meta.data} and \code{SingleCelMQC} slot in \code{misc}
#'
#' @export
CalculateMetricsPerCell <- function(
    object,
    mt.pattern = "^MT-",
    mt.genes   = NULL,
    rb.pattern = "^RP[SL]",
    rb.genes   = NULL,
    hb.pattern = NULL,
    hb.genes   = c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBZ"),
    dis.pattern = "^FOS|^JUN|^DNAJ|^HSP",
    dis.genes   = c("DUSP1", "IER2", "IER3", "BTG1", "BTG2", "ATF3", "EGR1", "CEBPB",
                    "CEBPB-AS1", "CEBPD", "CXCL1", "FOS","FOSB", "FOSL1", "FOSL1P1", "FOSL2",
                    "ID3","JUN", "JUNB", "JUND", "MT1A", "MT1B", "MT1E", "MT1F", "MT1G", "MT1H",
                    "MT1L", "MT1M", "MT1X", "MT2A", "NFKBIA", "NR4A1", "PPP1R15A", "SOCS3",
                    "UBC", "ZFP36"),
    exclude.genes = c("HSPG2")
) {

  object <- calculateMetricsPerCell(object, add.Seurat=TRUE,
                               mt.pattern    = mt.pattern,
                               mt.genes      = mt.genes,
                               rb.pattern    = rb.pattern,
                               rb.genes      = rb.genes,
                               hb.pattern    = hb.pattern,
                               hb.genes      = hb.genes,
                               dis.pattern   = dis.pattern,
                               dis.genes     = dis.genes,
                               exclude.genes =exclude.genes)

  return(object)
}

.checkMetricsInput <- function(x) {
  if (inherits(x, "Seurat") || inherits(x, "SingleCellExperiment") || "Gene Expression" %in% names(x)) {
    return("GEX")
  } else if (is.list(x) && ("chain" %in% colnames(x[[1]]) || "chain" %in% colnames(x))) {
    return("VDJ")
  } else {
    return("Other")
  }
}
calculateMetricsPerCell <- function(object, add.Seurat, ...) {
  UseMethod("calculateMetricsPerCell")
}

#' @method calculateMetricsPerCell Seurat
calculateMetricsPerCell.Seurat <- function(object,add.Seurat, ...) {
  # ... 会自动包含你传入的 mt.pattern, rb.pattern 等
  out_GEX <- CalculateMetricsPerCell.GEX(object = object,add.Seurat, ...)
  return(out_GEX)
}

#' @method calculateMetricsPerCell list
calculateMetricsPerCell.list <- function(object,add.Seurat, ...) {

  # 1. 识别列表中每个元素的类型
  # .checkMetricsInput 一般返回 "GEX", "VDJ", "Other" 等
  types <- sapply(object, .checkMetricsInput)

  # 分类索引
  gex_idx <- which(types == "GEX")
  other_idx <- which(types == "VDJ")

  # 2. 使用 lapply 遍历处理所有非 GEX 模态 (VDJ, ADT 等)
  # 这样无论列表中有 1 个还是 10 个 VDJ 模态，都会被逐一处理
  all_other_metrics_list <- lapply(other_idx, function(i) {
    item <- object[[i]]
    type <- types[i]
    # 针对不同模态进行分发
    if (type == "VDJ") {
      # 显式调用 VDJ 处理函数
      return(CalculateMetricsPerCell.VDJ(item))
    }
  })

  # 3. 过滤掉可能产生的空结果并合并 DataFrames
  all_other_metrics_list <- Filter(function(x) !is.null(x), all_other_metrics_list)

  merged_other_metrics <- NULL
  if (length(all_other_metrics_list) > 0) {
    if (length(all_other_metrics_list) == 1) {
      merged_other_metrics <- all_other_metrics_list[[1]]
    } else {
      # 使用 Reduce 和内部合并函数 .merge.qc 将所有模态的指标合并为一张大表
      merged_other_metrics <- Reduce(function(x, y) .merge.qc(x, y), all_other_metrics_list)
    }
  }

  # 4. 场景 A：如果列表中没有 GEX 模态，直接返回非 GEX 的合并指标表
  if (length(gex_idx) == 0) {
    if (!is.null(merged_other_metrics)) {
      return(merged_other_metrics)
    } else {
      stop("Input list contains no valid GEX or VDJ modalities for metric calculation.")
    }
  }

  # 5. 场景 B：有 GEX 模态 (通常将其作为主分析对象)
  # 取第一个发现的 GEX 模态作为容器
  res_obj <- toSeuratObject(object[[gex_idx[1]]])

  # 计算 GEX 指标 (如线粒体、核糖体、细胞周期等)
  # 此时会通过泛函数自动匹配到 calculateMetricsPerCell.Seurat
  metrics_GEX <- calculateMetricsPerCell(res_obj, add.Seurat =F, ...)
  res_obj <- Seurat::AddMetaData(res_obj, metadata = metrics_GEX)
  metrics_name = colnames(metrics_GEX)

  # 6. 将合并后的其他模态指标 (VDJ 等) 同步到 Seurat 对象中
  if (!is.null(merged_other_metrics)) {
    # 注入 Metadata
    res_obj <- Seurat::AddMetaData(res_obj, metadata = merged_other_metrics)
    metrics_name = unique(c(metrics_name, colnames(merged_other_metrics)))
  }
  full_meta <- getMetaData(res_obj)
  metrics_out <- full_meta[, match(metrics_name, colnames(full_meta)), drop=F]
  suppressWarnings(Seurat::Misc(res_obj, slot = "SingleCellMQC")$perQCMetrics$perCell  <- metrics_out)

  return(res_obj)
}



#' @method calculateMetricsPerCell default
calculateMetricsPerCell.default <- function(object, add.Seurat, ...) {
  stop(paste("Object system does not know how to handle class:", class(object)[1]))
}



CalculateMetricsPerCell.GEX <- function(object, add.Seurat=TRUE,
                                        mt.pattern,
                                        mt.genes,
                                        rb.pattern,
                                        rb.genes,
                                        hb.pattern,
                                        hb.genes,
                                        dis.pattern,
                                        dis.genes,
                                        exclude.genes){

  noise_list <- GetNoiseGeneList(
    object        = object,
    mt.pattern    = mt.pattern,
    mt.genes      = mt.genes,
    rb.pattern    = rb.pattern,
    rb.genes      = rb.genes,
    hb.pattern    = hb.pattern,
    hb.genes      = hb.genes,
    dis.pattern   = dis.pattern,
    dis.genes     = dis.genes,
    exclude.genes =exclude.genes
  )
  name_mapping <- c(
    "mt"           = "percent.mt",
    "rb"           = "percent.rb",
    "hb"           = "percent.hb",
    "dissociation" = "percent.dissociation"
  )
  genelist <- list()
  for (old_name in names(name_mapping)) {
    new_name <- name_mapping[old_name]
    genes <- noise_list[[old_name]]

    # Only assign if the list element exists and contains at least one gene
    if (!is.null(genes) && length(genes) > 0) {
      genelist[[new_name]] <- genes
    }
  }


  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------CalculateMetricsPerCell.GEX"))
  GEX <- toSeuratObject(object)
  metadata_temp <- getMetaData(GEX)
  # 自动计算 RNA 的 nCount 和 nFeature
  if (!("nCount_RNA" %in% colnames(metadata_temp)) ) {
    GEX@meta.data$nCount_RNA <- Matrix::colSums(getMatrix(GEX, assay = "RNA", slot = "counts"))

  }
  if (!("nFeature_RNA" %in% colnames(metadata_temp)) ) {
    GEX@meta.data$nFeature_RNA <- Matrix::colSums(getMatrix(GEX, assay = "RNA", slot = "counts") > 0)
  }


  #GEX <- Seurat::CellCycleScoring(GEX, s.features = genelist$CellCycle_s, g2m.features = genelist$CellCycle_g2m)
  metric_name <- names(genelist)
  metric_list <- mapply(function(x,y){

    message(paste0(length(y)," features for calculating ", x) )

    out <- Seurat::PercentageFeatureSet(GEX, features = y, assay ="RNA")
    out <- as.data.frame(out)
    colnames(out) <- x
    return(out)
  }, x = metric_name, y = genelist, SIMPLIFY = F)

  metric_data <- do.call(cbind, metric_list)

  metadata = getMetaData(GEX)
  #RNA
  metric_data$per_feature_count_RNA <- metadata$nCount_RNA / metadata$nFeature_RNA

  if("ADT" %in% names(GEX@assays)){
    # 自动计算 RNA 的 nCount 和 nFeature
    if (!("nCount_RNA" %in% colnames(metadata)) ) {
      metadata$nCount_ADT <- Matrix::colSums(getMatrix(GEX, assay = "ADT", slot = "counts"))
    }
    if (!("nFeature_RNA" %in% colnames(metadata)) ) {
      metadata$nFeature_ADT <- Matrix::colSums(getMatrix(GEX, assay = "ADT", slot = "counts") > 0)
    }
    metric_data$per_feature_count_ADT <- metadata$nCount_ADT / metadata$nFeature_ADT
    iso <- Seurat::PercentageFeatureSet(GEX, pattern = "(?i).*isotype.*", assay ="ADT")
    metric_data$percent.isotype <- NA
    metric_data$percent.isotype[match(names(iso), rownames(metric_data))] <- iso
    metric_data <- data.frame(metadata[, c("nCount_ADT", "nFeature_ADT")], metric_data)
  }

  metric_data <- data.frame(metadata[,c("orig.ident","nCount_RNA", "nFeature_RNA")], metric_data)
  GEX  <- Seurat::AddMetaData(GEX, metadata = metric_data)
  suppressWarnings(Seurat::Misc(GEX, slot = "SingleCellMQC")$perQCMetrics$perCell  <- metric_data)
  if(!add.Seurat){
    return(metric_data)
  }
  return(GEX)
}


CalculateMetricsPerCell.VDJ <- function(object){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------CalculateMetricsPerCell.VDJ"))
  if( !("list" %in% class(object)) ){
    object <- list(object)
  }

  if(sum(grepl("TR",unique(object[[1]]$chain))) > 0){
    type = "TCR"
  }else{
    type = "BCR"
  }

  vdj_data <- lapply(object, function(x){
    x <- subset(x, chain != "Multi" & productive %in% c(TRUE, "TRUE", "True", "true"))
    #chain count
    x <- data.table::as.data.table(x)
    out <- x[, list(unique_values = paste(unique(chain), collapse = "/"), count=.N, unique_count = length(unique(chain))), by = id]
    ##chain type
    out$receptor_type <- type

    if (type == "TCR") {
      out[, receptor_subtype := dplyr::case_when(
        count > 2 ~ "multichain_T",
        count == 2 &
          unique_count == 1 &
          (
            grepl("TRA", unique_values) |
              grepl("TRB", unique_values) |
              grepl("TRG", unique_values) |
              grepl("TRD", unique_values)
          ) ~ "multichain_T",
        # grepl("Multi", unique_values) ~ "multichain",
        unique_count == 1 &
          count == 1 & grepl("TRA", unique_values) ~ "single TRA",
        unique_count == 1 &
          count == 1 & grepl("TRB", unique_values) ~ "single TRB",
        unique_count == 1 &
          count == 1 & grepl("TRG", unique_values) ~ "single TRG",
        unique_count == 1 &
          count == 1 & grepl("TRD", unique_values) ~ "single TRD",
        unique_count == 2 &
          grepl("TRA", unique_values) &
          grepl("TRB", unique_values) ~ "TRA + TRB",
        unique_count == 2 &
          grepl("TRG", unique_values) &
          grepl("TRD", unique_values) ~ "TRG + TRD",
        TRUE ~ "ambiguous_T"
      )]
    } else{
      out[, receptor_subtype := dplyr::case_when(
        count > 2 ~ "multichain_B",
        count == 2 &
          unique_count == 1 &
          (
            grepl("IGH", unique_values) |
              grepl("IGK", unique_values) |
              grepl("IGL", unique_values)
          ) ~ "multichain_B",
        # grepl("Multi", unique_values) ~ "multichain",
        unique_count == 1 &
          count == 1 & grepl("IGH", unique_values) ~ "single IGH",
        unique_count == 1 &
          count == 1 & grepl("IGK", unique_values) ~ "single IGK",
        unique_count == 1 &
          count == 1 & grepl("IGL", unique_values) ~ "single IGL",
        unique_count == 2 &
          grepl("IGH", unique_values) &
          grepl("IGK", unique_values) ~ "IGH + IGK",
        unique_count == 2 &
          grepl("IGH", unique_values) &
          grepl("IGL", unique_values) ~ "IGH + IGL",
        TRUE ~ "ambiguous_B"
      )]
    }


    out[, chain_pair := dplyr::case_when(
      count > 2 ~ "multichain",
      count == 2 & unique_count==1 ~ "multichain",
      # grepl("Multi", unique_values) ~ "multichain",
      unique_count == 1 & count ==1  ~ "single chain",
      unique_count == 2  ~ "pair",
      TRUE ~ "ambiguous"
    )]
    out$chain_pair[out$receptor_subtype=="ambiguous_B" | out$receptor_subtype=="ambiguous_T"] <- "ambiguous"
    colnames(out)[3] <- paste0("nChain_", type)
    out <- out[, -c("unique_values", "unique_count")]
    ##chain split
    chain_split <- data.table::dcast( x[,c("id","chain")] , id~chain, value.var="chain", length)
    colnames(chain_split)[-1] <- paste0("nChain_",colnames(chain_split)[-1])
    ##merge data
    out <- merge(out, chain_split, by="id")

    out <- as.data.frame(out)
    row.names(out) <- out$id
    out <- out[, -1]
    return(out)
  })

  names(vdj_data) <- NULL
  vdj_data <- do.call(dplyr::bind_rows, lapply(vdj_data, function(x) x))
  return(vdj_data)
}

.merge.qc <- function(x=NULL, y=NULL){
  index <- intersect(rownames(x)[!is.na(x$chain_pair)], rownames(y)[!is.na(y$chain_pair)] )
  x$cell <- rownames(x)
  y$cell <- rownames(y)

  x <- data.table::data.table(x)
  y <- data.table::data.table(y)
  if ("receptor_type" %in% colnames(x) & "receptor_type" %in% colnames(y)){

    merge_data <- merge(x, y[, .SD, .SDcols = !c("receptor_type", "receptor_subtype", "chain_pair")],
                        by = "cell", all = TRUE)
    merge_data[merge_data$cell %in% index, `:=` (
      receptor_type = "TCR&BCR",
      receptor_subtype = "ambiguous_TB",
      chain_pair = "ambiguous"
    )]

    diff_idx <- setdiff(y$cell, index)
    merge_data[ match(diff_idx, cell)  , `:=` (
      receptor_type = y[match(diff_idx, cell), receptor_type],
      receptor_subtype = y[match(diff_idx, cell), receptor_subtype],
      chain_pair = y[match(diff_idx, cell), chain_pair]
    )]

  } else {
    merge_data <- merge(x, y, by = "cell", all = TRUE)
  }

  setkey(merge_data, cell)
  cell_name <- merge_data$cell
  merge_data[, cell := NULL]

  merge_data <- data.frame(merge_data)
  rownames(merge_data) <- cell_name
  return(merge_data)
}



#' Calculate Sample QC Metrics
#'
#' Calculate sample qc metrics, including nCell, nGene/Protein, nChain, etc.
#'
#' @param object A Seurat object.
#' @param sample.by The column name by which to group samples for metrics computation.
#' @param add.Seurat Logical indicating whether to add metrics directly to a Seurat object. If TRUE, adds the computed metrics to the Seurat object's metadata and stored in \code{SingleCelMQC} slot in \code{misc}. Otherwise, returns a list with the metrics.
#' @param metrics A vector of qc metric names to be computed. Default: "nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT", "percent.mt", "percent.rb", "percent.hb", "percent.dissociation", "per_feature_count_RNA", "per_feature_count_ADT".
#' @return The Seurat object with sample metrics added, or a list of metrics if `add.Seurat` is FALSE.
#' @seealso \code{\link{CalculateMetricsPerSample.count}}, \code{\link{CalculateMetricsPerSample.summary}}
#' @export
#' @examples
#' \dontrun{
#' data <- CalculateMetricsPerSample(data, sample.by = "orig.ident", add.Seurat = T)
#' }
CalculateMetricsPerSample <- function(object, #seurat object or metadata list
                                      sample.by = "orig.ident",
                                      add.Seurat = T,
                                      metrics=c("nCount_RNA","nFeature_RNA","nCount_ADT","nFeature_ADT",
                                                "percent.mt","percent.rb","percent.hb","percent.dissociation",
                                                "per_feature_count_RNA", "per_feature_count_ADT")
){
  count_table <- CalculateMetricsPerSample.count(object=object,sample.by=sample.by)
  metrics_summary <- CalculateMetricsPerSample.summary(object=object,metrics=metrics,sample.by=sample.by )
  result <- list(count=count_table, summary = metrics_summary)
  if(add.Seurat){
    suppressWarnings(Seurat::Misc(object, slot = "SingleCellMQC")$perQCMetrics$perSample$count  <- result$count)
    suppressWarnings(Seurat::Misc(object, slot = "SingleCellMQC")$perQCMetrics$perSample$summary  <- result$summary)
    return(object)
  }
  return(result)
}


#' Calculate Count Metrics for Samples
#'
#' Calculate such as nCell, nGene/Protein for each sample in a Seurat object.
#'
#' @param object A Seurat object.
#' @param sample.by The column name by which to group samples for metrics computation.
#'
#' @return A data frame with samples as rows and count types as columns.
#' @seealso \code{\link{CalculateMetricsPerSample}}, \code{\link{CalculateMetricsPerSample.summary}}
#' @examples
#' \dontrun{
#' # Assuming `seuratObj` is a Seurat object
#' seuratObj <- CalculateMetricsPerSample.count(seuratObj, sample.by = "orig.ident")
#' }
#' @export
#'
#' @import data.table
CalculateMetricsPerSample.count <- function(object, #seurat object
                                            sample.by = "orig.ident"
){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------CalculateMetricsPerSample.count "))
  if (!("Seurat" %in% class(object))) {
    stop("Error: Seurat object must be as input!!")
  }

  metadata <- getMetaData(object)
  split_name <- split( as.character(rownames(metadata)), as.character(metadata[[sample.by]]) )

  if( "RNA" %in% names(object@assays)){
    SeuratObject::DefaultAssay(object) <- "RNA"
    data <- RunPseudobulkData(object, assay = "RNA", slot = "counts", sample.by = sample.by)
    gene_num <- Matrix::colSums(data$pseudobulk_mat>0)
    out_rna <- data.frame(data$ncells, nGene_RNA= gene_num)
    colnames(out_rna) <- c("sample", "nCell", "nGene_RNA")
  }

  if( "ADT" %in% names(object@assays)){
    SeuratObject::DefaultAssay(object) <- "ADT"
    data <- RunPseudobulkData(object, assay = "ADT", slot = "counts", sample.by = sample.by)
    pro_num <- Matrix::colSums(data$pseudobulk_mat>0)
    out_adt <- data.frame(sample= names(pro_num), nPro_ADT= pro_num)
    out <- merge(out_rna, out_adt, by = "sample", all = TRUE)
  }else{
    out <- out_rna
  }

  if( "nChain_TCR" %in% colnames(metadata) ){
    ##
    TCR_data <- as.data.table(metadata)
    TCR_table <- TCR_data[, list("nChain_TCR"=sum(nChain_TCR, na.rm = T),
                                 "nCell_TCR" = sum(!is.na(nChain_TCR)),
                                 "TCR%" = round(sum(!is.na(nChain_TCR)) / length(nChain_TCR) * 100,3) ),
                          by=sample.by]
    TCR_table <- as.data.frame(TCR_table)
    colnames(TCR_table)[1] <- "sample"

    ##subtype
    temp <- intersect( paste0("nChain_",c("TRA","TRB","TRD","TRG")),colnames(TCR_data) )
    TCR_subtype <- lapply(temp, function(x){
      TCR_subtype <- TCR_data[, list(nChain=sum(get(x), na.rm = T),
                                     nCell = sum( get(x) !=0 , na.rm = T),
                                     pct = round(sum(get(x) !=0, na.rm = T) / length(get(x)) *100,3)
      ),
      by = sample.by]
      name <- stringr::str_split(x[1],"_")[[1]][2]
      colnames(TCR_subtype) <- c("sample", paste0(c("nChain","nCell"), "_",name), paste0(name,"%") )
      return(TCR_subtype)
    })
    TCR_subtype <- Reduce(function(x, y) {
      merge_data <- merge(x, y, by="sample")
    } , TCR_subtype)
    TCR_table <- merge(TCR_table, TCR_subtype, by = "sample", all = TRUE)
    out <- merge(out, TCR_table, by = "sample", all = TRUE)
  }

  if( "nChain_BCR" %in% colnames(metadata) ){
    ##
    BCR_data <- as.data.table(metadata)
    BCR_table <- BCR_data[, list("nChain_BCR"=sum(nChain_BCR, na.rm = T),
                                 "nCell_BCR" = sum(!is.na(nChain_BCR)),
                                 "BCR%" = round(sum(!is.na(nChain_BCR)) / length(nChain_BCR) *100,3) ),
                          by=sample.by]
    BCR_table <- as.data.frame(BCR_table)
    colnames(BCR_table)[1] <- "sample"

    ##subtype
    temp <- intersect( paste0("nChain_",c("IGH","IGK","IGL")),colnames(BCR_data) )
    BCR_subtype <- lapply(temp, function(x){
      BCR_subtype <- BCR_data[, list(nChain=sum(get(x), na.rm = T),
                                     nCell = sum( get(x) !=0 , na.rm = T),
                                     pct = round(sum(get(x) !=0, na.rm = T) / length(get(x))*100,3)
      ),
      by = sample.by]
      name <- stringr::str_split(x[1],"_")[[1]][2]
      colnames(BCR_subtype) <- c("sample", paste0(c("nChain","nCell"), "_",name), paste0(name,"%") )
      return(BCR_subtype)
    })
    BCR_subtype <- Reduce(function(x, y) {
      merge_data <- merge(x, y, by="sample")
    } , BCR_subtype)
    BCR_table <- merge(BCR_table, BCR_subtype, by = "sample", all = TRUE)
    out <- merge(out, BCR_table, by = "sample", all = TRUE)
  }

  if("receptor_subtype" %in% colnames(metadata)){
    metadata <- as.data.table(metadata)
    receptor_table <- data.table::dcast( metadata[,list(sample.by=get(sample.by),receptor_subtype)] ,sample.by~receptor_subtype, value.var="receptor_subtype", length)
    pct <- round(receptor_table[,-1]/ apply(receptor_table[,-1],1,sum) * 100,3 )
    pct <- pct[, -"NA"]
    receptor_table <- receptor_table[,-"NA"]
    colnames(pct) <- paste0(colnames(pct),"%")
    receptor_table <- cbind(receptor_table, pct)
    colnames(receptor_table)[1] <- "sample"
    receptor_table <- as.data.frame(receptor_table)
    out <- merge(out, receptor_table, by = "sample", all = TRUE)
  }
  return(out)
}


#' Calculate Summary Metrics for Samples
#'
#' Calculate summary statistics for specified metrics for each sample in a Seurat object or a metadata list.
#'
#' @param object A Seurat object or a metadata list to compute the summaries on.
#' @param sample.by The column name by which to group samples for metrics computation.
#' @param metrics A character vector of metric names for which to calculate summaries. Default:"nCount_RNA","nFeature_RNA","nCount_ADT","nFeature_ADT","percent.mt","percent.rb","percent.hb","percent.dissociation", "per_feature_count_RNA", "per_feature_count_ADT".
#' @return A list of data frames, each containing summary statistics for a metric across samples.
#' @seealso \code{\link{CalculateMetricsPerSample}}, \code{\link{CalculateMetricsPerSample.count}}
#' @export
CalculateMetricsPerSample.summary <- function(object,
                                              sample.by="orig.ident",
                                              metrics=c("nCount_RNA","nFeature_RNA","nCount_ADT","nFeature_ADT",
                                                        "percent.mt","percent.rb","percent.hb","percent.dissociation",
                                                        "per_feature_count_RNA", "per_feature_count_ADT")){
  message(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "------CalculateMetricsPerSample.summary "))
  if("Seurat" %in% class(object)){
    metadata <- object@meta.data
  }else{
    metadata <- object
  }
  metadata <- as.data.table(metadata)
  metadata[[sample.by]] <- as.character( metadata[[sample.by]])

  metrics <- intersect(metrics, colnames(metadata))
  # metrics=c("per_feature_count_RNA", "per_feature_count_ADT")
  stat_table <- metadata[, lapply(.SD, function(x) {
    if(sum(is.na(x)) == length(x)  ){
      return( as.numeric(rep(NA,6) ) )
    }else{
      x=stats::na.omit(x)
      return(as.numeric(summary(x)[1:6]))
    }
  }) , .SDcols = metrics, by = sample.by]


  stat_table <- stat_table[, stat_method := rep( c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max"), length(unique(get(sample.by)))) ]
  setcolorder(stat_table, c("stat_method", names(stat_table)[-length(names(stat_table))]))

  stat_list <- lapply(metrics, function(x){
    out <- as.data.frame(data.table::dcast(stat_table, paste0(sample.by," ~ stat_method"), value.var=x))
    out[, -1] <- lapply(out[, -1], function(x) round(x, 3))
    colnames(out)[1] <- "sample"
    return(out)
  })
  names(stat_list) <- metrics
  return(stat_list)
}


#' Calculate Feature QC Metrics
#'
#' Calculate metrics at the feature level, such as variable features and expression statistics.
#'
#' @param object A Seurat object.
#' @param add.Seurat Logical, if TRUE, adds the computed metrics to the Seurat object. Otherwise, returns a list or data.frame with the metrics.
#' @param verbose A logical value. If `TRUE` (default), the message will be printed.
#'              If `FALSE` , the message will be suppressed (not printed).
#' @param sample.by The column name by which to group samples for metrics computation.
#'
#'
#' @return The Seurat object with feature metrics added if `add.Seurat` is TRUE, otherwise a list or data frame of feature qc metrics.
#' @export
#' @examples
#' \dontrun{
#' data <- CalculateMetricsPerFeature(data, add.Seurat= T)
#' }
CalculateMetricsPerFeature <- function(object, add.Seurat= T, verbose=TRUE, sample.by= "orig.ident"){
  if (!("Seurat" %in% class(object))) {
    stop("Error: Seurat object must be as input!!")
  }

  index <- names(object@assays)
  out <- lapply(index, function(x){

    .log_msg(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   "------CalculateMetricsPerFeature assay:", x), verbose = verbose)

    exp <- Seurat::GetAssayData(object, slot = "counts", assay = x)
    row_indices <- split(1:ncol(exp), object@meta.data[[sample.by]])
    split_matrices <- lapply(row_indices, function(idx) exp[, idx, drop = FALSE])
    out_list <- lapply( names(split_matrices), function(y){

      .log_msg(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     "------CalculateMetricsPerFeature ", y), verbose = verbose)

      return(.calPerFeature(split_matrices[[y]],x))
    })
    names(out_list) <- names(split_matrices)
    return(out_list)
  })
  names(out) <- index
  if(add.Seurat){
    suppressWarnings(Seurat::Misc(object, slot = "SingleCellMQC")$perQCMetrics$perFeature  <- out)
    return(object)
  }
  return(out)
}

.calPerFeature <- function(exp, assay){
  ##
  vst_exp <- Seurat::FindVariableFeatures(exp, verbose=F)
  index <- intersect(c("mean", "variance", "variance.standardized", "vst.mean",
                       "vst.variance", "vst.variance.standardized"), colnames(vst_exp))
  vst_exp <- vst_exp[,index]
  colnames(vst_exp) <- gsub("^vst\\.", "",  colnames(vst_exp) )
  # exp@x <- log(exp@x+1)
  nCell <- Matrix::rowSums(exp>0)

  if(assay=="RNA"){
    ## mean log1p-expression
    if(is(exp, "dgCMatrix") ){
      exp <- as(exp, "IterableMatrix")
    }
    exp <- BPCells::multiply_cols(exp, 1/Matrix::colSums(exp))
    exp <- log1p(exp * 10000) # Log normalization
    stats <- BPCells::matrix_stats(exp, row_stats="variance")
    exp_bind <- data.frame(Feature=rownames(exp), nCell=nCell, pct=nCell/dim(exp)[2], vst_exp,
                           mean_lognorm = stats$row_stats["mean", ],
                           variance_lognorm = stats$row_stats["variance", ],
                           cv_lognorm = sqrt(stats$row_stats["variance", ]) /stats$row_stats["mean", ] )
  }else if(assay=="ADT"){
    ## mean CLR
    exp <- .clrBP(exp)
    stats <- BPCells::matrix_stats(exp, row_stats="variance")
    exp_bind <- data.frame(Feature=rownames(exp), nCell=nCell, pct=nCell/dim(exp)[2], vst_exp,
                           mean_clr = stats$row_stats["mean", ],
                           variance_clr = stats$row_stats["variance", ],
                           cv_clr = sqrt(stats$row_stats["variance", ]) /stats$row_stats["mean", ] )
  } else {
    stop("Unsupported assay type provided. Only 'RNA' and 'ADT' are supported.")
  }
  return(exp_bind)
}
