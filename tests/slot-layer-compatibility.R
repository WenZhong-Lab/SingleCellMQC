library(SingleCellMQC)

# Turn the deprecated assay slot argument into an error even with older v5
# installations, where it would otherwise only emit an occasional warning.
local({
  old_options <- options(lifecycle_verbosity = "warning")
  on.exit(options(old_options))
  withCallingHandlers({
    set.seed(42)
    counts <- Matrix::Matrix(matrix(rpois(60 * 80, 3), nrow = 60), sparse = TRUE)
    dimnames(counts) <- list(paste0("gene", seq_len(60)), paste0("cell", seq_len(80)))
    object <- SeuratObject::CreateSeuratObject(counts)
    object$orig.ident <- rep(c("sample1", "sample2"), each = 40)
    adt <- counts
    rownames(adt) <- paste0("protein", seq_len(nrow(adt)))
    object[["ADT"]] <- if (inherits(object[["RNA"]], "Assay5")) {
      SeuratObject::CreateAssay5Object(adt)
    } else {
      SeuratObject::CreateAssayObject(adt)
    }

    # Public slot remains supported, and layer takes precedence in both helpers.
    get_matrix <- SingleCellMQC:::getMatrix.Seurat
    set_matrix <- SingleCellMQC:::setAssayData
    changed <- set_matrix(object, slot = "counts", layer = "data", new.data = counts * 2)
    stopifnot(
      isTRUE(all.equal(as.matrix(get_matrix(changed, slot = "counts")), as.matrix(counts))),
      isTRUE(all.equal(as.matrix(get_matrix(changed, slot = "counts", layer = "data")),
                       as.matrix(counts * 2)))
    )

    # Exercise the PCA and UMAP dispatch that RunScType uses after splitting.
    for (assay in c("RNA", "ADT")) {
      prefix <- tolower(assay)
      result <- RunPipeline(object, preprocess = paste0(prefix, ".umap"),
                            nfeatures = 55, exclude_patterns = NULL, dims = 1:5)
      stopifnot(
        nrow(SeuratObject::Embeddings(result, paste0(prefix, ".pca"))) == ncol(object),
        nrow(SeuratObject::Embeddings(result, paste0(prefix, ".umap"))) == ncol(object),
        "seurat_clusters" %in% colnames(result[[]])
      )
    }

    # Force both conversion branches and read back normalized/scaled BPCells layers.
    root <- tempfile("SingleCellMQC-slot-layer-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    for (assay in c("RNA", "ADT")) {
      result <- RunPipeline(object, preprocess = paste0(tolower(assay), ".pca"),
                            nfeatures = 55, exclude_patterns = NULL, BP.cutoff = 0,
                            BPdir = file.path(root, assay), BPtmpdir = root,
                            add.BPscale_data = TRUE)
      stopifnot(
        inherits(get_matrix(result, assay = assay, slot = "counts"), "IterableMatrix"),
        all(is.finite(Matrix::colSums(get_matrix(result, assay = assay, slot = "data")))),
        all(is.finite(Matrix::colSums(get_matrix(result, assay = assay, slot = "scale.data"))))
      )
    }

    metrics <- CalculateMetricsPerFeature(object, add.Seurat = FALSE, verbose = FALSE)
    stopifnot(identical(names(metrics), c("RNA", "ADT")))
    for (assay in names(metrics)) {
      stopifnot(identical(names(metrics[[assay]]), c("sample1", "sample2")))
      for (sample in names(metrics[[assay]])) {
        expected <- Matrix::rowSums(get_matrix(object, assay = assay)[,
                                      object$orig.ident == sample] > 0)
        stopifnot(isTRUE(all.equal(unname(metrics[[assay]][[sample]]$nCell),
                                  unname(expected))))
      }
    }

    cutoffs <- data.frame(Sample = "sample1", feature1 = "protein1", cutoff1 = 1,
                          feature2 = "protein2", cutoff2 = 1)
    plot <- PlotADTCutoff(object, cutoffs, feature1 = "protein1", feature2 = "protein2",
                          Sample = "sample1")
    stopifnot(inherits(plot, "ggplot"), nrow(plot$data) == 40)
  }, lifecycle_warning_deprecated = function(w) {
    msg <- conditionMessage(w)
    if (utils::packageVersion("SeuratObject") >= "5.0.0" &&
        grepl("slot", msg, fixed = TRUE) && grepl("(Get|Set)AssayData", msg)) {
      stop(msg, call. = FALSE)
    }
  })
})
