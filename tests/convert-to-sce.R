library(SingleCellMQC)

local({
  counts <- Matrix::Matrix(matrix(seq_len(24), nrow = 4), sparse = TRUE)
  dimnames(counts) <- list(paste0("gene", 1:4), paste0("cell", 1:6))
  object <- SeuratObject::CreateSeuratObject(counts)
  object$sample <- rep(c("a", "b"), each = 3)
  object[["RNA"]][["feature_label"]] <- setNames(paste0("label", 1:4), rownames(counts))
  object <- SingleCellMQC:::setAssayData(object, slot = "data", new.data = counts * 2)
  scaled <- as.matrix(counts[1:2, , drop = FALSE]) - 10
  object <- SingleCellMQC:::setAssayData(object, slot = "scale.data", new.data = scaled)
  embeddings <- matrix(seq_len(12), nrow = 6,
                       dimnames = list(colnames(counts), c("PC_1", "PC_2")))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings, key = "PC_", assay = "RNA")

  # The default still delegates to Seurat, retaining its normal conversion output.
  expected <- suppressWarnings(Seurat::as.SingleCellExperiment(object, assay = "RNA"))
  default <- ConvertToSCE(object)
  stopifnot(isTRUE(all.equal(default, expected)))

  for (args in list(list(slot = "data"), list(layer = "data"),
                    list(slot = "scale.data", layer = "data"))) {
    selected <- do.call(ConvertToSCE, c(list(object = object), args))
    stopifnot(
      isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(selected, "counts")),
                       as.matrix(counts * 2))),
      identical(as.character(selected$sample), unname(object$sample)),
      identical(as.character(SummarizedExperiment::rowData(selected)$feature_label),
                object[["RNA"]][[]]$feature_label)
    )
  }
  selected <- ConvertToSCE(object, layer = "scale.data")
  stopifnot(isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(selected, "counts")), scaled)))
  error <- tryCatch(ConvertToSCE(object, layer = "missing"), error = identity)
  stopifnot(inherits(error, "error"), grepl("missing", conditionMessage(error)))
  raw <- SeuratObject::CreateSeuratObject(counts)
  error <- tryCatch(ConvertToSCE(raw, layer = "data"), error = identity)
  stopifnot(inherits(error, "error"))

  legacy <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssayObject(counts))
  legacy <- SingleCellMQC:::setAssayData(legacy, slot = "data", new.data = counts * 4)
  selected <- ConvertToSCE(legacy, slot = "data")
  stopifnot(isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(selected, "counts")),
                             as.matrix(counts * 4))))

  # A custom v5 layer may contain only a subset of cells and features.
  if (inherits(object[["RNA"]], "Assay5")) {
    corrected <- counts[3:4, c(6, 2, 4), drop = FALSE] * 3
    object <- SingleCellMQC:::setAssayData(object, layer = "corrected", new.data = corrected)
    # Seurat orders layer cells by its own cell map when assigning the matrix.
    corrected <- SingleCellMQC:::getMatrix.Seurat(object, layer = "corrected")
    selected <- ConvertToSCE(object, layer = "corrected")
    stopifnot(
      isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(selected, "counts")),
                       as.matrix(corrected))),
      identical(rownames(SummarizedExperiment::colData(selected)), colnames(corrected)),
      identical(as.character(selected$sample), object[[]][colnames(corrected), "sample"]),
      identical(as.character(SummarizedExperiment::rowData(selected)$feature_label),
                object[["RNA"]][[]][rownames(corrected), "feature_label"])
    )

    root <- tempfile("sce-layer-")
    on.exit(unlink(root, recursive = TRUE))
    disk <- BPCells::write_matrix_dir(corrected, root)
    object <- SingleCellMQC:::setAssayData(object, layer = "corrected", new.data = disk)
    selected <- ConvertToSCE(object, layer = "corrected")
    stopifnot(
      inherits(SummarizedExperiment::assay(selected, "counts"), "IterableMatrix"),
      isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(selected, "counts")),
                       as.matrix(corrected)))
    )
  }

  stopifnot(identical(ConvertToSCE(default), default))
  stopifnot(isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(ConvertToSCE(counts), "counts")),
                             as.matrix(counts))))
})
