library(SingleCellMQC)

local({
  metadata <- data.frame(cluster = factor(rep(c("18", "g18", "T_cell", "T-cell"), 2),
                                          levels = c("unused", "18", "g18", "T_cell", "T-cell")))
  mapped <- SingleCellMQC:::.scCDC_cluster_idents(metadata, "cluster")
  stopifnot(
    all(grepl("^cluster[0-9]+$", as.character(mapped))),
    nlevels(mapped) == 4L,
    identical(outer(as.character(mapped), as.character(mapped), "=="),
              outer(as.character(metadata$cluster), as.character(metadata$cluster), "=="))
  )
  reversed <- metadata
  reversed$cluster <- factor(as.character(metadata$cluster),
                             levels = c("unused", "T-cell", "T_cell", "g18", "18"))
  reordered <- SingleCellMQC:::.scCDC_cluster_idents(reversed, "cluster")
  stopifnot(identical(levels(reordered), paste0("cluster", 1:4)),
            identical(as.integer(reordered), as.integer(droplevels(reversed$cluster))))
  reversed$cluster <- as.character(reversed$cluster)
  stopifnot(identical(as.integer(SingleCellMQC:::.scCDC_cluster_idents(reversed, "cluster")),
                      match(reversed$cluster, unique(reversed$cluster))))
  for (bad in list(data.frame(cluster = NA_character_), data.frame(cluster = ""),
                   data.frame(cluster = " "), data.frame(other = "18"))) {
    error <- tryCatch(SingleCellMQC:::.scCDC_cluster_idents(bad, "cluster"), error = identity)
    stopifnot(inherits(error, "error"), grepl("group.by", conditionMessage(error), fixed = TRUE))
  }

  # Override dependency probes on a local function copy, leaving namespaces intact.
  check_version <- function(run, installed = TRUE, sccdc = "1.4", seurat = "5.0.0", expected) {
    env <- new.env(parent = environment(run))
    env$requireNamespace <- function(...) installed
    env$packageVersion <- function(pkg) {
      package_version(if (pkg == "scCDC") sccdc else seurat)
    }
    env$splitObject <- function(...) stop("Data splitting must not start")
    environment(run) <- env
    error <- tryCatch(run(NULL), error = identity)
    stopifnot(inherits(error, "error"), grepl(expected, conditionMessage(error), fixed = TRUE))
  }
  for (run in list(FindContaminationFeature, RunCorrection_scCDC)) {
    check_version(run, installed = FALSE, expected = "scCDC >= 1.4")
    check_version(run, sccdc = "1.3", expected = "scCDC >= 1.4")
    check_version(run, seurat = "4.4.0", expected = "Seurat >= 5.0.0")
  }

  if (!requireNamespace("scCDC", quietly = TRUE) ||
      packageVersion("scCDC") < "1.4" || packageVersion("Seurat") < "5.0.0") {
    message("Skipping scCDC integration checks: scCDC >= 1.4 and Seurat >= 5 are required.")
  } else {
    # scCDC 1.4 calls Seurat helpers that it does not import into its namespace.
    library(Seurat)
    root <- tempfile("sccdc-compatibility-")
    dir.create(root)
    old_wd <- setwd(root)
    old_options <- options(Seurat.object.assay.version = "v5")
    on.exit({
      options(old_options)
      setwd(old_wd)
      unlink(root, recursive = TRUE)
    })

    # Counts-only Assay5 must reach detection without legacy-assay coercion.
    set.seed(42)
    detection_counts <- matrix(rnbinom(200 * 240,
      mu = rep(exp(seq(log(0.2), log(50), length.out = 200)), 240), size = 1), 200)
    detection_counts[1, ] <- 5
    detection_counts[2, ] <- 10
    dimnames(detection_counts) <- list(paste0("protein", 1:200), paste0("cell", 1:240))
    detection_counts <- Matrix::Matrix(detection_counts, sparse = TRUE)
    detection_object <- SeuratObject::CreateSeuratObject(detection_counts, assay = "ADT")
    detection_object$orig.ident <- "sample"
    detection_object$seurat_clusters <- factor(rep(c("18", "g18"), each = 120),
                                               levels = c("g18", "18"))
    detection_metadata <- detection_object[[]]
    detection_idents <- SeuratObject::Idents(detection_object)
    stopifnot(inherits(detection_object[["ADT"]], "Assay5"),
              identical(SeuratObject::Layers(detection_object[["ADT"]]), "counts"))
    set.seed(100)
    features <- FindContaminationFeature(detection_object, assay = "ADT", restriction_factor = 0.1)
    stopifnot(setequal(features, c("protein1", "protein2")),
              identical(detection_object[[]], detection_metadata),
              identical(SeuratObject::Idents(detection_object), detection_idents),
              isTRUE(all.equal(SingleCellMQC:::getMatrix.Seurat(detection_object, assay = "ADT"),
                               detection_counts)))

    set.seed(42)
    counts <- matrix(rpois(20 * 240, 10), nrow = 20,
                     dimnames = list(paste0("protein", 1:20), paste0("cell", 1:240)))
    clusters <- rep(rep(0:1, each = 60), 2)
    counts[1, ] <- rpois(240, ifelse(clusters == 0, 5, 50))
    counts <- Matrix::Matrix(counts, sparse = TRUE)
    object <- SeuratObject::CreateSeuratObject(counts, assay = "ADT")
    object$orig.ident <- rep(c("sample1", "sample2"), each = 120)
    object$seurat_clusters <- clusters
    object$custom_clusters <- paste0("cluster", clusters)
    get_counts <- function(x, assay) SingleCellMQC:::getMatrix.Seurat(x, assay = assay)

    # Real correction exercises numeric cluster labels and per-sample features.
    corrected <- RunCorrection_scCDC(object, assay = "ADT",
                                     features = list(sample1 = "protein1", sample2 = "protein1"))
    actual <- get_counts(corrected, "scCDC_ADT")[rownames(counts), colnames(counts)]
    stopifnot(
      inherits(corrected[["scCDC_ADT"]], "Assay5"),
      isTRUE(all.equal(get_counts(corrected, "ADT"), counts)),
      isTRUE(all.equal(actual[-1, ], counts[-1, ])),
      all(actual[1, ] >= 0), all(actual[1, ] <= counts[1, ]),
      any(actual[1, ] < counts[1, ])
    )

    # The selected metadata column, rather than the original Idents, drives correction.
    single <- subset(object, cells = colnames(object)[1:120])
    single$seurat_clusters <- NULL
    corrected <- RunCorrection_scCDC(single, assay = "ADT", features = "protein1",
                                     group.by = "custom_clusters")
    stopifnot(isTRUE(all.equal(get_counts(corrected, "scCDC_ADT"),
                               actual[, colnames(single)])))

    # Labels that Seurat/scCDC rewrite must still describe the same partition.
    for (labels in list(c("0", "18"), c("g18", "18"),
                        c("T_cell", "T-cell"), c("all", "B cell"))) {
      single$custom_clusters <- factor(rep(labels, each = 60),
                                       levels = c("unused", labels))
      original_metadata <- single[[]]
      original_idents <- SeuratObject::Idents(single)
      corrected <- RunCorrection_scCDC(single, assay = "ADT", features = "protein1",
                                       group.by = "custom_clusters")
      stopifnot(
        isTRUE(all.equal(get_counts(corrected, "scCDC_ADT"), actual[, colnames(single)])),
        identical(corrected[[]][, names(original_metadata)], original_metadata),
        identical(SeuratObject::Idents(corrected), original_idents)
      )
    }
  }
})
