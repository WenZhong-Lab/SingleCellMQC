library(SingleCellMQC)

local({
  # Override dependency probes on a local function copy, leaving namespaces intact.
  check_version <- function(installed = TRUE, sccdc = "1.4", seurat = "5.0.0", expected) {
    run <- RunCorrection_scCDC
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
  check_version(installed = FALSE, expected = "scCDC >= 1.4")
  check_version(sccdc = "1.3", expected = "scCDC >= 1.4")
  check_version(seurat = "4.4.0", expected = "Seurat >= 5.0.0")

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
  }
})
