


#' @importFrom grDevices colorRampPalette
#' @importFrom stats cmdscale dist sd aggregate
#' @importFrom utils packageVersion
#' @importFrom utils write.csv
#' @importFrom graphics par lines axis
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom methods as is
#' @importFrom utils head



utils::globalVariables(c("."))


fake_function <- function(){
  MatrixGenerics::colCounts
  sparkline::sparkline
}


if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "BCR%", "CellType", "Cell_Marker", "Cell_standard", "Filtered", "Group",
      "IGH%", "IGK%", "IGL%", "Library Type", "Metric Name", "Metric Value",
      "N", "PCT", "PC1", "PC2", "PCs", "Post", "Pre", "Pvalue", "Sample",
      "TCR%", "TRA%", "TRB%", "Variable", "X1", "X2", "binomial", "cell",
      "chain", "chain_pair", "cluster", "cosine_similarity", "count", "feature",
      "freq", "gene", "group_combo", "i.count", "id", "isOutlier", "metrics_name",
      "nChain_BCR", "nChain_TCR", "new_name", "num.all", "num.select",
      "orig.ident", "parallel_paramSweep", "productive", "proportion",
      "receptor_subtype", "receptor_type", "row_id", "sample", "scores",
      "stat_method", "total", "value", "values", "variable", "x", "top", "name", "object_temp","method",
      "overall_rank", "asw", "ch", "db", "ari", "nmi", "overall","Feature", "Proportion","variance"
    )
  )
}

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers(progressr::handler_progress(
      format = ":spin :bar :percent :elapsed :eta",
      clear = FALSE
    ))
  }
}

.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("progressr", quietly = TRUE)) {
    packageStartupMessage("Package 'progressr' is not installed. Progress bars will not be displayed.")
  }
}
