#' Pollen single-cell RNA-seq data
#'
#' Single-cell RNA-seq data from pollen grains, obtained from the scRNAseq package.
#' This dataset contains expression data for multiple cell types in pollen development.
#'
#' @format A SummarizedExperiment object with:
#' \describe{
#'   \item{assays}{Gene expression count matrix}
#'   \item{colData}{Cell metadata including cell type annotations}
#'   \item{rowData}{Gene metadata}
#' }
#' @source \code{scRNAseq::PollenGliaData()}
#' @references 
#' Li et al. (2017). Single-cell transcriptomes reveal characteristic features 
#' of human pancreatic islet cell types. EMBO Reports, 18(9), 1589-1600.
#' @examples
#' data(Pollen)
#' dim(Pollen)
#' \donttest{
#' # colData function requires SummarizedExperiment to be loaded
#' library(SummarizedExperiment)
#' table(colData(Pollen)$"Inferred Cell Type")
#' }
"Pollen"

#' Cell type annotations for Pollen data
#'
#' Character vector containing inferred cell type labels for each cell 
#' in the Pollen dataset.
#'
#' @format A character vector with cell type labels:
#' \describe{
#'   \item{length}{Number of cells in the dataset}
#'   \item{values}{Cell type names as character strings}
#' }
#' @source Extracted from \code{colData(PollenGliaData())[["Inferred Cell Type"]]}
#' @examples
#' data(ann)
#' table(ann)
#' length(ann)
"ann"
