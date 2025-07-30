#' IndividualResults Class
#'
#' An S4 class to store individual clustering results
#'
#' @slot clustering_results A list of clustering outputs
#' @slot ari_scores A list of ARI scores
#' @slot runtimes A list of runtime measurements
#' @slot errors A list of error messages
#' @slot embedding_data A list of embedding data
#' 
#' @importFrom methods setClass
#' @export
setClass("IndividualResults",
         slots = list(
           clustering_results = "list",
           ari_scores = "list",
           runtimes = "list",
           errors = "list",
           embedding_data = "list"
         ))

#' ValidationResults Class
#'
#' An S4 class to store validation results
#'
#' @slot validation_indices A list of validation indices
#' @slot normalized_indices A list of normalized indices
#' 
#' @importFrom methods setClass
#' @export
setClass("ValidationResults",
         slots = list(
           validation_indices = "list",
           normalized_indices = "list"
         ))

#' EnsembleResults Class
#'
#' An S4 class to store ensemble clustering results
#'
#' @slot ensemble_clusters A list of ensemble clustering results
#' @slot ensemble_ari A list of ensemble ARI scores
#' @slot ensemble_quality A list of ensemble quality metrics
#' 
#' @importFrom methods setClass
#' @export
setClass("EnsembleResults",
         slots = list(
           ensemble_clusters = "list",
           ensemble_ari = "list",
           ensemble_quality = "list"
         ))

#' ScEnsemble Class
#'
#' An S4 class to store data and results for ensemble clustering
#'
#' @slot sce A SingleCellExperiment object
#' @slot annotation A character vector of true labels
#' @slot individual_results A list of individual clustering outputs
#' @slot validation_metrics A list of internal validation scores
#' @slot hypergraphs A list of hypergraph representations
#' @slot ensemble_results A list of ensemble clustering results
#' 
#' @importFrom methods setClass setClassUnion
#' @export
setClassUnion("numericORNULL", c("numeric", "NULL"))
setClass("ScEnsemble",
         slots = list(
           sce = "SingleCellExperiment",
           annotation = "numericORNULL",
           individual_results = "IndividualResults",
           validation_metrics = "ValidationResults",
           hypergraphs = "list",
           ensemble_results = "EnsembleResults"
         )
)