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
#' @slot performance_metrics A list of performance evaluations
#' @export

setClassUnion("numericORNULL", c("numeric", "NULL"))

setClass("ScEnsemble",
         slots = list(
           sce = "SingleCellExperiment",
           annotation = "numericORNULL",
           individual_results = "IndividualResults",
           validation_metrics = "ValidationResults",
           hypergraphs = "list",
           hypergraph_lists = "list",
           ensemble_results = "list",
           performance_metrics = "list"
         )
)

setClass("IndividualResults",
         slots = list(
           clustering_results = "list",
           ari_scores = "list",
           runtimes = "list",
           errors = "list",
           embedding_data = "list"
         ))


setClass("ValidationResults",
         slots = list(
           validation_indices = "list",
           normalized_indices = "list",
           average_index = "numeric"
         ))
