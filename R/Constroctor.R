#' Create a ScEnsemble object
#'
#' Initializes a ScEnsemble object from a SingleCellExperiment and optional annotation vector.
#'
#' @param sce A SingleCellExperiment object
#' @param annotation A character vector of true labels (optional)
#'
#' @return A ScEnsemble S4 object
#' @export
ScEnsemble <- function(sce, annotation = NULL) {
  if (!is(sce, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object.")
  }
  
  if (!is.null(annotation)) {
    if (!is.character(annotation) || length(annotation) != ncol(sce)) {
      stop("annotation must be a character vector with length equal to number of cells.")
    }
    annotation <- as.numeric(factor(annotation))
  }
  
  new("ScEnsemble",
      sce = sce,
      annotation = annotation,
      individual_results = new("IndividualResults",
                               clustering_results = list(),
                               ari_scores = list(),
                               runtimes = list(),
                               errors = list(),
                               embedding_data = list()
      ),
      validation_metrics = list(),
      hypergraphs = list(),
      hypergraph_lists = list(),
      ensemble_results = list(),
      performance_metrics = list())
}
