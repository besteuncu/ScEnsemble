#' Create a ScEnsemble object
#'
#' Initializes a ScEnsemble object from a SingleCellExperiment and optional annotation vector.
#'
#' @param sce A SingleCellExperiment object
#' @param annotation A character vector of true labels (optional)
#' 
#' @importFrom methods is new
#' @import scRNAseq
#'
#' @return A ScEnsemble S4 object
#' 
#' @examples
#' # Load required packages
#' library(scRNAseq)
#' library(SingleCellExperiment)
#' 
#' # Load example data
#' Pollen <- PollenGliaData()
#' 
#' # Extract annotations if available
#' ann <- colData(Pollen)[["Inferred Cell Type"]]
#' 
#' # Create ScEnsemble object
#' scens <- CreateScEnsemble(Pollen, ann)
#' 
#' @export
CreateScEnsemble <- function(sce, annotation = NULL) {
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
      validation_metrics = new("ValidationResults",
                               validation_indices = list(),
                               normalized_indices = list()
                               ),
      hypergraphs = list(),
      ensemble_results = new("EnsembleResults",
                             ensemble_clusters = list() ,
                             ensemble_ari = list(),
                             ensemble_quality = list()))
}
