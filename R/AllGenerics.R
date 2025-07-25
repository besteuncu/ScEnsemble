#' Generic functions for ScEnsemble package
#' 
#' @name ScEnsemble-generics
#' @rdname ScEnsemble-generics

#' @title Run Individual Clustering Algorithms
#' @description Apply multiple individual clustering algorithms to single-cell data
#' 
#' @param object ScEnsemble object
#' @param algorithms Character vector of algorithms to run. Options: "SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID"
#' @param seed Random seed for reproducibility
#' @param verbose Logical, whether to print progress messages
#' @param n_cores Number of cores to use for parallel processing
#' @param ... Additional arguments
#' 
#' @return ScEnsemble object with individual_results slot filled
#' @export
setGeneric("run_individual_algorithms", 
           function(object, 
                    algorithms = c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID"),
                    seed = 42, 
                    verbose = TRUE, 
                    n_cores = 1, 
                    ...) 
             standardGeneric("run_individual_algorithms"))

#' @title Calculate Validation Indices
#' @description Calculate validation indices for clustering results
#' 
#' @param object ScEnsemble object with completed individual algorithms
#' @param verbose Logical, whether to print progress messages
#' @param ... Additional arguments
#' 
#' @return ScEnsemble object with validation_indices slot filled
#' @export
setGeneric("calculate_all_validation_indices", 
           function(object,
                    clustering_results,
                    embedding_data,
                    verbose = TRUE,
                    ...) 
             standardGeneric("calculate_all_validation_indices"))
