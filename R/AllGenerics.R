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
                    verbose = TRUE,
                    ...) 
             standardGeneric("calculate_all_validation_indices"))

#' @title Generate Hypergraphs
#' @description Generate hypergraph representations from clustering results
#' 
#' @param object ScEnsemble object with completed validation indices
#' @param algorithms Character vector of algorithms to include
#' @param ... Additional arguments
#' 
#' @return ScEnsemble object with hypergraphs slot filled
#' @export
setGeneric("generate_all_hypergraphs", function(object, ...) standardGeneric("generate_all_hypergraphs"))

#' @title Ensemble Clustering Algorithms
#' @description Apply ensemble clustering methods
#' 
#' @param object ScEnsemble object with completed hypergraph lists
#' @param k Number of clusters (if NULL, will be estimated)
#' @param ensemble_methods Character vector of ensemble methods to apply
#' @param ... Additional arguments
#' 
#' @return ScEnsemble object with ensemble_results slot filled
#' @export
setGeneric("ensemble_clustering", 
           function(object,
                    expression_data,
                    true_labels = NULL,
                    H_matrix,
                    k = NULL,
                    ensemble_methods = c("CSPA_Hc", "CSPA_Louvain", "CSPA_Leiden",
                                         "MCLA_Hc", "MCLA_Louvain", "MCLA_Leiden", "HGSC"),
                    ...) 
             standardGeneric("ensemble_clustering"))
