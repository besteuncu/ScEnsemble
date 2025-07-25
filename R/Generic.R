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
#' @param object SCEnsemble object with completed individual algorithms
#' @param verbose Logical, whether to print progress messages
#' @param ... Additional arguments
#' 
#' @return SCEnsemble object with validation_indices slot filled
#' @export
setGeneric("calculate_all_validation_indices", 
           function(object, verbose = TRUE, ...) 
             standardGeneric("calculate_all_validation_indices"))

#' @title Generate Hypergraphs
#' @description Generate hypergraph representations from clustering results
#' 
#' @param object SCEnsemble object with completed validation indices
#' @param algorithms Character vector of algorithms to include
#' @param ... Additional arguments
#' 
#' @return SCEnsemble object with hypergraphs slot filled
#' @export
setGeneric("generate_all_hypergraphs", 
           function(object, 
                    algorithms = c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID"),
                    ...) 
             standardGeneric("generate_all_hypergraphs"))

#' @title Extract Hypergraph Lists
#' @description Extract hypergraph information for ensemble methods
#' 
#' @param object SCEnsemble object with completed hypergraphs
#' @param metrics Character vector of metrics to use
#' @param ... Additional arguments
#' 
#' @return SCEnsemble object with hypergraph_lists slot filled
#' @export
setGeneric("extract_hypergraph_lists", 
           function(object, 
                    metrics = c("sil", "ch", "db", "dunn", "ave"),
                    ...) 
             standardGeneric("extract_hypergraph_lists"))

#' @title Ensemble Clustering Algorithms
#' @description Apply ensemble clustering methods
#' 
#' @param object SCEnsemble object with completed hypergraph lists
#' @param k Number of clusters (if NULL, will be estimated)
#' @param ensemble_methods Character vector of ensemble methods to apply
#' @param ... Additional arguments
#' 
#' @return SCEnsemble object with ensemble_results slot filled
#' @export
setGeneric("ensemble_clustering_algorithms", 
           function(object, 
                    k = NULL,
                    ensemble_methods = c("CSPA_Hc", "CSPA_Louvain", "CSPA_Leiden", 
                                         "MCLA_Hc", "MCLA_Louvain", "MCLA_Leiden", "HGSC"),
                    ...) 
             standardGeneric("ensemble_clustering_algorithms"))

#' @title Analyze Clustering Performance
#' @description Perform comprehensive performance analysis of clustering results
#' 
#' @param object SCEnsemble object with completed ensemble clustering
#' @param ensemble_algorithms Character vector of ensemble algorithms to analyze
#' @param individual_algorithms Character vector of individual algorithms to analyze
#' @param variants Character vector of validation metric variants
#' @param output_dir Directory to save analysis outputs
#' @param save_plots Logical, whether to save plots
#' @param plot_format Format for saved plots
#' @param plot_width Width of saved plots
#' @param plot_height Height of saved plots
#' @param seed Random seed for reproducibility
#' @param ... Additional arguments
#' 
#' @return SCEnsemble object with performance_analysis slot filled
#' @export
setGeneric("analyze_clustering_performance", 
           function(object,
                    ensemble_algorithms = c("cspa", "cspa_louvain", "cspa_leiden",
                                            "mcla", "mcla_louvain", "mcla_leiden", "HGSC"),
                    individual_algorithms = c("SC3", "CIDR", "Seurat", "SIMLR",
                                              "TSNE_Kmeans", "Monocle", "RaceID"),
                    variants = c("standard", "silhouette", "ch", "db", "average"),
                    output_dir = "clustering_analysis_output",
                    save_plots = TRUE,
                    plot_format = "png",
                    plot_width = 12,
                    plot_height = 8,
                    seed = 42,
                    ...) 
           standardGeneric("analyze_clustering_performance"))

#' @title Get True Labels
#' @description Extract true cell type labels from SCEnsemble object
#' 
#' @param object SCEnsemble object
#' @param ... Additional arguments
#' 
#' @return Character vector of true labels
#' @export
setGeneric("getTrueLabels", 
           function(object, ...) 
             standardGeneric("getTrueLabels"))

#' @title Get Individual Results
#' @description Extract individual clustering algorithm results
#' 
#' @param object SCEnsemble object
#' @param ... Additional arguments
#' 
#' @return List of individual clustering results
#' @export
setGeneric("getIndividualResults", 
           function(object, ...) 
             standardGeneric("getIndividualResults"))

#' @title Get Validation Indices
#' @description Extract validation indices results
#' 
#' @param object SCEnsemble object
#' @param ... Additional arguments
#' 
#' @return List of validation indices
#' @export
setGeneric("getValidationIndices", 
           function(object, ...) 
             standardGeneric("getValidationIndices"))

#' @title Get Ensemble Results
#' @description Extract ensemble clustering results
#' 
#' @param object SCEnsemble object
#' @param ... Additional arguments
#' 
#' @return List of ensemble clustering results
#' @export
setGeneric("getEnsembleResults", 
           function(object, ...) 
             standardGeneric("getEnsembleResults"))

