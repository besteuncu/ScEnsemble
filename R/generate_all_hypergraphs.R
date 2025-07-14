#' Generate All Hypergraphs from Clustering Results
#'
#' This function creates hypergraph matrices from clustering results and applies
#' various weighting schemes based on validation indices. It generates both
#' standard and weighted hypergraph representations for ensemble clustering.
#'
#' @param clustering_results A named list containing clustering results from
#'   different algorithms. Each element should be a vector of cluster assignments.
#' @param weight_list A named list of normalized validation indices used for
#'   weighting the hypergraphs. Each element should be a named list with
#'   algorithm names as keys and weights as values.
#' @param algorithms A character vector specifying which algorithms to include
#'   in the hypergraph generation. Default includes common single-cell
#'   clustering algorithms.
#'
#' @return A list containing:
#' \describe{
#'   \item{H}{Standard combined hypergraph matrix}
#'   \item{HH}{Standard hypergraph similarity matrix}
#'   \item{H_[metric]}{Weighted hypergraph matrices for each validation metric}
#'   \item{HH_[metric]}{Weighted hypergraph similarity matrices for each validation metric}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates binary hypergraph matrices for each clustering algorithm
#'   \item Combines matrices to form standard hypergraph representation
#'   \item Applies validation-based weights to create weighted hypergraphs
#'   \item Computes hypergraph similarity matrices (H * H^T)
#'   \item Normalizes matrices appropriately
#' }
#'
#' The hypergraph matrix H has dimensions n_samples x n_hyperedges, where
#' each hyperedge represents a cluster from a specific algorithm.
#'
#' @examples
#' \dontrun{
#' # Assuming clustering_results and weight_list are available
#' hypergraphs <- generate_all_hypergraphs(
#'   clustering_results = clustering_results,
#'   weight_list = validation_weights,
#'   algorithms = c("SC3", "CIDR", "Seurat")
#' )
#'
#' # Access standard hypergraph
#' standard_H <- hypergraphs$H
#'
#' # Access weighted hypergraph for silhouette metric
#' weighted_H_sil <- hypergraphs$H_sil
#' }
#'
#' @seealso
#' \code{\link{run_individual_algorithms}} for generating clustering results
#' \code{\link{calculate_all_validation_indices}} for computing validation weights
#' \code{\link{ensemble_clustering_algorithms}} for ensemble clustering
#'
#' @export
generate_all_hypergraphs <- function(clustering_results,
                                     weight_list,
                                     algorithms = c("SC3", "CIDR", "Seurat",
                                                    "SIMLR", "TSNE_Kmeans",
                                                    "Monocle", "RaceID")) {

  # Input validation
  if (!is.list(clustering_results)) {
    stop("clustering_results must be a list")
  }

  if (!is.list(weight_list)) {
    stop("weight_list must be a list")
  }

  if (!is.character(algorithms)) {
    stop("algorithms must be a character vector")
  }

  # Check if all algorithms are present in clustering_results
  missing_algorithms <- setdiff(algorithms, names(clustering_results))
  if (length(missing_algorithms) > 0) {
    stop("Missing algorithms in clustering_results: ",
         paste(missing_algorithms, collapse = ", "))
  }

  message("Generating hypergraph matrices...")

  # 1. Generate the clustering matrix for each algorithm
  create_hypergraph_matrix <- function(clustering_vector) {
    n_samples <- length(clustering_vector)
    n_clusters <- length(unique(clustering_vector))
    H <- matrix(0, nrow = n_samples, ncol = n_clusters)
    for (i in 1:n_samples) {
      H[i, clustering_vector[i]] <- 1
    }
    return(H)
  }

  # 2. Generate all H matrices
  H_list <- lapply(algorithms, function(algo) {
    create_hypergraph_matrix(clustering_results[[algo]])
  })
  names(H_list) <- algorithms

  # 3. Calculate standard H and HH
  H_combined <- do.call(cbind, H_list)
  HH_combined <- H_combined %*% t(H_combined)
  HH_combined <- HH_combined / length(algorithms)

  # 4. Helper function that handles weighted metrics
  apply_weight_and_combine <- function(metric_name) {
    weighted_H_list <- lapply(algorithms, function(algo) {
      weight <- weight_list[[metric_name]][[algo]]
      H_list[[algo]] * sqrt(weight)
    })
    H_weighted <- do.call(cbind, weighted_H_list)
    HH_weighted <- H_weighted %*% t(H_weighted)
    HH_weighted <- HH_weighted / sum(unlist(weight_list[[metric_name]][algorithms]))
    return(list(H = H_weighted, HH = HH_weighted))
  }

  # 5. Calculate weighted H and HH for each metric
  metrics <- names(weight_list)
  hypergraph_results <- list(
    H = H_combined,
    HH = HH_combined
  )

  for (metric in metrics) {
    result <- apply_weight_and_combine(metric)
    hypergraph_results[[paste0("H_", metric)]] <- result$H
    hypergraph_results[[paste0("HH_", metric)]] <- result$HH
  }

  return(hypergraph_results)
}

#' Extract Hypergraph Lists from Results
#'
#' Helper function to extract and organize hypergraph matrices into separate
#' lists for H and HH matrices, replicating the original workflow.
#'
#' @param hypergraph_results Output from \code{\link{generate_all_hypergraphs}}
#' @param metrics Character vector of metric names to extract. Default extracts
#'   standard metrics: silhouette, ch, db, dunn, average.
#'
#' @return A list containing:
#' \describe{
#'   \item{H_list}{Named list of hypergraph matrices}
#'   \item{HH_list}{Named list of hypergraph similarity matrices}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate hypergraphs
#' hypergraphs <- generate_all_hypergraphs(
#'   clustering_results = results_ind$clustering_results,
#'   weight_list = results_val$normalized_indices,
#'   algorithms = c("SC3", "CIDR", "Seurat")
#' )
#'
#' # Extract organized lists
#' extracted <- extract_hypergraph_lists(hypergraphs)
#' H_list <- extracted$H_list
#' HH_list <- extracted$HH_list
#' }
#'
#' @export
extract_hypergraph_lists <- function(hypergraph_results,
                                     metrics = c("sil", "ch", "db", "dunn", "ave")) {

  if (!is.list(hypergraph_results)) {
    stop("hypergraph_results must be a list")
  }

  # Create H_list following original pattern
  H_list <- list(
    standard = hypergraph_results$H
  )

  # Add weighted matrices for each metric
  for (metric in metrics) {
    key <- paste0("H_", metric)
    if (key %in% names(hypergraph_results)) {
      H_list[[metric]] <- hypergraph_results[[key]]
    }
  }

  # Create HH_list following original pattern
  HH_list <- list(
    standard = hypergraph_results$HH
  )

  # Add weighted matrices for each metric
  for (metric in metrics) {
    key <- paste0("HH_", metric)
    if (key %in% names(hypergraph_results)) {
      HH_list[[metric]] <- hypergraph_results[[key]]
    }
  }

  return(list(H_list = H_list, HH_list = HH_list))
}
