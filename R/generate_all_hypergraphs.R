#' Generate Hypergraphs from Clustering Results
#'
#' This function computes hypergraph incidence matrices (H) and pairwise co-association matrices (HH)
#' from individual clustering algorithm results. It also computes weighted versions of these matrices
#' based on internal validation metrics.
#'
#' @rdname generate_all_hypergraphs
#' 
#' @param object A \code{ScEnsemble} object containing individual clustering results and validation metrics
#' @param verbose Logical; whether to print progress messages (default is TRUE)
#'
#' @return An updated \code{ScEnsemble} object with added slot \code{@hypergraphs} storing the H and HH matrices.
#' @export
setMethod("generate_all_hypergraphs", "ScEnsemble", 
          function(object, verbose = TRUE) {
            
  if (verbose) message("Generating hypergraph matrices...")
  
  clustering_results <- object@individual_results@clustering_results
  weight_list <- object@validation_metrics@normalized_indices
  algorithms <- names(clustering_results)
  
  create_hypergraph_matrix <- function(clustering_vector) {
    n_samples <- length(clustering_vector)
    n_clusters <- length(unique(clustering_vector))
    H <- matrix(0, nrow = n_samples, ncol = n_clusters)
    for (i in seq_len(n_samples)) {
      cluster <- clustering_vector[i]
      if (!is.na(cluster) && cluster > 0) {
        H[i, cluster] <- 1
      }
    }
    return(H)
  }
  
  H_list <- lapply(algorithms, function(algo) {
    create_hypergraph_matrix(clustering_results[[algo]])
  })
  names(H_list) <- algorithms
  
  H_combined <- do.call(cbind, H_list)
  HH_combined <- H_combined %*% t(H_combined)
  HH_combined <- HH_combined / length(algorithms)
  
  apply_weight_and_combine_fixed <- function(metric_name, weight_list, H_list, algorithms) {
    weighted_H_list <- lapply(algorithms, function(algo) {

      weight <- weight_list[[algo]][[metric_name]] 
      
      if (is.null(weight) || is.na(weight) || !is.numeric(weight)) {
        warning(sprintf("Invalid weight for metric '%s' and algorithm '%s'. Using weight = 0.", 
                        metric_name, algo))
        weight <- 0
      }
      H_list[[algo]] * sqrt(weight)
    })
    
    H_weighted <- do.call(cbind, weighted_H_list)
    HH_weighted <- H_weighted %*% t(H_weighted)
    
    
    total_weight <- sum(sapply(algorithms, function(algo) {
      w <- weight_list[[algo]][[metric_name]]
      if (is.null(w) || is.na(w) || !is.numeric(w)) 0 else w
    }), na.rm = TRUE)
    
    if (total_weight > 0) {
      HH_weighted <- HH_weighted / total_weight
    } else {
      warning(sprintf("Total weight for metric '%s' is 0 or invalid.", metric_name))
    }
    
    return(list(H = H_weighted, HH = HH_weighted))
  }
  
  all_metrics <- unique(unlist(lapply(weight_list, names)))
  hypergraph_results <- list(
    H = H_combined,
    HH = HH_combined
  )
  
  for (metric in all_metrics) {
    result <- apply_weight_and_combine_fixed(metric, weight_list, H_list, algorithms)
    hypergraph_results[[paste0("H_", metric)]] <- result$H
    hypergraph_results[[paste0("HH_", metric)]] <- result$HH
  }
  
  object@hypergraphs <- hypergraph_results
  return(object)
})
