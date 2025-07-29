#' Calculate All Validation Indices - Refactored Version
#'
#' This function calculates multiple clustering validation indices for each algorithm's
#' clustering results and returns both raw and normalized scores using a simplified,
#' more robust approach.
#'
#' @param object A ScEnsemble object containing clustering results and embedding data
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#'
#' @return A ScEnsemble object with updated validation_metrics slot containing:
#' \describe{
#'   \item{validation_indices}{Named list of raw validation indices for each algorithm}
#'   \item{normalized_indices}{Named list of normalized validation indices}
#'   \item{average_index}{Named vector of average normalized scores for each algorithm}
#' }
#'
#' @details
#' The function calculates the following validation indices:
#' \itemize{
#'   \item \strong{Silhouette}: Measures how similar each point is to its own cluster compared to other clusters.
#'   \item \strong{Calinski-Harabasz (CH)}: Ratio of between-cluster dispersion to within-cluster dispersion.
#'   \item \strong{Davies-Bouldin (DB)}: Average similarity between clusters.
#'   \item \strong{Dunn}: Ratio of minimum inter-cluster distance to maximum intra-cluster distance.
#' }
#'
#' @importFrom cluster silhouette
#' @importFrom fpc calinhara
#' @importFrom clusterSim index.DB
#' @importFrom clValid dunn
#' @importFrom stats dist
#' 
#' 
#' @export
setMethod("calculate_all_validation_indices", "ScEnsemble", 
          function(object,
                   verbose = TRUE,
                   ...) {
            
            # Retrieve data 
            clustering_results <- object@individual_results@clustering_results
            embedding_data <- object@individual_results@embedding_data
            
            if (length(clustering_results) == 0) {
              stop("clustering_results must be a non-empty list")
            }
            
            if (verbose) {
              cat("Starting validation index calculations for", length(clustering_results), "algorithms\n")
            }
            
            # Prepare result containers
            all_indices <- list()          # For raw metric values
            all_indices_norm <- list()     # For normalized metric values
            avg_scores <- numeric()        # For average scores
            
            # Temporarily store all raw values for global normalization
            all_raw_values <- list(
              silhouette = numeric(),
              ch = numeric(),
              db = numeric(),
              dunn = numeric()
            )
            
            # Calculate raw metrics for each algorithm
            for (method in names(clustering_results)) {
              if (verbose) {
                cat(sprintf("Calculating raw indices for method: %s\n", method))
              }
              
              clusters <- clustering_results[[method]]
              embed <- embedding_data[[method]]
              
              if (is.null(embed)) {
                if (verbose) {
                  cat(sprintf("Warning: No embedding data found for method: %s. Skipping.\n", method))
                }
                next
              }
              
              raw_metrics <- calculate_single_algorithm_metrics(embed, clusters, method, verbose)
              all_indices[[method]] <- raw_metrics
              
              for (metric_name in names(raw_metrics)) {
                if (!is.na(raw_metrics[[metric_name]])) {
                  all_raw_values[[metric_name]] <- c(all_raw_values[[metric_name]], raw_metrics[[metric_name]])
                }
              }
            }
            
            # Global normalization
            if (verbose) {
              cat("Performing global normalization...\n")
            }
            
            normalized_metrics <- perform_global_normalization(all_raw_values, verbose)
            
            # Assign normalized values back to each method
            for (method in names(all_indices)) {
              method_normalized <- list()
              raw_values <- all_indices[[method]]
              
              for (metric_name in names(raw_values)) {
                raw_val <- raw_values[[metric_name]]
                if (!is.na(raw_val)) {
                  method_normalized[[metric_name]] <- calculate_normalized_value(
                    raw_val, metric_name, normalized_metrics
                  )
                } else {
                  method_normalized[[metric_name]] <- 0  # Assign 0 to NA values
                }
              }
              
              valid_scores <- unlist(method_normalized)
              method_normalized[["average"]] <- mean(valid_scores, na.rm = TRUE)
              all_indices_norm[[method]] <- method_normalized
            }
           
            # Create result object and return
            validation_metrics <- new("ValidationResults",
                                      validation_indices = all_indices,
                                      normalized_indices = all_indices_norm)
            
            object@validation_metrics <- validation_metrics
            
            if (verbose) {
              cat("Validation index calculation completed successfully!\n")
            }
            
            return(object)
          })

