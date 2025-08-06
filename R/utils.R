#' Calculate Validation Metrics for a Single Clustering Algorithm
#'
#' This function calculates multiple clustering validation indices (Silhouette, 
#' Calinski-Harabasz, Davies-Bouldin, and Dunn) for a single clustering algorithm.
#' It handles edge cases and provides detailed error reporting.
#'
#' @param data A data frame or matrix containing the data points to be clustered.
#' @param clustering A numeric vector containing cluster assignments for each data point.
#' @param method_name A character string specifying the name of the clustering method 
#'   (used for verbose output and error reporting).
#' @param verbose A logical value indicating whether to print detailed progress 
#'   information. Default is \code{FALSE}.
#'
#' @return A list containing four validation indices:
#' \describe{
#'   \item{silhouette}{Average silhouette width (higher is better, range: -1 to 1)}
#'   \item{ch}{Calinski-Harabasz index (higher is better)}
#'   \item{db}{Davies-Bouldin index (lower is better)}
#'   \item{dunn}{Dunn index (higher is better)}
#' }
#' Returns \code{NA} values for indices that cannot be calculated.
#'
#' @details 
#' The function requires at least 2 unique clusters to calculate meaningful indices.
#' If fewer clusters are found, all indices return \code{NA}.
#' 
#' Individual index calculations are wrapped in \code{tryCatch} to handle errors
#' gracefully and continue processing even if one index fails.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' data <- matrix(rnorm(100), ncol = 2)
#' 
#' # Create cluster assignments
#' clustering <- c(rep(1, 25), rep(2, 25))
#' 
#' # Calculate metrics
#' metrics <- calculate_single_algorithm_metrics(
#'   data = data, 
#'   clustering = clustering, 
#'   method_name = "kmeans",
#'   verbose = TRUE
#' )
#' 
#' # View results
#' print(metrics)
#'
#' @importFrom cluster silhouette
#' @importFrom fpc calinhara
#' @importFrom clusterSim index.DB
#' @importFrom stats dist
#' @importFrom clusterCrit intCriteria
#'
#' @keywords internal
#' @export
calculate_single_algorithm_metrics <- function(data, clustering, method_name, verbose = FALSE) {
  
  unique_clusters <- unique(clustering)
  if (length(unique_clusters) < 2) {
    if (verbose) {
      message(sprintf("  %s: Only %d unique cluster(s), returning NA values\n", 
                  method_name, length(unique_clusters)))
    }
    return(list(silhouette = NA, ch = NA, db = NA, dunn = NA))
  }
  
  data_matrix <- as.matrix(data)
  
  if (verbose) {
    message(sprintf("  %s: Processing %d points with %d clusters\n", 
                method_name, nrow(data_matrix), length(unique_clusters)))
  }
  
  sil_index <- tryCatch({
    sil <- silhouette(clustering, dist(data_matrix))
    mean(sil[, 3], na.rm = TRUE)
  }, error = function(e) {
    if (verbose) message(sprintf("    Silhouette calculation failed: %s\n", e$message))
    NA
  })
  
  ch_index <- tryCatch({
    calinhara(data_matrix, clustering, cn = length(unique(clustering)))
  }, error = function(e) {
    if (verbose) message(sprintf("    CH calculation failed: %s\n", e$message))
    NA
  })
  
  db_index <- tryCatch({
    result <- clusterSim::index.DB(data_matrix, clustering, centrotypes = "centroids")
    db_val <- result$DB
  }, error = function(e) {
    if (verbose) message(sprintf("    DB calculation failed: %s\n", e$message))
    NA
  })
  
  dunn_index <- tryCatch({
    dunn(dist(data_matrix), clustering)
  }, error = function(e) {
    if (verbose) message(sprintf("    Dunn calculation failed: %s\n", e$message))
    NA
  })
  
  return(list(
    silhouette = sil_index,
    ch = ch_index,
    db = db_index,
    dunn = dunn_index
  ))
}

#' Perform Global Normalization of Clustering Validation Indices
#'
#' This function determines appropriate normalization parameters for different 
#' clustering validation indices based on their mathematical properties and 
#' the distribution of raw values across all clustering methods.
#'
#' @param all_raw_values A list containing vectors of raw validation index values
#'   for all clustering methods. Expected components are:
#'   \describe{
#'     \item{silhouette}{Vector of silhouette index values}
#'     \item{ch}{Vector of Calinski-Harabasz index values}
#'     \item{db}{Vector of Davies-Bouldin index values}
#'     \item{dunn}{Vector of Dunn index values}
#'   }
#' @param verbose A logical value indicating whether to print normalization 
#'   details. Default is \code{FALSE}.
#'
#' @return A list containing normalization parameters for each validation index.
#'   Each component contains:
#'   \describe{
#'     \item{type}{The normalization method used}
#'     \item{min_val, max_val}{For min-max normalization (Silhouette)}
#'     \item{sum_val}{For sum-based or inverse-sum normalization (CH, DB, Dunn)}
#'   }
#'
#' @details 
#' Different normalization strategies are applied based on index properties:
#' \itemize{
#'   \item \strong{Silhouette}: Min-max normalization with shift (range: -1 to 1)
#'   \item \strong{Calinski-Harabasz}: Sum-based normalization (higher is better)
#'   \item \strong{Davies-Bouldin}: Inverse sum normalization (lower is better)
#'   \item \strong{Dunn}: Sum-based normalization (higher is better)
#' }
#'
#' @examples
#' # Create sample raw validation values from multiple methods
#' raw_values <- list(
#'   silhouette = c(0.7, 0.5, 0.8, 0.6),
#'   ch = c(120.5, 95.2, 150.3, 80.1),
#'   db = c(0.8, 1.2, 0.6, 1.5),
#'   dunn = c(0.25, 0.18, 0.32, 0.15)
#' )
#' 
#' # Calculate normalization parameters
#' norm_params <- perform_global_normalization(raw_values, verbose = TRUE)
#' 
#' # View normalization parameters
#' str(norm_params)
#'
#' @keywords internal
#' @export
perform_global_normalization <- function(all_raw_values, verbose = FALSE) {
  
  normalization_params <- list()
  
  sil_values <- all_raw_values$silhouette
  if (length(sil_values) > 0) {
    normalization_params$silhouette <- list(
      type = "minmax_shifted",
      min_val = -1,
      max_val = 1
    )
    if (verbose) message("  Silhouette: Using min-max normalization with shift\n")
  }
  
  ch_values <- all_raw_values$ch[!is.na(all_raw_values$ch)]
  if (length(ch_values) > 0) {
    ch_sum <- sum(ch_values)
    normalization_params$ch <- list(
      type = "sum_based",
      sum_val = ch_sum
    )
    if (verbose) message("  CH: Using sum-based normalization\n")
  }
  
  db_values <- all_raw_values$db[!is.na(all_raw_values$db) & all_raw_values$db > 0]
  if (length(db_values) > 0) {
    inv_db_values <- 1 / db_values
    inv_db_sum <- sum(inv_db_values)
    normalization_params$db <- list(
      type = "inverse_sum",
      sum_val = inv_db_sum
    )
    if (verbose) message("  DB: Using inverse sum normalization\n")
  }
  
  dunn_values <- all_raw_values$dunn[!is.na(all_raw_values$dunn)]
  if (length(dunn_values) > 0) {
    dunn_sum <- sum(dunn_values)
    normalization_params$dunn <- list(
      type = "sum_based",
      sum_val = dunn_sum
    )
    if (verbose) message("  Dunn: Using sum-based normalization\n")
  }
  
  return(normalization_params)
}

#' Calculate Normalized Value for a Clustering Validation Index
#'
#' This function applies the appropriate normalization transformation to a raw 
#' validation index value based on pre-computed normalization parameters.
#'
#' @param raw_value A numeric value representing the raw validation index score.
#'   Can be \code{NA}.
#' @param metric_name A character string specifying the name of the validation 
#'   index (e.g., "silhouette", "ch", "db", "dunn").
#' @param normalization_params A list containing normalization parameters as 
#'   returned by \code{\link{perform_global_normalization}}.
#'
#' @return A numeric value between 0 and 1 representing the normalized score.
#'   Returns 0 for \code{NA} input values or when normalization parameters 
#'   are unavailable.
#'
#' @details 
#' The normalization is applied according to the method specified in the 
#' normalization parameters:
#' \itemize{
#'   \item \strong{minmax_shifted}: \code{(raw_value + 1) / 2} for Silhouette
#'   \item \strong{sum_based}: \code{raw_value / sum_val} for CH and Dunn
#'   \item \strong{inverse_sum}: \code{(1 / raw_value) / sum_val} for DB
#' }
#'
#' All normalized values are scaled to (0, 1) where higher values indicate 
#' better clustering performance.
#'
#' @examples
#' # First create normalization parameters
#' raw_values <- list(
#'   silhouette = c(0.7, 0.5, 0.8),
#'   ch = c(120.5, 95.2, 150.3),
#'   db = c(0.8, 1.2, 0.6),
#'   dunn = c(0.25, 0.18, 0.32)
#' )
#' norm_params <- perform_global_normalization(raw_values)
#' 
#' # Normalize individual values
#' norm_sil <- calculate_normalized_value(0.7, "silhouette", norm_params)
#' norm_ch <- calculate_normalized_value(120.5, "ch", norm_params)
#' norm_db <- calculate_normalized_value(0.8, "db", norm_params)
#' norm_dunn <- calculate_normalized_value(0.25, "dunn", norm_params)
#' 
#' # View results
#' cat("Normalized values:\n")
#' cat("Silhouette:", norm_sil, "\n")
#' cat("CH:", norm_ch, "\n") 
#' cat("DB:", norm_db, "\n")
#' cat("Dunn:", norm_dunn, "\n")
#' 
#' # Handle NA values
#' norm_na <- calculate_normalized_value(NA, "silhouette", norm_params)
#' cat("NA value result:", norm_na, "\n")
#'
#' @seealso \code{\link{perform_global_normalization}}
#'
#' @keywords internal
#' @export
calculate_normalized_value <- function(raw_value, metric_name, normalization_params) {
  
  if (is.na(raw_value)) {
    return(0)
  }
  
  params <- normalization_params[[metric_name]]
  if (is.null(params)) {
    return(0)
  }
  
  switch(params$type,
         "minmax_shifted" = {
           (raw_value + 1) / 2
         },
         "sum_based" = {
           if (params$sum_val > 0) {
             raw_value / params$sum_val
           } else {
             0
           }
         },
         "inverse_sum" = {
           if (raw_value > 0 && params$sum_val > 0) {
             (1 / raw_value) / params$sum_val
           } else {
             0
           }
         },
         0
  )
}
