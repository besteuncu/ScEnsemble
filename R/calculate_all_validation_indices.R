#' Calculate All Validation Indices
#'
#' This function calculates multiple clustering validation indices for each algorithm's
#' clustering results and returns both raw and normalized scores.
#'
#' @param clustering_results A named list containing clustering labels for each algorithm.
#'   Each element should be a vector of cluster assignments.
#' @param embedding_data A named list containing embedding or distance matrices for each algorithm.
#'   Each element should be a matrix or data frame with the same number of rows as clustering labels.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#'
#' @return A list containing:
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
#' @examples
#' # Calculate validation indices
#' validation_results <- calculate_all_validation_indices(
#'   clustering_results = clustering_results,
#'   embedding_data = embedding_data
#' )
#'
#' @importFrom cluster silhouette
#' @importFrom fpc calinhara
#' @importFrom clusterSim index.DB
#' @importFrom clValid dunn
#' @importFrom stats dist
#'
#' @export
setMethod("calculate_all_validation_indices", "ScEnsemble", 
          function(object,
                   clustering_results,
                   embedding_data,
                   verbose = TRUE,
                   ...) {

  clustering_results <- object@individual_results@clustering_results
  embedding_data     <- object@individual_results@embedding_data
  
  
  # Input validation
  if (!is.list(clustering_results) || length(clustering_results) == 0) {
    stop("clustering_results must be a non-empty list")
  }

  if (!is.list(embedding_data) || length(embedding_data) == 0) {
    stop("embedding_data must be a non-empty list")
  }
  
  # Internal function that calculates the indexes for each algorithm
  calculate_indices <- function(data, clustering) {
    if (length(unique(clustering)) < 2) {
      return(list(silhouette = NA, ch = NA, db = NA, dunn = NA))
    }
    data_matrix <- as.matrix(data)
    
    sil_index <- tryCatch({
      sil <- silhouette(clustering, dist(data_matrix))
      mean(sil[, 3], na.rm = TRUE)
    }, error = function(e) NA)
    
    ch_index <- tryCatch({
      calinhara(data_matrix, clustering, cn = length(unique(clustering)))
    }, error = function(e) NA)
    
    db_index <- tryCatch({
      result <- clusterSim::index.DB(data_matrix, clustering, centrotypes = "centroids")
      db_val <- result$DB
      return(db_val)
    }, error = function(e) NA)
    
    dunn_index <- tryCatch({
      dunn(dist(data_matrix), clustering)
    }, error = function(e) NA)
    
    return(list(
      silhouette = sil_index,
      ch = ch_index,
      db = db_index,
      dunn = dunn_index
    ))
  }

  # Function to normalize and calculate average scores
  normalize_all <- function(indices_list) {
    extract_metric <- function(metric) sapply(indices_list, function(x) x[[metric]])
    sil_values  <- extract_metric("silhouette")
    ch_values   <- extract_metric("ch")
    db_values   <- extract_metric("db")
    dunn_values <- extract_metric("dunn")

    replace_na <- function(x) ifelse(is.na(x), 0, x)

    sil_norm  <- replace_na((sil_values + 1) / 2)
    ch_norm   <- replace_na(ch_values / sum(ch_values, na.rm = TRUE))
    inv_db    <- replace_na(1 / db_values)
    db_norm   <- inv_db / sum(inv_db, na.rm = TRUE)
    dunn_norm <- replace_na(dunn_values / sum(dunn_values, na.rm = TRUE))

    # Average score (average of 4 normalized indices)
    avg_index <- rowMeans(cbind(sil_norm, ch_norm, db_norm, dunn_norm), na.rm = TRUE)

    return(list(
      silhouette = sil_norm,
      ch = ch_norm,
      db = db_norm,
      dunn = dunn_norm,
      average_index = avg_index
    ))
  }

  algorithm_names <- intersect(names(clustering_results), names(embedding_data))
  validation_indices <- list()

  for (alg in algorithm_names) {
    if (verbose) {
      message(paste("Calculating indices for:", alg))
    }
    labels <- clustering_results[[alg]]
    data   <- embedding_data[[alg]]
    validation_indices[[alg]] <- calculate_indices(data, labels)
  }

  normalized_indices <- normalize_all(validation_indices)

  result_list <- new("ValidationResults",
                     validation_indices = validation_indices,
                     normalized_indices = normalized_indices,
                     average_index = normalized_indices$average_index 
  )
    

  object@validation_metrics <- result_list
  return(object)
}
)

# Fonksiyonun içindeki tüm değişkenleri listeleyin
ls(environment(calculate_all_validation_indices))

# Fonksiyonun parent environment'ını kontrol edin
parent.env(environment(calculate_all_validation_indices))
