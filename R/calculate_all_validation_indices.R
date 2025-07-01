
calculate_all_validation_indices <- function(clustering_results, embedding_data) {
  library(cluster)       # silhouette
  library(fpc)           # calinhara
  library(clusterSim)    # index.DB
  library(clValid)       # dunn

  # Her algoritma i??in indeksleri hesaplayan i?? fonksiyon
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
      index.DB(data_matrix, clustering, centrotypes = "centroids")$DB
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

  # Normalize edip ortalama puanlar?? hesaplayan fonksiyon
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

    # Ortalama skor (normalize edilmi?? 4 indeksin ortalamas??)
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
    cat("Hesaplan??yor:", alg, "\n")
    labels <- clustering_results[[alg]]
    data   <- embedding_data[[alg]]
    validation_indices[[alg]] <- calculate_indices(data, labels)
  }

  normalized_indices <- normalize_all(validation_indices)

  return(list(
    validation_indices = validation_indices,
    normalized_indices = normalized_indices,
    average_index = normalized_indices$average_index
  ))
}


