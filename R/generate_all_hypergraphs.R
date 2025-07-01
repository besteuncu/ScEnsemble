generate_all_hypergraphs <- function(clustering_results, weight_list, algorithms = c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID")) {

  message("Hypergraph matrisleri olusturuluyor...")

  # 1. Her algoritma i??in clustering matrisini ??ret
  create_hypergraph_matrix <- function(clustering_vector) {
    n_samples <- length(clustering_vector)
    n_clusters <- length(unique(clustering_vector))
    H <- matrix(0, nrow = n_samples, ncol = n_clusters)
    for (i in 1:n_samples) {
      H[i, clustering_vector[i]] <- 1
    }
    return(H)
  }

  # 2. T??m H matrislerini ??ret
  H_list <- lapply(algorithms, function(algo) {
    create_hypergraph_matrix(clustering_results[[algo]])
  })
  names(H_list) <- algorithms

  # 3. Standart H ve HH hesapla
  H_combined <- do.call(cbind, H_list)
  HH_combined <- H_combined %*% t(H_combined)
  HH_combined <- HH_combined / length(algorithms)

  # 4. A????rl??kl?? metrikler i??in i??lem yapan yard??mc?? fonksiyon
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

  # 5. Her metrik i??in a????rl??kl?? H ve HH hesapla
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

Hypergraphs <- generate_all_hypergraphs(
  clustering_results = results_ind_pollen3$clustering_results,
  weight_list = results_all_val_indices$normalized_indices,
  algorithms = c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID")
)

H_list <- list(
  standard = Hypergraphs$H,
  silhouette = Hypergraphs$H_sil,
  ch = Hypergraphs$H_ch,
  db = Hypergraphs$H_db,
  dunn = Hypergraphs$H_dunn,
  average = Hypergraphs$H_ave
)

HH_list <- list(
  standard = Hypergraphs$HH,
  silhouette = Hypergraphs$HH_sil,
  ch = Hypergraphs$HH_ch,
  db = Hypergraphs$HH_db,
  dunn = Hypergraphs$HH_dunn,
  average = Hypergraphs$HH_ave
)
