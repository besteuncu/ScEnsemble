run_individual_algorithms <- function(data, true_labels = NULL) {
  library(SingleCellExperiment)
  library(SC3)
  library(cidr)
  library(Seurat)
  library(SIMLR)
  library(Rtsne)
  library(cluster)
  library(monocle3)
  library(RaceID)
  library(SHARP)
  library(mclust)

  clustering_results <- list()
  ari_scores <- list()
  runtimes <- list()
  errors <- list()
  embedding_data <- list()

  algorithms <- list(
    SC3 = function(data) {
      set.seed(42)
      sce <- SingleCellExperiment(
        assays = list(counts = as.matrix(data),
                      logcounts = log2(as.matrix(data) + 1))
      )
      rowData(sce)$feature_symbol <- rownames(data)
      sce <- sc3_estimate_k(sce)
      estimated_k <- sce@metadata[["sc3"]][["k_estimation"]]
      sce <- sc3(sce, ks = estimated_k, biology = TRUE, n_cores = 1)
      sc3_data <- sce@metadata$sc3$consensus[[as.character(estimated_k)]]$consensus
      sc3_labels <- as.integer(sce@colData@listData[[paste0("sc3_", estimated_k, "_clusters")]])
      return(list(labels = sc3_labels, data = sc3_data))
    },

    CIDR = function(data) {
      scData <- scDataConstructor(as.matrix(data))
      scData <- determineDropoutCandidates(scData)
      scData <- wThreshold(scData)
      scData <- scDissim(scData)
      scData <- scPCA(scData)
      scData <- nPC(scData)
      cidr_data <- scData@dissim
      max_k <- min(15, ncol(data) - 1)
      sil_scores <- numeric(max_k - 1)
      for (k_val in 2:max_k) {
        scData_temp <- scCluster(scData, nCluster = k_val)
        cluster_labels <- scData_temp@clusters
        sil <- silhouette(cluster_labels, cidr_data)
        sil_scores[k_val - 1] <- mean(sil[, "sil_width"])
      }
      estimated_k <- which.max(sil_scores) + 1
      scData <- scCluster(scData, nCluster = estimated_k)
      cidr_labels <- scData@clusters
      return(list(labels = cidr_labels, data = cidr_data))
    },

    Seurat = function(data) {
      seu <- CreateSeuratObject(counts = as.matrix(data))
      seu <- NormalizeData(seu, verbose = FALSE)
      seu <- FindVariableFeatures(seu, verbose = FALSE)
      seu <- ScaleData(seu, verbose = FALSE)
      seu <- RunPCA(seu, verbose = FALSE)
      seu <- FindNeighbors(seu, dims = 1:10, verbose = FALSE)
      seu <- FindClusters(seu, resolution = 0.7, verbose = FALSE)
      seurat_data <- Embeddings(seu, reduction = "pca")[, 1:10]
      seurat_labels <- as.numeric(Idents(seu))
      return(list(labels = seurat_labels, data = seurat_data))
    },

    SIMLR = function(data) {
      k_estimate <- SIMLR_Estimate_Number_of_Clusters(data, NUMC = 2:10, cores.ratio = 0.5)
      best_k1 <- which.min(k_estimate$K1)
      best_k2 <- which.min(k_estimate$K2)
      best_k <- if (best_k1 <= best_k2) best_k1 else best_k2
      res_simlr <- SIMLR(X = data, c = best_k, normalize = TRUE)
      simlr_data <- res_simlr[["S"]]
      simlr_labels <- res_simlr$y$cluster
      return(list(labels = simlr_labels, data = simlr_data))
    },

    TSNE_Kmeans = function(data) {
      pca_result <- prcomp(t(data))
      pca_data <- pca_result$x[, 1:50]
      set.seed(123)
      tsne_result <- Rtsne(pca_data, dims = 2, perplexity = 30, check_duplicates = FALSE)
      tsne_data <- tsne_result$Y
      gap_stat <- clusGap(tsne_data, FUN = kmeans, nstart = 25, K.max = 15, B = 50)
      optimal_k <- with(gap_stat, maxSE(Tab[, "gap"], Tab[, "SE.sim"]))
      kmeans_result <- kmeans(tsne_data, centers = optimal_k, nstart = 25)
      tsne_kmeans_labels <- kmeans_result$cluster
      return(list(labels = tsne_kmeans_labels, data = tsne_data))
    },

    Monocle = function(data) {
      gene_metadata <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
      cell_metadata <- data.frame(row.names = colnames(data))
      cds <- new_cell_data_set(as.matrix(data), cell_metadata = cell_metadata, gene_metadata = gene_metadata)
      cds <- preprocess_cds(cds, num_dim = 50)
      cds <- reduce_dimension(cds, preprocess_method = "PCA")
      cds <- cluster_cells(cds)
      monocle_clusters <- clusters(cds)
      monocle_data <- cds@int_colData@listData[["reducedDims"]][["PCA"]]
      monocle_labels <- as.integer(monocle_clusters)
      return(list(labels = monocle_labels, data = monocle_data))
    },

    RaceID = function(data) {
      # Veri kontrol?? ve haz??rl??????
      if(is.data.frame(data)) {
        data <- as.matrix(data)
      }
      sc <- SCseq(data)
      sc <- filterdata(sc, mintotal = 1000, minexpr = 5, minnumber = 1)
      fdata <- getfdata(sc)
      sc <- compdist(sc, metric = "pearson")
      sc <- clustexp(sc, clustnr = 30, bootnr = 50, rseed = 12345)
      raceID_data <- sc@distances
      raceid_labels <- as.integer(sc@cluster$kpart)
      return(list(labels = raceid_labels, data = raceID_data))
    }

  )

  for (alg_name in names(algorithms)) {
    message("------------------------------------------------")
    message(paste("Algorithm Running:", alg_name))
    func <- algorithms[[alg_name]]
    result <- tryCatch({
      start_time <- Sys.time()
      res <- func(data)
      end_time <- Sys.time()
      runtime <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
      ari_val <- if (!is.null(true_labels)) adjustedRandIndex(true_labels, res$labels) else NA

      clustering_results[[alg_name]] <- res$labels
      embedding_data[[alg_name]] <- res$data
      runtimes[[alg_name]] <- runtime
      ari_scores[[alg_name]] <- ari_val
      errors[[alg_name]] <- NA

      message(paste(alg_name, "completed successfully. Time:", runtime, "s"))
      if (!is.null(true_labels)) message(paste("ARI:", round(ari_val, 4)))
    }, error = function(e) {
      message(paste("ERROR:", alg_name, "running:", e$message))

      clustering_results[[alg_name]] <- list(
        labels = NA, data = NA, runtime = NA, ARI = NA, error = e$message
      )

      embedding_data[[alg_name]] <- NA
      runtimes[[alg_name]] <- NA
      ari_scores[[alg_name]] <- NA
      errors[[alg_name]] <- e$message
    })
  }

  return(list(
    clustering_results = clustering_results,
    ari_scores = ari_scores,
    runtimes = runtimes,
    errors = errors,
    embedding_data = embedding_data
  ))
}


