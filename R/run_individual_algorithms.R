#' Run Individual Clustering Algorithms
#'
#' This function runs multiple single-cell clustering algorithms on the provided data
#' and returns clustering results, performance metrics, and runtime information.
#'
#' @param data A numeric matrix or data frame containing single-cell expression data.
#'   Rows represent genes and columns represent cells.
#' @param true_labels Optional vector of true cluster labels for validation.
#'   If provided, Adjusted Rand Index (ARI) will be calculated.
#' @param algorithms Character vector specifying which algorithms to run.
#'   Default is c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID").
#' @param seed Integer value for random seed to ensure reproducibility. Default is 42.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @param n_cores Integer specifying number of cores to use for parallel processing where applicable.
#'   Default is 1.
#'
#' @return A list containing:
#' \describe{
#'   \item{clustering_results}{Named list of clustering labels for each algorithm}
#'   \item{ari_scores}{Named list of ARI scores (if true_labels provided)}
#'   \item{runtimes}{Named list of algorithm runtimes in seconds}
#'   \item{errors}{Named list of error messages for failed algorithms}
#'   \item{embedding_data}{Named list of embedding/distance matrices for each algorithm}
#' }
#'
#' @details
#' The function supports the following clustering algorithms:
#' \itemize{
#'   \item SC3: Single-Cell Consensus Clustering
#'   \item CIDR: Clustering through Imputation and Dimensionality Reduction
#'   \item Seurat: Graph-based clustering with Louvain algorithm
#'   \item SIMLR: Single-cell Interpretation via Multi-kernel Learning
#'   \item TSNE_Kmeans: t-SNE followed by K-means clustering
#'   \item Monocle: Trajectory-based clustering
#'   \item RaceID: Rare Cell Identification
#' }
#'
#' @examples
#' # Load example data
#' data(example_scdata)
#'
#' # Run all algorithms
#' results <- run_individual_algorithms(example_scdata)
#'
#' # Run specific algorithms with true labels
#' results <- run_individual_algorithms(
#'   data = example_scdata,
#'   true_labels = example_labels,
#'   algorithms = c("SC3", "Seurat", "SIMLR")
#' )
#'
#' @importFrom SingleCellExperiment SingleCellExperiment rowData colData
#' @importFrom SC3 sc3_estimate_k sc3
#' @importFrom cidr scDataConstructor determineDropoutCandidates wThreshold scDissim scPCA nPC scCluster
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters Embeddings Idents
#' @importFrom SIMLR SIMLR_Estimate_Number_of_Clusters SIMLR
#' @importFrom Rtsne Rtsne
#' @importFrom cluster silhouette clusGap
#' @importFrom monocle3 new_cell_data_set preprocess_cds reduce_dimension cluster_cells clusters
#' @importFrom RaceID SCseq filterdata getfdata compdist clustexp
#' @importFrom mclust adjustedRandIndex
#' @importFrom stats prcomp kmeans
#'
#' @export
setMethod("run_individual_algorithms", "ScEnsemble", 
          function(object, 
                   algorithms = c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID"),
                   seed = 42,
                   verbose = TRUE,
                   n_cores = 1,
                   ...) {
  
  data <- assay(object@sce, "counts")
  true_labels <- object@annotation
  
  validate_input <- function(object) {
    if (!is(object, "ScEnsemble")) {
      stop("object must be a ScEnsemble instance")
    }
    
    if (!"counts" %in% assayNames(object@sce)) {
      stop("No 'counts' assay found in SCE object")
    }
    
    data <- assay(object@sce, "counts")
    if (nrow(data) == 0 || ncol(data) == 0) {
      stop("Empty data matrix")
    }
    
    return(TRUE)
  }

  available_algorithms <- c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID")
  if (!all(algorithms %in% available_algorithms)) {
    stop("Invalid algorithm(s) specified. Available algorithms: ", paste(available_algorithms, collapse = ", "))
  }


  # Set seed for reproducibility
  set.seed(seed)

  # Initialize result lists
  clustering_results <- list()
  ari_scores <- list()
  runtimes <- list()
  errors <- list()
  embedding_data <- list()

  # Define algorithm functions
  algorithm_functions <- list(
    SC3 = function(data) {
      set.seed(seed)
      sce <- SingleCellExperiment(
        assays = list(counts = as.matrix(data),
                      logcounts = log2(as.matrix(data) + 1))
      )
      rowData(sce)$feature_symbol <- rownames(data)
      sce <- sc3_estimate_k(sce)
      estimated_k <- sce@metadata[["sc3"]][["k_estimation"]]
      sce <- sc3(sce, ks = estimated_k, biology = TRUE, n_cores = n_cores)
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
      if (is(data, "dgCMatrix") || is(data, "Matrix")) {
        data <- as.matrix(data)
      }
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
      pca_data <- pca_result$x[, 1:min(50, ncol(pca_result$x))]
      set.seed(seed)
      tsne_result <- Rtsne(pca_data, dims = 2, perplexity = 30,
                           check_duplicates = FALSE)
      tsne_data <- tsne_result$Y
      gap_stat <- clusGap(tsne_data, FUN = kmeans, nstart = 25, K.max = min(15, nrow(tsne_data)-1), B = 50)
      optimal_k <- with(gap_stat, maxSE(Tab[, "gap"], Tab[, "SE.sim"]))
      set.seed(seed)
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
      if(is.data.frame(data)) {
        data <- as.matrix(data)
      }
      sc <- SCseq(data)
      sc <- filterdata(sc, mintotal = 1000, minexpr = 5, minnumber = 1)
      fdata <- getfdata(sc)
      sc <- compdist(sc, metric = "pearson")
      sc <- clustexp(sc, clustnr = 30, bootnr = 50, rseed = seed)
      raceID_data <- sc@distances
      raceid_labels <- as.integer(sc@cluster$kpart)
      return(list(labels = raceid_labels, data = raceID_data))
    }
  )

  # Run selected algorithms
  for (alg_name in algorithms) {
    if (verbose) {
      message("------------------------------------------------")
      message(paste("Running algorithm:", alg_name))
    }

    func <- algorithm_functions[[alg_name]]
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

      if (verbose) {
        message(paste(alg_name, "completed successfully. Runtime:", runtime, "seconds"))
        if (!is.null(true_labels)) {
          message(paste("ARI:", round(ari_val, 4)))
        }
      }

    }, error = function(e) {
      if (verbose) {
        message(paste("ERROR in", alg_name, ":", e$message))
      }

      clustering_results[[alg_name]] <- NA
      embedding_data[[alg_name]] <- NA
      runtimes[[alg_name]] <- NA
      ari_scores[[alg_name]] <- NA
      errors[[alg_name]] <- e$message
    })
  }

  # Return results
  result_list <- new("IndividualResults",
    clustering_results = clustering_results,
    ari_scores = ari_scores,
    runtimes = runtimes,
    errors = errors,
    embedding_data = embedding_data
  )

  object@individual_results <- result_list
  return(object)
}
)
