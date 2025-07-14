#' Visualize Dimension Reduction with Clustering
#'
#' Generates 2D or 3D visualizations of clustering results using t-SNE or UMAP after PCA.
#'
#' @param expr_matrix A numeric matrix of gene expression (cells x genes).
#' @param cluster_labels A vector of clustering labels from individual or ensemble clustering.
#' @param method Character. Dimension reduction method: either `"tsne"` or `"umap"`.
#' @param dims Integer. Number of dimensions for visualization (2 or 3). Default is 2.
#' @param title Character. Plot title. Default is auto-generated.
#' @return A `ggplot` (2D) or `plotly` (3D) object.
#'
#' @import ggplot2
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
#' @importFrom plotly plot_ly layout
#'
#' @export
visualize_tsne_umap <- function(expr_matrix, cluster_labels,
                               method = c("tsne", "umap"),
                               dims = 2,
                               title = NULL) {
  method <- match.arg(method)
  stopifnot(dims %in% c(2, 3))

  set.seed(123)
  pca_result <- prcomp(t(expr_matrix), center = TRUE)
  pca_data <- pca_result$x[, 1:50]

  if (method == "tsne") {
    tsne_result <- Rtsne(pca_data, dims = dims, perplexity = 30, verbose = FALSE)
    embedding <- tsne_result$Y
  } else {
    embedding <- umap(pca_data, n_components = dims, n_neighbors = 15, min_dist = 0.1)
  }

  cluster_labels <- as.factor(cluster_labels)
  colnames(embedding) <- paste0("Dim", seq_len(ncol(embedding)))
  df <- as.data.frame(embedding)
  df$Cluster <- cluster_labels

  if (is.null(title)) {
    title <- paste(toupper(method), paste0("(", dims, "D)"), "Clustering Visualization")
  }

  if (dims == 2) {
    ggplot(df, aes_string(x = "Dim1", y = "Dim2", color = "Cluster")) +
      geom_point(alpha = 0.7) +
      labs(title = title) +
      theme_minimal()
  } else {
    plot_ly(df, x = ~Dim1, y = ~Dim2, z = ~Dim3, color = ~Cluster, colors = "Set1",
            type = "scatter3d", mode = "markers", marker = list(size = 3)) %>%
      layout(title = title)
  }
}
