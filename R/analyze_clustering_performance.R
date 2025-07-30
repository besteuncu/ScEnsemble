#' Enhanced Clustering Algorithm Performance Analysis and Visualization
#'
#' This function provides comprehensive analysis and visualization of clustering
#' algorithm performance using multiple evaluation metrics. It generates detailed
#' barplots, normalized performance comparisons, rankings, and correlation analyses.
#'
#' @param ensemble_results A list containing ensemble clustering results with
#'   quality_indices component
#' @param individual_results A list containing individual algorithm results with
#'   validation_indices component
#' @param ensemble_algorithms Character vector of ensemble algorithm names
#'   (default: c("cspa", "cspa_louvain", "cspa_leiden"))
#' @param individual_algorithms Character vector of individual algorithm names
#'   (default: c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID"))
#' @param variants Character vector of algorithm variants
#'   (default: c("standard", "silhouette", "ch", "db", "average"))
#' @param output_dir Character string specifying output directory for plots
#'   (default: "clustering_analysis_output")
#' @param save_plots Logical indicating whether to save plots to files
#'   (default: TRUE)
#' @param plot_format Character string specifying plot format ("png", "pdf", "svg")
#'   (default: "png")
#' @param plot_width Numeric value for plot width in inches (default: 12)
#' @param plot_height Numeric value for plot height in inches (default: 8)
#' @param seed Integer for random seed (default: 42)
#'
#' @return A list containing:
#' \describe{
#'   \item{plots}{List of all generated ggplot objects}
#'   \item{data}{Combined clustering performance data}
#'   \item{normalized_data}{Normalized performance metrics}
#'   \item{top_performers}{Top 5 performing algorithms}
#'   \item{correlation_matrix}{Correlation matrix of evaluation metrics}
#'   \item{summary_stats}{Summary statistics of performance metrics}
#' }
#'
#' @details
#' The function performs the following analyses:
#' \itemize{
#'   \item Individual metric barplot analysis (Silhouette, CH, DB, Dunn)
#'   \item Normalized performance comparison across all metrics
#'   \item Overall performance ranking with weighted scoring
#'   \item Correlation analysis between evaluation metrics
#'   \item Performance heatmap visualization
#'   \item Sankey diagram for ranking flow analysis
#' }
#'
#' Evaluation metrics interpretation:
#' \itemize{
#'   \item Silhouette: Higher values indicate better clustering (range: -1 to 1)
#'   \item CH (Calinski-Harabasz): Higher values indicate better clustering
#'   \item DB (Davies-Bouldin): Lower values indicate better clustering
#'   \item Dunn: Higher values indicate better clustering
#' }
#'

#' @import ggplot2
#' @importFrom dplyr select filter mutate group_by summarise arrange bind_rows across desc
#' @import tidyr
#' @importFrom gridExtra grid.arrange
#' @import RColorBrewer
#' @import viridis
#' @import corrplot
#' @import ggalluvial
#' @importFrom grDevices png dev.off pdf
#' @importFrom stats cor complete.cases reorder
#' @importFrom utils head
#'
#' @export
analyze_clustering_performance <- function(ensemble_results = NULL,
                                           individual_results = NULL,
                                           ensemble_algorithms = c("cspa", "cspa_louvain", "cspa_leiden",
                                                                   "mcla", "mcla_louvain", "mcla_leiden","HGSC"),
                                           individual_algorithms = c("SC3", "CIDR", "Seurat", "SIMLR",
                                                                     "TSNE_Kmeans", "Monocle", "RaceID"),
                                           variants = c("standard", "silhouette", "ch", "db", "average"),
                                           output_dir = "clustering_analysis_output",
                                           save_plots = TRUE,
                                           plot_format = "png",
                                           plot_width = 12,
                                           plot_height = 8,
                                           seed = 42) {

  # Input validation
  if (is.null(ensemble_results) && is.null(individual_results)) {
    stop("At least one of ensemble_results or individual_results must be provided")
  }

  # Set random seed for reproducibility
  set.seed(seed)

  # Create output directory if saving plots
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Initialize storage for plots and results
  all_plots <- list()
  plot_counter <- 1

  # Helper function to add plots to collection
  add_plot <- function(plot_obj, title = "") {
    all_plots[[plot_counter]] <<- list(plot = plot_obj, title = title)
    plot_counter <<- plot_counter + 1
    return(plot_obj)
  }

  # Helper function to extract ensemble data
  extract_ensemble_data <- function(internal_indices, algorithms, variants) {
    ensemble_data <- data.frame()

    for (alg in algorithms) {
      if (alg %in% names(internal_indices)) {
        for (var in variants) {
          if (var %in% names(internal_indices[[alg]])) {
            indices <- internal_indices[[alg]][[var]]

            # Extract metrics with error handling
            silhouette <- ifelse("silhouette" %in% names(indices), indices$silhouette, NA)
            ch <- ifelse("ch" %in% names(indices), indices$ch, NA)
            db <- ifelse("db" %in% names(indices), indices$db, NA)
            dunn <- ifelse("dunn" %in% names(indices), indices$dunn, NA)

            ensemble_data <- rbind(ensemble_data, data.frame(
              Algorithm = alg,
              Variant = var,
              Silhouette = silhouette,
              CH = ch,
              DB = db,
              Dunn = dunn,
              Type = "Ensemble",
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }

    return(ensemble_data)
  }

  # Helper function to extract individual algorithm data
  extract_individual_data <- function(validation_indices) {
    individual_data <- data.frame()

    for (alg in names(validation_indices)) {
      indices <- validation_indices[[alg]]

      # Extract metrics with error handling
      silhouette <- ifelse("silhouette" %in% names(indices), indices$silhouette, NA)
      ch <- ifelse("ch" %in% names(indices), indices$ch, NA)
      db <- ifelse("db" %in% names(indices), indices$db, NA)
      dunn <- ifelse("dunn" %in% names(indices), indices$dunn, NA)

      individual_data <- rbind(individual_data, data.frame(
        Algorithm = alg,
        Variant = "base",
        Silhouette = silhouette,
        CH = ch,
        DB = db,
        Dunn = dunn,
        Type = "Base",
        stringsAsFactors = FALSE
      ))
    }

    return(individual_data)
  }

  # Extract data from results
  clustering_data <- data.frame()

  if (!is.null(ensemble_results) && "quality_indices" %in% names(ensemble_results)) {
    ensemble_data <- extract_ensemble_data(
      internal_indices = ensemble_results$quality_indices,
      algorithms = ensemble_algorithms,
      variants = variants
    )
    clustering_data <- rbind(clustering_data, ensemble_data)
  }

  if (!is.null(individual_results) && "validation_indices" %in% names(individual_results)) {
    individual_data <- extract_individual_data(individual_results$validation_indices)
    clustering_data <- rbind(clustering_data, individual_data)
  }

  # Check if data was extracted successfully
  if (nrow(clustering_data) == 0) {
    stop("No valid clustering data could be extracted from the provided results")
  }

  # Create full algorithm names
  clustering_data$Full_Name <- ifelse(clustering_data$Type == "Base",
                                      clustering_data$Algorithm,
                                      paste(clustering_data$Algorithm, clustering_data$Variant, sep = "_"))

  # Enhanced color palettes
  type_colors <- c("Ensemble" = "#E31A1C", "Base" = "#1F78B4")

  message("=== COMPREHENSIVE CLUSTERING ALGORITHM ANALYSIS ===")
  message(paste("Total algorithms analyzed:", nrow(clustering_data)))
  message("Analysis includes performance metrics, barplots, and comprehensive reporting")

  # ============================================================================
  # 1. DETAILED BARPLOT ANALYSIS FOR EACH METRIC
  # ============================================================================

  message("=== STEP 1: CREATING DETAILED BARPLOT ANALYSIS FOR EACH METRIC ===")

  # Transform data to long format for individual metric analysis
  data_long <- clustering_data %>%
    select(-Full_Name) %>%
    pivot_longer(cols = c(Silhouette, CH, DB, Dunn),
                 names_to = "Metric",
                 values_to = "Value") %>%
    filter(!is.na(Value))

  # Create individual barplots for each metric
  metrics_list <- unique(data_long$Metric)
  barplot_list <- list()

  for(i in seq_along(metrics_list)) {
    current_metric <- metrics_list[i]
    current_data <- data_long[data_long$Metric == current_metric, ]

    # Metric evaluation criteria descriptions
    metric_info <- switch(current_metric,
                          "Silhouette" = "Higher value = Better clustering",
                          "CH" = "Higher value = Better clustering",
                          "DB" = "Lower value = Better clustering",
                          "Dunn" = "Higher value = Better clustering"
    )

    # Create enhanced barplot for current metric
    p <- ggplot(current_data, aes(x = reorder(paste(Algorithm, Variant, sep="_"), Value),
                                  y = Value, fill = Type)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      coord_flip() +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 9),
            plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 10),
            legend.position = "bottom") +
      labs(title = paste(current_metric, "Score Comparison"),
           subtitle = metric_info,
           x = "Algorithm_Variant",
           y = current_metric,
           fill = "Algorithm Type") +
      scale_fill_manual(values = type_colors) +
      geom_text(aes(label = round(Value, 3)),
                hjust = ifelse(current_data$Value > 0, -0.1, 1.1),
                size = 3)

    barplot_list[[i]] <- p
    add_plot(p, paste("Barplot Analysis:", current_metric, "Metric"))

    # Save individual plots if requested
    if (save_plots) {
      filename <- file.path(output_dir, paste0("barplot_", current_metric, ".", plot_format))
      ggsave(filename, plot = p, width = plot_width, height = plot_height)
    }
  }

  # ============================================================================
  # 2. NORMALIZATION AND OVERALL PERFORMANCE CALCULATION
  # ============================================================================

  message("=== STEP 2: CALCULATING NORMALIZED PERFORMANCE METRICS ===")

  # Advanced normalization function
  normalize_metric <- function(x, higher_is_better = TRUE, method = "minmax") {
    if(method == "minmax") {
      if(higher_is_better) {
        return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
      } else {
        return((max(x, na.rm = TRUE) - x) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
      }
    }
  }

  # Calculate normalized metrics with weighted importance
  normalized_data <- clustering_data %>%
    mutate(
      Silhouette_norm = normalize_metric(Silhouette, TRUE),
      CH_norm = normalize_metric(CH, TRUE),
      DB_norm = normalize_metric(DB, FALSE),  # Lower DB values are better
      Dunn_norm = normalize_metric(Dunn, TRUE)
    ) %>%
    mutate(
      # Weighted overall score (equal weights for all metrics)
      Overall_Score = (Silhouette_norm * 0.25 + CH_norm * 0.25 + DB_norm * 0.25 + Dunn_norm * 0.25)
    )

  message("Normalization complete. All metrics now scaled 0-1 for fair comparison.")

  # ============================================================================
  # 3. NORMALIZED PERFORMANCE BARPLOT
  # ============================================================================

  message("=== STEP 3: CREATING NORMALIZED PERFORMANCE COMPARISON ===")

  # Transform normalized data to long format
  normalized_long <- normalized_data %>%
    select(Algorithm, Variant, Type, ends_with("_norm")) %>%
    pivot_longer(cols = ends_with("_norm"),
                 names_to = "Metric",
                 values_to = "Normalized_Value") %>%
    mutate(Metric = gsub("_norm", "", Metric)) %>%
    filter(!is.na(Normalized_Value))

  # Create normalized performance barplot
  p_normalized <- ggplot(normalized_long,
                         aes(x = paste(Algorithm, Variant, sep="_"),
                             y = Normalized_Value,
                             fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "bottom",
          plot.title = element_text(size = 12, face = "bold")) +
    labs(title = "Normalized Performance Comparison",
         subtitle = "All metrics normalized to 0-1 scale (1 = best performance)",
         x = "Algorithm_Variant",
         y = "Normalized Value",
         fill = "Metric") +
    scale_fill_brewer(palette = "Set3") +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5)

  add_plot(p_normalized, "Normalized Performance Comparison")

  if (save_plots) {
    filename <- file.path(output_dir, paste0("normalized_performance.", plot_format))
    ggsave(filename, plot = p_normalized, width = plot_width, height = plot_height)
  }

  # ============================================================================
  # 4. OVERALL PERFORMANCE RANKING
  # ============================================================================

  message("=== STEP 4: CREATING OVERALL PERFORMANCE RANKING ===")

  # Create overall performance ranking barplot
  p_overall <- ggplot(normalized_data,
                      aes(x = reorder(paste(Algorithm, Variant, sep="_"), Overall_Score),
                          y = Overall_Score,
                          fill = Type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(size = 12, face = "bold")) +
    labs(title = "Overall Performance Score Ranking",
         subtitle = "Ranking by weighted average of all metrics",
         x = "Algorithm_Variant",
         y = "Overall Performance Score (0-1)",
         fill = "Algorithm Type") +
    scale_fill_manual(values = type_colors) +
    geom_text(aes(label = round(Overall_Score, 3)),
              hjust = -0.1, size = 3) +
    geom_hline(yintercept = mean(normalized_data$Overall_Score, na.rm = TRUE),
               linetype = "dashed", color = "red", alpha = 0.7)

  add_plot(p_overall, "Overall Performance Ranking")

  if (save_plots) {
    filename <- file.path(output_dir, paste0("overall_ranking.", plot_format))
    ggsave(filename, plot = p_overall, width = plot_width, height = plot_height)
  }

  # ============================================================================
  # 5. IDENTIFY TOP PERFORMERS
  # ============================================================================

  message("=== STEP 5: IDENTIFYING TOP PERFORMERS ===")

  # Get top 5 performers overall
  top_performers <- normalized_data %>%
    arrange(desc(Overall_Score)) %>%
    head(5)

  message("Top 5 Algorithms by Overall Performance:")
  for(i in 1:nrow(top_performers)) {
    alg <- top_performers[i, ]
    message(paste(i, ".", alg$Full_Name, "- Overall Score:", round(alg$Overall_Score, 3)))
  }

  # ============================================================================
  # 6. CORRELATION ANALYSIS
  # ============================================================================

  message("=== STEP 6: ANALYZING METRIC CORRELATIONS ===")

  # Understanding how different evaluation metrics relate to each other
  metrics_for_correlation <- normalized_data %>%
    select(Silhouette, CH, DB, Dunn) %>%
    filter(complete.cases(.))

  correlation_matrix <- cor(metrics_for_correlation, use = "complete.obs")

  # Save correlation plot if requested
  if (save_plots) {
    filename <- file.path(output_dir, paste0("correlation_plot.", plot_format))
    if (plot_format == "png") {
      png(filename, width = 600, height = 600)
    } else if (plot_format == "pdf") {
      pdf(filename, width = 8, height = 8)
    }

    corrplot(correlation_matrix, method = "color", type = "upper",
             addCoef.col = "black", tl.col = "black", tl.srt = 45,
             title = "Correlation Between Clustering Evaluation Metrics",
             mar = c(0,0,2,0))
    dev.off()
  }

  # ============================================================================
  # 7. PERFORMANCE HEATMAP
  # ============================================================================

  message("=== STEP 7: CREATING COMPREHENSIVE PERFORMANCE HEATMAP ===")

  # Multi-metric heatmap for visual performance comparison
  heatmap_data <- normalized_data %>%
    select(Full_Name, Type, Silhouette_norm, CH_norm, DB_norm, Dunn_norm, Overall_Score) %>%
    pivot_longer(cols = c(Silhouette_norm, CH_norm, DB_norm, Dunn_norm),
                 names_to = "Metric", values_to = "Normalized_Value") %>%
    mutate(Metric = gsub("_norm", "", Metric)) %>%
    filter(!is.na(Normalized_Value))

  p_heatmap <- ggplot(heatmap_data,
                      aes(x = Metric, y = reorder(Full_Name, Overall_Score), fill = Normalized_Value)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_viridis_c(name = "Performance\nScore", option = "plasma") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)) +
    labs(title = "Algorithm Performance Heatmap",
         subtitle = "Comprehensive view of normalized performance across all metrics",
         x = "Evaluation Metrics",
         y = "Algorithms (ordered by overall performance)")

  add_plot(p_heatmap, "Performance Heatmap")

  if (save_plots) {
    filename <- file.path(output_dir, paste0("performance_heatmap.", plot_format))
    ggsave(filename, plot = p_heatmap, width = plot_width, height = plot_height)
  }

  # ============================================================================
  # 8. SUMMARY STATISTICS
  # ============================================================================

  message("=== STEP 8: CALCULATING SUMMARY STATISTICS ===")

  summary_stats <- normalized_data %>%
    select(Type, Silhouette, CH, DB, Dunn, Overall_Score) %>%
    group_by(Type) %>%
    summarise(
      across(everything(), list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE)
      )),
      .groups = "drop"
    )

  # ============================================================================
  # 9. DISPLAY GRID OF BARPLOTS
  # ============================================================================

  message("=== DISPLAYING COMPREHENSIVE ANALYSIS RESULTS ===")

  # Display all metric barplots in a grid
  if (length(barplot_list) > 0) {
    grid_plot <- do.call(grid.arrange, c(barplot_list, ncol = 2))

    if (save_plots) {
      filename <- file.path(output_dir, paste0("barplot_grid.", plot_format))
      ggsave(filename, plot = grid_plot, width = plot_width * 1.5, height = plot_height * 1.5)
    }
  }

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  message("=== ANALYSIS COMPLETE ===")
  message(paste("Generated", length(all_plots), "plots"))
  if (save_plots) {
    message(paste("Plots saved to:", output_dir))
  }

  return(list(
    plots = all_plots,
    data = clustering_data,
    normalized_data = normalized_data,
    top_performers = top_performers,
    correlation_matrix = correlation_matrix,
    summary_stats = summary_stats,
    barplot_list = barplot_list,
    normalized_plot = p_normalized,
    overall_ranking_plot = p_overall,
    heatmap_plot = p_heatmap
  ))
}
