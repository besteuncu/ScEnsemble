% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_clustering_performance.R
\name{analyze_clustering_performance}
\alias{analyze_clustering_performance}
\title{Enhanced Clustering Algorithm Performance Analysis and Visualization}
\usage{
analyze_clustering_performance(
  ensemble_results = NULL,
  individual_results = NULL,
  ensemble_algorithms = c("cspa", "cspa_louvain", "cspa_leiden", "mcla", "mcla_louvain",
    "mcla_leiden", "HGSC"),
  individual_algorithms = c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle",
    "RaceID"),
  variants = c("standard", "silhouette", "ch", "db", "average"),
  output_dir = "clustering_analysis_output",
  save_plots = TRUE,
  plot_format = "png",
  plot_width = 12,
  plot_height = 8,
  seed = 42
)
}
\arguments{
\item{ensemble_results}{A list containing ensemble clustering results with
quality_indices component}

\item{individual_results}{A list containing individual algorithm results with
validation_indices component}

\item{ensemble_algorithms}{Character vector of ensemble algorithm names
(default: c("cspa", "cspa_louvain", "cspa_leiden"))}

\item{individual_algorithms}{Character vector of individual algorithm names
(default: c("SC3", "CIDR", "Seurat", "SIMLR", "TSNE_Kmeans", "Monocle", "RaceID"))}

\item{variants}{Character vector of algorithm variants
(default: c("standard", "silhouette", "ch", "db", "average"))}

\item{output_dir}{Character string specifying output directory for plots
(default: "clustering_analysis_output")}

\item{save_plots}{Logical indicating whether to save plots to files
(default: TRUE)}

\item{plot_format}{Character string specifying plot format ("png", "pdf", "svg")
(default: "png")}

\item{plot_width}{Numeric value for plot width in inches (default: 12)}

\item{plot_height}{Numeric value for plot height in inches (default: 8)}

\item{seed}{Integer for random seed (default: 42)}
}
\value{
A list containing:
\describe{
\item{plots}{List of all generated ggplot objects}
\item{data}{Combined clustering performance data}
\item{normalized_data}{Normalized performance metrics}
\item{top_performers}{Top 5 performing algorithms}
\item{correlation_matrix}{Correlation matrix of evaluation metrics}
\item{summary_stats}{Summary statistics of performance metrics}
}
}
\description{
This function provides comprehensive analysis and visualization of clustering
algorithm performance using multiple evaluation metrics. It generates detailed
barplots, normalized performance comparisons, rankings, and correlation analyses.
}
\details{
The function performs the following analyses:
\itemize{
\item Individual metric barplot analysis (Silhouette, CH, DB, Dunn)
\item Normalized performance comparison across all metrics
\item Overall performance ranking with weighted scoring
\item Correlation analysis between evaluation metrics
\item Performance heatmap visualization
\item Sankey diagram for ranking flow analysis
}

Evaluation metrics interpretation:
\itemize{
\item Silhouette: Higher values indicate better clustering (range: -1 to 1)
\item CH (Calinski-Harabasz): Higher values indicate better clustering
\item DB (Davies-Bouldin): Lower values indicate better clustering
\item Dunn: Higher values indicate better clustering
}
}
