% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{perform_global_normalization}
\alias{perform_global_normalization}
\title{Perform Global Normalization of Clustering Validation Indices}
\usage{
perform_global_normalization(all_raw_values, verbose = FALSE)
}
\arguments{
\item{all_raw_values}{A list containing vectors of raw validation index values
for all clustering methods. Expected components are:
\describe{
\item{silhouette}{Vector of silhouette index values}
\item{ch}{Vector of Calinski-Harabasz index values}
\item{db}{Vector of Davies-Bouldin index values}
\item{dunn}{Vector of Dunn index values}
}}

\item{verbose}{A logical value indicating whether to print normalization
details. Default is \code{FALSE}.}
}
\value{
A list containing normalization parameters for each validation index.
Each component contains:
\describe{
\item{type}{The normalization method used}
\item{min_val, max_val}{For min-max normalization (Silhouette)}
\item{sum_val}{For sum-based or inverse-sum normalization (CH, DB, Dunn)}
}
}
\description{
This function determines appropriate normalization parameters for different
clustering validation indices based on their mathematical properties and
the distribution of raw values across all clustering methods.
}
\details{
Different normalization strategies are applied based on index properties:
\itemize{
\item \strong{Silhouette}: Min-max normalization with shift (range: -1 to 1)
\item \strong{Calinski-Harabasz}: Sum-based normalization (higher is better)
\item \strong{Davies-Bouldin}: Inverse sum normalization (lower is better)
\item \strong{Dunn}: Sum-based normalization (higher is better)
}
}
\keyword{internal}
