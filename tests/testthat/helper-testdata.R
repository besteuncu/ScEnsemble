library(SingleCellExperiment)
library(scRNAseq)
library(Matrix)

suppressMessages({
  Pollen <- PollenGliaData()
})

# Take a subset of 100 cells and 500 most variable genes
cell_subset <- colnames(Pollen)[1:100]
Pollen_sub <- Pollen[, cell_subset]


rowVars_base <- function(x) {
  x <- as.matrix(x)
  rowMeans((x - rowMeans(x))^2) * (ncol(x) / (ncol(x) - 1))
}
gene_vars <- rowVars_base(assay(Pollen_sub))
top_genes <- order(gene_vars, decreasing = TRUE)[1:500]
Pollen_sub <- Pollen_sub[top_genes, ]

# Assign globally for tests
test_pollen <- Pollen_sub
test_labels <- colData(Pollen_sub)[["Inferred Cell Type"]]
