create_test_data <- function() {
  library(scRNAseq)
  sce <- PollenGliaData()
  ann <- colData(sce)[["Inferred Cell Type"]]
  colData(sce)$cell_type <- ann
  return(list(sce = sce, ann = ann))
}
