
library(testthat)
library(ScEnsemble)

test_check("ScEnsemble")

library(SingleCellExperiment)
library(scRNAseq)

create_test_data <- function(n_genes = 50, n_cells = 30) {
  set.seed(123)
  
  counts <- matrix(0, nrow = n_genes, ncol = n_cells)
  
  cell_types <- rep(c("TypeA", "TypeB", "TypeC"), length.out = n_cells)
  
  for(i in 1:n_cells) {
    for(j in 1:n_genes) {
      if(runif(1) > 0.7) {  
        counts[j, i] <- rpois(1, lambda = 2)
      }
    }
  }
  
  rownames(counts) <- paste0("Gene", 1:n_genes)  
  colnames(counts) <- paste0("Cell", 1:n_cells)
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = DataFrame(cell_type = cell_types)
  )
  
  return(list(sce = sce, ann = cell_types))
}

# Test 1: CreateScEnsemble
test_that("CreateScEnsemble works", {
  data <- create_test_data()
  
  expect_no_error({
    scens <- CreateScEnsemble(data$sce, data$ann)
  })
  
  scens <- CreateScEnsemble(data$sce, data$ann)
  expect_true(is(scens, "S4"))
})

# Test 2: run_individual_algorithms  
test_that("run_individual_algorithms works", {
  data <- create_test_data()
  scens <- CreateScEnsemble(data$sce, data$ann)
  
  expect_no_error({
    scens <- run_individual_algorithms(scens)
  })
})

# Test 3: calculate_all_validation_indices
test_that("calculate_all_validation_indices works", {
  data <- create_test_data()
  scens <- CreateScEnsemble(data$sce, data$ann)
  scens <- run_individual_algorithms(scens)
  
  expect_no_error({
    scens <- calculate_all_validation_indices(scens)
  })
})

# Test 4: generate_all_hypergraphs
test_that("generate_all_hypergraphs works", {
  data <- create_test_data()
  scens <- CreateScEnsemble(data$sce, data$ann)
  scens <- run_individual_algorithms(scens)
  scens <- calculate_all_validation_indices(scens)
  
  expect_no_error({
    scens <- generate_all_hypergraphs(scens)
  })
})

# Test 5: ensemble_clustering (pipeline)
test_that("full pipeline works", {
  data <- create_test_data()
  
  expect_no_error({
    scens <- CreateScEnsemble(data$sce, data$ann)
    scens <- run_individual_algorithms(scens)
    scens <- calculate_all_validation_indices(scens)
    scens <- generate_all_hypergraphs(scens)
    scens <- ensemble_clustering(scens)
  })
})

# Test 6: Error handling
test_that("error handling works", {
  data <- create_test_data()
  
  # Wrong annotation length
  wrong_ann <- data$ann[1:10]
  expect_error(CreateScEnsemble(data$sce, wrong_ann))
  
  # NA in annotation
  na_ann <- data$ann
  na_ann[1] <- NA
  expect_error(CreateScEnsemble(data$sce, na_ann))
})
