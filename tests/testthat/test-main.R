library(testthat)
library(ScEnsemble)

test_check("ScEnsemble")

library(SingleCellExperiment)
library(scRNAseq)

create_test_data <- function() {
  library(scRNAseq)
  sce <- PollenGliaData()
  ann <- colData(sce)[["Inferred Cell Type"]]
  colData(sce)$cell_type <- ann
  return(list(sce = sce, ann = ann))
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
