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
