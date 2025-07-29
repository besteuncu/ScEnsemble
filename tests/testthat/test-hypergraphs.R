test_that("generate_all_hypergraphs works", {
  data <- create_test_data()
  scens <- CreateScEnsemble(data$sce, data$ann)
  scens <- run_individual_algorithms(scens)
  scens <- calculate_all_validation_indices(scens)
  
  expect_no_error({
    scens <- generate_all_hypergraphs(scens)
  })
})
