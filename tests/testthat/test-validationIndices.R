
test_that("calculate_all_validation_indices works", {
  data <- create_test_data()
  scens <- CreateScEnsemble(data$sce, data$ann)
  scens <- run_individual_algorithms(scens)
  
  expect_no_error({
    scens <- calculate_all_validation_indices(scens)
  })
})
