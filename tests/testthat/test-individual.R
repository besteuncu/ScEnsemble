test_that("run_individual_algorithms works", {
  data <- create_test_data()
  scens <- CreateScEnsemble(data$sce, data$ann)
  
  expect_no_error({
    scens <- run_individual_algorithms(scens)
  })
})