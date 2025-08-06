test_that("calculate_all_validation_indices works", {

  scens <- CreateScEnsemble(test_pollen, test_labels)
  scens <- run_individual_algorithms(scens)
  
  expect_no_error({
    scens <- calculate_all_validation_indices(scens)
  })
})
