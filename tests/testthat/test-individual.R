test_that("run_individual_algorithms works", {

  scens <- CreateScEnsemble(test_pollen, test_labels)
  
  expect_no_error({
    scens <- run_individual_algorithms(scens)
  })
})
