test_that("calculate_all_validation_indices works", {
  Pollen <- PollenGliaData()
  ann <- colData(Pollen)[["Inferred Cell Type"]]
  scens <- CreateScEnsemble(Pollen, ann)
  scens <- run_individual_algorithms(scens)
  
  expect_no_error({
    scens <- calculate_all_validation_indices(scens)
  })
})
