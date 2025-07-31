test_that("run_individual_algorithms works", {
  Pollen <- PollenGliaData()
  ann <- colData(Pollen)[["Inferred Cell Type"]]
  scens <- CreateScEnsemble(Pollen, ann)
  
  expect_no_error({
    scens <- run_individual_algorithms(scens)
  })
})
