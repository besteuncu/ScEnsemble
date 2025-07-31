test_that("CreateScEnsemble works with Pollen", {
  Pollen <- PollenGliaData()
  ann <- colData(Pollen)[["Inferred Cell Type"]]
  
  expect_no_error({
    scens <- CreateScEnsemble(Pollen, ann)
  })
  
  expect_s4_class(scens, "ScEnsemble")
})
