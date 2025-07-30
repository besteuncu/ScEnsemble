test_that("CreateScEnsemble works with Pollen subset", {
  Pollen <- PollenGliaData()
  ann <- colData(Pollen)[["Inferred Cell Type"]]
  
  expect_no_error({
    scens <- CreateScEnsemble(data, ann)
  })
  
  expect_s4_class(scens, "ScEnsemble")
})
