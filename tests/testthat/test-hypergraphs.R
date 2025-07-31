test_that("generate_all_hypergraphs works", {
  Pollen <- PollenGliaData()
  ann <- colData(Pollen)[["Inferred Cell Type"]]
  scens <- CreateScEnsemble(Pollen, ann)
  scens <- run_individual_algorithms(scens)
  scens <- calculate_all_validation_indices(scens)
  
  expect_no_error({
    scens <- generate_all_hypergraphs(scens)
  })
})
