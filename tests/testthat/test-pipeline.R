test_that("full pipeline works", {
  Pollen <- PollenGliaData()
  ann <- colData(Pollen)[["Inferred Cell Type"]]
  
  expect_no_error({
    scens <- CreateScEnsemble(Pollen, ann)
    scens <- run_individual_algorithms(scens)
    scens <- calculate_all_validation_indices(scens)
    scens <- generate_all_hypergraphs(scens)
    scens <- ensemble_clustering(scens)
  })
})
