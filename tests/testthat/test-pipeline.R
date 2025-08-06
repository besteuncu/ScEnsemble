test_that("full pipeline works", {

  
  expect_no_error({
    scens <- CreateScEnsemble(test_pollen, test_labels)
    scens <- run_individual_algorithms(scens)
    scens <- calculate_all_validation_indices(scens)
    scens <- generate_all_hypergraphs(scens)
    scens <- ensemble_clustering(scens)
  })
})
