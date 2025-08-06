test_that("generate_all_hypergraphs works", {

  scens <- CreateScEnsemble(test_pollen, test_labels)
  scens <- run_individual_algorithms(scens)
  scens <- calculate_all_validation_indices(scens)
  
  expect_no_error({
    scens <- generate_all_hypergraphs(scens)
  })
})
