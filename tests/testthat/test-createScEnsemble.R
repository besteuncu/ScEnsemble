test_that("CreateScEnsemble works with Pollen subset", {
  data <- create_test_data()
  
  expect_no_error({
    scens <- CreateScEnsemble(data$sce, data$ann)
  })
  
  expect_s4_class(scens, "ScEnsemble")
})
