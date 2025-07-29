test_that("error handling works", {
  data <- create_test_data()
  
  # Wrong annotation length
  wrong_ann <- data$ann[1:10]
  expect_error(CreateScEnsemble(data$sce, wrong_ann))
  
  # NA in annotation
  na_ann <- data$ann
  na_ann[1] <- NA
  expect_error(CreateScEnsemble(data$sce, na_ann))
})