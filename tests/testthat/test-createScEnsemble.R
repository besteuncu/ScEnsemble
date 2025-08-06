test_that("CreateScEnsemble works with Pollen subset", {

  expect_no_error({
    scens <- CreateScEnsemble(test_pollen, test_labels)
  })
  
  expect_s4_class(scens, "ScEnsemble")
})
