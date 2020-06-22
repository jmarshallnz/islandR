library(islandR)

context("Sampling distribution")

test_that("iterations works", {
  expect_equal(iterations(st), 100)
})

test_that("st_fit gives correct types", {
  expect_true(all(st$types %in% unique(manawatu$ST)))
  expect_equal(length(st$types), nrow(st$sequences))
  expect_equal(ncol(st$sequences), 7)

  expect_equal(class(st$sampling_distribution), "list")
  expect_equal(length(st$sampling_distribution), 100)
  expect_equal(dim(st$sampling_distribution[[1]]), c(348, 4))
  expect_equal(st$sources, dimnames(st$sampling_distribution[[1]])[[2]])
})

test_that("st_fit plots correctly", {
  p = plot(st)
  expect_equal(length(p), length(st$types))
})

test_that("st_fit summary", {
  s = summary(st)
  expect_equal(nrow(s), length(st$types))
  expect_equal(colnames(s), st$sources)
})

test_that("st_fit print", {
  expect_output(print(st), "Sampling distribution of genotypes")
  expect_output(print(st), paste0("Genotypes:[ ]+", nrow(st$types)))
  expect_output(print(st), "Model:[ ]+island")
  expect_output(print(st), "Sources:[ ]+4")
  expect_output(print(st), "Iterations:[ ]+100")
})
