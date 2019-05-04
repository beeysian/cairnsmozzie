library(testthat)
library(mozzie)

context("Correct age conversion")

test_that("An in-range age returns the correct stage",{
  expect_equal(init_juv_stage(0),1)
  expect_equal(init_juv_stage(5),1)
  expect_equal(init_juv_stage(6),2)
  expect_equal(init_juv_stage(13),3)
})
