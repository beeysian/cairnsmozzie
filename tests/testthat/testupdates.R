library(testthat)
library(mozzie)

context("Mozzies fly appropriate distances")
#' We're setting lambda to be 0.012 here as it roughly corresponds to the peak of a gamma(3,22) distribution

test_that("Mozzie don't fly too far",{
  expect_equal(as.numeric(random_dispersal(-16.918, 145.760, 0.012)),
               as.numeric(list(-16.918, 145.760)), tolerance = 1e-03)
})
