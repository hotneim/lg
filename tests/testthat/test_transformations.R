library(lg)
library(lgde)
context("Transformations")

n <- 100
set.seed(1)

x <- mvtnorm::rmvt(n, df = 10, sigma = diag(3))
grid <- cbind(c(1,2,3), c(1,2,3), c(1,2,3))

new <- trans_normal(x)
old <- transLocal(x, grid = grid, return.normalizing.constants = TRUE)

test_that("New equals old", {
    expect_equal(new$transformed_data, old$transformed.data)
    expect_equal(new$trans_new(grid)$trans, old$transformed.grid)
    expect_equal(new$trans_new(grid)$normalizing_constants,
                 old$normalizing.constants)
})
