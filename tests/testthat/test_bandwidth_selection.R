library(lg)
context("Bandwidth selection")

n <- 50
x <- rt(n, df = 10)

result <- bw_select_cv_univariate(x)

test_that("Univariate bw selection works", {
    expect_true(is.numeric(result$bw))
    expect_equal(result$convergence, 0)
})
