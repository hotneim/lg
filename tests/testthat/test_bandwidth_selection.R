library(lg)
context("Bandwidth selection")

set.seed(1)
n <- 50
x <- rt(n, df = 10)

result <- bw_select_cv_univariate(x)

test_that("Univariate bw selection works", {
    expect_true(is.numeric(result$bw))
    expect_equal(result$convergence, 0)
})

n <- 30
x <- mvtnorm::rmvt(n, df = 10, sigma <- diag(2))
x_3col <- cbind(x, rnorm(n))
est_method_wrong <- "BLA"

result1 <- bw_select_cv_bivariate(x, est_method = "1par")
result2 <- bw_select_cv_bivariate(x, est_method = "5par")
result3 <- bw_select_cv_bivariate(x, est_method = "5par_marginals_fixed", bw_marginal = c(1, 1))

test_that("Bivariate bw selection works", {
    expect_error(bw_select_cv_bivariate(x_3col),
                 "The data can only have 2 variables")
    expect_error(bw_select_cv_bivariate(x, est_method = est_method_wrong),
                 "Estimation method must be either '1par', '5par' or '5par_marginals_fixed'")
    expect_error(bw_select_cv_bivariate(x, est_method = "5par_marginals_fixed"),
                 "If estimation method is '5par_marginals_fixed', then bw_marginal must be supplied.")
    expect_true(is.numeric(result1$bw))
    expect_true(is.numeric(result2$bw))
    expect_true(is.numeric(result3$bw))
    expect_equal(result1$convergence, 0)
    expect_equal(result2$convergence, 0)
    expect_equal(result3$convergence, 0)
})


