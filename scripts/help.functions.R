



# The likelihood function is now defined outside this function
lik.1par <- function(rho,m1,m2,m3,m4,x1_0,x2_0,h1,h2) {
  - log(2*pi*sqrt(1 - rho^2))*m1 - m2/(2*(1 - rho^2)) - m3/(2*(1 - rho^2)) +
    rho*m4/(1 - rho^2) - 1/2*exp(-1/2*(x2_0^2*h1^2 + x1_0^2 + x1_0^2*h2^2 -
                                         2*x1_0*rho*x2_0 + x2_0^2)/(-rho^2 + h2^2 + 1 + h1^2 + h1^2*h2^2))/
    (pi*(-rho^2 + h2^2 + 1 + h1^2 + h1^2*h2^2)^(1/2))
}
