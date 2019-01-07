library(DEoptimR)
# NOTE: Examples were excluded from testing
#       to reduce package check time.

# Use a preset seed so test values are reproducible.
set.seed(1234)

# Bound-constrained optimization

#   Griewank function
#
#   -600 <= xi <= 600, i = {1, 2, ..., n}
#   The function has a global minimum located at
#   x* = (0, 0, ..., 0) with f(x*) = 0. Number of local minima
#   for arbitrary n is unknown, but in the two dimensional case
#   there are some 500 local minima.
#
#   Source:
#     Ali, M. Montaz, Khompatraporn, Charoenchai, and
#     Zabinsky, Zelda B. (2005).
#     A numerical evaluation of several stochastic algorithms
#     on selected continuous global optimization test problems.
#     Journal of Global Optimization 31, 635-672.
griewank <- function(x) {
  1 + crossprod(x)/4000 - prod( cos(x/sqrt(seq_along(x))) )
}
  
init <- JDEoptim(rep(-600, 10), rep(600, 10), griewank,
         tol = 1e-7, trace = FALSE, triter = 100, maxiter = 2000)

write.csv(rep(-600, 10), "Submodule=JDEOPTIM.Data=lower.csv")
write.csv(rep(600, 10), "Submodule=JDEOPTIM.Data=upper.csv")
write.csv(init$par, "Submodule=JDEOPTIM.Data=par.csv")
write.csv(init$value, "Submodule=JDEOPTIM.Data=value.csv")

