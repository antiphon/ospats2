# Test the c code

devtools::load_all()
library(dplyr)

d0 <- readRDS("test-case-2mc.rds")

# Needs data frame with x, y, value, variance
dat <- d0$pop |>
  select(x, y, pred = pred_mean, var = pred_var)
Z <- d0$Z

see0 <- 1

set.seed(see0 <- 1 + see0)

t2 <- system.time(
  s2 <- ospatsMC(Z,
                 nstrata = 3,
                 verbose = 2, coolingrate = 0.995, temperature = 1,
                 niter = 30, niter_outer = 2,
                 v = 2  )
)
print(t2)

set.seed(see0)
t1 <- system.time(
  s <- ospatsMC(Z,
                 nstrata = 3,
                 verbose = 2, coolingrate = 0.995, temperature = 1,
                 niter = 30, niter_outer = 2,
                 v = 1)
)#
#
#
#
 print(rbind(c1=t1, c2=t2))
#
 print(rbind(OA_c1=c(s$OA, Obar=s$Obar),
             OA_c2=c(s2$OA, Obar=s2$Obar)))
#
#
#
a <- s$stratification
b <- s2$stratification
print(table(a,b) )
#
# # Check obar
# of <- function(s) {
#   N <- table(s)
#   sum((split(dat$pred, s) |> sapply(\(y) sqrt(mean((y-mean(y))^2)) )) * N ) /sum(N)
# }
# print(c(of(a), of(b)))
# #
#
# par(mfrow=c(2,1))
# image(z = matrix(a, length(unique(dat$x))))
image(z = matrix(b, length(unique(dat$x))))
#
#
