# Test the c code

devtools::load_all()
library(dplyr)

d0 <- readRDS("test-case-1.rds")

# Needs data frame with x, y, value, variance
dat <- d0 |>
  filter(obs) |>
  select(x, y, pred = var1.pred, var = var1.var)

set.seed(1)
t1 <- system.time(
  s <- ospatsF_ref(dat,
                   covmodel_range= 5000,
                   nstrata = 5,
                   verbose = 0, coolingrate = 0.995, temperature = 1,
                   niter = 5000, niter_outer = 100 # does not seem to do anything
  ))

set.seed(1)
t2 <- system.time(
  s2 <- ospatsF_refc(dat,
                   covmodel_range= 5000,
                   nstrata = 5,
                   verbose = 0, coolingrate = 0.995, temperature = 1,
                   niter = 5000, niter_outer = 100 # does not seem to do anything
  ))

print(rbind(r=t1, c=t2))

print(rbind(OA_r=c(s$OA, Obar=s$Obar),
            OA_c=c(s2$OA, Obar=s2$Obar)))



a <- s$stratification
b <- s2$stratification
print(table(a,b) )

# Check obar
of <- function(s) {
  N <- table(s)
  sum((split(dat$pred, s) |> sapply(\(y) sqrt(mean((y-mean(y))^2)) )) * N ) /sum(N)
}
print(c(of(a), of(b)))
#

par(mfrow=c(2,1))
image(z = matrix(a, 11), asp = 1)
image(z = matrix(b, 11), asp = 1)


