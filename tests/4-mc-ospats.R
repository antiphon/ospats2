# Test the MC target

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
                   nstrata = 3,
                   verbose = FALSE, coolingrate = 0.95, temperature = 1,
                   niter = 3000, niter_outer = 30 # does not seem to do anything
  ))
