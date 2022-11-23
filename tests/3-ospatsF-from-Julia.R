# Original Julian translation code testing

devtools::load_all()
library(dplyr)

#data("Nowley")
#dat <- Nowley[,-3]

# Needs data frame with x, y, value, variance
d0 <- readRDS("test-case-1.rds")
d <- d0 |> filter(obs) |>
  select(x, y, pred=var1.pred, var = var1.var)
dat <- d |>
  as.data.frame()

set.seed(1234)
s <- ospatsF_julia(data = dat,
             dRange = 3000,
             dRSquare = 1,
             dStrata = 3,
             dMaxrun = 30,
             nCycles = 3000,
)

# Refractor
set.seed(1234)
devtools::load_all()
tb <- system.time(
  sb <- ospatsF_julia_ref(x = d,
                          covmodel_range = 3000,
                          rsquared = 1,
                          nstrata = 3,
                          niter_outer = 30,
                          niter = 3000,
))

a  <- s$stratification
ab <- sb$stratification

print(table(a, ab) )
# looks the same to me.

s$OA
sb$OA
