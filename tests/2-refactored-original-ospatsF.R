# Refactored original code testing

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

# Compare to original
set.seed(1)
library(Ospats)
to <- system.time( so <- ospatsF(dat |> as.data.frame(),
             nCycles = 3000, dMaxrun = 30,dRSquare = 1, dStart = 0,
             dRange= 5000, verbose = FALSE, debug = FALSE,
             dStrata = 3,coolingRate = .95, initialTemperature = 1)
)


# Compare

a <- s$stratFrame
b <- so$stratFrame$ospats_strata
print(table(a,b) )

# ok, apart from arbitrary label order we got the same!

print(rbind(t1, to))

par(mfrow=c(2,1))
plot(d0$x, d0$y, col = d0$var1.pred, asp = 1, pch = 19)
plot(d0$x, d0$y, col = a, asp = 1, pch = 19)
