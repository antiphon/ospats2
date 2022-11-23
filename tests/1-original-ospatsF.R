# Original code testing

library(Ospats)
library(dplyr)


d0 <- readRDS("test-case-1.rds")

# Needs data frame with x, y, value, variance
dat <- d0 |> filter(obs) |>
  select(x, y, var1.pred, var1.var) |>
  as.data.frame()
s <- ospatsF(data = dat,
             dRange = 5000,
             dStart = 0,
             dRSquare = 1,
             dStrata = 3,
             dMaxrun = 300,
             coolingRate = 0.9,
             initialTemperature = 1,
             verbose =FALSE,
             nCycles = 100, # does not seem to do anything
             )



## Illustrate
