# Debugging specific data

devtools::load_all()

dl <- readRDS("D:/temp/debug_ospats.rda")

set.seed(1)
r <- ospatsF_ref( x = dl$dat,
                  nstrata = dl$nstrata,
                  covmodel_range = dl$covrange)
