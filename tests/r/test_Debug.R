# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))
set.seed(32421)

np = 100
seed = 14313

data = Db_createFillRandom(np, 2, 0)
model = Model_createFromParam(ECov_EXPONENTIAL(), 0.5, 1.)
model$setMeans(0)
model$setDriftIRF(-1)

simtub(NULL, data, model, NULL, 1, seed)

nburn = 10
niter = 30
nbsimu = 5

gibbs_sampler(
    data, model, nbsimu, seed, nburn, niter, TRUE, TRUE, FALSE, FALSE, FALSE, 0,
    5., FALSE, FALSE, TRUE)

data$getContentsAsTable()$display()
