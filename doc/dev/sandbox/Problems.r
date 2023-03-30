library(gstlearn)

# Here are collected all the misfunctioning codes in R

nech = 10
tab = cbind(runif(nech),runif(nech),rnorm(nech))
dat = Db_createFromSamples(nech, ELoadBy_SAMPLE(), as.numeric(tab))
dat$setLocators(c("New.1","New.2"),ELoc_X())
dat$setName("New.3","Var")
dat$setLocator("Var",ELoc_Z())  
dat 

# Problem with dumpNF (function of the mother class ASerializable)
# dat$dumpNF("Test_Dump_Db.NF")
# Only possibility:
ASerializable_dumpToNF(dat, "Test_Dump_Db")
 
# Problem using the basic assessors on Db
# This one crahes R!
# dat[1,2]

varioParamOmni = VarioParam_createOmniDirection(5, 1)
varioexp = Vario(varioParamOmni, dat)
err = varioexp$compute()

# Assessors on variogram
varioexp[1] # This one is OK and returns the second value of gg
# But the following fail
# varioexp[1]$nlag
# varioexp[1]$gg

# List of ENUM
model = Model()
# err = model$fit(varioexp, types=c(ECov_CUBIC(), ECov_EXPONENTIAL()))
err = model$fit(varioexp, types=ECov_fromValues(c(1,2,3,4)))
# Swig is waiting for: _p_std__vectorT_ECov_std__allocatorT_ECov_t_t
# but we transmit a list

