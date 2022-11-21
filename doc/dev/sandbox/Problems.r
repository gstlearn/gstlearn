library(gstlearn)

# Here are collected all the misfunctioning codes in R

nech = 10
tab = cbind(runif(nech),runif(nech),rnorm(nech))
dat = Db_createFromSamples(nech, ELoadBy_SAMPLE(), as.numeric(tab))
dat$setLocators(c("New.1","New.2"),ELoc_X())
dat$setName("New.3","Var")
dat$setLocator("Var",ELoc_Z())  
dat 

# Problem with dumpNF (function of the Mother class)
# dat$dumpNF("Test_Dump_Db.NF")
 
# Problem using the basic assessors on Db
# dat[1,2]
# dat[]

varioParamOmni = VarioParam_createOmniDirection(5, 1)
varioexp = Vario(varioParamOmni, dat)
err = varioexp$compute()

# Assessors sur le variogramme
# varioexp[1]
# varioexp[1]$nlag
# varioexp[1]$gg

# List of ENUM
model = Model()
# err = model$fit(varioexp, types=c(ECov_CUBIC(), ECov_EXPONENTIAL()))
err = model$fit(varioexp, types=ECov_fromValues(c(1,2,3,4)))
# Le swig attend: _p_std__vectorT_ECov_std__allocatorT_ECov_t_t"
# et on lui transmet un: list

grid = DbGrid::create([10,10], [0.1,0.1])
neigh = NeighUnique()
# Probleme lie a l'impossibilite de nommer les arguments
err = kriging(dat, grid, model, neigh)

# Ce qui ne fonctionne pas
# err = kriging(dat, grid, model, neigh, namconv=NamingConvention("OK"))
#
# Correlativement, placer tous les arguments en utilisant les valeurs par defaut
# err = kriging(dat, grid, model, neigh, TRUE, TRUE, FALSE,
# VectorInt(), VectorInt(), VectorVectorDouble(), NamingConvention("OK"))

