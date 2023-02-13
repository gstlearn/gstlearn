    rm(list=ls())
    library(gstlearn)
    library(ggplot2)
    library(ggpubr)

    fileNF = paste(Sys.getenv('GSTLEARN_DATA'),'Scotland',
      'Scotland_Temperatures.NF',sep="/")
    dat = Db_createFromNF(fileNF)
    fileNF = paste(Sys.getenv('GSTLEARN_DATA'),'Scotland',
      'Scotland_Elevations.NF',sep="/")
    grid = DbGrid_createFromNF(fileNF)

------------------------------------------------------------------------

    p = plot.hist(dat,"January*")
    plot.decoration(p, title="Temperatures")


------------------------------------------------------------------------

    anam = AnamHermite(30)
    err = anam$fitFromLocator(dat)
    err = anam$rawToGaussian(dat, "January_temp")
    anam

    ## 
    ## Hermitian Anamorphosis
    ## ----------------------
    ## Minimum absolute value for Y  = -2.8
    ## Maximum absolute value for Y  = 2.7
    ## Minimum absolute value for Z  = 0.62599
    ## Maximum absolute value for Z  = 5.24756
    ## Minimum practical value for Y = -2.8
    ## Maximum practical value for Y = 2.7
    ## Minimum practical value for Z = 0.62599
    ## Maximum practical value for Z = 5.24756
    ## Mean                          = 2.81457
    ## Variance                      = 1.01677
    ## Number of Hermite polynomials = 30
    ## Normalized coefficients for Hermite polynomials (punctual variable)
    ##                [,  0]    [,  1]    [,  2]    [,  3]    [,  4]    [,  5]    [,  6]
    ##      [  0,]     2.815    -1.003     0.010     0.067     0.005     0.030    -0.007
    ##      [  7,]    -0.035     0.009     0.027    -0.011    -0.019     0.014     0.013
    ##      [ 14,]    -0.017    -0.008     0.019     0.004    -0.020    -0.001     0.020
    ##      [ 21,]    -0.002    -0.018     0.004     0.016    -0.005    -0.014     0.006
    ##      [ 28,]     0.011    -0.005

------------------------------------------------------------------------

    p = plot.XY(dat["Y.January_temp"], dat["January_temp"])
    plot.decoration(p, xlab="Gaussian", ylab="Raw")


------------------------------------------------------------------------

    p = plot.hist(dat,"Y.January*")
    plot.decoration(p, title="Temperatures (Gaussian scale)")


------------------------------------------------------------------------

We calculate the experimental directional variogram of the gaussian
scores and fit the Model (with the constraints that sill should be 1)

    varioparam = VarioParam_createMultiple(ndir=2, npas=40, dpas=10)
    vario_gauss2dir = Vario_create(varioparam, dat)
    err = vario_gauss2dir$compute()

    fitmodgauss = Model()
    err = fitmodgauss$fit(vario_gauss2dir, 
                          types=ECov_fromKeys(c("NUGGET", "SPHERICAL","CUBIC")),
                          constraints = Constraints(1))

    plot.varmod(vario_gauss2dir, fitmodgauss)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


------------------------------------------------------------------------

    neighU = NeighUnique_create()

    err = kriging(dat, grid, fitmodgauss, neighU)

    plot.setDefault(1, dims=c(8,8))
    p = plot(grid,"*estim")
    p = plot(dat, padd=p)
    plot.decoration(p, title="Kriging of Gaussian scores")


------------------------------------------------------------------------

    p = plot(grid,"*stdev")
    p = plot(dat, padd=p, flagCst=TRUE)
    plot.decoration(p, title="St. Dev. of Gaussian scores")


------------------------------------------------------------------------

Use the Turning Bands method with 1000 simulations

    selectivity = Selectivity_createByKeys(c("Z"), flag_est=TRUE, flag_std=TRUE)
    err = ConditionalExpectation(grid, anam, selectivity, "K*.estim", "K*.stdev",
                                 nbsimu=100,
                                 namconv=NamingConvention("CE",FALSE,TRUE,FALSE))

------------------------------------------------------------------------

    p = plot(grid,"CE*estim")
    p = plot(dat, padd=p)
    plot.decoration(p, title = "Conditional Expectation")


------------------------------------------------------------------------

    p = plot(grid, "CE*stdev")
    p = plot(dat, padd=p, flagCst=TRUE)
    plot.decoration(p, title="Conditional Standard Deviation")


------------------------------------------------------------------------

    selectivity = Selectivity_createByKeys(c("PROP"), zcuts=c(0),
                                           flag_est=TRUE, flag_std=TRUE)
    err = ConditionalExpectation(grid, anam, selectivity, 
                                 "K*.estim", "K*.stdev",
                                 namconv=NamingConvention("CE",FALSE,TRUE,FALSE))

    p = plot(grid,"CE.Proba*estim")
    p = plot(dat, padd=p)
    plot.decoration(p, title = "Conditional Probability below 0")


------------------------------------------------------------------------

    selectivity = Selectivity_createByKeys(c("T"), zcuts=c(1),
                                           flag_est=TRUE, flag_std=TRUE)
    err = ConditionalExpectation(grid, anam, selectivity, "K*.estim", "K*.stdev",
                                 namconv=NamingConvention("CE",FALSE,TRUE,FALSE))

------------------------------------------------------------------------

    p = plot(grid,"CE.T*estim-1")
    p = plot(dat, padd=p)
    plot.decoration(p, title = "Conditional Probability above 1")


------------------------------------------------------------------------

    p = plot(grid, "CE.T*stdev-1")
    p = plot(dat, padd=p, flagCst=TRUE)
    plot.decoration(p,title = "Conditional probability (Standard Deviation)")

