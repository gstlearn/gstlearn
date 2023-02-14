In this preamble, we load the **gstlearn** library, clean the workspace.

    rm(list=ls())
    library(gstlearn)
    library(ggplot2)
    library(ggpubr)

Then we download the data base **dat**.

    fileNF = paste(Sys.getenv('GSTLEARN_DATA'),'Scotland',
      'Scotland_Temperatures.NF',sep="/")
    dat = Db_createFromNF(fileNF)

Calculate the experimental variogram **vario2dir** (in 2 directions)

    varioParamMulti = VarioParam_createMultiple(ndir=2, npas=15, dpas=15.)
    vario2dir = Vario(varioParamMulti, dat)
    err = vario2dir$compute()

------------------------------------------------------------------------

Calculate the fitted model **fitmodOK** (add the Universality Condition)

    fitmodOK = Model()
    types = ECov_fromKeys(c("NUGGET","EXPONENTIAL","GAUSSIAN"))
    err = fitmodOK$fit(vario2dir,types=types)
    err = fitmodOK$addDrift(Drift1())
    plot.varmod(vario2dir, fitmodOK)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


------------------------------------------------------------------------

Define the Unique Neighborhood **unique.neigh**:

    unique.neigh = NeighUnique()

Get the extension of the Data:

    dat$getExtremas()

    ## [[1]]
    ## [1]  78.2 460.7
    ## 
    ## [[2]]
    ## [1]  530.4 1208.9

------------------------------------------------------------------------

Create the Target file **grid**:

    grid = DbGrid_create(x0=c(65,535),dx=c(4.94, 4.96),nx=c(81,137))
    dbfmt = DbStringFormat_createFromFlags(flag_resume=FALSE, flag_vars=TRUE,
                                           flag_extend=TRUE, flag_stats=FALSE,
                                           flag_array=FALSE)
    grid$display(dbfmt)

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Extension
    ## -------------------
    ## Coor #1 - Min =     65.000 - Max =    460.200 - Ext = 395.2
    ## Coor #2 - Min =    535.000 - Max =   1209.560 - Ext = 674.56
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2

    ## NULL

------------------------------------------------------------------------

    err = kriging(dbin=dat, dbout=grid, model=fitmodOK, neighparam=unique.neigh,
                  flag_est=TRUE, flag_std=TRUE, flag_varz=FALSE,
                  namconv=NamingConvention("OK"))

------------------------------------------------------------------------

    p = plot.grid(grid, legend.name.raster="°C", zlim=c(0,8))
    p = plot.decoration(p, title="Ordinary Kriging over whole Grid")
    p


------------------------------------------------------------------------

    fileNF = paste(Sys.getenv('GSTLEARN_DATA'),'Scotland','Scotland_Elevations.NF', sep="/")
    grid = DbGrid_createFromNF(fileNF)
    grid$display(dbfmt)

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Extension
    ## -------------------
    ## Coor #1 - Min =     65.000 - Max =    455.123 - Ext = 390.123
    ## Coor #2 - Min =    535.000 - Max =   1200.109 - Ext = 665.109
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = Longitude - Locator = x1
    ## Column = 1 - Name = Latitude - Locator = x2
    ## Column = 2 - Name = Elevation - Locator = f1
    ## Column = 3 - Name = inshore - Locator = sel

    ## NULL

------------------------------------------------------------------------

The output grid now contains the selection **inshore**. Estimation is
restricted to the active cells only.

    err = kriging(dbin=dat, dbout=grid, model=fitmodOK, neighparam=unique.neigh,
                  flag_est=TRUE, flag_std=TRUE, flag_varz=FALSE,
                  namconv=NamingConvention("OK"))

------------------------------------------------------------------------

    p = plot.grid(grid, "OK*estim", legend.name.raster="°C", zlim=c(0.,8.))
    p = plot.point(dat,name_size="January_temp",padd=p,sizmax=300,color='black')
    p = plot.decoration(p, title="Estimation by Ordinary Kriging")
    p


------------------------------------------------------------------------

    p = plot.grid(grid,"OK*stdev", legend.name.raster="°C", zlim=c(0.,1.))
    p = plot.point(dat,name_size="January_temp",padd=p,sizmax=300,color='black')
    p = plot.decoration(p, title="St. dev. by Ordinary Kriging")
    p


------------------------------------------------------------------------

The Model **fitmodOK** is first duplicated into **fitmodSK**. Then the
Universality Condition is deleted.

    fitmodSK = fitmodOK$clone()
    err = fitmodSK$delDrift(rank=0)
    err = fitmodSK$setMean(ivar=0, mean=20.)

Simple Kriging is performed

    err = kriging(dbin=dat, dbout=grid, model=fitmodSK, neighparam=unique.neigh,
                  flag_est=TRUE, flag_std=TRUE, 
                  namconv=NamingConvention("SK"))

------------------------------------------------------------------------

    p = plot.grid(grid,"SK*estim",legend.name.raster="°C", zlim=c(0.,8.))
    p = plot.point(dat,name_size="January_temp",padd=p,sizmax=300,color='black')
    p = plot.decoration(p, title="Estimation by Simple Kriging")
    p


------------------------------------------------------------------------

    p = plot.grid(grid,"SK*stdev",legend.name.raster="°C", zlim=c(0.,1.))
    p = plot.point(dat,name_size="January_temp",padd=p,sizmax=300,color='black')
    p = plot.decoration(p, title="St. dev. by Simple Kriging")
    p


------------------------------------------------------------------------

    p = plot.correlation(grid,name1="OK*estim",name2="SK*estim", flagDiag=TRUE)
    p = plot.decoration(p, title="Estimation Simple vs. Ordinary", 
                        xlab="Ordinary Kriging", ylab="Simple Kriging")
    p


------------------------------------------------------------------------

    p = plot.correlation(grid,name1 = "OK*stdev",name2="SK*stdev", flagDiag=TRUE)
    p = plot.decoration(p, title="St. dev. Simple vs. Ordinary", 
                        xlab="Ordinary Kriging", ylab="Simple Kriging")
    p


------------------------------------------------------------------------

    err = xvalid(db=dat, model=fitmodOK, neighparam=unique.neigh, 
                 flag_xvalid_est=1, flag_xvalid_std=1,  
                 namconv=NamingConvention_create("Xvalid", flag_locator = FALSE))

------------------------------------------------------------------------

    p = plot.hist(dat,name="*esterr*",nbins=30,fill="blue")
    p = plot.decoration(p, xlab="Estimation Errors", title="Cross-Validation")
    p


------------------------------------------------------------------------

    p = plot.hist(dat,name="*stderr*",nbins=30,fill="blue")
    p = plot.decoration(p, xlab="Standardized Errors", title="Cross-Validation")
    p


------------------------------------------------------------------------

    mean(dat$getColumn("*esterr*"),na.rm=TRUE)

    ## [1] -0.004167872

    mean(dat$getColumn("*esterr*")^2,na.rm=TRUE)

    ## [1] 0.2393726

    mean(dat$getColumn("*stderr*")^2,na.rm=TRUE)

    ## [1] 75.16447

------------------------------------------------------------------------

    p = plot.grid(grid,"inshore")
    p = plot.point(dat,name_size="*esterr",padd=p,sizmax=300)
    p = plot.decoration(p, title="Cross-Validation scores")
    p


------------------------------------------------------------------------

    p = plot.grid(grid,"inshore")
    p = plot.point(dat,name_size="*esterr",padd=p,sizmax=300,flagAbsSize=TRUE)
    p = plot.decoration(p, title="Cross-Validation scores (abs. value)")
    p


------------------------------------------------------------------------

We design a small Moving Neighborhood **small.neigh** with only 1 sample
per neighborhood.

    small.neigh = NeighMoving_create(nmini=1, nmaxi=1, radius=1000000)

We perform Ordinary Kriging

    err = kriging(dbin=dat, dbout=grid, model=fitmodOK, neighparam=small.neigh,
                  flag_est=TRUE, flag_std=TRUE, 
                  namconv=NamingConvention("Small"))

------------------------------------------------------------------------

    p = plot.grid(grid,"Small*estim",zlim=c(0.,8.), legend.name.raster="°C")
    p = plot.point(dat,name_size="January_temp",padd=p)
    p = plot.decoration(p, title="Estimation by Ordinary Kriging (Small Neigh.)")
    p


------------------------------------------------------------------------

Building a reasonable Moving Neighborhood, although with a limited
extension (*radius*)

    moving.neigh = NeighMoving_create(nmini=1, nmaxi=10, radius=20)

Running the Ordinary Kriging

    err = kriging(dat,grid,fitmodOK,moving.neigh,
                  flag_est=TRUE, flag_std=TRUE, 
                  namconv=NamingConvention("Reduced"))

------------------------------------------------------------------------

    p = plot.grid(grid,"Reduced*estim", zlim=c(0.,8.), legend.name.raster="°C")
    p = plot.point(dat,name_size="January_temp",padd=p)
    p = plot.decoration(p, title="Estimation by Ordinary Kriging (Reduced Moving Neigh.)")
    p


Lots of target sites are not estimated as no sample is found within the
neighborhood.

------------------------------------------------------------------------

Building a reasonable Moving Neighborhood correctly tuned: 10 samples
(maximum) selected in a radius of 150 around the target site.

    moving.neigh = NeighMoving_create(nmini=1, nmaxi=10, radius=150)

Running the Ordinary Kriging

    err = kriging(dat,grid,fitmodOK,moving.neigh,
                  flag_est=TRUE, flag_std=TRUE, 
                  namconv=NamingConvention("Moving"))

------------------------------------------------------------------------

    p = plot.grid(grid,"Moving*estim",zlim=c(0.,8.), legend.name.raster="°C")
    p = plot.point(dat,name_size="January_temp",padd=p)
    p = plot.decoration(p, title="Estimation by Ordinary Kriging (Moving Neigh.)")
    p


------------------------------------------------------------------------

    p = plot.grid(grid,"Moving*stdev", zlim=c(0.,1.), legend.name.raster="°C")
    p = plot.point(dat,name_size="January_temp",padd=p)
    p = plot.decoration(p, title="St. dev. by Ordinary Kriging (Moving Neigh.)")
    p


------------------------------------------------------------------------

    p = plot.correlation(grid,name1 = "OK*estim",name2="Moving*estim", flagDiag=TRUE)
    p = plot.decoration(p, title="Ordinary Kriging Estimation", 
                        xlab="Unique", ylab="Moving")
    p


------------------------------------------------------------------------

    p = plot.correlation(grid,name1 = "OK*stdev",name2="Moving*stdev", flagDiag=TRUE)
    p = plot.decoration(p, title="Ordinary Kriging St. Dev.", 
                        xlab="Unique", ylab="Moving")
    p

