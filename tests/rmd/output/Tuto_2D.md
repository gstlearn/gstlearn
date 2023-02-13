Global variables

    verbose = TRUE
    graphics = TRUE
    OptCst_define(ECst_NTCOL(),6)

    ## NULL

Reading data
============

The data are stored in a CSV format in the file called Pollution.dat

** At this stage: - it is not possible to pass the list of variable
names *c("X","Y")*. - to pass the FLAGS for *DbStringFormat* (they are
not an ENUM) **

    dataDir = "../../data/Pollution"
    filepath = paste(dataDir,"Pollution.dat",sep="/")
    mydb = Db_createFromCSV(filepath,CSVformat())
    err = mydb$setLocator("X",ELoc_X(),0)
    err = mydb$setLocator("Y",ELoc_X(),1)
    err = mydb$setLocator("Zn",ELoc_Z())
    if (verbose)
    {
      dbfmt = DbStringFormat()
      dbfmt$setFlags(flag_extend = TRUE)
      mydb$display(dbfmt)
    }

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 5
    ## Maximum Number of UIDs       = 5
    ## Total number of samples      = 102
    ## 
    ## Data Base Extension
    ## -------------------
    ## Coor #1 - Min =    109.850 - Max =    143.010 - Ext = 33.16
    ## Coor #2 - Min =    483.660 - Max =    513.040 - Ext = 29.38
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = X - Locator = x1
    ## Column = 2 - Name = Y - Locator = x2
    ## Column = 3 - Name = Zn - Locator = z1
    ## Column = 4 - Name = Pb - Locator = p1

    ## NULL

Accessing to the variable names

    cat("List of all variable names =",mydb$getAllNames())

    ## List of all variable names = rank X Y Zn Pb

Extracting the vector containing the Zn variable in order to perform a
selection

    tabZn = mydb$getColumn("Zn")
    selZn = as.numeric(tabZn < 20)
    mydb$addSelection(selZn,"sel")

    ## [1] 5

    mydb$setLocator('Pb',ELoc_Z())

    ## NULL

    if (verbose)
        mydb$display()

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 6
    ## Maximum Number of UIDs       = 6
    ## Total number of samples      = 102
    ## Number of active samples     = 99
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = X - Locator = x1
    ## Column = 2 - Name = Y - Locator = x2
    ## Column = 3 - Name = Zn - Locator = NA
    ## Column = 4 - Name = Pb - Locator = z1
    ## Column = 5 - Name = sel - Locator = sel

    ## NULL

Display my Data (with samples represented by color and size)

    if (graphics)
        plot.point(mydb,color_name="Pb",title="Data Set")

    ## Warning in geom_point(data = df, mapping = aes(x = x, y = y, color = colval, :
    ## Ignoring unknown parameters: `colour_name` and `title`


Variograms
==========

We first define the geometry of the variogram calculations

    myVarioParamOmni = VarioParam()
    mydir = DirParam_create(npas=10,dpas=1.)
    myVarioParamOmni$addDir(mydir)

    ## NULL

We use the variogram definition in order to calculate the variogram
cloud.

    dbcloud = db_variogram_cloud(db=mydb, varioparam=myVarioParamOmni)

We recall that the Variogram cloud is calculated by filling an
underlying grid where each cell is painted according to the number of
pairs at the given distance and given variability. Representing the
variogram cloud.

    if (graphics)
        plot.grid(dbcloud,"Cloud*",title="Variogram Cloud")

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `title`


Calculating the experimental omni-directional variogram

    myVarioOmni = Vario(myVarioParamOmni,mydb)
    err = myVarioOmni$compute(ECalcVario_VARIOGRAM())
    if (verbose)
        myVarioOmni$display()

    ## 
    ## Variogram characteristics
    ## =========================
    ## Number of variable(s)       = 1
    ## Number of direction(s)      = 1
    ## Space dimension             = 2
    ## Variance-Covariance Matrix     2.881
    ## 
    ## Direction #1
    ## ------------
    ## Number of lags              = 10
    ## Direction coefficients      =      1.000     0.000
    ## Direction angles (degrees)  =      0.000     0.000
    ## Tolerance on direction      =     90.000 (degrees)
    ## Calculation lag             =      1.000
    ## Tolerance on distance       =     50.000 (Percent of the lag value)
    ## 
    ## For variable 1
    ##       Rank    Npairs  Distance     Value
    ##          0     3.000     0.389     0.462
    ##          1   123.000     1.081     1.495
    ##          2   183.000     2.038     1.620
    ##          3   205.000     3.006     2.526
    ##          4   231.000     4.013     2.240
    ##          5   229.000     5.036     2.524
    ##          6   198.000     5.962     2.396
    ##          7   187.000     7.000     2.708
    ##          8   204.000     7.996     2.772
    ##          9   184.000     8.990     2.868

    ## NULL

The variogram is represented graphically for a quick check

    if (graphics)
        plot.varmod(myVarioOmni,
                    title="Omni-directional Variogram for Pb")

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), ...): Ignoring
    ## unknown parameters: `title`


Calculate a variogram in several directions

    myvarioParam = VarioParam()
    mydirs = DirParam_createMultiple(ndir=4, npas=10, dpas=1.)
    myvarioParam$addMultiDirs(mydirs)

    ## NULL

    myvario = Vario(myvarioParam,mydb)
    myvario$compute(ECalcVario_VARIOGRAM())

    ## [1] 0

    if (verbose)
        myvario$display()

    ## 
    ## Variogram characteristics
    ## =========================
    ## Number of variable(s)       = 1
    ## Number of direction(s)      = 4
    ## Space dimension             = 2
    ## Variance-Covariance Matrix     2.881
    ## 
    ## Direction #1
    ## ------------
    ## Number of lags              = 10
    ## Direction coefficients      =      1.000     0.000
    ## Direction angles (degrees)  =      0.000     0.000
    ## Tolerance on direction      =     22.500 (degrees)
    ## Calculation lag             =      1.000
    ## Tolerance on distance       =     50.000 (Percent of the lag value)
    ## 
    ## For variable 1
    ##       Rank    Npairs  Distance     Value
    ##          0     1.000     0.410     0.180
    ##          1    29.000     1.094     1.634
    ##          2    47.000     2.079     1.415
    ##          3    53.000     3.003     2.824
    ##          4    63.000     3.999     2.348
    ##          5    66.000     5.035     2.319
    ##          6    60.000     5.978     3.115
    ##          7    52.000     7.045     2.746
    ##          8    52.000     8.020     3.927
    ##          9    37.000     8.980     2.554
    ## 
    ## Direction #2
    ## ------------
    ## Number of lags              = 10
    ## Direction coefficients      =      0.707     0.707
    ## Direction angles (degrees)  =     45.000     0.000
    ## Tolerance on direction      =     22.500 (degrees)
    ## Calculation lag             =      1.000
    ## Tolerance on distance       =     50.000 (Percent of the lag value)
    ## 
    ## For variable 1
    ##       Rank    Npairs  Distance     Value
    ##          0     1.000     0.344     0.080
    ##          1    31.000     1.051     1.113
    ##          2    50.000     1.960     1.890
    ##          3    62.000     2.999     2.443
    ##          4    58.000     4.014     2.701
    ##          5    51.000     5.016     2.702
    ##          6    36.000     5.999     1.833
    ##          7    37.000     7.015     2.130
    ##          8    50.000     7.997     2.060
    ##          9    53.000     8.995     2.381
    ## 
    ## Direction #3
    ## ------------
    ## Number of lags              = 10
    ## Direction coefficients      =      0.000     1.000
    ## Direction angles (degrees)  =     90.000     0.000
    ## Tolerance on direction      =     22.500 (degrees)
    ## Calculation lag             =      1.000
    ## Tolerance on distance       =     50.000 (Percent of the lag value)
    ## 
    ## For variable 1
    ##       Rank    Npairs  Distance     Value
    ##          1    32.000     1.149     1.631
    ##          2    39.000     2.080     1.670
    ##          3    39.000     2.979     2.511
    ##          4    48.000     4.012     2.120
    ##          5    51.000     5.029     3.055
    ##          6    47.000     5.939     2.856
    ##          7    49.000     6.965     2.386
    ##          8    42.000     7.952     2.708
    ##          9    41.000     9.018     2.320
    ## 
    ## Direction #4
    ## ------------
    ## Number of lags              = 10
    ## Direction coefficients      =     -0.707     0.707
    ## Direction angles (degrees)  =    135.000     0.000
    ## Tolerance on direction      =     22.500 (degrees)
    ## Calculation lag             =      1.000
    ## Tolerance on distance       =     50.000 (Percent of the lag value)
    ## 
    ## For variable 1
    ##       Rank    Npairs  Distance     Value
    ##          0     1.000     0.411     1.125
    ##          1    31.000     1.028     1.606
    ##          2    47.000     2.044     1.496
    ##          3    51.000     3.040     2.330
    ##          4    62.000     4.028     1.791
    ##          5    61.000     5.058     2.155
    ##          6    55.000     5.939     1.587
    ##          7    49.000     6.975     3.425
    ##          8    60.000     8.004     2.408
    ##          9    53.000     8.972     3.996

    ## NULL

    if (graphics)
        plot.varmod(myvario,title="Multi-Directional Variogram of Pb")

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), ...): Ignoring unknown parameters: `title`
    ## Ignoring unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), ...): Ignoring
    ## unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), ...): Ignoring
    ## unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


Calculating the Variogram Map

    myvmap = db_vmap_compute(db=mydb,calcul_type=ECalcVario_VARIOGRAM(),nxx=c(20,20))
    if (verbose)
        myvmap$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 5
    ## Maximum Number of UIDs       = 5
    ## Total number of samples      = 1681
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    -33.160   -29.380
    ## Mesh   :      1.658     1.469
    ## Number :         41        41
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2
    ## Column = 3 - Name = VMAP.Pb.Var - Locator = z1
    ## Column = 4 - Name = VMAP.Pb.Nb - Locator = NA

    ## NULL

    if (graphics)
        plot.grid(myvmap,"*Var",title="Variogram Map")

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `title`


Model
=====

Fitting a Model. We call the Automatic Fitting procedure providing the
list of covariance functions to be tested.

    mymodel = Model_createFromDb(mydb)
    err = mymodel$fit(vario=myvario,types=ECov_fromKeys(c("EXPONENTIAL","SPHERICAL")))

Visualizing the resulting model, overlaid on the experimental variogram

    if (graphics)
        plot.varmod(myvario,mymodel,title="Model for Pb")

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), ...): Ignoring
    ## unknown parameters: `title`

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), na.rm = TRUE, :
    ## Ignoring unknown parameters: `title`

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), ...): Ignoring
    ## unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), na.rm = TRUE, : Ignoring unknown parameters: `title`
    ## Ignoring unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), na.rm = TRUE, : Ignoring unknown parameters: `title`
    ## Ignoring unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), na.rm = TRUE, :
    ## Ignoring unknown parameters: `title`


Model with equality constraints
-------------------------------

We can impose some constraints on the parameters during the fit. For
instance here, we impose an equality constraint on the range (range =
1).

    myModelConstrained = Model_createFromDb(mydb)
    constr = Constraints()
    paramid = CovParamId(0,0,EConsElem_RANGE(),0,0)
    constr$addItem(ConsItem(paramid,EConsType_EQUAL(),1.))

    ## NULL

    err = myModelConstrained$fit(vario=myVarioOmni,types=ECov_fromKeys(c("EXPONENTIAL","SPHERICAL")),constraints=constr)
    myModelConstrained

    ## 
    ## Model characteristics
    ## =====================
    ## Space dimension              = 2
    ## Number of variable(s)        = 1
    ## Number of basic structure(s) = 2
    ## Number of drift function(s)  = 0
    ## Number of drift equation(s)  = 0
    ## 
    ## Covariance Part
    ## ---------------
    ## Exponential
    ## - Sill         =      1.032
    ## - Range        =      1.000
    ## - Theo. Range  =      0.334
    ## Spherical
    ## - Sill         =      1.605
    ## - Range        =      5.880
    ## Total Sill     =      2.638

We can impose inequality constraints by using EConsType.LOWER or
EConsType.UPPER.

Adding a drift
--------------

    mymodel$addDrift(Drift1(mymodel$getContext()))

    ## NULL

    if (verbose)
        mymodel$display()

    ## 
    ## Model characteristics
    ## =====================
    ## Space dimension              = 2
    ## Number of variable(s)        = 1
    ## Number of basic structure(s) = 2
    ## Number of drift function(s)  = 1
    ## Number of drift equation(s)  = 1
    ## 
    ## Covariance Part
    ## ---------------
    ## Exponential
    ## - Sill         =      1.035
    ## - Ranges       =      1.786     0.366
    ## - Theo. Ranges =      0.596     0.122
    ## - Angles       =    405.023     0.000
    ## - Rotation Matrix
    ##                [,  0]    [,  1]
    ##      [  0,]     0.707    -0.707
    ##      [  1,]     0.707     0.707
    ## Spherical
    ## - Sill         =      1.621
    ## - Ranges       =      7.051     5.132
    ## - Angles       =   1576.897     0.000
    ## - Rotation Matrix
    ##                [,  0]    [,  1]
    ##      [  0,]    -0.730    -0.683
    ##      [  1,]     0.683    -0.730
    ## Total Sill     =      2.656
    ## 
    ## Drift Part
    ## ----------
    ## Universality Condition

    ## NULL

Defining the Neighborhood
=========================

We initiate a Neigbourhood (Moving with a small number of samples for
Demonstration)

    myneigh = NeighMoving_create(flag_xvalid=FALSE,nmaxi=6,radius=10)
    if (verbose)
        myneigh$display()

    ## 
    ## Moving Neighborhood
    ## ===================
    ## Space dimension = 2
    ## Minimum number of samples           = 1
    ## Maximum number of samples           = 6
    ## Maximum horizontal distance         = 10

    ## NULL

Checking the Moving Neighborhood
--------------------------------

We must first create a Grid which covers the area of interest

    mygrid = DbGrid_createCoveringDb(dbin=mydb,nodes=c(80,72),dcell=c(0.5,0.5),
                                     origin=c(107.,481.),margin=c(2,2))
    if (verbose)
        mygrid$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 2
    ## Maximum Number of UIDs       = 2
    ## Total number of samples      = 5913
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    105.000   479.000
    ## Mesh   :      0.500     0.500
    ## Number :         81        73
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2

    ## NULL

We can now test the neighborhood characteristics for each node of the
previously defined grid.

    err = test_neigh(mydb,mygrid,mymodel,myneigh)
    if (verbose)
        mygrid$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 7
    ## Maximum Number of UIDs       = 7
    ## Total number of samples      = 5913
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    105.000   479.000
    ## Mesh   :      0.500     0.500
    ## Number :         81        73
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = Neigh.Pb.Number - Locator = NA
    ## Column = 3 - Name = Neigh.Pb.MaxDist - Locator = NA
    ## Column = 4 - Name = Neigh.Pb.MinDist - Locator = NA
    ## Column = 5 - Name = Neigh.Pb.NbNESect - Locator = NA
    ## Column = 6 - Name = Neigh.Pb.NbCESect - Locator = z1

    ## NULL

We can visualize some of the newly created variables, such as:

-   the number of points per neighborhood

<!-- -->

    if (graphics)
        plot.grid(mygrid,"Neigh*Number",
                  title="Number of Samples per Neighborhood")

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `title`


-   the one giving the maximum distance per neighborhood

<!-- -->

    if (graphics)
        plot.grid(mygrid,"Neigh*MaxDist",
                  title="Maximum Distance per Neighborhood")

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `title`


Cross-validation
================

We can now process the cross-validation step

    err = xvalid(mydb,mymodel,myneigh)
    if (verbose)
        mydb$display()

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 8
    ## Maximum Number of UIDs       = 8
    ## Total number of samples      = 102
    ## Number of active samples     = 99
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = X - Locator = x1
    ## Column = 2 - Name = Y - Locator = x2
    ## Column = 3 - Name = Zn - Locator = NA
    ## Column = 4 - Name = Pb - Locator = NA
    ## Column = 5 - Name = sel - Locator = sel
    ## Column = 6 - Name = Xvalid.Pb.esterr - Locator = z1
    ## Column = 7 - Name = Xvalid.Pb.stderr - Locator = NA

    ## NULL

    if (graphics)
        plot.hist(mydb,"Xvalid.Pb.stderr")


Estimating by Kriging
=====================

We now perform the Estimation by Ordinary Kriging. The Neighborhood is
changed into a Unique Neighborhood.

    mydb$setLocator("Pb",ELoc_Z())

    ## NULL

    myneigh = NeighUnique_create()
    err = kriging(mydb,mygrid,mymodel,myneigh)
    if (verbose)
        mygrid$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 9
    ## Maximum Number of UIDs       = 9
    ## Total number of samples      = 5913
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    105.000   479.000
    ## Mesh   :      0.500     0.500
    ## Number :         81        73
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = Neigh.Pb.Number - Locator = NA
    ## Column = 3 - Name = Neigh.Pb.MaxDist - Locator = NA
    ## Column = 4 - Name = Neigh.Pb.MinDist - Locator = NA
    ## Column = 5 - Name = Neigh.Pb.NbNESect - Locator = NA
    ## Column = 6 - Name = Neigh.Pb.NbCESect - Locator = NA
    ## Column = 7 - Name = Kriging.Pb.estim - Locator = z1
    ## Column = 8 - Name = Kriging.Pb.stdev - Locator = NA

    ## NULL

Visualizing the results

    if (graphics)
        ax = plot.grid(mygrid,"Kriging.Pb.estim",end.plot = FALSE)

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `end.plot`

        plot.point(mydb,"Pb",title="Estimate of Pb",padd=ax)

    ## Warning in geom_point(data = df, mapping = aes(x = x, y = y, color = colval, :
    ## Ignoring unknown parameters: `title`


    if (graphics)
        ax = plot.grid(mygrid,"Kriging.Pb.stdev", end.plot=FALSE)

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `end.plot`

        plot.point(mydb,"Pb",title="St. Deviation of Pb",padd=ax)

    ## Warning in geom_point(data = df, mapping = aes(x = x, y = y, color = colval, :
    ## Ignoring unknown parameters: `title`


Simulations
===========

We must first transform the Data into Gaussian

    myanamPb = AnamHermite_create(nbpoly=30)
    err = myanamPb$fitFromLocator(mydb)
    if (verbose)
        myanamPb

    ## 
    ## Hermitian Anamorphosis
    ## ----------------------
    ## Minimum absolute value for Y  = -2.7
    ## Maximum absolute value for Y  = 2.6
    ## Minimum absolute value for Z  = 3.0029
    ## Maximum absolute value for Z  = 12.9777
    ## Minimum practical value for Y = -2.7
    ## Maximum practical value for Y = 2.6
    ## Minimum practical value for Z = 3.0029
    ## Maximum practical value for Z = 12.9777
    ## Mean                          = 5.65758
    ## Variance                      = 2.86296
    ## Number of Hermite polynomials = 30
    ## Normalized coefficients for Hermite polynomials (punctual variable)
    ##                [,  0]    [,  1]    [,  2]    [,  3]    [,  4]    [,  5]    [,  6]
    ##      [  0,]     5.658    -1.625     0.440    -0.069    -0.017     0.082    -0.061
    ##      [  7,]     0.001     0.036    -0.044     0.004     0.047    -0.030    -0.029
    ##      [ 14,]     0.037     0.007    -0.031     0.010     0.018    -0.019    -0.003
    ##      [ 21,]     0.019    -0.010    -0.014     0.019     0.006    -0.023     0.004
    ##      [ 28,]     0.022    -0.013

We can produce the Gaussian Anamorphosis graphically within its
definition domain.

    if (graphics)
        plot.anam(myanamPb)

The next step consists in translating the target variable ('Pb') into
its Gaussian transform

    mydb$setLocator("Pb",ELoc_Z())

    ## NULL

    err = myanamPb$rawToGaussianByLocator(mydb)
    if (verbose)
        mydb$display()

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 9
    ## Maximum Number of UIDs       = 9
    ## Total number of samples      = 102
    ## Number of active samples     = 99
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = X - Locator = x1
    ## Column = 2 - Name = Y - Locator = x2
    ## Column = 3 - Name = Zn - Locator = NA
    ## Column = 4 - Name = Pb - Locator = NA
    ## Column = 5 - Name = sel - Locator = sel
    ## Column = 6 - Name = Xvalid.Pb.esterr - Locator = NA
    ## Column = 7 - Name = Xvalid.Pb.stderr - Locator = NA
    ## Column = 8 - Name = Y.Pb - Locator = z1

    ## NULL

We quickly calculate experimental (omni-directional) variograms using
the already defined directions.

    myvarioParam = VarioParam()
    mydir = DirParam_create(npas=10,dpas=1.)
    myvarioParam$addDir(mydir)

    ## NULL

    myVario = Vario(myvarioParam,mydb)
    err = myvario$compute(ECalcVario_VARIOGRAM())

We fit the model by automatic fit (with the constraints that the total
sill be equal to 1).

    mymodelG = Model_createFromDb(mydb)
    err = mymodelG$fit(myvario,types=ECov_fromKeys(c("EXPONENTIAL")))
    if (graphics)
        plot.varmod(myvario,mymodelG,title="Model for Gaussian Pb")

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), ...): Ignoring
    ## unknown parameters: `title`

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), na.rm = TRUE, :
    ## Ignoring unknown parameters: `title`

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), ...): Ignoring
    ## unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), na.rm = TRUE, : Ignoring unknown parameters: `title`
    ## Ignoring unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), na.rm = TRUE, : Ignoring unknown parameters: `title`
    ## Ignoring unknown parameters: `title`

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning in geom_line(data = df, mapping = aes(x = hh, y = gg), na.rm = TRUE, :
    ## Ignoring unknown parameters: `title`


We perform a set of 10 conditional simulations using the Turning Bands
Method.

    err = simtub(mydb,mygrid,mymodel,myneigh,nbsimu=10)
    if (verbose)
        mygrid$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 19
    ## Maximum Number of UIDs       = 19
    ## Total number of samples      = 5913
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    105.000   479.000
    ## Mesh   :      0.500     0.500
    ## Number :         81        73
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = Neigh.Pb.Number - Locator = NA
    ## Column = 3 - Name = Neigh.Pb.MaxDist - Locator = NA
    ## Column = 4 - Name = Neigh.Pb.MinDist - Locator = NA
    ## Column = 5 - Name = Neigh.Pb.NbNESect - Locator = NA
    ## Column = 6 - Name = Neigh.Pb.NbCESect - Locator = NA
    ## Column = 7 - Name = Kriging.Pb.estim - Locator = NA
    ## Column = 8 - Name = Kriging.Pb.stdev - Locator = NA
    ## Column = 9 - Name = Simu.Y.Pb.1 - Locator = z1
    ## Column = 10 - Name = Simu.Y.Pb.2 - Locator = z2
    ## Column = 11 - Name = Simu.Y.Pb.3 - Locator = z3
    ## Column = 12 - Name = Simu.Y.Pb.4 - Locator = z4
    ## Column = 13 - Name = Simu.Y.Pb.5 - Locator = z5
    ## Column = 14 - Name = Simu.Y.Pb.6 - Locator = z6
    ## Column = 15 - Name = Simu.Y.Pb.7 - Locator = z7
    ## Column = 16 - Name = Simu.Y.Pb.8 - Locator = z8
    ## Column = 17 - Name = Simu.Y.Pb.9 - Locator = z9
    ## Column = 18 - Name = Simu.Y.Pb.10 - Locator = z10

    ## NULL

Some statistics on the Conditional simulations in Gaussian scale

\*\* Still impossible due to use of VectorString \*\*

    if (verbose)
    {
      err = mygrid$statisticsByLocator(locatorType = ELoc_Z(),
                                       opers = EStatOption_fromKeys(c("MEAN","STDV","MINI","MAXI")))
    }

    ##                  MEAN      STDV      MINI      MAXI
    ## Simu.Y.Pb.1      -0.052     1.582    -6.234     6.963
    ## Simu.Y.Pb.2       0.081     1.451    -5.167     6.113
    ## Simu.Y.Pb.3       0.133     1.529    -5.363     5.639
    ## Simu.Y.Pb.4       0.453     1.478    -6.322     5.836
    ## Simu.Y.Pb.5      -0.084     1.577    -6.403     5.578
    ## Simu.Y.Pb.6       0.316     1.551    -4.816     6.297
    ## Simu.Y.Pb.7       0.185     1.598    -6.441     6.127
    ## Simu.Y.Pb.8       0.705     1.599    -4.713     7.175
    ## Simu.Y.Pb.9       0.159     1.639    -6.660     7.135
    ## Simu.Y.Pb.10      0.310     1.553    -5.955     4.807

We visualize a conditional simulation in Gaussian scale

    if (graphics)
    {
        ax = plot.grid(mygrid,"Simu.Y.Pb.1",end.plot=FALSE)
        plot.point(mydb,"Pb",title="One Simulation of Pb in Gaussian Scale",padd=ax)
    }

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `end.plot`

    ## Warning in geom_point(data = df, mapping = aes(x = x, y = y, color = colval, :
    ## Ignoring unknown parameters: `title`


We turn the Gaussian conditional simulations into Raw scale (using the
Anamorphosis back transform) and get rid of the Gaussian conditional
simulations.

    err = myanamPb$gaussianToRaw(mygrid,name="Simu.Y.*")
    mygrid$deleteColumn("Simu.Y.*")

    ## NULL

    if (verbose)
        mygrid$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 19
    ## Maximum Number of UIDs       = 29
    ## Total number of samples      = 5913
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    105.000   479.000
    ## Mesh   :      0.500     0.500
    ## Number :         81        73
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = Neigh.Pb.Number - Locator = NA
    ## Column = 3 - Name = Neigh.Pb.MaxDist - Locator = NA
    ## Column = 4 - Name = Neigh.Pb.MinDist - Locator = NA
    ## Column = 5 - Name = Neigh.Pb.NbNESect - Locator = NA
    ## Column = 6 - Name = Neigh.Pb.NbCESect - Locator = NA
    ## Column = 7 - Name = Kriging.Pb.estim - Locator = NA
    ## Column = 8 - Name = Kriging.Pb.stdev - Locator = NA
    ## Column = 9 - Name = Z.Simu.Y.Pb.1 - Locator = z1
    ## Column = 10 - Name = Z.Simu.Y.Pb.2 - Locator = z2
    ## Column = 11 - Name = Z.Simu.Y.Pb.3 - Locator = z3
    ## Column = 12 - Name = Z.Simu.Y.Pb.4 - Locator = z4
    ## Column = 13 - Name = Z.Simu.Y.Pb.5 - Locator = z5
    ## Column = 14 - Name = Z.Simu.Y.Pb.6 - Locator = z6
    ## Column = 15 - Name = Z.Simu.Y.Pb.7 - Locator = z7
    ## Column = 16 - Name = Z.Simu.Y.Pb.8 - Locator = z8
    ## Column = 17 - Name = Z.Simu.Y.Pb.9 - Locator = z9
    ## Column = 18 - Name = Z.Simu.Y.Pb.10 - Locator = z10

    ## NULL

We calculate some statistics on the Conditional Simulations in Raw
scale.

    if (verbose)
      err = mygrid$statisticsByLocator(locatorType = ELoc_Z(),
                                       opers = EStatOption_fromKeys(c("MEAN","STDV","MINI","MAXI")))

    ##                  MEAN      STDV      MINI      MAXI
    ## Z.Simu.Y.Pb.1      -2.666     0.224    -2.700     0.883
    ## Z.Simu.Y.Pb.2      -2.680     0.158    -2.700     0.453
    ## Z.Simu.Y.Pb.3      -2.664     0.225    -2.700     0.175
    ## Z.Simu.Y.Pb.4      -2.655     0.246    -2.700     0.294
    ## Z.Simu.Y.Pb.5      -2.677     0.172    -2.700     0.137
    ## Z.Simu.Y.Pb.6      -2.649     0.267    -2.700     0.552
    ## Z.Simu.Y.Pb.7      -2.657     0.246    -2.700     0.460
    ## Z.Simu.Y.Pb.8      -2.600     0.390    -2.700     0.985
    ## Z.Simu.Y.Pb.9      -2.652     0.261    -2.700     0.966
    ## Z.Simu.Y.Pb.10     -2.655     0.242    -2.700    -0.381

We visualize a Conditional Simulation in Raw Scale

    if (graphics)
    {
      ax = plot.grid(mygrid,"Z.Simu.Y.Pb.1",end.plot=FALSE)
      plot.point(mydb,"Pb",title="One simulation of Pb in Raw Scale", padd=ax)
    }

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `end.plot`

    ## Warning in geom_point(data = df, mapping = aes(x = x, y = y, color = colval, :
    ## Ignoring unknown parameters: `title`


Let us now average the conditional simulations in order to have a
comparison with the estimation by kriging.

    err = mygrid$statisticsByLocator(locatorType = ELoc_Z(),
                                     opers = EStatOption_fromKeys(c("MEAN")),
                                     flagStoreInDb = TRUE,
                                     verbose = FALSE)
    if (verbose)
        mygrid$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 20
    ## Maximum Number of UIDs       = 30
    ## Total number of samples      = 5913
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    105.000   479.000
    ## Mesh   :      0.500     0.500
    ## Number :         81        73
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = Neigh.Pb.Number - Locator = NA
    ## Column = 3 - Name = Neigh.Pb.MaxDist - Locator = NA
    ## Column = 4 - Name = Neigh.Pb.MinDist - Locator = NA
    ## Column = 5 - Name = Neigh.Pb.NbNESect - Locator = NA
    ## Column = 6 - Name = Neigh.Pb.NbCESect - Locator = NA
    ## Column = 7 - Name = Kriging.Pb.estim - Locator = NA
    ## Column = 8 - Name = Kriging.Pb.stdev - Locator = NA
    ## Column = 9 - Name = Z.Simu.Y.Pb.1 - Locator = NA
    ## Column = 10 - Name = Z.Simu.Y.Pb.2 - Locator = NA
    ## Column = 11 - Name = Z.Simu.Y.Pb.3 - Locator = NA
    ## Column = 12 - Name = Z.Simu.Y.Pb.4 - Locator = NA
    ## Column = 13 - Name = Z.Simu.Y.Pb.5 - Locator = NA
    ## Column = 14 - Name = Z.Simu.Y.Pb.6 - Locator = NA
    ## Column = 15 - Name = Z.Simu.Y.Pb.7 - Locator = NA
    ## Column = 16 - Name = Z.Simu.Y.Pb.8 - Locator = NA
    ## Column = 17 - Name = Z.Simu.Y.Pb.9 - Locator = NA
    ## Column = 18 - Name = Z.Simu.Y.Pb.10 - Locator = NA
    ## Column = 19 - Name = Stats.MEAN - Locator = z1

    ## NULL

Displaying the average of the Conditional Simulations

    if (graphics)
    {
      ax = plot.grid(mygrid,"Stats*MEAN", end.plot=FALSE)
      plot.point(mydb,"Pb",title="Mean of Pb simulations",padd=ax)
    }

    ## Warning in geom_tile(data = df, mapping = aes(x = x, y = y, fill = data), :
    ## Ignoring unknown parameters: `end.plot`

    ## Warning in geom_point(data = df, mapping = aes(x = x, y = y, color = colval, :
    ## Ignoring unknown parameters: `title`


Multivariate case
=================

The Gaussian transform of the Pb variable has already been calculated.
It suffices to perform the Gaussian transform of the Zn variable.

    mydb$setLocator("Zn",ELoc_Z())

    ## NULL

    myanamZn = AnamHermite(nbpoly=30)
    myanamZn$fit(mydb, "Zn")

    ## [1] 0

    if (verbose)
        myanamZn

    ## 
    ## Hermitian Anamorphosis
    ## ----------------------
    ## Minimum absolute value for Y  = -2.5
    ## Maximum absolute value for Y  = 2.6
    ## Minimum absolute value for Z  = 1.1469
    ## Maximum absolute value for Z  = 12.1276
    ## Minimum practical value for Y = -2.5
    ## Maximum practical value for Y = 2.6
    ## Minimum practical value for Z = 1.1469
    ## Maximum practical value for Z = 12.1276
    ## Mean                          = 2.88061
    ## Variance                      = 2.76263
    ## Number of Hermite polynomials = 30
    ## Normalized coefficients for Hermite polynomials (punctual variable)
    ##                [,  0]    [,  1]    [,  2]    [,  3]    [,  4]    [,  5]    [,  6]
    ##      [  0,]     2.881    -1.277     0.877    -0.447    -0.095     0.294    -0.121
    ##      [  7,]    -0.087     0.134    -0.029    -0.087     0.069     0.034    -0.065
    ##      [ 14,]     0.005     0.044    -0.026    -0.020     0.034     0.001    -0.033
    ##      [ 21,]     0.010     0.027    -0.016    -0.019     0.016     0.012    -0.014
    ##      [ 28,]    -0.005     0.011

    if (graphics)
    {
      p = plot.anam(myanamZn)
      plot.decoration(p, title="Gaussian Anamorphosis for Zn")
    }


We convert the raw data into its Gaussian equivalent

    mydb$setLocator("Zn",ELoc_Z())

    ## NULL

    err = myanamZn$rawToGaussianByLocator(mydb)
    if (verbose)
        mydb$display()

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 10
    ## Maximum Number of UIDs       = 20
    ## Total number of samples      = 102
    ## Number of active samples     = 99
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = X - Locator = x1
    ## Column = 2 - Name = Y - Locator = x2
    ## Column = 3 - Name = Zn - Locator = NA
    ## Column = 4 - Name = Pb - Locator = NA
    ## Column = 5 - Name = sel - Locator = sel
    ## Column = 6 - Name = Xvalid.Pb.esterr - Locator = NA
    ## Column = 7 - Name = Xvalid.Pb.stderr - Locator = NA
    ## Column = 8 - Name = Y.Pb - Locator = NA
    ## Column = 9 - Name = Y.Zn - Locator = z1

    ## NULL

We now perform the multivariate variogram caculation

    mydb$setLocator("Y.Pb",ELoc_Z(),0)

    ## NULL

    mydb$setLocator("Y.Zn",ELoc_Z(),1)

    ## NULL

    myvario = Vario(myvarioParam,mydb)
    err = myvario$compute(ECalcVario_VARIOGRAM())
    mymodelM = Model_createFromDb(mydb)
    err = mymodelM$fit(myvario,ECov_fromKeys(c("EXPONENTIAL")))
    if (graphics)
    {
      p = plot.varmod(myvario,mymodelM)
      plot.decoration(p, title="Multivariate Model")
    }


We perform 10 bivariate conditional simulations (deleting the previous
monovariate simulation outcomes first for better legibility).

    mygrid$deleteColumn("Z.Simu*")

    ## NULL

    err = simtub(mydb,mygrid,mymodelM,myneigh,nbsimu=10)
    if (verbose)
        mygrid$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 30
    ## Maximum Number of UIDs       = 50
    ## Total number of samples      = 5913
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    105.000   479.000
    ## Mesh   :      0.500     0.500
    ## Number :         81        73
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x1 - Locator = x1
    ## Column = 1 - Name = x2 - Locator = x2
    ## Column = 2 - Name = Neigh.Pb.Number - Locator = NA
    ## Column = 3 - Name = Neigh.Pb.MaxDist - Locator = NA
    ## Column = 4 - Name = Neigh.Pb.MinDist - Locator = NA
    ## Column = 5 - Name = Neigh.Pb.NbNESect - Locator = NA
    ## Column = 6 - Name = Neigh.Pb.NbCESect - Locator = NA
    ## Column = 7 - Name = Kriging.Pb.estim - Locator = NA
    ## Column = 8 - Name = Kriging.Pb.stdev - Locator = NA
    ## Column = 9 - Name = Stats.MEAN - Locator = NA
    ## Column = 10 - Name = Simu.Y.Pb.1 - Locator = z1
    ## Column = 11 - Name = Simu.Y.Pb.2 - Locator = z2
    ## Column = 12 - Name = Simu.Y.Pb.3 - Locator = z3
    ## Column = 13 - Name = Simu.Y.Pb.4 - Locator = z4
    ## Column = 14 - Name = Simu.Y.Pb.5 - Locator = z5
    ## Column = 15 - Name = Simu.Y.Pb.6 - Locator = z6
    ## Column = 16 - Name = Simu.Y.Pb.7 - Locator = z7
    ## Column = 17 - Name = Simu.Y.Pb.8 - Locator = z8
    ## Column = 18 - Name = Simu.Y.Pb.9 - Locator = z9
    ## Column = 19 - Name = Simu.Y.Pb.10 - Locator = z10
    ## Column = 20 - Name = Simu.Y.Zn.1 - Locator = z11
    ## Column = 21 - Name = Simu.Y.Zn.2 - Locator = z12
    ## Column = 22 - Name = Simu.Y.Zn.3 - Locator = z13
    ## Column = 23 - Name = Simu.Y.Zn.4 - Locator = z14
    ## Column = 24 - Name = Simu.Y.Zn.5 - Locator = z15
    ## Column = 25 - Name = Simu.Y.Zn.6 - Locator = z16
    ## Column = 26 - Name = Simu.Y.Zn.7 - Locator = z17
    ## Column = 27 - Name = Simu.Y.Zn.8 - Locator = z18
    ## Column = 28 - Name = Simu.Y.Zn.9 - Locator = z19
    ## Column = 29 - Name = Simu.Y.Zn.10 - Locator = z20

    ## NULL

We back-transform each set of simulation outcomes using its own Gaussian
Anamorphosis function. Finally we delete the Gaussian variables and ask
for the statistics on the simulated variables in the Raw Scale.

    err = myanamZn$gaussianToRaw(mygrid,"Simu.Y.Zn*")
    err = myanamPb$gaussianToRaw(mygrid,"Simu.Y.Pb*")
    mygrid$deleteColumn("Simu.Y*")

    ## NULL

    if (verbose)
        err = mygrid$statisticsByLocator(locatorType = ELoc_Z(),
                          opers = EStatOption_fromKeys(c("MEAN","STDV","MINI","MAXI")))

    ##                  MEAN      STDV      MINI      MAXI
    ## Z.Simu.Y.Pb.1      -2.697     0.060    -2.700    -1.221
    ## Z.Simu.Y.Pb.2      -2.699     0.029    -2.700    -1.525
    ## Z.Simu.Y.Pb.3      -2.699     0.031    -2.700    -1.476
    ## Z.Simu.Y.Pb.4      -2.698     0.035    -2.700    -1.626
    ## Z.Simu.Y.Pb.5      -2.699     0.029    -2.700    -1.443
    ## Z.Simu.Y.Pb.6      -2.699     0.036    -2.700    -1.164
    ## Z.Simu.Y.Pb.7      -2.699     0.037    -2.700    -1.488
    ## Z.Simu.Y.Pb.8      -2.698     0.043    -2.700    -1.191
    ## Z.Simu.Y.Pb.9      -2.696     0.065    -2.700    -0.621
    ## Z.Simu.Y.Pb.10     -2.698     0.047    -2.700    -1.384

Categorical Variable¶
=====================

We compare the initial variable 'Pb' with a set of disjoint intervals.
The 'Pb' values varying from 3 to 12.7, we consider three classes:

-   values below 4
-   values between 4 and 6
-   values above 6

We first build the indicators for each class.

    limits = Limits(c(NA, 4., 6., NA))
    if (verbose)
        limits$display()

    ## Bound( 1 ) : ] -Inf ; 4 [
    ## Bound( 2 ) : [ 4 ; 6 [
    ## Bound( 3 ) : [ 6 ;  +Inf [

    ## NULL

We apply the set of limits previously defined in order to transform the
input variable into Indicators of the different classes.

    err = limits$toIndicator(mydb,name="Pb")
    if (verbose)
        mydb$display()

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 13
    ## Maximum Number of UIDs       = 43
    ## Total number of samples      = 102
    ## Number of active samples     = 99
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = X - Locator = x1
    ## Column = 2 - Name = Y - Locator = x2
    ## Column = 3 - Name = Zn - Locator = NA
    ## Column = 4 - Name = Pb - Locator = NA
    ## Column = 5 - Name = sel - Locator = sel
    ## Column = 6 - Name = Xvalid.Pb.esterr - Locator = NA
    ## Column = 7 - Name = Xvalid.Pb.stderr - Locator = NA
    ## Column = 8 - Name = Y.Pb - Locator = NA
    ## Column = 9 - Name = Y.Zn - Locator = NA
    ## Column = 10 - Name = Indicator.Pb.Class.1 - Locator = z1
    ## Column = 11 - Name = Indicator.Pb.Class.2 - Locator = z2
    ## Column = 12 - Name = Indicator.Pb.Class.3 - Locator = z3

    ## NULL

We calculate the variogram of the Indicators for future use

    myvarioindParam = VarioParam()
    myvarioindParam$addDir(mydir)

    ## NULL

    myvarioInd = Vario(myvarioindParam,mydb)
    err = myvarioInd$compute(ECalcVario_VARIOGRAM())
    if (verbose)
        myvarioInd$display()

    ## 
    ## Variogram characteristics
    ## =========================
    ## Number of variable(s)       = 3
    ## Number of direction(s)      = 1
    ## Space dimension             = 2
    ## Variance-Covariance Matrix
    ##                [,  0]    [,  1]    [,  2]
    ##      [  0,]     0.107    -0.062    -0.044
    ##      [  1,]    -0.062     0.250    -0.187
    ##      [  2,]    -0.044    -0.187     0.231
    ## 
    ## Direction #1
    ## ------------
    ## Number of lags              = 10
    ## Direction coefficients      =      1.000     0.000
    ## Direction angles (degrees)  =      0.000     0.000
    ## Tolerance on direction      =     90.000 (degrees)
    ## Calculation lag             =      1.000
    ## Tolerance on distance       =     50.000 (Percent of the lag value)
    ## 
    ## For variable 1
    ##       Rank    Npairs  Distance     Value
    ##          0     3.000     0.389     0.000
    ##          1   123.000     1.081     0.081
    ##          2   183.000     2.038     0.126
    ##          3   205.000     3.006     0.156
    ##          4   231.000     4.013     0.132
    ##          5   229.000     5.036     0.159
    ##          6   198.000     5.962     0.152
    ##          7   187.000     7.000     0.107
    ##          8   204.000     7.996     0.096
    ##          9   184.000     8.990     0.068
    ## 
    ## For variables 2 and 1
    ##       Rank    Npairs  Distance     Value
    ##          0     3.000     0.389     0.000
    ##          1   123.000     1.081    -0.065
    ##          2   183.000     2.038    -0.077
    ##          3   205.000     3.006    -0.085
    ##          4   231.000     4.013    -0.093
    ##          5   229.000     5.036    -0.085
    ##          6   198.000     5.962    -0.061
    ##          7   187.000     7.000    -0.045
    ##          8   204.000     7.996    -0.042
    ##          9   184.000     8.990    -0.038
    ## 
    ## For variable 2
    ##       Rank    Npairs  Distance     Value
    ##          0     3.000     0.389     0.167
    ##          1   123.000     1.081     0.199
    ##          2   183.000     2.038     0.221
    ##          3   205.000     3.006     0.251
    ##          4   231.000     4.013     0.292
    ##          5   229.000     5.036     0.258
    ##          6   198.000     5.962     0.237
    ##          7   187.000     7.000     0.254
    ##          8   204.000     7.996     0.228
    ##          9   184.000     8.990     0.234
    ## 
    ## For variables 3 and 1
    ##       Rank    Npairs  Distance     Value
    ##          0     3.000     0.389     0.000
    ##          1   123.000     1.081    -0.016
    ##          2   183.000     2.038    -0.049
    ##          3   205.000     3.006    -0.071
    ##          4   231.000     4.013    -0.039
    ##          5   229.000     5.036    -0.074
    ##          6   198.000     5.962    -0.091
    ##          7   187.000     7.000    -0.061
    ##          8   204.000     7.996    -0.054
    ##          9   184.000     8.990    -0.030
    ## 
    ## For variables 3 and 2
    ##       Rank    Npairs  Distance     Value
    ##          0     3.000     0.389    -0.167
    ##          1   123.000     1.081    -0.134
    ##          2   183.000     2.038    -0.145
    ##          3   205.000     3.006    -0.166
    ##          4   231.000     4.013    -0.199
    ##          5   229.000     5.036    -0.172
    ##          6   198.000     5.962    -0.177
    ##          7   187.000     7.000    -0.209
    ##          8   204.000     7.996    -0.186
    ##          9   184.000     8.990    -0.196
    ## 
    ## For variable 3
    ##       Rank    Npairs  Distance     Value
    ##          0     3.000     0.389     0.167
    ##          1   123.000     1.081     0.150
    ##          2   183.000     2.038     0.194
    ##          3   205.000     3.006     0.237
    ##          4   231.000     4.013     0.238
    ##          5   229.000     5.036     0.247
    ##          6   198.000     5.962     0.268
    ##          7   187.000     7.000     0.270
    ##          8   204.000     7.996     0.240
    ##          9   184.000     8.990     0.226

    ## NULL

    if (graphics)
        plot.varmod(myvarioInd)


Then we build a categorical variable which gives the index of the class
to which each sample belongs

    err = limits$toCategory(mydb,"Pb")
    if (verbose)
        dbfmt = DbStringFormat()
        dbfmt$setFlags(flag_resume = FALSE,
                       flag_vars = FALSE,
                       flag_stats = TRUE)

    ## NULL

        dbfmt$setNames("Category*")

    ## NULL

        dbfmt$setMode(mode=2) # Consider the variable categorical

    ## NULL

        mydb$display(dbfmt)

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Statistics
    ## --------------------
    ## 14 - Name Category.Pb - Locator z1
    ##  Nb of data          =        102
    ##  Nb of active values =         99
    ##  Class         1 =         12 (    12.121%)
    ##  Class         2 =         51 (    51.515%)
    ##  Class         3 =         36 (    36.364%)

    ## NULL
