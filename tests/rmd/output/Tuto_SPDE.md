General Introduction
====================

This notebook presents sevral possibilities of the SPDE procedure

Creating the Data Bases
-----------------------

Some global constants

    ndat = 1000
    rangev = 0.2
    sill = 1.
    nugget = 0.1
    law_set_random_seed(123)

    ## NULL

Creating the Data Base for conditioning information

    dat = Db_create()
    dat["x"] = VectorHelper_simulateUniform(ndat)
    dat["y"] = VectorHelper_simulateUniform(ndat)
    dat$setLocators(c("x","y"),ELoc_X())

    ## NULL

    dat

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 2
    ## Maximum Number of UIDs       = 2
    ## Total number of samples      = 1000
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = x - Locator = x1
    ## Column = 1 - Name = y - Locator = x2

Creating the Data Base for the output grid

    grid = DbGrid_create(nx=c(50,50),dx=c(0.02,0.02))

Creating the Meshing (Turbo) based on extended grid

    gridExt = DbGrid_create(nx=c(75,75),dx=c(0.02,0.02),x0=c(-0.25,-0.25))
    mesh = MeshETurbo(gridExt)
    mesh

    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :     -0.250    -0.250
    ## Mesh   :      0.020     0.020
    ## Number :         75        75
    ## 
    ## Turbo Meshing
    ## =============
    ## Euclidean Geometry
    ## Space Dimension           = 2
    ## Number of Apices per Mesh = 3
    ## Number of Meshes          = 10952
    ## Number of Apices          = 5625
    ## 
    ## Bounding Box Extension
    ## ----------------------
    ## Dim #1 - Min:-0.25 - Max:1.23
    ## Dim #2 - Min:-0.25 - Max:1.23

Creating the Model

    model = Model_createFromParam(type=ECov_BESSEL_K(),param=1,range=rangev,sill=sill)
    model$addCovFromParam(ECov_NUGGET(),sill=nugget)

    ## NULL

    model

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
    ## K-Bessel (Third Parameter = 1)
    ## - Sill         =      1.000
    ## - Range        =      0.200
    ## - Theo. Range  =      0.058
    ## Nugget Effect
    ## - Sill         =      0.100
    ## Total Sill     =      1.100

SPDE processing
===============

Non-conditional simulation
--------------------------

Building the SPDE environment for Non-conditional environments

    spdeS = SPDE()
    spdeS$init(model=model, field=grid, calc=ESPDECalcMode_SIMUNONCOND(), mesh=mesh)

    ## NULL

    spdeS$compute()

    ## NULL

Apply the SPDE system in order to obtain a non-conditional simulation on
the input and output Data Bases

    iuid = spdeS$query(dat)
    iuid = spdeS$query(grid)

Representing the resulting non-conditional simulation on the Grid

    p = plot.grid(grid, "spde.simu")
    p = plot.decoration(p, title="Non Conditional Simulation")
    p


Estimation by Kriging
---------------------

Preparing the SPDE environment

    spdeK = SPDE()
    spdeK$init(model=model, field=grid, data=dat, calc=ESPDECalcMode_KRIGING(), mesh=mesh)

    ## NULL

    spdeK$compute()

    ## NULL

Apply the SPDE system in order to obtain an estimation on the output
Data Base

    iuid = spdeK$query(grid, namconv=NamingConvention("spde",FALSE))

Representing the resulting estimation on the Grid

    p = plot.grid(grid, "spde.kriging")
    p = plot.decoration(p, title="Estimation")
    p


Producing the internal elements of the SPDE environment
-------------------------------------------------------

Extracting the internal elements from 'spdeK' in order to perform
calculations by hand (if necessary).

    Qtr = csToTriplet(PrecisionOpCs(mesh, model)$getQ(), flag_from_1=TRUE)
    Atr = csToTriplet(ProjMatrix(dat,mesh)$getAproj(), flag_from_1=TRUE)

    Q = sparseMatrix(i=Qtr$rows,j=Qtr$cols,x=Qtr$values)
    Aproj = sparseMatrix(i=Atr$rows, j=Atr$cols, x=Atr$values,
                         dims=c(Atr$nrows,Atr$ncols))

Posterior calculations

    size = dim(Q)[1]
    cholQ = Cholesky(Q)
    u = VectorHelper_simulateGaussian(size)

    Qop = Q + 1/nugget * t(Aproj) %*% Aproj
