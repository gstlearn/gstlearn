In this preamble, we load the **gstlearn** library and clean the
workspace.

Then the necessary data set is downloaded and named **dat**: the target
variable is **January\_temp**

    fileNF = paste(Sys.getenv('GSTLEARN_DATA'),'Scotland',
      'Scotland_Temperatures.NF',sep="/")
    dat = Db_createFromNF(fileNF)

------------------------------------------------------------------------

------------------------------------------------------------------------

    varioParamOmni = VarioParam_createOmniDirection(npas = 100)
    grid.cloud = db_variogram_cloud(dat, varioParamOmni)
    grid.cloud$display()

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 4
    ## Maximum Number of UIDs       = 4
    ## Total number of samples      = 10000
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      7.789     0.068
    ## Number :        100       100
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2
    ## Column = 3 - Name = Cloud.January_temp - Locator = NA

    ## NULL

    p = plot.grid(grid.cloud, "Cloud.January*")
    p = plot.geometry(p, asp=0)
    p


------------------------------------------------------------------------

We calculate the omni-directional variogram of the temperatures.

    varioParamOmni = VarioParam_createOmniDirection(npas=40, dpas=10)
    varioexp = Vario(varioParamOmni, dat)
    err = varioexp$compute()

Print the variogram contents

    varioexp

------------------------------------------------------------------------

Plot the omni-directional variogram

    plot.varmod(varioexp)


------------------------------------------------------------------------

    plot.varmod(varioexp,draw_psize=TRUE)


------------------------------------------------------------------------

    fitmod = Model()
    err = fitmod$fit(varioexp)
    plot.varmod(varioexp, fitmod)


------------------------------------------------------------------------

    fitmod

    ## 
    ## Model characteristics
    ## =====================
    ## Space dimension              = 2
    ## Number of variable(s)        = 1
    ## Number of basic structure(s) = 1
    ## Number of drift function(s)  = 0
    ## Number of drift equation(s)  = 0
    ## 
    ## Covariance Part
    ## ---------------
    ## Spherical
    ## - Sill         =      1.156
    ## - Range        =    135.130
    ## Total Sill     =      1.156

------------------------------------------------------------------------

    ECov_printAll()

    ##   -2 -     UNKNOWN : Unknown covariance
    ##   -1 -    FUNCTION : External covariance function
    ##    0 -      NUGGET : Nugget effect
    ##    1 - EXPONENTIAL : Exponential
    ##    2 -   SPHERICAL : Spherical
    ##    3 -    GAUSSIAN : Gaussian
    ##    4 -       CUBIC : Cubic
    ##    5 -     SINCARD : Sine Cardinal
    ##    6 -    BESSEL_J : Bessel J
    ##    7 -    BESSEL_K : Bessel K
    ##    8 -       GAMMA : Gamma
    ##    9 -      CAUCHY : Cauchy
    ##   10 -      STABLE : Stable
    ##   11 -      LINEAR : Linear
    ##   12 -       POWER : Power
    ##   13 -   ORDER1_GC : First Order Generalized covariance
    ##   14 -   SPLINE_GC : Spline Generalized covariance
    ##   15 -   ORDER3_GC : Third Order Generalized covariance
    ##   16 -   ORDER5_GC : Fifth Order Generalized covariance
    ##   17 -     COSINUS : Cosine
    ##   18 -    TRIANGLE : Triangle
    ##   19 -      COSEXP : Cosine Exponential
    ##   20 -       REG1D : 1-D Regular
    ##   21 -       PENTA : Pentamodel
    ##   22 -  SPLINE2_GC : Order-2 Spline
    ##   23 -     STORKEY : Storkey covariance in 1-D
    ##   24 -   WENDLAND0 : Wendland covariance (2,0)
    ##   25 -   WENDLAND1 : Wendland covariance (3,1)
    ##   26 -   WENDLAND2 : Wendland covariance (4,2)
    ##   27 -      MARKOV : Markovian covariances

    ## NULL

------------------------------------------------------------------------

    types = ECov_fromKeys(c("NUGGET","CUBIC","SPHERICAL"))
    err = fitmod$fit(varioexp, types=types)
    plot.varmod(varioexp, fitmod)


------------------------------------------------------------------------

    fitmod

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
    ## Cubic
    ## - Sill         =      0.414
    ## - Range        =     76.040
    ## Spherical
    ## - Sill         =      0.892
    ## - Range        =    241.113
    ## Total Sill     =      1.306

------------------------------------------------------------------------

    constraints = Constraints()
    err = constraints$addItemFromParamId(EConsElem_RANGE(),icov=1,type=EConsType_UPPER(),value=20.)
    err = constraints$addItemFromParamId(EConsElem_SILL(),icov=1,type=EConsType_LOWER(),value=0.03)
    err = fitmod$fit(varioexp, types=types, constraints, Option_VarioFit(TRUE))
    plot.varmod(varioexp, fitmod)


------------------------------------------------------------------------

    fitmod

    ## 
    ## Model characteristics
    ## =====================
    ## Space dimension              = 2
    ## Number of variable(s)        = 1
    ## Number of basic structure(s) = 3
    ## Number of drift function(s)  = 0
    ## Number of drift equation(s)  = 0
    ## 
    ## Covariance Part
    ## ---------------
    ## Nugget Effect
    ## - Sill         =      0.000
    ## Cubic
    ## - Sill         =      0.109
    ## - Range        =     20.000
    ## Spherical
    ## - Sill         =      1.056
    ## - Range        =    155.566
    ## Total Sill     =      1.166

------------------------------------------------------------------------

    constraints = Constraints()
    err = constraints$addItemFromParamId(EConsElem_RANGE(),icov=1,type=EConsType_EQUAL(),value=1000.)
    err = constraints$addItemFromParamId(EConsElem_SILL(),icov=1,type=EConsType_EQUAL(),value=0.4)
    err = fitmod$fit(varioexp, types=types, constraints, Option_VarioFit(flag_noreduce=TRUE))
    plot.varmod(varioexp, fitmod)


------------------------------------------------------------------------

    fitmod

    ## 
    ## Model characteristics
    ## =====================
    ## Space dimension              = 2
    ## Number of variable(s)        = 1
    ## Number of basic structure(s) = 3
    ## Number of drift function(s)  = 0
    ## Number of drift equation(s)  = 0
    ## 
    ## Covariance Part
    ## ---------------
    ## Nugget Effect
    ## - Sill         =      0.053
    ## Cubic
    ## - Sill         =      0.400
    ## - Range        =   1000.000
    ## Spherical
    ## - Sill         =      1.003
    ## - Range        =    130.042
    ## Total Sill     =      1.457

------------------------------------------------------------------------

    varioParamMulti = VarioParam_createMultiple(ndir=4, npas=15, dpas=15.)
    vario.4dir = Vario(varioParamMulti, dat)
    err = vario.4dir$compute()
    plot.varmod(vario.4dir)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


------------------------------------------------------------------------

    model.4dir = Model()
    err = model.4dir$fit(vario.4dir,types=types)
    plot.varmod(vario.4dir, model.4dir)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.


------------------------------------------------------------------------

    grid.vmap = db_vmap_compute(dat, ECalcVario_VARIOGRAM())
    plot.grid(grid.vmap)


------------------------------------------------------------------------

    modelVM = Model()
    err = modelVM$fitFromVMap(grid.vmap, types=types)
    modelVM

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
    ## Nugget Effect
    ## - Sill         =      0.251
    ## Cubic
    ## - Sill         =      0.949
    ## - Ranges       =    154.795   215.122
    ## - Angles       =    -24.923     0.000
    ## - Rotation Matrix
    ##                [,  0]    [,  1]
    ##      [  0,]     0.907     0.421
    ##      [  1,]    -0.421     0.907
    ## Total Sill     =      1.200

------------------------------------------------------------------------

    err = dbgrid_model(grid.vmap, modelVM)
    plot.grid(grid.vmap)


------------------------------------------------------------------------

    plot.varmod(vario.4dir, modelVM)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

