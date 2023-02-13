We define the Space dimension

    defineDefaultSpace(ESpaceType_RN(), 2)

    ## NULL

Creating a dummy Model used for simulating a random field

    mymodel = Model_createFromParam(ECov_CUBIC(), range=10, sill=1)

Grid representations
====================

Standard Grid
-------------

We construct a standard non rotated 2-D grid with two simulated
variables

    nx = c(70,25)
    dx = c(1,2)
    x0 = c(-40, 20)
    mygrid = DbGrid_create(nx,dx,x0)

    err = simtub(NULL,mygrid,mymodel,nbsimu=2)
    mygrid$display()

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
    ## Total number of samples      = 1750
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    -40.000    20.000
    ## Mesh   :      1.000     2.000
    ## Number :         70        25
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2
    ## Column = 3 - Name = Simu.1 - Locator = z1
    ## Column = 4 - Name = Simu.2 - Locator = z2

    ## NULL

We construct a Selection

    mygrid["sel"] = 1. - 
      (mygrid["x1"] > 0) * (mygrid["x1"] < 15) * (mygrid["x2"] > 40) * (mygrid["x2"] < 50)
    mygrid$setLocator("sel",ELoc_SEL())

    ## NULL

    mygrid

    ## 
    ## Data Base Grid Characteristics
    ## ==============================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a regular grid
    ## Space dimension              = 2
    ## Number of Columns            = 6
    ## Maximum Number of UIDs       = 6
    ## Total number of samples      = 1750
    ## Number of active samples     = 1694
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    -40.000    20.000
    ## Mesh   :      1.000     2.000
    ## Number :         70        25
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2
    ## Column = 3 - Name = Simu.1 - Locator = z1
    ## Column = 4 - Name = Simu.2 - Locator = z2
    ## Column = 5 - Name = sel - Locator = sel

Displaying the grid cells

    p = plot.grid(mygrid, name_raster="Simu.1", name_contour="Simu.2", 
                  show.legend.raster=TRUE, legend.name.raster="ma Legende")
    p = plot.decoration(p,title="Display of Grid Cells")
    p

    ## Warning: Removed 56 rows containing non-finite values (`stat_contour()`).


Displaying the grid nodes only

    p = plot.point(mygrid, name_color="Simu.1",name_size="Simu.2", show.legend.symbol = TRUE,
                   legend.name.color="myColor")
    plot.decoration(p,title="Display of Grid Nodes")


Rotated grid
------------

    mygrid = DbGrid_create(nx,dx,x0,angles=c(10,0))
    err = simtub(NULL,mygrid,mymodel,nbsimu=2)
    mygrid

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
    ## Total number of samples      = 1750
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :    -40.000    20.000
    ## Mesh   :      1.000     2.000
    ## Number :         70        25
    ## Rotation Angles        =     10.000     0.000
    ## Direct Rotation Matrix
    ##                [,  0]    [,  1]
    ##      [  0,]     0.985     0.174
    ##      [  1,]    -0.174     0.985
    ## Inverse Rotation Matrix
    ##                [,  0]    [,  1]
    ##      [  0,]     0.985    -0.174
    ##      [  1,]     0.174     0.985
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2
    ## Column = 3 - Name = Simu.1 - Locator = z1
    ## Column = 4 - Name = Simu.2 - Locator = z2

Displaying the cell contents in the rotated grid

    p = plot.grid(mygrid, name_raster="Simu.1", show.legend.raster=FALSE)
    plot.decoration(p,title="Display of Rotated Grid")


As a set of grid nodes

    p = plot.point(mygrid,name_color="Simu.1", show.legend.symbol = TRUE)
    plot.decoration(p, title="Display of Rotated Grid Nodes")


Points and Polygon
==================

A set of points is sampled from the previous Grid and stored in a new
Point Db. The number of samples if fixed to 1% of the number of grid
nodes.

    mypoint = Db_createSamplingDb(mygrid,0.01)
    mypoint$display()

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
    ## Total number of samples      = 17
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2
    ## Column = 3 - Name = Simu.1 - Locator = z1
    ## Column = 4 - Name = Simu.2 - Locator = z2

    ## NULL

We create a polygon as the convex hull of the samples

    mypoly = Polygons_createFromDb(mypoint)

We now display the points and the polygon on top of the grid: the
overlay is ensured by using the argument 'p'.

    p = plot(mygrid)
    p = plot(mypoly, fill=NA, colour='yellow', linewidth=1, padd=p)
    p = plot(mypoint,color="black", padd=p)
    p = plot.decoration(p, title="mon titre", xlab="mon axe des X", ylab="Mon axe des Y")
    p


Several plots on the same figure
================================

We create two layers containing the Point and the Grid representation.
We then use *ggarrange* to display them side-by-side.

    p1 = plot(mygrid)
    p2 = plot(mypoint)
    ggarrange(p1, p2, labels = c("Plot #1", "Plot #2"), ncol = 2, nrow = 1)


Variograms and Models
=====================

We calculate the variogram along the two main directions of the 2-D
grid. We compute it for the 2 variables currently defined in the Grid
file.

    varioparam = VarioParam_createMultipleFromGrid(npas=10)
    vario = Vario(varioparam, mygrid)
    err = vario$compute(ECalcVario_VARIOGRAM())

We fit a Model (automatic procedure)

    model = Model()
    err = model$fit(vario,type=ECov_fromKeys(c("SPHERICAL", "CUBIC")))

We display the experimental variogram for the first variable in the
first direction.

    p = plot(vario, draw_psize=F, linewidth=1, color="red")

    ## Warning: Duplicated aesthetics after name standardisation: colour

    p = plot.decoration(p, title="First Variable - First Direction")
    p

    ## Warning: Removed 9 rows containing missing values (`geom_line()`).


We display the experimental variogram for the first variable in the
second direction.

    p = plot(vario, idir=1)
    p = plot.decoration(p, title="First Variable - Second Direction")
    p


We simply overlay the contents of the two previous plots.

    p = plot(vario, color="black")

    ## Warning: Duplicated aesthetics after name standardisation: colour

    p = plot(vario, idir=1, color='red', padd=p)

    ## Warning: Duplicated aesthetics after name standardisation: colour

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    p = plot.decoration(p, title="First Variable - Two directions")
    p

    ## Warning: Removed 9 rows containing missing values (`geom_line()`).

    ## Warning: Removed 9 rows containing missing values (`geom_line()`).


We display the cross-variogram between both variables along both
directions.

    p = plot(vario, ivar=1, jvar=0, idir=-1, linewidth = 1)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    p = plot.decoration(p, title="Cross-variogram - Two directions")
    p


We display the simple and cross variograms in both calculation
directions.

    p = plot(vario, ivar=-1, jvar=-1, idir=-1, linewidth=1)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    p = plot.decoration(p,title="Simple and Cross Variograms in all directions")
    p


We display the Model as calculated in the first direction of the
variogram

    p = plot(model, ivar=1, jvar=1, vario=vario, idir=0)
    p = plot.decoration(p, title="Model for the Second Variable in First Direction")
    p


We now represent all the simple and cross variograms together with the
fitted models.

    p = plot.varmod(vario, model, linewidth=1)

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    p = plot.decoration(p, title="All variograms in First Direction")
    p

