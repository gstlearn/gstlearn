Standard Meshing based on Irregular data set
============================================

We construct a Meshing Standard based on a set of Data Points

    nech = 40
    extendmin = c(0,0)
    extendmax = c(150,100)
    data = Db_createFromBox(nech,extendmin, extendmax)
    data

    ## 
    ## Data Base Characteristics
    ## =========================
    ## 
    ## Data Base Summary
    ## -----------------
    ## File is organized as a set of isolated points
    ## Space dimension              = 2
    ## Number of Columns            = 3
    ## Maximum Number of UIDs       = 3
    ## Total number of samples      = 40
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x-1 - Locator = x1
    ## Column = 2 - Name = x-2 - Locator = x2

    p = plot.point(data)
    plot.decoration(p, title="Display of Data Set")


Creating the Meshing

    mesh1 = MeshEStandardExt()
    err = mesh1$resetFromDb(data)
    mesh1$display()

    ## 
    ## Standard Meshing
    ## ================
    ## Euclidean Geometry
    ## Space Dimension           = 2
    ## Number of Apices per Mesh = 3
    ## Number of Meshes          = 128
    ## Number of Apices          = 79
    ## 
    ## Bounding Box Extension
    ## ----------------------
    ## Dim #1 - Min:0.648392 - Max:142.336
    ## Dim #2 - Min:3.38503 - Max:97.8622

    ## NULL

    p = plot.mesh(mesh1,flagApex=TRUE)
    plot.decoration(p , title="Standard Meshing")

