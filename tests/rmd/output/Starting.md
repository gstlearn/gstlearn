Introduction
============

This document gives the template preliminary statements that need to be
specified in order to use gstlearn package under R.

    library(gstlearn)

    ## 
    ## Attachement du package : 'gstlearn'

    ## Les objets suivants sont masqués depuis 'package:base':
    ## 
    ##     message, toString

    library(ggplot2)

Demonstration
=============

This script: - constructs a regular grid - load values generated
randomly

    # Grid size
    nx = 60
    ny = 30
    mygrid = DbGrid_create(c(nx,ny))

    # Add a uniform random field
    var = VectorHelper_simulateUniform(nx * ny)
    uid = mygrid$addColumns(var, "var1")

    # Display the current contents of the Data Base
    mygrid$display()

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
    ## Total number of samples      = 1800
    ## 
    ## Grid characteristics:
    ## ---------------------
    ## Origin :      0.000     0.000
    ## Mesh   :      1.000     1.000
    ## Number :         60        30
    ## 
    ## Variables
    ## ---------
    ## Column = 0 - Name = rank - Locator = NA
    ## Column = 1 - Name = x1 - Locator = x1
    ## Column = 2 - Name = x2 - Locator = x2
    ## Column = 3 - Name = var1 - Locator = NA

    ## NULL

The grid is plotted (this uses a function provided in plot.r)

    # Display the field
    plot.grid(mygrid, "var1")

