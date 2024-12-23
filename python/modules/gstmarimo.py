################################################################################
#                                                                              #
#                         gstlearn Python package                              #
#                                                                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://gstlearn.org                                                #
# License: BSD 3-clause                                                        #
#                                                                              #
################################################################################
# This part is meant to distribute (as is) a set of functions written in Python
# to be included in Marimo interface
# Reminder: all methods staring by "W" are dedicated to UI

import gstlearn as gl
import gstlearn.plot as gp

import numpy as np
import matplotlib.pyplot as plt
import copy

import marimo as mo

def getCovarianceDict():
    keys = gl.ECov.getAllKeys(0)
    names = gl.ECov.getAllDescr(0)
    options = {}
    for k in np.arange(len(names)):
        options[names[k]] = keys[k]
    return options

def WCovariance(ic = 0, ncovmax = 1, distmax = 100, varmax = 100):
    '''
    Returns the widget for inquiring the parameters for a single Basic structure
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    distmax: Maximum distance
    varmax: Maximum Variance value
    '''
    typeRef = "Spherical"
    distRef = distmax * (ic+1) / (ncovmax + 1)
    varRef  = varmax / ncovmax

    WUsed   = mo.ui.switch(True, label="Basic Structure Used")
    WType   = mo.ui.dropdown(options=getCovarianceDict(), 
                            value=typeRef, label="Structure")
    WRange  = mo.ui.slider(1, distmax, value = distRef, label="Range")
    WSill   = mo.ui.slider(0, varmax, value=varRef, label="Sill")
    WAniso  = mo.ui.switch(label="Anisotropy")
    WRange2 = mo.ui.slider(1, distmax, value = distRef, label="Range Aux.")
    WAngle  = mo.ui.slider(0, 180, value = 0, label="Angle")

    return mo.ui.array([WUsed, WType, WRange, WSill, WAniso, WRange2, WAngle])

def WModel(ncovmax=1, distmax=100, varmax=100):
    '''
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    '''
    return mo.ui.array([WCovariance(ic, ncovmax, distmax, varmax) 
                         for ic in range(ncovmax)])

def getWModel(WAlls):
    '''
    Create a gstlearn Model from the WModel widget
    '''
    model = gl.Model()
    for WAll in WAlls:
        used   = WAll.value[0]
        type   = gl.ECov.fromKey(WAll.value[1])
        range  = WAll.value[2]
        sill   = WAll.value[3]
        aniso  = WAll.value[4]
        range2 = WAll.value[5]
        angle  = WAll.value[6]
        if used:
            if aniso:
                model.addCovFromParam(type=type, sill=sill,
                                      angles = [angle, 0],
                                      ranges = [range, range2])
            else:
                model.addCovFromParam(type=type, range=range, sill=sill)
    return model

def WGrid(nxdef = 50):
    '''
    Widget to inquire the parameters for constructing a Grid
    '''
    WNX = mo.ui.slider(start=1, stop=200, value = nxdef, label="NX")
    WNY = mo.ui.slider(start=1, stop=200, value = nxdef, label="NY")
    WDX = mo.ui.number(start=1, stop=None, value = 1,  label="DX")
    WDY = mo.ui.number(start=1, stop=None, value = 1,  label="DY")
    WX0 = mo.ui.number(start=0, stop=None, value = 0,  label="X0")
    WY0 = mo.ui.number(start=0, stop=None, value = 0,  label="Y0")

    return mo.ui.array([WNX, WNY, WDX, WDY, WX0, WY0])

def getWGrid(WAll):
    '''
    Create the gstlearn Grid from the widget WGrid
    '''
    grid = gl.DbGrid.create(nx = [WAll[0].value, WAll[1].value],
                            dx = [WAll[2].value, WAll[3].value],
                            x0 = [WAll[4].value, WAll[5].value])
    return grid

def WSimtub(seed = 13134):
    '''
    Inquire for performing a Simulation using Turning Bands Method
    '''

    WNbtuba = mo.ui.number(start=1, stop=None, value = 100, 
                          label = "Number of Turning Bands")
    WDefineSeed = mo.ui.switch(True, label="Define Seed")
    WSeed = mo.ui.number(start=0, stop=None, value = seed, 
                         label = "Seed")
    WNewSimu = mo.ui.button(label="Perform a New Simulation")

    return mo.ui.array([WNbtuba, WDefineSeed, WSeed, WNewSimu])

def getWSimtub(WAll):
    '''
    Returns the parameters for simulation using Turning Bands Method
    '''
    nbtuba = WAll[0].value
    defineSeed = WAll[1].value
    seed = 0
    if defineSeed:
        seed = WAll[2].value
    return nbtuba, seed

def WDbRandom(nech = 100, nvar = 1, xmin = 0, ymin = 0, xmax = 100, ymax = 100, seed = 14543):
    '''
    Inquiry the parameters for generating a Random Data Base
    '''
    WNech = mo.ui.number(start=1, stop=None, value = nech, 
                         label = "Number of Samples")
    WNvar = mo.ui.number(start=1, stop=None, value = nvar, 
                         label = "Number of Variables")
    WXmin = mo.ui.number(start=None, stop=None, value = xmin, 
                         label = "Minimum along X")
    WYmin = mo.ui.number(start=None, stop=None, value = ymin, 
                         label = "Minimum along Y")
    WXmax = mo.ui.number(start=None, stop=None, value = xmax, 
                         label = "Maximum along X")
    WYmax = mo.ui.number(start=None, stop=None, value = ymax, 
                         label = "Maximum along Y")
    WSeed = mo.ui.number(start=None, stop=None, value=seed,
                         label="Seed")
    
    return mo.ui.array([WNech, WNvar, WXmin, WYmin, WXmax, WYmax, WSeed])

def getWDbRandom(WAll):
    '''
    Create a Db with Radom Samples 
    '''
    coormin = [WAll[2].value, WAll[3].value]
    coormax = [WAll[4].value, WAll[5].value]
    db = gl.Db.createFillRandom(ndat = WAll[0].value, 
                                ndim = 2,  
                                nvar = WAll[1].value, 
                                nfex = 0,
                                ncode = 0,
                                varmax = 0.,
                                selRatio = 0.,
                                heteroRatio = gl.VectorDouble(),
                                coormin = coormin,
                                coormax = coormax,
                                seed = WAll[6].value)
    return db

def WDbNFFile():
    '''
    Inquiry to load a file from a Neutral File
    '''
    WFile = mo.ui.file_browser(label="Select a Db Neutral File")

    return mo.ui.array([WFile])

def getWDbNFFile(WAll):
    '''
    Create the gstlearn Db from the widget WDbNFFile
    '''
    filename = WAll[0].name()
    if filename is None:
        return None
    db = gl.Db.createFromNF(filename)

    return db

def WDbGridRandom(nxdef = 10):
    '''
    Widget to inquire the parameters for constructing a Db as a randomized grid
    '''
    WNX = mo.ui.number(start=1, stop=100, value = nxdef, label="NX")
    WNY = mo.ui.number(start=1, stop=100, value = nxdef, label="NY")
    WDX = mo.ui.number(start=1, stop=None, value = 1,  label="DX")
    WDY = mo.ui.number(start=1, stop=None, value = 1,  label="DY")
    WX0 = mo.ui.number(start=0, stop=None, value = 0,  label="X0")
    WY0 = mo.ui.number(start=0, stop=None, value = 0,  label="Y0")
    WPerc = mo.ui.number(start=0, stop=100, value=10, 
                         label="Randomization percentage")

    return mo.ui.array([WNX, WNY, WDX, WDY, WX0, WY0, WPerc])

def getWDbGridRandom(WAll):
    '''
    Create the gstlearn Grid from the widget WDbRandomGrid
    '''
    grid = gl.DbGrid.create(nx = [WAll[0].value, WAll[1].value],
                            dx = [WAll[2].value, WAll[3].value],
                            x0 = [WAll[4].value, WAll[5].value])
    db = gl.Db.createFromGridRandomized(grid, WAll[6].value);
    return db
