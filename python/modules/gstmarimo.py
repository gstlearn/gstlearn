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
    Returns the widget for enquiring the parameters for a single Basic structure
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
    WSeed = mo.ui.number(start=0, stop=None, value = seed, 
                         label = "Seed")
    
    return mo.ui.array([WNbtuba, WSeed])

def getWSimtub(WAll):
    '''
    Returns the parameters for simulation using Turning Bands Method
    '''
    nbtuba = WAll[0].value
    seed = WAll[1].value
    return nbtuba, seed