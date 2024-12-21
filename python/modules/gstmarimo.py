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
    Returns the widgets for enquiring the parameters for a single Basic structure
    '''
    typeRef = "Spherical"
    distRef = distmax / (ncovmax + 1)
    varRef  = varmax / ncovmax

    WRange = mo.ui.slider(1, 100, value = (ic+1) * distRef, label="Range")
    WSill  = mo.ui.slider(1, 100, value=varRef, label="Sill")
    WType  = mo.ui.dropdown(options=getCovarianceDict(), 
                            value=typeRef, label="Structure")

    return WRange, WSill, WType

def WCovariances(ncovmax, distmax, varmax):
    '''
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    '''

    TWRanges = [None] * ncovmax
    TWSills = [None] * ncovmax
    TWTypes = [None] * ncovmax
    for ic in range(ncovmax):
        TWRanges[ic], TWSills[ic], TWTypes[ic] = WCovariance(ic, ncovmax, distmax, varmax)
    WRanges = mo.ui.array(TWRanges)
    WSills  = mo.ui.array(TWSills)
    WTypes  = mo.ui.array(TWTypes)

    return WRanges, WSills, WTypes

def WModel(ncovmax, distmax, varmax):

    WAll = WCovariances(ncovmax, distmax, varmax)

    UI = mo.accordion({"Covariance #"+str(i+1): mo.md(            
            f"""
            {WAll[0][i]} # Range

            {WAll[1][i]} # Sill

            {WAll[2][i]} # Type
            """) for i in range(ncovmax)
        }
    )
    return UI, WAll

def getWModel(WAll):
    model = gl.Model()
    ncovmax = len(WAll[0])
    for ic in range(ncovmax):
         model.addCovFromParam(type=gl.ECov.fromKey(WAll[2][ic].selected_key), 
                               range= WAll[0][ic].value, 
                               sill = WAll[1][ic].value)
    return model

def WGrid():
    WNX = mo.ui.slider(1, 100, value = 10, label="NX")
    WNY = mo.ui.slider(1, 100, value = 10, label="NY")
    WDX = mo.ui.number(1, 100, value = 1,  label="DX")
    WDY = mo.ui.number(1, 100, value = 1,  label="DY")
    WX0 = mo.ui.number(0, 100, value = 0,  label="X0")
    WY0 = mo.ui.number(0, 100, value = 0,  label="Y0")

    UI = mo.accordion({"Grid Parameters": mo.md(            
            f"""
            {WNX} {WNY}

            {WDX} {WDY}

            {WX0} {WY0}
            """)})
    WAll = mo.ui.array([WNX, WNY, WDX, WDY, WX0, WY0])
    return UI, WAll

def getWGrid(WAll):
    print(WAll[0].value, WAll[1].value)
    grid = gl.DbGrid.create(nx = [WAll[0].value, WAll[1].value],
                            dx = [WAll[2].value, WAll[3].value],
                            x0 = [WAll[4].value, WAll[5].value])
    return grid
