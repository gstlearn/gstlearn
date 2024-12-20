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

def getCovarianceOptions():
    keys = gl.ECov.getAllKeys(0)
    names = gl.ECov.getAllDescr(0)
    options = {}
    for k in np.arange(len(names)):
        options[names[k]] = keys[k]
    return options

def WgetOneStructure(ic, ncovmax, distmax, varmax):
    '''
    Returns the widgets for enquiring the parameters for a single Basic structure
    '''

    # Calculate the default values
    typeRef = "Spherical"
    distRef = distmax / (ncovmax + 1)
    varRef  = varmax / ncovmax

    # Instantiate the widgets
    WRange = mo.ui.slider(1, 100, value = (ic+1) * distRef, label="Range")
    WSill  = mo.ui.slider(1,50, value=varRef, label="Sill")
    WType  = mo.ui.dropdown(options=getCovarianceOptions(), 
                            value=typeRef, label="Structure")

    return WRange, WSill, WType

def WgetStructures(ncovmax, distmax, varmax):
    '''
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    '''

    TWRange = [None] * ncovmax
    TWSill = [None] * ncovmax
    TWType = [None] * ncovmax
    for ic in range(ncovmax):
        TWRange[ic], TWSill[ic], TWType[ic] = WgetOneStructure(ic, ncovmax, distmax, varmax)
    WRange = mo.ui.array(TWRange)
    WSill = mo.ui.array(TWSill)
    WType = mo.ui.array(TWType)

    return WRange, WSill, WType

def WStructures(ncovmax, distmax, varmax):

    WAll = WgetStructures(ncovmax, distmax, varmax)

    UI = mo.accordion({"Covariance #"+str(i+1): mo.md(            
            f"""
            {WAll[0][i]} # Range

            {WAll[1][i]} # Sill

            {WAll[2][i]} # Type
            """) for i in range(ncovmax)
        }
    )
    return UI, WAll

def getModelFromW(WAll):
    model = gl.Model()
    ncovmax = len(WAll[0])
    for ic in range(ncovmax):
         model.addCovFromParam(type=gl.ECov.fromKey(WAll[2][ic].selected_key), 
                               range= WAll[0][ic].value, 
                               sill = WAll[1][ic].value)
    return model
