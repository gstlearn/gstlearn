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

import marimo as mo

def _getCovarianceDict():
    '''
    Returns the list of covariances available as a Dictionary
    '''
    keys = gl.ECov.getAllKeys(0)
    names = gl.ECov.getAllDescr(0)
    options = {}
    for k in np.arange(len(names)):
        options[names[k]] = keys[k]
    return options


def _WLock(WTest, condition, colorBackground = "white", colorText = "black"):
    '''
    Turns the Widget to grey (as if it was locked if 'condition' is fulfilled)
    '''
    if not condition:
        newWTest = WTest.style({"backgroundColor": colorBackground, 
            "color": colorText
        })
    else:
        newWTest = WTest.style({"backgroundColor": "#f0f0f0",
            "color": "#a0a0a0"
        })
    return newWTest

def WdefineCovariance(ic = 0, ncovmax = 1, distmax = 100, varmax = 100):
    '''
    Returns the widget for inquiring the parameters for a single Basic structure
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    distmax: Maximum distance
    varmax:  Maximum Variance value
    '''
    typeRef = "Spherical"
    distRef = distmax * (ic+1) / (ncovmax + 1)
    varRef  = varmax / ncovmax

    WUsed   = mo.ui.switch(True, label="Basic Structure Used")
    WType   = mo.ui.dropdown(options=_getCovarianceDict(), value=typeRef, label="Structure")
    WRange  = mo.ui.slider(1, distmax, value = distRef, label="Range")
    WSill   = mo.ui.slider(0, varmax, value=varRef, label="Sill")
    WAniso  = mo.ui.switch(label="Anisotropy")
    WRange2 = mo.ui.slider(1, distmax, value = distRef, label="Range Aux.")
    WAngle  = mo.ui.slider(0, 180, value = 0, label="Angle")

    return mo.ui.array([WUsed, WType, WRange, WSill, WAniso, WRange2, WAngle])

def WshowCovariance(WAll):
    '''
    Returns the contents of the Covariance Widget as HTML, ready to be displayed
    '''
    [WUsed, WType, WRange, WSill, WAniso, WRange2, WAngle] = WAll

    WTypeupd   = _WLock(WType, not WUsed.value)
    WRangeupd  = _WLock(WRange, not WUsed.value)
    WSillupd   = _WLock(WSill, not WUsed.value)
    WAnisoupd  = _WLock(WAniso, not WUsed.value)
    WRange2upd = _WLock(WRange2, not WUsed.value or not WAniso.value)
    WAngleupd  = _WLock(WAngle, not WUsed.value or not WAniso.value)

    return mo.vstack([WUsed, WTypeupd, WRangeupd, WSillupd, WAnisoupd, WRange2upd, WAngleupd])

def WgetCovariance(WAll):
    '''
    Uses the contents of the Covariance Widget to produce a CovAniso item
    '''
    [WUsed, WType, WRange, WSill, WAniso, WRange2, WAngle] = WAll

    type = gl.ECov.fromKey(WType.value)
    cova = None
    if WUsed.value:
        ctxt = gl.CovContext(1, 2)
        if not WAniso.value:
            cova = gl.CovAniso.createIsotropic(ctxt, type = type,
                                 range = WRange.value, sill = WSill.value,
                                 param = 1., flagRange= True)
        else:
            cova = gl.CovAniso.createAnisotropic(ctxt, type = type,
                                   ranges = [WRange.value, WRange2.value], sill = WSill.value,
                                   param = 1.,
                                   angles = [WAngle.value, 0.],
                                   flagRange = True)
    return cova

def WdefineModel(ncovmax=1, distmax=100, varmax=100):
    '''
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    '''
    return mo.ui.array([WdefineCovariance(ic, ncovmax, distmax, varmax) 
                         for ic in range(ncovmax)])

def WshowModel(WAlls):
    ncov = len(WAlls)
    UI = mo.accordion({
        "Covariance "+str(ic+1) : WshowCovariance(WAlls[ic]) for ic in range(ncov)
    })
    return UI

def WgetModel(WAlls):
    '''
    Create a gstlearn Model
    '''
    model = gl.Model()
    for WAll in WAlls:
        cova = WgetCovariance(WAll)
        model.addCov(cova)
    return model

def WdefineGrid(nxdef = 50):
    '''
    Widget to inquire the parameters for constructing a Grid
    '''
    WNX = mo.ui.slider(start=1, stop=200, value = nxdef)
    WNY = mo.ui.slider(start=1, stop=200, value = nxdef)
    WDX = mo.ui.number(start=1, stop=None, value = 1)
    WDY = mo.ui.number(start=1, stop=None, value = 1)
    WX0 = mo.ui.number(start=0, stop=None, value = 0)
    WY0 = mo.ui.number(start=0, stop=None, value = 0)

    return mo.ui.array([WNX, WNY, WDX, WDY, WX0, WY0])

def WshowGrid(WAll):
    [WNX, WNY, WDX, WDY, WX0, WY0] = WAll

    return mo.hstack([mo.vstack([mo.md("Grid"), mo.md("Number"), mo.md("Mesh"), mo.md("Origin")]),
                      mo.vstack([mo.md("along X"), WNX, WDX, WX0],align='end'),
                      mo.vstack([mo.md("along Y"), WNY, WDY, WY0],align='end')],justify='start')
    
def WgetGrid(WAll):
    '''
    Create the gstlearn Grid
    '''
    [WNX, WNY, WDX, WDY, WX0, WY0] = WAll
    grid = gl.DbGrid.create(nx = [WNX.value, WNY.value],
                            dx = [WDX.value, WDY.value],
                            x0 = [WX0.value, WY0.value])
    return grid

def WdefineSimtub(seed = 13134):
    '''
    Inquire for performing a Simulation using Turning Bands Method
    '''
    WNbtuba = mo.ui.number(start=1, stop=None, value = 100, 
                          label = "Number of Turning Bands")
    WSeed = mo.ui.number(start=0, stop=None, value = seed, 
                         label = "Seed")


    return mo.ui.array([WNbtuba, WSeed])

def WshowSimtub(WAll):
    return mo.vstack(WAll)

def WgetSimtub(WAll):
    '''
    Returns the parameters for simulation using Turning Bands Method
    '''
    [WNbtuba, WSeed] = WAll
    return WNbtuba.value, WSeed.value

def WdefineDbFromBox(nech = 100, nvar = 1, xmin = 0, ymin = 0, xmax = 100, ymax = 100, seed = 14543):
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

def WgetDbFromBox(WAll):
    '''
    Create a Db with Radom Samples 
    '''
    [WNech, WNvar, WXmin, WYmin, WXmax, WYmax, WSeed] = WAll
    coormin = [WXmin.value, WYmin.value]
    coormax = [WXmax.value, WYmax.value]
    return gl.Db.createFillRandom(ndat = WNech.value, 
                                ndim = 2,  
                                nvar = WNvar.value, 
                                nfex = 0,
                                ncode = 0,
                                varmax = 0.,
                                selRatio = 0.,
                                heteroRatio = gl.VectorDouble(),
                                coormin = coormin,
                                coormax = coormax,
                                seed = WSeed.value)

def WdefineDbFromNF():
    '''
    Inquiry to load a file from a Neutral File
    '''
    WFile = mo.ui.file_browser(label="Select a Db Neutral File", multiple=False)

    return mo.ui.array([WFile])

def WgetDbFromNF(WAll):
    '''
    Create the gstlearn Db 
    '''
    [WFile] = WAll
    filename = WFile.name()
    if filename is None:
        return None
    return gl.Db.createFromNF(filename)

def WdefineDbFromGrid(nxdef = 10):
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

def WgetDbFromGrid(WAll):
    '''
    Create the gstlearn Grid from the widget WDbRandomGrid
    '''
    [WNX, WNY, WDX, WDY, WX0, WY0, WPerc] = WAll
    grid = gl.DbGrid.create(nx = [WNX.value, WNY.value],
                            dx = [WDX.value, WDY.value],
                            x0 = [WX0.value, WY0.value])
    return gl.Db.createFromGridRandomized(grid, WPerc.value);

def WdefineVarioParamOmni(nlag = 10, dlag = 1):
    '''
    Widget to define the Omnidirectional variogram
    '''
    WNlag = mo.ui.number(start=1, stop=100, value = nlag, 
                         label="Number of Lags")
    WDlag = mo.ui.number(start=1, stop=100, value = dlag, 
                         label="Lag Value")
    WToldis = mo.ui.number(start=0, stop=1, value = 0.5, 
                           label="Tolerance on Distance")
    WCylrad = mo.ui.number(start=0, stop=None, value = 0, 
                           label="Cylinder Radius")
    return mo.ui.array([WNlag, WDlag, WToldis, WCylrad])

def WshowVarioParamOmni(WAll):
    mo.vstack(WAll)

def WgetVarioParamOmni(WAll):
    [WNlag, WDlag, WToldis, WCylrad] = WAll
    if WCylrad.value > 0:
        varioparam = gl.VarioParam.createOmniDirection(npas = WNlag.value,
                                                       dpas = WDlag.value,
                                                       toldis = WToldis.value,
                                                       cylrad = WCylrad.value)
    else:
        varioparam = gl.VarioParam.createOmniDirection(npas = WNlag.value,
                                                       dpas = WDlag.value,
                                                       toldis = WToldis.value)
    return varioparam

def WshowVarioParamOmni(WAll):
    return mo.vstack(WAll)

def WdefineVarioParamMulti(ndir = 4, nlag = 10, dlag = 1):
    '''
    Widget to define the Multi-directional (regular) variogram
    '''
    WNdir = mo.ui.number(start=1, stop=10, value = ndir, 
                         label="Number of Directions")
    WNlag = mo.ui.number(start=1, stop=100, value = nlag, 
                         label="Number of Lags")
    WDlag = mo.ui.number(start=1, stop=100, value = dlag, 
                         label="Lag Value")
    WAngref = mo.ui.number(start=0, stop=180, value = 0., 
                           label="Reference angle (degree)")
    WToldis = mo.ui.number(start=0, stop=1, value = 0.5, 
                           label="Tolerance on Distance")
    return mo.ui.array([WNdir, WNlag, WDlag, WAngref, WToldis])

def WshowVarioParamMulti(WAll):
    return mo.vstack(WAll)

def WgetVarioParamMulti(WAll):
    [WNdir, WNlag, WDlag, WAngref, WToldis] = WAll
    varioparam = gl.VarioParam.createMultiple(ndir = WNdir.value,
                                              npas = WNlag.value,
                                              dpas = WDlag.value,
                                              toldis = WToldis.value,
                                              angref = WAngref.value)
    return varioparam

def WdefineDb(nech = 100, nvar = 1, xmin = 0, ymin = 0, xmax = 100, ymax = 100, 
              nxdef = 10, seed = 14543):
     
    WidgetDbFromBox = WdefineDbFromBox(nech = nech, nvar = nvar, 
                                  xmin = xmin, ymin = ymin, 
                                  xmax = xmax, ymax = ymax, seed = seed)

    WidgetDbFromNF = WdefineDbFromNF()

    WidgetDbFromGrid = WdefineDbFromGrid(nxdef = nxdef)

    return mo.ui.array([WidgetDbFromBox, WidgetDbFromNF, WidgetDbFromGrid])

def WshowDb(WAll):
    [WidgetDbFromBox, WidgetDbFromNF, WidgetDbFromGrid] = WAll

    return mo.ui.tabs(
        tabs = {
            "Random in Box":  mo.vstack(WidgetDbFromBox),
            "Random on Grid": mo.vstack(WidgetDbFromGrid),
            "From File":      mo.vstack(WidgetDbFromNF)
        })

def WgetDb(WAll, DbLayout):
    [WidgetDbFromBox, WidgetDbFromNF, WidgetDbFromGrid] = WAll

    db = None
    if DbLayout == "From File":
        db = WgetDbFromNF(WidgetDbFromNF)
    if DbLayout == "Random in Box":
        db = WgetDbFromBox(WidgetDbFromBox)
    if DbLayout == "Random on Grid":
        db = WgetDbFromGrid(WidgetDbFromGrid)

    if db is None:
        print("You must define a valid Db")

    return db

def WdefineVarioFromNF():
    '''
    Inquiry to load a Variogram File from a Neutral File
    '''
    WFile = mo.ui.file_browser(label="Select a Vario Neutral File", multiple=False)

    return mo.ui.array([WFile])

def WgetVarioFromNF(WAll):
    '''
    Create the gstlearn Vario 
    '''
    [WFile] = WAll
    filename = WFile.name()
    if filename is None:
        return None
    return gl.Vario.createFromNF(filename)

def WdefineVario(nlag = 10, dlag = 1, ndir = 4):
     
    WidgetVarioParamOmni = WdefineVarioParamOmni(nlag = nlag, dlag = dlag)

    WidgetVarioParamMulti = WdefineVarioParamMulti(ndir = ndir, nlag = nlag, dlag = dlag)

    WidgetVarioFromNF = WdefineVarioFromNF()

    return mo.ui.array([WidgetVarioParamOmni, WidgetVarioParamMulti, WidgetVarioFromNF])

def WshowVario(WAll):
    [WidgetVarioParamOmni, WidgetVarioParamMulti, WidgetVarioFromNF] = WAll

    return mo.ui.tabs(
        tabs = {
            "Omni-directional":    mo.vstack(WidgetVarioParamOmni),
            "Multiple Directions": mo.vstack(WidgetVarioParamMulti),
            "Read from NF":        mo.vstack(WidgetVarioFromNF)
        })

def WgetVario(WAll, VarioLayout, db):
    [WidgetVarioParamOmni, WidgetVarioParamMulti, WidgetVarioFromNF] = WAll

    varioparam = None
    if VarioLayout == "Omni-directional":
        varioparam = WgetVarioParamOmni(WidgetVarioParamOmni)
    if VarioLayout == "Multiple Directions":
        varioparam = WgetVarioParamMulti(WidgetVarioParamMulti)

    if VarioLayout == "Read from NF":
        vario = WgetVarioFromNF(WidgetVarioFromNF)
    else:
        if varioparam is None:
            print("You must define a valid VarioParam")
            return None

        if db is None:
            print("You must define a valid Db")
            return None
    
        vario = gl.Vario.computeFromDb(varioparam, db, 
                                       calcul = gl.ECalcVario.VARIOGRAM, 
                                       verbose = True)
    return vario

def WdefineCovList():
    '''
    Returns the widget for inquiring the list of basic structures to be used
    for fitting a Model to an Experimental variogram
    '''
    WTypes = mo.ui.multiselect(options=_getCovarianceDict(), value = ["Nugget effect"])

    return mo.ui.array([WTypes])

def WshowCovList(WAll):
    [WTypes] = WAll

    title = mo.md("Select the Basic Structures used for Fitting:")
    return mo.vstack([title, WTypes])

def WgetCovList(WAll, vario):
    [WTypes] = WAll

    if vario is None:
        print("You must define a valid Vario")
        return None
    
    # Create a list of ECov from 'options'
    model = None
    types = WTypes.value
    if types:
        model = gl.Model.createFromVario(vario, gl.ECov.fromKeys(types))
    return model
