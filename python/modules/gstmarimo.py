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
import pandas as pd
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

def WdefineCovariance(ic = 0, ncovmax = 1, distmax = 100, varmax = 100, defmodel=None):
    '''
    Returns the widget for inquiring the parameters for a single Basic structure
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    distmax: Maximum distance
    varmax:  Maximum Variance value
    defmodel: Model used for providing default values (optional)
    '''
    if  defmodel is None or ic >= defmodel.getCovaNumber():
        typeRef   = "Spherical"
        distRef   = distmax * (ic+1) / (ncovmax + 1)
        distAux   = distRef
        varRef    = varmax / ncovmax
        angRef    = 0
        flagAniso = False
    else:
        typeRef   = defmodel.getCovaType(ic).getDescr()
        distRef   = defmodel.getRange(ic)
        varRef    = defmodel.getSill(ic,0,0)
        distAux   = defmodel.getRanges(ic)[1]
        angRef    = defmodel.getCovAniso(ic).getAnisoAngles(1)
        flagAniso = defmodel.getCovAniso(ic).getFlagAniso()

    WUsed   = mo.ui.switch(True, label="Basic Structure Used")
    WType   = mo.ui.dropdown(options=_getCovarianceDict(), value = typeRef, label="Structure")
    WRange  = mo.ui.number(start=None, stop=None, value = distRef, label="Range")
    WSill   = mo.ui.number(start=0, stop=None, value = varRef, label="Sill")
    WAniso  = mo.ui.switch(value = False, label="Anisotropy")
    WRange2 = mo.ui.number(start=0, stop=None, value = distAux, label="Range Aux.")
    WAngle  = mo.ui.number(start=0, stop=None, value = angRef, label="Angle")

    return mo.ui.array([WUsed, WType, WRange, WSill, WAniso, WRange2, WAngle])

def WshowCovariance(WAll, flagTitle=True):
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

    return mo.vstack([WgetTitle("Covariance Definition", flagTitle),
                      WUsed, WTypeupd, WRangeupd, WSillupd, WAnisoupd, WRange2upd, WAngleupd])

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

def WdefineModel(ncovmax=1, distmax=100, varmax=100, defmodel=None):
    '''
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    distmax: Maximum distance
    varmax:  Maximum Variance value
    defmodel: Model used for providing default values
    '''
    if defmodel is not None:
        ncovmax = defmodel.getCovaNumber()
        distmax = defmodel.getMaximumDistance()
        varmax  = defmodel.getTotalSill()
    return mo.ui.array([WdefineCovariance(ic, ncovmax, distmax, varmax, defmodel) 
                         for ic in range(ncovmax)])

def WshowModel(WAlls, flagTitle=True):
    ncov = len(WAlls)
    UI = mo.accordion({
        "Covariance "+str(ic+1) : WshowCovariance(WAlls[ic], False) for ic in range(ncov)
    })
    return mo.vstack([WgetTitle("Model Definition", flagTitle), UI], justify='start')

def WgetModel(WAlls):
    '''
    Create a gstlearn Model
    '''
    model = gl.Model()
    for WAll in WAlls:
        cova = WgetCovariance(WAll)
        model.addCov(cova)
    return model

def WgetTitle(string, flagTitle=True):
    WTitle = mo.md("")
    if flagTitle:
        WTitle = mo.md("##" + string)
    return WTitle

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

def WshowGrid(WAll, flagTitle = True):
    [WNX, WNY, WDX, WDY, WX0, WY0] = WAll
    Wgrid = mo.hstack([
        mo.vstack([mo.md("Parameters"), mo.md("Number"), mo.md("Mesh"), mo.md("Origin")]),
        mo.vstack([mo.md("along X"), WNX, WDX, WX0],align='end'),
        mo.vstack([mo.md("along Y"), WNY, WDY, WY0],align='end')],justify='start')

    return mo.vstack([WgetTitle("Grid Definition", flagTitle), Wgrid])
   
def WgetGrid(WAll):
    '''
    Create the gstlearn Grid
    '''
    [WNX, WNY, WDX, WDY, WX0, WY0] = WAll
    grid = gl.DbGrid.create(nx = [WNX.value, WNY.value],
                            dx = [WDX.value, WDY.value],
                            x0 = [WX0.value, WY0.value])
    return grid

def WdefineSimtub(nbtuba=100, seed = 13134):
    '''
    Inquire for performing a Simulation using Turning Bands Method
    '''
    WNbtuba = mo.ui.number(start=1, stop=None, value = nbtuba, 
                          label = "Number of Turning Bands")
    WSeed = mo.ui.number(start=0, stop=None, value = seed, 
                         label = "Seed")

    return mo.ui.array([WNbtuba, WSeed])

def WshowSimtub(WAll, flagTitle=True):
    [WNbtuba, WSeed] = WAll
    return mo.vstack([WgetTitle("Parameters for Turning Bands Simulations", flagTitle), WNbtuba, WSeed])

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

def WdefineDbFromCSV():
    '''
    Inquiry to load a file from a CSV File
    '''
    WFile = mo.ui.file_browser(label="Select a CSV File", multiple=False)

    return mo.ui.array([WFile])

def WgetDbFromCSV(WAll, flagHeader=True, charSep=";", charDec=","):
    '''
    Create the gstlearn Db 
    '''
    [WFile] = WAll
    filename = WFile.name()
    if filename is None:
        return None
    
    csvformat = gl.CSVformat.create(flagHeader=flagHeader, charSep=charSep, charDec=charDec)
    return gl.Db.createFromCSV(filename, csvformat)

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
    WPerc = mo.ui.number(start=0, stop=100, value=10,  label="Rand. Percent.")
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
    WDlag = mo.ui.number(start=0, stop=100, value = dlag, 
                         label="Lag Value")
    WToldis = mo.ui.number(start=0, stop=1, value = 0.5, 
                           label="Tolerance on Distance")
    WCylrad = mo.ui.number(start=0, stop=None, value = 0, 
                           label="Cylinder Radius")
    return mo.ui.array([WNlag, WDlag, WToldis, WCylrad])

def WshowVarioParamOmni(WAll, flagTitle=True):
    [WNlag, WDlag, WToldis, WCylrad] = WAll
    mo.vstack(WgetTitle("Variogram Parameters", flagTitle), WNlag, WDlag, WToldis, WCylrad)

def WgetVarioParamOmni(WAll):
    [WNlag, WDlag, WToldis, WCylrad] = WAll
    if WCylrad.value > 0:
        varioparam = gl.VarioParam.createOmniDirection(nlag = WNlag.value,
                                                       dlag = WDlag.value,
                                                       toldis = WToldis.value,
                                                       cylrad = WCylrad.value)
    else:
        varioparam = gl.VarioParam.createOmniDirection(nlag = WNlag.value,
                                                       dlag = WDlag.value,
                                                       toldis = WToldis.value)
    return varioparam

def WdefineVarioParamMulti(ndir = 4, nlag = 10, dlag = 1):
    '''
    Widget to define the Multi-directional (regular) variogram
    '''
    WNdir = mo.ui.number(start=1, stop=10, value = ndir, 
                         label="Number of Directions")
    WNlag = mo.ui.number(start=1, stop=100, value = nlag, 
                         label="Number of Lags")
    WDlag = mo.ui.number(start=0, stop=100, value = dlag, 
                         label="Lag Value")
    WAngref = mo.ui.number(start=0, stop=180, value = 0., 
                           label="Reference angle (degree)")
    WToldis = mo.ui.number(start=0, stop=1, value = 0.5, 
                           label="Tolerance on Distance")
    return mo.ui.array([WNdir, WNlag, WDlag, WAngref, WToldis])

def WshowVarioParamMulti(WAll, flagTitle=True):
    [WNdir, WNlag, WDlag, WAngref, WToldis] = WAll
    return mo.vstack(WgetTitle("Variogram Definition", flagTitle), WNdir, WNlag, WDlag, WAngref, WToldis)

def WgetVarioParamMulti(WAll):
    [WNdir, WNlag, WDlag, WAngref, WToldis] = WAll
    varioparam = gl.VarioParam.createMultiple(ndir = WNdir.value,
                                              nlag = WNlag.value,
                                              dlag = WDlag.value,
                                              toldis = WToldis.value,
                                              angref = WAngref.value)
    return varioparam

def WdefineDb(nech = 100, nvar = 1, xmin = 0, ymin = 0, xmax = 100, ymax = 100, 
              nxdef = 10, seed = 14543, valdef = "From Box"):
     
    WidgetDbFromBox = WdefineDbFromBox(nech = nech, nvar = nvar, 
                                  xmin = xmin, ymin = ymin, 
                                  xmax = xmax, ymax = ymax, seed = seed)

    WidgetDbFromGrid = WdefineDbFromGrid(nxdef = nxdef)

    WidgetDbFromNF = WdefineDbFromNF()

    WidgetDbFromCSV = WdefineDbFromCSV()

    WidgetDbChoice = mo.ui.radio(
        options = {
            "From Box":  1,
            "From Grid": 2,
            "From NF":   3,
            "From CSV":  4
        },
        value = valdef)
    
    return mo.ui.array([WidgetDbChoice, WidgetDbFromBox, WidgetDbFromGrid, 
                        WidgetDbFromNF, WidgetDbFromCSV])

def WshowDb(WAll, flagTitle=True):
    [WidgetDbChoice, WidgetDbFromBox, WidgetDbFromGrid, WidgetDbFromNF, WidgetDbFromCSV] = WAll

    WTitle = WgetTitle("Data Base Parameters", flagTitle)
    option = WidgetDbChoice.value
    if option == 1:
        return mo.vstack([WTitle, WidgetDbChoice, *WidgetDbFromBox])
    elif option == 2:
        return mo.vstack([WTitle, WidgetDbChoice, *WidgetDbFromGrid])
    elif option == 3:
        return mo.vstack([WTitle, WidgetDbChoice, *WidgetDbFromNF])
    elif option == 4:
        return mo.vstack([WTitle, WidgetDbChoice, *WidgetDbFromCSV])
    else:
        return None

def WgetDb(WAll):
    [WidgetDbChoice, WidgetDbFromBox, WidgetDbFromGrid, WidgetDbFromNF, WidgetDbFromCSV] = WAll

    option = WidgetDbChoice.value
    db = None
    if option == 1:
        db = WgetDbFromBox(WidgetDbFromBox)
    elif option == 2:
        db = WgetDbFromGrid(WidgetDbFromGrid)
    elif option == 3:
        db = WgetDbFromNF(WidgetDbFromNF)
    elif option == 4:
        db = WgetDbFromCSV(WidgetDbFromCSV)
    else:
        db = None

    if db is None:
        print("You must define a valid Db beforehand")

    return db

def WdefineSaveNF(filename = "file.ascii"):
    '''
    Save the contents into a Neutral File
    '''
    Wbutton = mo.ui.run_button(label="Save")

    Wfilename = mo.ui.text(value=filename, label="Saving in Neutral File")

    return mo.ui.array([Wfilename, Wbutton])

def WshowSaveNF(WAll):
    return mo.hstack(WAll, justify='start')

def WperformSaveNF(Wall, contents):
    '''
    Save the contents into a Neutral File
    '''
    [Wfilename, Wbutton] = Wall

    if Wbutton.value:
        filename = Wfilename.value
        if filename is not None:
            contents.dumpToNF(filename)
        else:
            print("You must define a valid filename")

def WdefineVarioFromNF():
    '''
    Inquiry to load a Variogram File from a Neutral File
    '''
    WFile = mo.ui.file_browser(label="Select a Variogram Neutral File", multiple=False)

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

def WdefineVario(nlag = 10, dlag = 1, ndir = 4, valdef = "Omni"):
    
    WidgetVarioParamOmni = WdefineVarioParamOmni(nlag = nlag, dlag = dlag)

    WidgetVarioParamMulti = WdefineVarioParamMulti(ndir = ndir, nlag = nlag, dlag = dlag)

    WidgetVarioFromNF = WdefineVarioFromNF()

    WidgetVarioChoice = mo.ui.radio(
        options = {
            "Omni":  1,
            "Multi": 2,
            "From NF": 3
        },
        value = valdef)
    
    WidgetVarioSaveNF = WdefineSaveNF("vario.ascii")

    return mo.ui.array([WidgetVarioChoice, WidgetVarioParamOmni, WidgetVarioParamMulti, WidgetVarioFromNF, WidgetVarioSaveNF])

def WshowVario(WAll, flagTitle=True):
    [WidgetVarioChoice, WidgetVarioParamOmni, WidgetVarioParamMulti, WidgetVarioFromNF, WidgetVarioSaveNF] = WAll

    WTitle = WgetTitle("Variogram Parameters", flagTitle)
    option = WidgetVarioChoice.value
    if option == 1:
        return mo.vstack([WTitle, WidgetVarioChoice, *WidgetVarioParamOmni,
                          WshowSaveNF(WidgetVarioSaveNF)])
    elif option == 2:
        return mo.vstack([WTitle, WidgetVarioChoice, *WidgetVarioParamMulti, 
                          WshowSaveNF(WidgetVarioSaveNF)])
    elif option == 3:
        return mo.vstack([WTitle, WidgetVarioChoice, *WidgetVarioFromNF])
    else:
        return None

def WgetVario(WAll, db):
    [WidgetVarioChoice, WidgetVarioParamOmni, WidgetVarioParamMulti, WidgetVarioFromNF, WidgetVarioSaveNF] = WAll

    varioparam = None
    option = WidgetVarioChoice.value
    
    if option == 1:
        varioparam = WgetVarioParamOmni(WidgetVarioParamOmni)
    elif option == 2:
        varioparam = WgetVarioParamMulti(WidgetVarioParamMulti)
    elif option == 3:
        vario = WgetVarioFromNF(WidgetVarioFromNF)
        return vario
    else:
        varioparam = None
        
    if varioparam is None:
        print("You must define a valid VarioParam")
        return None

    if db is None:
        print("You must define a valid Db")
        return None
    
    vario = gl.Vario.computeFromDb(varioparam, db, 
                                   calcul = gl.ECalcVario.VARIOGRAM, 
                                   verbose = True)
    
    WperformSaveNF(WidgetVarioSaveNF, vario)

    return vario

def WdefineCovList(deftypes = ["Spherical"]):
    '''
    Returns the widget for inquiring the list of basic structures to be used
    for fitting a Model to an Experimental variogram

    deftypes: List containined the types of the defaulted basic structures 
    '''
    WTypes = mo.ui.multiselect(options=_getCovarianceDict(), value = deftypes)

    return mo.ui.array([WTypes])

def WshowCovList(WAll, flagTitle=True):
    [WTypes] = WAll

    return mo.vstack([WgetTitle("Select the Basic Structures used for Fitting", flagTitle), WTypes])

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

def WdefineBox(db = None):
    if db is not None and db.getNLoc(gl.ELoc.Z) == 2:
        box = db.getExtremas()
        longmin = box[0][0]
        longmax = box[1][0]
        latmin  = box[0][1]
        latmax  = box[1][1]
    else:
        longmin = -180
        longmax =  180
        latmin  =  -90
        latmax  =   90

    WLongMin = mo.ui.number(start=None, stop=None, value=longmin)
    WLongMax = mo.ui.number(start=None, stop=None, value=longmax)
    WLatMin  = mo.ui.number(start=None, stop=None, value=latmin)
    WLatMax  = mo.ui.number(start=None, stop=None, value=latmax)

    return mo.ui.array([WLongMin, WLongMax, WLatMin, WLatMax])

def WshowBox(WAll, flagTitle=True):
    [WLongMin, WLongMax, WLatMin, WLatMax] = WAll

    grid = mo.vstack([mo.hstack([mo.md("Longitude"), WLongMin, WLongMax]),
                      mo.hstack([mo.md("Latitude"),  WLatMin,  WLatMax])
                    ])
    
    Wgrid = mo.hstack([
        mo.vstack([mo.md("Parameters"), mo.md("Minimum"), mo.md("Maximum")]),
        mo.vstack([mo.md("Longitude"), WLongMin, WLongMax],align='end'),
        mo.vstack([mo.md("Latitude"),  WLatMin,  WLatMax],align='end')],justify='start')

    return mo.vstack([WgetTitle("Box Definition", flagTitle), Wgrid])

def WgetBox(WAll):
    [WLongMin, WLongMax, WLatMin, WLatMax] = WAll

    box = np.ndarray(shape=[2,2])
    box[0,0] = WLongMin.value
    box[0,1] = WLongMax.value
    box[1,0] = WLatMin.value
    box[1,1] = WLatMax.value
    return box

def WdefineGridN(nxdef = 50):
    WNX = mo.ui.number(start=1, stop=None, value = nxdef)
    WNY = mo.ui.number(start=1, stop=None, value = nxdef)
    return mo.ui.array([WNX, WNY])

def WshowGridN(WAll, flagTitle = True):
    [WNX, WNY] = WAll
    return mo.vstack([WgetTitle("Grid Discretization", flagTitle), WNX, WNY])
   
def WgetGridN(WAll, box):
    [WNX, WNY] = WAll

    nx = WNX.value
    ny = WNY.value

    deltax = box[0,1] - box[0,0]
    deltay = box[1,1] - box[1,0]
    dx = deltax / (nx-1)
    dy = deltay / (ny-1)
    x0 = box[0,0]
    y0 = box[1,0]
    return gl.DbGrid.create(nx = [nx,ny], dx = [dx,dy], x0 = [x0, y0])

def WdefineTable(db):
    names = db.getAllNames()
    df = pd.DataFrame(db[:], columns=names)
    WDF = mo.ui.table(df)
    return mo.ui.array([WDF])

def WshowTable(WAll, flagTitle = True):
    [WDF] = WAll
    return mo.vstack([WgetTitle("Grid Discretization", flagTitle), WDF])

