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

from distro import name
import gstlearn as gl
import numpy as np
import pandas as pd
import marimo as mo
import os

nxdef = 100
globalPath = "./"
debugOption = False
os.environ["GSTLEARN_OUTPUT_DIR"] = ""

def _getCovarianceDict():
    """
    Returns the list of covariances available as a Dictionary
    """
    keys = gl.ECov.getAllKeys(0)
    names = gl.ECov.getAllDescr(0)
    options = {}
    for k in np.arange(len(names)):
        options[names[k]] = keys[k]
    return options

def _WLock(WTest, condition, colorBackground="white", colorText="black"):
    """
    Turns to grey (as if it was locked if 'condition' is fulfilled)
    """
    if not condition:
        newWTest = WTest.style({"backgroundColor": colorBackground, "color": colorText})
    else:
        newWTest = WTest.style({"backgroundColor": "#f0f0f0", "color": "#a0a0a0"})
    return newWTest


def _WgetTitle(title, flagTitle=True):
    if flagTitle:
        return mo.md(f"## {title}")
    else:
        return mo.md("")

def _displayItem(contents = None):
    if contents is not None and debugOption:
        contents.display()

def _saveNF(contents = None, filename = "myFile.NF"):
    if contents is not None:
        contents.dumpToNF(filename)

#================================================================
# Widget to manage one Covariance
# The covariance is ranked within a list of 'ncovmax' covariances
#================================================================

def WdefineOneCovariance(ic=0, ncovmax=1, distmax=100, varmax=100, typeRef="Spherical"):
    """
    Returns the widget for inquiring the parameters for a single Basic structure
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    distmax: Maximum distance
    varmax:  Maximum Variance value
    typeRef: Type of covariance used as default
    """
    distRef = distmax * (ic + 1) / (ncovmax + 1)
    distAux = distRef
    varRef  = varmax / ncovmax
    angRef  = 0

    return {
        "WUsed": mo.ui.switch(True, label="Basic Structure Used"),
        "WType": mo.ui.dropdown(options=_getCovarianceDict(), value=typeRef, label="Structure"),
        "WRange": mo.ui.number(start=None, stop=None, value=distRef, label="Range"),
        "WSill": mo.ui.number(start=0, stop=None, value=varRef, label="Sill"),
        "WAniso": mo.ui.switch(value=False, label="Anisotropy"),
        "WRange2": mo.ui.number(start=0, stop=None, value=distAux, label="Range Aux."),
        "WAngle": mo.ui.number(start=0, stop=None, value=angRef, label="Angle")
    }

def WshowOneCovariance(Wall, flagTitle=True):
    WUsed   = Wall["WUsed"]
    WType   = Wall["WType"]
    WRange  = Wall["WRange"]
    WSill   = Wall["WSill"]
    WAniso  = Wall["WAniso"]
    WRange2 = Wall["WRange2"]
    WAngle  = Wall["WAngle"]

    WTypeupd   = _WLock(WType,   not WUsed.value)
    WRangeupd  = _WLock(WRange,  not WUsed.value)
    WSillupd   = _WLock(WSill,   not WUsed.value)
    WAnisoupd  = _WLock(WAniso,  not WUsed.value)
    WRange2upd = _WLock(WRange2, not WUsed.value or not WAniso.value)
    WAngleupd  = _WLock(WAngle,  not WUsed.value or not WAniso.value)

    items = []

    if flagTitle:
        items.append(_WgetTitle("Covariance Definition", True))

    items.extend([
        WUsed,
        WTypeupd,
        WRangeupd,
        WSillupd,
        WAnisoupd,
        WRange2upd,
        WAngleupd,
    ])

    return mo.vstack(items, gap=4)

def WgetOneCovariance(WAll):
    type_cov = gl.ECov.fromKey(WAll["WType"].value)
    cova = None

    if WAll["WUsed"].value:
        ctxt = gl.CovContext(1, 2)
        if not WAll["WAniso"].value:
            # isotropic covariance
            cova = gl.CovAniso.createIsotropic(
                ctxt,
                type=type_cov,
                range=WAll["WRange"].value,
                sill=WAll["WSill"].value,
                param=1.0,
                flagRange=True,
            )
        else:
            # anisotropic covariance
            cova = gl.CovAniso.createAnisotropic(
                ctxt,
                type=type_cov,
                ranges=[WAll["WRange"].value, WAll["WRange2"].value],
                sill=WAll["WSill"].value,
                param=1.0,
                angles=[WAll["WAngle"].value, 0.0],
                flagRange=True,
            )

    return cova

#=========================================
# Widget to manage the list of Covariances
#=========================================

def WdefineCovariances(ncovmax=1, distmax=100, varmax=100):
    """
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    distmax: Maximum distance
    varmax:  Maximum Variance value
    """
    return 
    [
        WdefineOneCovariance(ic, ncovmax, distmax, varmax) for ic in range(ncovmax)
    ]

def WshowCovariances(Wall, flagTitle=True):
    items = []

    if flagTitle:
        items.append(_WgetTitle("Model Definition", True))

    acc = mo.accordion(
        {
            f"Covariance {i + 1}": WshowOneCovariance(cov, flagTitle=False)
            for i, cov in enumerate(Wall)
        }
    )
    items.append(acc)

    return mo.vstack(items, justify="start", gap=6)

def WgetCovariances(WAll):
    model = gl.Model()
    for cov in WAll:
        cova = WgetOneCovariance(cov)
        if cova is not None:
            model.addCov(cova)
    return model

#===============================================================
# Widget to manage the list of Basic Structures used for Fitting
#===============================================================

def WshowBasicList(basic_list, flagTitle=True):
    items = []

    if flagTitle:
        items.append(_WgetTitle("Basic Structures for Fitting", True))

    items.append(basic_list["types"])

    return mo.vstack(items, gap=4)

#=========================
# Widget to manage a Model
#=========================

def WdefineModel(ncovmax=1, distmax=100, varmax=100, vario=None, 
                 deftypes=["Spherical"], valdef="Fit"):
    """
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    distmax: Maximum distance
    varmax:  Maximum Variance value
    vario: Vario used for providing default values (if provided)
    valdef: Defaulted option for Model definition
    """
    if vario is not None:
        distmax = vario.getMaximumDistance()
        varmax = vario.getVar()

    return {
        "WChoice": mo.ui.radio(options={"Define": 1, "Fit": 2, "From NF": 3}, value=valdef),
        "WDefine": WdefineCovariances(ncovmax=ncovmax, distmax=distmax, varmax=varmax),
        "WFitVario": WdefineModelFitVario(deftypes=deftypes),
        "WFromNF": WdefineModelFromNF(),
    }

def WshowModel(WAll, flagTitle=True):
    WTitle = _WgetTitle("Model Definition", flagTitle)
    choice = WAll["WChoice"]    
    option = choice.value

    # Contenu à afficher selon le choix
    if option == 1:
        content = mo.vstack([WshowOneCovariance(cov, False) for cov in WAll["WDefine"]], gap=6)
    elif option == 2:
        content = mo.vstack(WAll["WFitVario"]) if isinstance(WAll["WFitVario"], list) else WAll["WFitVario"]
    elif option == 3:
        content = WAll["WFromNF"] if hasattr(WAll["WFromNF"], "style") else mo.vstack([WAll["WFromNF"]])
    else:
        content = mo.md("Invalid option selected")

    return mo.vstack([WTitle, choice, content], gap=6)

def WgetModel(WAll, vario=None):
    option = WAll["WChoice"].value
    model = None

    if option == 1:
        model = WgetCovariances(WAll["WDefine"])
    elif option == 2:
        model = WgetModelFitVario(WAll["WFitVario"], vario)
    elif option == 3:
        model = WgetModelFromNF(WAll["WFromNF"])
    else:
        print("You must define a valid Model")
        return None

    _saveNF(model, "myModel.NF")
    _displayItem(model)

    return model

def WdefineModelFromNF():
    return {
        "WFile": mo.ui.file_browser(label="Select a Model Neutral File",
                               multiple=False, filetypes=[".NF", ".ascii"])
    }

def WgetModelFromNF(WAll):
    WFile = WAll["WFile"]

    filename = WFile.name()
    if filename is None:
        return None

    return gl.Model.createFromNF(str(WFile.path(index=0)))

def WdefineModelFitVario(deftypes=["Spherical"]):
    return {
        "WTypes": mo.ui.multiselect(options=_getCovarianceDict(), value=deftypes)
    }

def WgetModelFitVario(WAll, vario):
    WTypes = WAll["WTypes"]

    if vario is None:
        print("You must define a valid Vario")
        return None

    types = WTypes.value
    if not types:
        return None

    model = gl.Model.createFromVario(vario, gl.ECov.fromKeys(types))
    return model

#========================
# Widget to manage a Grid
#========================

def WdefineGrid(nxdef=50):
    """
    Returns parameters for a regular 2-D grid
    nxdef: Number of grid meshes (same along X and Y)
    """
    return {
        "WNX": mo.ui.slider(start=1, stop=200, value=nxdef),
        "WNY": mo.ui.slider(start=1, stop=200, value=nxdef),
        "WDX": mo.ui.number(start=1, stop=None, value=1),
        "WDY": mo.ui.number(start=1, stop=None, value=1),
        "WX0": mo.ui.number(start=0, stop=None, value=0),
        "WY0": mo.ui.number(start=0, stop=None, value=0),
    }

def WshowGrid(WAll, flagTitle=True):
    WNX = WAll["WNX"]
    WNY = WAll["WNY"]
    WDX = WAll["WDX"]
    WDY = WAll["WDY"]
    WX0 = WAll["WX0"]
    WY0 = WAll["WY0"]

    Wgrid = mo.hstack(
        [
            mo.vstack([mo.md("Parameters"), mo.md("Number"), mo.md("Mesh"), mo.md("Origin")]),
            mo.vstack([mo.md("along X"), WNX, WDX, WX0], align="end"),
            mo.vstack([mo.md("along Y"), WNY, WDY, WY0], align="end"),
        ],
        justify="start",
        gap=6,
    )

    return mo.vstack([_WgetTitle("Grid Definition", flagTitle), Wgrid], gap=6)

def WgetGrid(WAll):
    grid = gl.DbGrid.create(
        nx=[WAll["WNX"].value, WAll["WNY"].value],
        dx=[WAll["WDX"].value, WAll["WDY"].value],
        x0=[WAll["WX0"].value, WAll["WY0"].value],
    )
    return grid

#=============================
# Widget to manage Simulations
#=============================

def WdefineSimtub(nbtuba=100, seed=13134):
    """
    Returns parameters for performing Turning Bands simulations
    nbtuba: Number of Turning Bands
    seed: Seed for random number generator
    """
    return {
        "WNbtuba": mo.ui.number(start=1, stop=None, value=nbtuba, label="Number of Turning Bands"),
        "WSeed": mo.ui.number(start=0, stop=None, value=seed, label="Seed"),
    }   

def WshowSimtub(WAll, flagTitle=True):
    items = []

    if flagTitle:
        items.append(_WgetTitle("Parameters for Turning Bands Simulations", True))

    items.append(WAll["WNbtuba"])
    items.append(WAll["WSeed"])

    return mo.vstack(items, gap=4)

def WgetSimtub(WAll):
    return WAll["WNbtuba"].value, WAll["WSeed"].value

#=========================
# Widget to manage a Vario
#=========================

def WdefineVario(nlag=10, ndir=4, dlag=None, db=None, valdef="Omni"):
    """
    Returns parameters for calculating experimental variograms
    nlag: Number of Lags
    ndir: Number of Directions
    dlag: Lag Distance
    db: Database for calculating max distance
    valdef: Defaulted option for Vario definition
    """

    # Calculate the lag distance if not provided
    if dlag is None and db is not None:
        maxdist = db.getExtensionDiagonal()
        dlag = maxdist / nlag / 2.0
    elif dlag is None:
        dlag = 1.0

    return {
        "WChoice": mo.ui.radio(options={"Omni": 1, "Multi": 2, "From NF": 3}, value=valdef),
        "WOmni": WdefineVarioParamOmni(nlag=nlag, dlag=dlag),
        "WMulti": WdefineVarioParamMulti(ndir=ndir, nlag=nlag, dlag=dlag),
        "WFromNF": WdefineVarioFromNF(),
    }

def WshowVario(WAll, flagTitle=True):
    WTitle = _WgetTitle("Variogram Parameters", flagTitle)
    choice = WAll["WChoice"]
    option = choice.value

    # Sélection du contenu selon le choix
    if option == 1:
        content = mo.vstack(WAll["WOmni"]) if isinstance(WAll["WOmni"], list) else WAll["WOmni"]
    elif option == 2:
        content = mo.vstack(WAll["WMulti"]) if isinstance(WAll["WMulti"], list) else WAll["WMulti"]
    elif option == 3:
        content = WAll["WFromNF"] if hasattr(WAll["WFromNF"], "style") else mo.vstack([WAll["WFromNF"]])
    else:
        content = mo.md("Invalid selection")

    return mo.vstack([WTitle, choice, content], gap=6)

def WgetVario(WAll, db=None):
    option = WAll["WChoice"].value
    varioparam = None

    if option == 1:
        varioparam = WgetVarioParamOmni(WAll["WOmni"])
    elif option == 2:
        varioparam = WgetVarioParamMulti(WAll["WMulti"])
    elif option == 3:
        return WgetVarioFromNF(WAll["WFromNF"])
    else:
        print("You must define a valid VarioParam")
        return None

    vario = None
    if db is None:
        print("To calculate a Variogram, you must define a valid Db")
    else:
        vario = gl.Vario.computeFromDb(
            varioparam, db, calculType=gl.ECalcVario.VARIOGRAM, verbose=True
        )

    _saveNF(vario, "myVario.NF")
    _displayItem(vario)

    return vario

def WdefineVarioParamOmni(nlag=10, dlag=1):
    return {
        "WNlag": mo.ui.number(start=1, stop=100, value=nlag, label="Number of Lags"),
        "WDlag": mo.ui.number(start=0, stop=100, value=dlag, label="Lag Value"),
        "WToldis": mo.ui.number(start=0, stop=1, value=0.5, label="Tolerance on Distance"),
        "WCylrad": mo.ui.number(start=0, stop=None, value=0, label="Cylinder Radius")
    }

def WshowVarioParamOmni(WAll, flagTitle=True):
    items = []

    if flagTitle:
        items.append(_WgetTitle("Variogram Parameters", True))

    items.extend([
        WAll["WNLag"],
        WAll["WDlag"],
        WAll["WToldis"],
        WAll["WCylrad"],
    ])

    return mo.vstack(items, gap=4)

def WgetVarioParamOmni(WAll):
    nlag = WAll["WNlag"].value
    dlag = WAll["WDlag"].value
    toldis = WAll["WToldis"].value
    cylrad = WAll["WCylrad"].value
    if cylrad > 0:
        varioparam = gl.VarioParam.createOmniDirection(
            nlag=nlag, dlag=dlag, toldis=toldis, cylrad=cylrad
        )
    else:
        varioparam = gl.VarioParam.createOmniDirection(
            nlag=nlag, dlag=dlag, toldis=toldis
        )

    return varioparam

def WdefineVarioParamMulti(ndir=4, nlag=10, dlag=1):
    return {
        "WNdir": mo.ui.number(start=1, stop=10, value=ndir, label="Number of Directions"),
        "WNlag": mo.ui.number(start=1, stop=100, value=nlag, label="Number of Lags"),
        "WDlag": mo.ui.number(start=0, stop=100, value=dlag, label="Lag Value"),
        "WAngref": mo.ui.number(start=0, stop=180, value=0.0, label="Reference angle (degree)"),
        "WToldis": mo.ui.number(start=0, stop=1, value=0.5, label="Tolerance on Distance"),
    }


def WshowVarioParamMulti(WAll, flagTitle=True):
    items = []

    if flagTitle:
        items.append(_WgetTitle("Variogram Definition", True))

    items.extend([
        WAll["WNdir"],
        WAll["WNlag"],
        WAll["WDlag"],
        WAll["WAngref"],
        WAll["WToldis"],
    ])

    return mo.vstack(items, gap=4)

def WgetVarioParamMulti(WAll):
    varioparam = gl.VarioParam.createMultiple(
        ndir=WAll["WNdir"].value,
        nlag=WAll["WNlag"].value,
        dlag=WAll["WDlag"].value,
        toldis=WAll["WToldis"].value,
        angref=WAll["WAngref"].value,
    )
    return varioparam

def WdefineVarioFromNF():
    return {
        "WFile": mo.ui.file_browser(label="Select a Variogram Neutral File", 
                               multiple=False, filetypes=[".NF", ".ascii"])
    }

def WgetVarioFromNF(WAll):
    WFile = WAll["WFile"]

    filename = WFile.name()
    if filename is None:
        return None

    return gl.Vario.createFromNF(str(WFile.path(index=0)))

#======================
# Widget to manage a Db
#======================

def WdefineDb(nech=100, nvar=1, xmin=0, ymin=0, xmax=100, ymax=100, nxdef=10, seed=145234, valdef="From Box"):
    """
    Returns parameters for constructing a Db
    nech: Number of Samples
    nvar: Number of Variables
    xmin: Minimum along X
    ymin: Minimum along Y
    xmax: Maximum along X
    ymax: Maximum along Y
    nxdef: Number of grid meshes (same along X and Y)   
    nbtuba: Number of Turning Bands
    seed: Seed for random number generator
    """
    return {
        "WChoice": mo.ui.radio(
            options={"From Box": 1, "From Grid": 2, "From NF": 3, "From CSV": 4}, 
            value=valdef
        ),
        "WFromBox": WdefineDbFromBox(nech, nvar, xmin, ymin, xmax, ymax, seed),
        "WFromGrid": WdefineDbFromGrid(nvar, nxdef),
        "WFromNF": WdefineDbFromNF(),
        "WFromCSV": WdefineDbFromCSV(),
    }

def WshowDb(WAll, flagTitle=True):
    WTitle = _WgetTitle("Data Base Parameters", flagTitle)
    choice = WAll["WChoice"]
    option = choice.value

    if option == 1:
        content = WAll["WFromBox"]
    elif option == 2:
        content = WAll["WFromGrid"]
    elif option == 3:
        content = WAll["WFromNF"]
    elif option == 4:
        content = WAll["WFromCSV"]
    else:
        content = mo.md("Invalid selection")

    return mo.vstack([WTitle, choice, content], gap=6)

def WgetDb(WAll):
    option = WAll["WChoice"].value
    db = None

    if option == 1:
        db = WgetDbFromBox(WAll["WFromBox"])
    elif option == 2:
        db = WgetDbFromGrid(WAll["WFromGrid"])
    elif option == 3:
        db = WgetDbFromNF(WAll["WFromNF"])
    elif option == 4:
        db = WgetDbFromCSV(WAll["WFromCSV"])
    else:
        print("No valid database option selected")
        db = None

    _displayItem(db)
    return db

def WdefineDbFromGrid(nvar=1, nxdef=10):
    return {
        "WNX": mo.ui.number(start=1, stop=100, value=nxdef, label="NX"),
        "WNY": mo.ui.number(start=1, stop=100, value=nxdef, label="NY"),
        "WDX": mo.ui.number(start=1, stop=None, value=1, label="DX"),
        "WDY": mo.ui.number(start=1, stop=None, value=1, label="DY"),
        "WX0": mo.ui.number(start=0, stop=None, value=0, label="X0"),
        "WY0": mo.ui.number(start=0, stop=None, value=0, label="Y0"),
        "WNvar": mo.ui.number(start=1, stop=None, value=nvar, label="Number of Variables"),
        "WPerc": mo.ui.number(start=0, stop=100, value=10, label="Rand. Percent.")
    }

def WgetDbFromGrid(WAll):
    grid = gl.DbGrid.create(
        nx=[WAll["WNX"].value, WAll["WNY"].value],
        dx=[WAll["WDX"].value, WAll["WDY"].value],
        x0=[WAll["WX0"].value, WAll["WY0"].value],
    )
    return gl.Db.createFromGridRandomized(
        grid,
        nvar=WAll["WNvar"].value,
        selRatio=WAll["WPerc"].value
    )

def WdefineDbFromBox(nech=100, nvar=1, xmin=0, ymin=0, xmax=100, ymax=100, seed=14543):
    return {
        "WNech": mo.ui.number(start=1, stop=None, value=nech, label="Number of Samples"),
        "WNvar": mo.ui.number(start=1, stop=None, value=nvar, label="Number of Variables"),
        "WXmin": mo.ui.number(start=None, stop=None, value=xmin, label="Minimum along X"),
        "WYmin": mo.ui.number(start=None, stop=None, value=ymin, label="Minimum along Y"),
        "WXmax": mo.ui.number(start=None, stop=None, value=xmax, label="Maximum along X"),
        "WYmax": mo.ui.number(start=None, stop=None, value=ymax, label="Maximum along Y"),
        "WSeed": mo.ui.number(start=None, stop=None, value=seed, label="Seed"),
    }

def WgetDbFromBox(WAll):
    coormin = [WAll["WXmin"].value, WAll["WYmin"].value]
    coormax = [WAll["WXmax"].value, WAll["WYmax"].value]

    return gl.Db.createFillRandom(
        ndat=WAll["WNech"].value,
        ndim=2,
        nvar=WAll["WNvar"].value,
        nfex=0,
        ncode=0,
        varmax=0.0,
        selRatio=0.0,
        heteroRatio=gl.VectorDouble(),
        coormin=coormin,
        coormax=coormax,
        seed=WAll["WSeed"].value,
    )

def WdefineDbFromNF():
    return {
        "WFile": mo.ui.file_browser(label="Select a Db Neutral File", 
                               multiple=False, filetypes=[".NF", ".ascii"])
    }

def WgetDbFromNF(WAll):
    WFile = WAll["WFile"]

    filename = WFile.name()
    if filename is None:
        return None

    return gl.Db.createFromNF(str(WFile.path(index=0)))

def WdefineDbFromCSV(nameX="Longitude", nameY="Latitude", nameVar="pH"):
    return {
        "WnameX": mo.ui.text(label="X Coordinate", value=nameX),
        "WnameY": mo.ui.text(label="Y Coordinate", value=nameY),
        "WnameVar": mo.ui.text(label="Variable Name", value=nameVar),
        "WFile": mo.ui.file_browser(label="Select a CSV File", multiple=False, filetypes=[".csv"])
    }

def WgetDbFromCSV(WAll, flagHeader=True, charSep=";", charDec=","):
    WFile = WAll["WFile"]
    WnameX = WAll["WnameX"]
    WnameY = WAll["WnameY"]
    WnameVar = WAll["WnameVar"]

    filename = WFile.name()
    if filename is None:
        return None

    path = WFile.path(index=0)
    dataframe = pd.read_csv(path, sep=charSep, decimal=charDec, header=0 if flagHeader else None)
    db = gl.Db_fromPandas(dataframe)

    db.setLocators([WnameX.value, WnameY.value], gl.ELoc.X)
    db.setLocator(WnameVar.value, gl.ELoc.Z)

    return db



#===============================================
# Widget to manage a Box based on an optional Db
#===============================================

def WdefineBox(db=None):
    """
    Returns parameters for defining a Box (by meshes only)
    db: Database (optional) for providing default values
    """
    if db is not None and db.getNLoc(gl.ELoc.Z) == 2:
        box = db.getExtremas()
        longmin = box[0][0]
        longmax = box[1][0]
        latmin = box[0][1]
        latmax = box[1][1]
    else:
        longmin = -180
        longmax = 180
        latmin = -90
        latmax = 90

    return {
        "WLongMin": mo.ui.number(start=None, stop=None, value=longmin),
        "WLongMax": mo.ui.number(start=None, stop=None, value=longmax),
        "WLatMin": mo.ui.number(start=None, stop=None, value=latmin),
        "WLatMax": mo.ui.number(start=None, stop=None, value=latmax),
    }

def WshowBox(WAll, flagTitle=True):
    WLongMin = WAll["WLongMin"]
    WLongMax = WAll["WLongMax"]
    WLatMin  = WAll["WLatMin"]
    WLatMax  = WAll["WLatMax"]
    Wgrid = mo.hstack(
        [
            mo.vstack([mo.md("Parameters"), mo.md("Minimum"), mo.md("Maximum")]),
            mo.vstack([mo.md("Longitude"), WLongMin, WLongMax], align="end"),
            mo.vstack([mo.md("Latitude"), WLatMin, WLatMax], align="end"),
        ],
        justify="start",
        gap=6,
    )

    return mo.vstack([_WgetTitle("Box Definition", flagTitle), Wgrid], gap=6)

def WgetBox(WAll):
    box = np.ndarray(shape=(2, 2))
    box[0, 0] = WAll["WLongMin"].value
    box[0, 1] = WAll["WLongMax"].value
    box[1, 0] = WAll["WLatMin"].value
    box[1, 1] = WAll["WLatMax"].value
    return box

#=======================================
# Widget to manage a discretization Grid
#=======================================

def WdefineGridN(nxdef=50):
    """
    Returns parameters for defining a discretization Grid
    nxdef: Number of grid meshes (same along X and Y)
    """
    return {
        "WNX": mo.ui.number(start=1, stop=None, value=nxdef),
        "WNY": mo.ui.number(start=1, stop=None, value=nxdef),
    }

def WshowGridN(WAll, flagTitle=True):
    items = []

    if flagTitle:
        items.append(_WgetTitle("Grid Discretization", True))

    items.append(WAll["WNX"])
    items.append(WAll["WNY"])

    return mo.vstack(items, gap=4)

def WgetGridN(WAll, box):
    nx = WAll["WNX"].value
    ny = WAll["WNY"].value

    deltax = box[0, 1] - box[0, 0]
    deltay = box[1, 1] - box[1, 0]
    dx = deltax / (nx - 1)
    dy = deltay / (ny - 1)
    x0 = box[0, 0]
    y0 = box[1, 0]

    return gl.DbGrid.create(nx=[nx, ny], dx=[dx, dy], x0=[x0, y0])

