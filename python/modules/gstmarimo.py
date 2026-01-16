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

from click import option
from distro import name
import gstlearn as gl
import numpy as np
import pandas as pd
import marimo as mo
import os

nxdef = 100
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


def _displayItem(contents=None, flagForced=False):
    if contents is None:
        return
    if debugOption or flagForced:
        contents.display()


def _saveNF(contents=None, filename="myFile.NF"):
    if contents is not None:
        contents.dumpToNF(filename)


def _WdefineFromNF():
    WFile = mo.ui.file_browser(
        label="Select a Neutral File", multiple=False, filetypes=[".NF", ".ascii"]
    )
    return mo.ui.array([WFile])


# ================================================================
# Widget to manage one Covariance
# The covariance is ranked within a list of 'ncovmax' covariances
# ================================================================


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
    varRef = varmax / ncovmax
    angRef = 0

    WUsed = mo.ui.switch(True, label="Basic Structure Used")
    WType = mo.ui.dropdown(
        options=_getCovarianceDict(), value=typeRef, label="Structure"
    )
    WRange = mo.ui.number(start=None, stop=None, value=distRef, label="Range")
    WSill = mo.ui.number(start=0, stop=None, value=varRef, label="Sill")
    WAniso = mo.ui.switch(value=False, label="Anisotropy")
    WRange2 = mo.ui.number(start=0, stop=None, value=distAux, label="Range Aux.")
    WAngle = mo.ui.number(start=0, stop=None, value=angRef, label="Angle")

    return mo.ui.array([WUsed, WType, WRange, WSill, WAniso, WRange2, WAngle])


def WshowOneCovariance(WAll, flagTitle=True):
    [WUsed, WType, WRange, WSill, WAniso, WRange2, WAngle] = WAll

    WTitle = _WgetTitle("Covariance Definition", flagTitle)

    WTypeupd = _WLock(WType, not WUsed.value)
    WRangeupd = _WLock(WRange, not WUsed.value)
    WSillupd = _WLock(WSill, not WUsed.value)
    WAnisoupd = _WLock(WAniso, not WUsed.value)
    WRange2upd = _WLock(WRange2, not WUsed.value or not WAniso.value)
    WAngleupd = _WLock(WAngle, not WUsed.value or not WAniso.value)

    return mo.ui.array(
        [WTitle, WUsed, WTypeupd, WRangeupd, WSillupd, WAnisoupd, WRange2upd, WAngleupd]
    )


def WgetOneCovariance(WAll):
    [WUsed, WType, WRange, WSill, WAniso, WRange2, WAngle] = WAll

    if WUsed.value:
        if not WAniso.value:
            # isotropic covariance
            return gl.CovAniso.createIsotropic(
                ctxt=gl.CovContext(1, 2),
                type=gl.ECov.fromKey(WType.value),
                range=WRange.value,
                sill=WSill.value,
                param=1.0,
                flagRange=True,
            )
        else:
            # anisotropic covariance
            return gl.CovAniso.createAnisotropic(
                ctxt=gl.CovContext(1, 2),
                type=gl.ECov.fromKey(WType.value),
                ranges=[WRange.value, WRange2.value],
                sill=WSill.value,
                param=1.0,
                angles=[WAngle.value, 0.0],
                flagRange=True,
            )
    else:
        return None


# =========================================
# Widget to manage the list of Covariances
# =========================================


def WdefineCovariances(ncovmax=1, distmax=100, varmax=100):
    """
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    distmax: Maximum distance
    varmax:  Maximum Variance value
    """
    return mo.ui.array(
        [WdefineOneCovariance(ic, ncovmax, distmax, varmax) for ic in range(ncovmax)]
    )


def WshowCovariances(WAll, flagTitle=True):
    WTitle = _WgetTitle("Model Definition", flagTitle)
    UI = mo.accordion(
        {
            f"Covariance {ic + 1}": WshowOneCovariance(cov, False)
            for ic, cov in enumerate(WAll)
        }
    )
    return mo.ui.array([WTitle, UI])


def WgetCovariances(WAll):
    model = gl.Model()
    for cov in WAll:
        cova = WgetOneCovariance(cov)
        if cova is not None:
            model.addCov(cova)
    return model


# ===============================================================
# Widget to manage the list of Basic Structures used for Fitting
# ===============================================================


def WshowBasicList(basic_list, flagTitle=True):
    WTitle = _WgetTitle("Basic Structures for Fitting", flagTitle)
    WList = basic_list["types"]
    return mo.ui.array([WTitle, WList])


# =========================
# Widget to manage a Model
# =========================


def WdefineModel(
    ncovmax=1, distmax=100, varmax=100, vario=None, deftypes=["Spherical"], valdef="Fit"
):
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

    WChoice = mo.ui.radio(
        options={"Interactive": 1, "Fit": 2, "From NF": 3}, value=valdef
    )
    WInter = WdefineCovariances(ncovmax=ncovmax, distmax=distmax, varmax=varmax)
    WFitVario = WdefineModelFitVario(deftypes=deftypes)
    WFromNF = _WdefineFromNF()

    return mo.ui.array([WChoice, WInter, WFitVario, WFromNF])


def WshowModel(WAll, flagTitle=True, gapv=2):
    [WChoice, WInter, WFitVario, WFromNF] = WAll

    WTitle = _WgetTitle("Model Definition", flagTitle)
    option = WChoice.value

    # Contenu à afficher selon le choix
    if option == 1:
        UI = mo.accordion(
            {f"Covariance {ic + 1}": mo.vstack(WInter[ic]) for ic in range(len(WInter))}
        )
        return mo.vstack([WTitle, WChoice, UI], gap=gapv)
    elif option == 2:
        return mo.vstack([WTitle, WChoice, *WFitVario], gap=gapv)
    elif option == 3:
        return mo.vstack([WTitle, WChoice, *WFromNF], gap=gapv)
    else:
        return mo.md("Invalid option selected")


def WgetModel(WAll, vario=None):
    [WChoice, WInter, WFitVario, WFromNF] = WAll

    option = WChoice.value

    model = None
    if option == 1:
        model = WgetCovariances(WInter)
    elif option == 2:
        model = WgetModelFitVario(WFitVario, vario)
    elif option == 3:
        model = WgetModelFromNF(WFromNF)
    else:
        return None

    if model is not None:
        _saveNF(model, "myModel.NF")
        _displayItem(model)

    return model


def WgetModelFromNF(WAll):
    [WFile] = WAll
    filename = WFile.name()
    if filename is None:
        return None
    return gl.Model.createFromNF(str(WFile.path(index=0)))


def WdefineModelFitVario(deftypes=["Spherical"]):
    WTypes = mo.ui.multiselect(options=_getCovarianceDict(), value=deftypes)
    return mo.ui.array([WTypes])


def WgetModelFitVario(WAll, vario):
    [WTypes] = WAll
    if vario is None:
        return None

    types = WTypes.value
    if not types:
        return None

    model = gl.Model.createFromVario(vario, gl.ECov.fromKeys(types))
    return model


# ========================
# Widget to manage a Grid
# ========================


def WdefineGrid(nxdef=50):
    """
    Returns parameters for a regular 2-D grid
    nxdef: Number of grid meshes (same along X and Y)
    """
    WNX = mo.ui.slider(start=1, stop=200, value=nxdef)
    WNY = mo.ui.slider(start=1, stop=200, value=nxdef)
    WDX = mo.ui.number(start=1, stop=None, value=1)
    WDY = mo.ui.number(start=1, stop=None, value=1)
    WX0 = mo.ui.number(start=0, stop=None, value=0)
    WY0 = mo.ui.number(start=0, stop=None, value=0)
    return mo.ui.array([WNX, WNY, WDX, WDY, WX0, WY0])


def WshowGrid(WAll, flagTitle=True, gaph=2, gapv=2):
    [WNX, WNY, WDX, WDY, WX0, WY0] = WAll
    WTitle = _WgetTitle("Grid Definition", flagTitle)
    Wgrid = mo.hstack(
        [
            mo.vstack(
                [mo.md("Parameters"), mo.md("Number"), mo.md("Mesh"), mo.md("Origin")],
                gap=gapv,
            ),
            mo.vstack([mo.md("along X"), WNX, WDX, WX0], align="end", gap=gapv),
            mo.vstack([mo.md("along Y"), WNY, WDY, WY0], align="end", gap=gapv),
        ],
        gap=gaph,
    )
    return mo.vstack([WTitle, Wgrid], gap=gapv)


def WgetGrid(WAll):
    [WNX, WNY, WDX, WDY, WX0, WY0] = WAll
    grid = gl.DbGrid.create(
        nx=[WNX.value, WNY.value],
        dx=[WDX.value, WDY.value],
        x0=[WX0.value, WY0.value],
    )
    return grid


# =============================
# Widget to manage Simulations
# =============================


def WdefineSimtub(nbtuba=100, seed=13134):
    """
    Returns parameters for performing Turning Bands simulations
    nbtuba: Number of Turning Bands
    seed: Seed for random number generator
    """
    WNbtuba = mo.ui.number(
        start=1, stop=None, value=nbtuba, label="Number of Turning Bands"
    )
    WSeed = mo.ui.number(start=0, stop=None, value=seed, label="Seed")
    return mo.ui.array([WNbtuba, WSeed])


def WshowSimtub(WAll, flagTitle=True):
    [WNbtuba, WSeed] = WAll

    WTitle = _WgetTitle("Parameters for Turning Bands Simulations", flagTitle)

    return mo.ui.array([WTitle, WNbtuba, WSeed])


def WgetSimtub(WAll):
    [WNbtuba, WSeed] = WAll
    return WNbtuba.value, WSeed.value


# =========================
# Widget to manage a Vario
# =========================


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

    WChoice = mo.ui.radio(options={"Omni": 1, "Multi": 2, "From NF": 3}, value=valdef)
    WOmni = WdefineVarioParamOmni(nlag=nlag, dlag=dlag)
    WMulti = WdefineVarioParamMulti(ndir=ndir, nlag=nlag, dlag=dlag)
    WFromNF = _WdefineFromNF()
    return mo.ui.array([WChoice, WOmni, WMulti, WFromNF])


def WshowVario(WAll, flagTitle=True):
    [WChoice, WOmni, WMulti, WFromNF] = WAll

    WTitle = _WgetTitle("Variogram Parameters", flagTitle)
    option = WChoice.value

    # Sélection du contenu selon le choix
    if option == 1:
        return mo.vstack([WTitle, WChoice, *WOmni])
    elif option == 2:
        return mo.vstack([WTitle, WChoice, *WMulti])
    elif option == 3:
        return mo.vstack([WTitle, WChoice, *WFromNF])
    else:
        return mo.md("Invalid selection")


def WgetVario(WAll, db=None):
    [WChoice, WOmni, WMulti, WFromNF] = WAll
    option = WChoice.value

    varioparam = None
    if option == 1:
        varioparam = WgetVarioParamOmni(WOmni)
    elif option == 2:
        varioparam = WgetVarioParamMulti(WMulti)
    elif option == 3:
        return WgetVarioFromNF(WFromNF)
    else:
        return None

    vario = None
    if varioparam is not None and db is not None:
        vario = gl.Vario.computeFromDb(
            varioparam, db, calculType=gl.ECalcVario.VARIOGRAM, verbose=True
        )

    if vario is not None:
        _saveNF(vario, "myVario.NF")
        _displayItem(vario)

    return vario


def WdefineVarioParamOmni(nlag=10, dlag=1):
    WNlag = mo.ui.number(start=1, stop=100, value=nlag, label="Number of Lags")
    WDlag = mo.ui.number(start=0, stop=100, value=dlag, label="Lag Value")
    WToldis = mo.ui.number(start=0, stop=1, value=0.5, label="Tolerance on Distance")
    WCylrad = mo.ui.number(start=0, stop=None, value=0, label="Cylinder Radius")
    return mo.ui.array([WNlag, WDlag, WToldis, WCylrad])


def WshowVarioParamOmni(WAll, flagTitle=True):
    [WNlag, WDlag, WToldis, WCylrad] = WAll

    WTitle = _WgetTitle("Variogram Parameters", flagTitle)
    return mo.ui.array([WTitle, WNlag, WDlag, WToldis, WCylrad])


def WgetVarioParamOmni(WAll):
    [WNlag, WDlag, WToldis, WCylrad] = WAll
    if WCylrad.value > 0:
        varioparam = gl.VarioParam.createOmniDirection(
            nlag=WNlag.value,
            dlag=WDlag.value,
            toldis=WToldis.value,
            cylrad=WCylrad.value,
        )
    else:
        varioparam = gl.VarioParam.createOmniDirection(
            nlag=WNlag.value, dlag=WDlag.value, toldis=WToldis.value
        )

    return varioparam


def WdefineVarioParamMulti(ndir=4, nlag=10, dlag=1):
    WNdir = mo.ui.number(start=1, stop=10, value=ndir, label="Number of Directions")
    WNlag = mo.ui.number(start=1, stop=100, value=nlag, label="Number of Lags")
    WDlag = mo.ui.number(start=0, stop=100, value=dlag, label="Lag Value")
    WAngref = mo.ui.number(
        start=0, stop=180, value=0.0, label="Reference angle (degree)"
    )
    WToldis = mo.ui.number(start=0, stop=1, value=0.5, label="Tolerance on Distance")
    return mo.ui.array([WNdir, WNlag, WDlag, WAngref, WToldis])


def WshowVarioParamMulti(WAll, flagTitle=True):
    [WNdir, WNlag, WDlag, WAngref, WToldis] = WAll

    WTitle = _WgetTitle("Variogram Definition", flagTitle)

    return mo.ui.array([WTitle, WNdir, WNlag, WDlag, WAngref, WToldis])


def WgetVarioParamMulti(WAll):
    [WNdir, WNlag, WDlag, WAngref, WToldis] = WAll
    return gl.VarioParam.createMultiple(
        ndir=WNdir.value,
        nlag=WNlag.value,
        dlag=WDlag.value,
        toldis=WToldis.value,
        angref=WAngref.value,
    )


def WgetVarioFromNF(WAll):
    [WFile] = WAll
    filename = WFile.name()
    if filename is None:
        return None
    return gl.Vario.createFromNF(str(WFile.path(index=0)))


# ======================
# Widget to manage a Db
# ======================


def WdefineDb(
    nech=100,
    nvarmax=1,
    xmin=0,
    ymin=0,
    xmax=100,
    ymax=100,
    nxdef=10,
    seed=145234,
    valdef="From Box",
):
    """
    Returns parameters for constructing a Db
    nech: Number of Samples
    nvarmax: Number of Variables
    xmin: Minimum along X
    ymin: Minimum along Y
    xmax: Maximum along X
    ymax: Maximum along Y
    nxdef: Number of grid meshes (same along X and Y)
    nbtuba: Number of Turning Bands
    seed: Seed for random number generator
    """

    WChoice = mo.ui.radio(
        options={"From Box": 1, "From Grid": 2, "From CSV": 3, "From NF": 4},
        value=valdef,
    )
    WFromBox = WdefineDbFromBox(nech, nvarmax, xmin, ymin, xmax, ymax, seed)
    WFromGrid = WdefineDbFromGrid(nvarmax, nxdef, seed)
    WFromCSV = WdefineDbFromCSV()
    WFromNF = _WdefineFromNF()

    return mo.ui.array([WChoice, WFromBox, WFromGrid, WFromCSV, WFromNF])


def WshowDb(WAll, flagTitle=True):
    [WChoice, WFromBox, WFromGrid, WFromCSV, WFromNF] = WAll

    WTitle = _WgetTitle("Data Base Parameters", flagTitle)
    option = WChoice.value

    if option == 1:
        return mo.vstack([WTitle, WChoice, *WFromBox])
    elif option == 2:
        return mo.vstack([WTitle, WChoice, *WFromGrid])
    elif option == 3:
        return mo.vstack([WTitle, WChoice, *WFromCSV])
    elif option == 4:
        return mo.vstack([WTitle, WChoice, *WFromNF])
    else:
        return None


def WgetDb(WAll):
    [WChoice, WFromBox, WFromGrid, WFromCSV, WFromNF] = WAll
    option = WChoice.value

    if option == 1:
        db = WgetDbFromBox(WFromBox)
    elif option == 2:
        db = WgetDbFromGrid(WFromGrid)
    elif option == 3:
        db = WgetDbFromCSV(WFromCSV)
    elif option == 4:
        db = WgetDbFromNF(WFromNF)
    else:
        db = None

    if db is not None:
        _saveNF(db, "myDb.NF")
        _displayItem(db)

    return db


def WdefineDbFromGrid(nvarmax=1, nxdef=10, seed=14543):
    WNX = mo.ui.number(start=1, stop=100, value=nxdef, label="NX")
    WNY = mo.ui.number(start=1, stop=100, value=nxdef, label="NY")
    WDX = mo.ui.number(start=1, stop=None, value=1, label="DX")
    WDY = mo.ui.number(start=1, stop=None, value=1, label="DY")
    WX0 = mo.ui.number(start=0, stop=None, value=0, label="X0")
    WY0 = mo.ui.number(start=0, stop=None, value=0, label="Y0")
    WNvar = mo.ui.number(start=1, stop=nvarmax, value=1, label="Number of Variables")
    WPerc = mo.ui.number(start=0, stop=100, value=10, label="Random Displacement")
    WSeed = mo.ui.number(start=None, stop=None, value=seed, label="Seed")

    return mo.ui.array([WNX, WNY, WDX, WDY, WX0, WY0, WNvar, WPerc, WSeed])


def WgetDbFromGrid(WAll):
    [WNX, WNY, WDX, WDY, WX0, WY0, WNvar, WPerc, WSeed] = WAll
    grid = gl.DbGrid.create(
        nx=[WNX.value, WNY.value], dx=[WDX.value, WDY.value], x0=[WX0.value, WY0.value]
    )
    db = gl.Db.createFromGridRandomized(grid, randperc=WPerc.value)
    db.addColumnsRandom(WNvar.value, "z", seed=WSeed.value)
    return db


def WdefineDbFromBox(
    nech=100, nvarmax=1, xmin=0, ymin=0, xmax=100, ymax=100, seed=14543
):
    WNech = mo.ui.number(start=1, stop=None, value=nech, label="Number of Samples")
    WNvar = mo.ui.number(start=1, stop=nvarmax, value=1, label="Number of Variables")
    WXmin = mo.ui.number(start=None, stop=None, value=xmin, label="Minimum along X")
    WYmin = mo.ui.number(start=None, stop=None, value=ymin, label="Minimum along Y")
    WXmax = mo.ui.number(start=None, stop=None, value=xmax, label="Maximum along X")
    WYmax = mo.ui.number(start=None, stop=None, value=ymax, label="Maximum along Y")
    WSeed = mo.ui.number(start=None, stop=None, value=seed, label="Seed")

    return mo.ui.array([WNech, WNvar, WXmin, WYmin, WXmax, WYmax, WSeed])


def WgetDbFromBox(WAll):
    [WNech, WNvar, WXmin, WYmin, WXmax, WYmax, WSeed] = WAll
    return gl.Db.createFillRandom(
        ndat=WNech.value,
        ndim=2,
        nvar=WNvar.value,
        coormin=[WXmin.value, WYmin.value],
        coormax=[WXmax.value, WYmax.value],
        seed=WSeed.value,
    )


def WgetDbFromNF(WAll):
    [WFile] = WAll
    filename = WFile.name()
    if filename is None:
        return None
    return gl.Db.createFromNF(str(WFile.path(index=0)))


def WdefineDbFromCSV(
    nameX="Longitude", nameY="Latitude", nameVar="pH", flagEnglishStyle=True
):
    WnameX = mo.ui.text(label="X Coordinate", value=nameX)
    WnameY = mo.ui.text(label="Y Coordinate", value=nameY)
    WnameVar = mo.ui.text(label="Variable Name", value=nameVar)
    WengStyle = mo.ui.checkbox(label="English Style", value=flagEnglishStyle)
    WFile = mo.ui.file_browser(
        label="Select a CSV File", multiple=False, filetypes=[".csv"]
    )
    return mo.ui.array([WnameX, WnameY, WnameVar, WengStyle, WFile])


def WgetDbFromCSV(WAll, flagHeader=True):
    [WnameX, WnameY, WnameVar, WengStyle, WFile] = WAll

    filename = WFile.name()
    if filename is None:
        return None
    path = WFile.path(index=0)
    if WengStyle.value:
        charSep = ","
        charDec = "."
    else:
        charSep = ";"
        charDec = ","
    dataframe = pd.read_csv(
        path, sep=charSep, decimal=charDec, header=0 if flagHeader else None
    )
    db = gl.Db_fromPandas(dataframe)
    if db.getNSample() <= 0:
        print("Reading of CSV file failed: Check its Style")
    else:
        db.setLocators([WnameX.value, WnameY.value], gl.ELoc.X)
        db.setLocator(WnameVar.value, gl.ELoc.Z)
    return db


# ===============================================
# Widget to manage a Box based on an optional Db
# ===============================================


def WdefineBox(db=None):
    """
    Returns parameters for defining a Box (by meshes only)
    db: Database (optional) for providing default values
    """
    if db is not None:
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

    WLongMin = mo.ui.number(start=None, stop=None, value=longmin)
    WLongMax = mo.ui.number(start=None, stop=None, value=longmax)
    WLatMin = mo.ui.number(start=None, stop=None, value=latmin)
    WLatMax = mo.ui.number(start=None, stop=None, value=latmax)
    WFlagProj = mo.ui.checkbox(
        label="Background (if coordinates are Long/Lat)", value=False
    )

    return mo.ui.array([WLongMin, WLongMax, WLatMin, WLatMax, WFlagProj])


def WshowBox(WAll, flagTitle=True, gapv=0, gaph=1):
    [WLongMin, WLongMax, WLatMin, WLatMax, WFlagProj] = WAll

    WTitle = _WgetTitle("Box Definition", flagTitle)
    Wgrid = mo.hstack(
        [
            mo.vstack(
                [mo.md("Parameters"), mo.md("Minimum"), mo.md("Maximum")], gap=gapv
            ),
            mo.vstack([mo.md("Longitude"), WLongMin, WLongMax], align="end", gap=gapv),
            mo.vstack([mo.md("Latitude"), WLatMin, WLatMax], align="end", gap=gapv),
        ],
        gap=gaph,
    )
    return mo.vstack([WTitle, Wgrid, WFlagProj], gap=gapv)


def WgetBox(WAll):
    [WLongMin, WLongMax, WLatMin, WLatMax, WFlagProj] = WAll
    box = np.ndarray(shape=(2, 2))
    box[0, 0] = WLongMin.value
    box[1, 0] = WLongMax.value
    box[0, 1] = WLatMin.value
    box[1, 1] = WLatMax.value
    return box, WFlagProj.value


# =======================================
# Widget to manage a discretization Grid
# =======================================


def WdefineGridN(nxdef=50):
    """
    Returns parameters for defining a discretization Grid
    nxdef: Number of grid meshes (same along X and Y)
    """
    WNX = mo.ui.number(start=1, stop=None, value=nxdef, label="Nodes along X")
    WNY = mo.ui.number(start=1, stop=None, value=nxdef, label="Nodes along Y")
    return mo.ui.array([WNX, WNY])


def WshowGridN(WAll, flagTitle=True):
    [WNX, WNY] = WAll

    WTitle = _WgetTitle("Grid Discretization", flagTitle=flagTitle)

    return mo.vstack([WTitle, WNX, WNY])


def WgetGridN(WAll, box):
    [WNX, WNY] = WAll

    nx = WNX.value
    ny = WNY.value
    deltax = box[0, 1] - box[0, 0]
    deltay = box[1, 1] - box[1, 0]
    dx = deltax / (nx - 1)
    dy = deltay / (ny - 1)
    x0 = box[0, 0]
    y0 = box[1, 0]
    return gl.DbGrid.create(nx=[nx, ny], dx=[dx, dy], x0=[x0, y0])
