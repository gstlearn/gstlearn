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
# This part is intended to distribute (as is) a set of functions written in Python
# to be included in Marimo interface
# Reminder: all methods staring by "W" are dedicated to UI


from curses.panel import panel

from curses.panel import panel
from pyexpat import model

from shapely import node

import gstlearn as gl
import gstlearn.plot as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import marimo as mo
import contextily as ctx
import os
import pathlib

optionGlobalDisplay = False
optionGlobalBackup = True
optionGlobalRadix = "My"


def setEnvironment(optionBackup=True, optionDisplay=False, optionRadix="My"):
    """
    Use this function to set the options for the gstmarimo environment

    optionBackup: Triggers saving newly created objects as a Neutral File
    optionDisplay: Triggers displaying newly created objects
    optionRadix: Default prefix for the Neutral File name
    """
    global optionGlobalDisplay
    global optionGlobalBackup
    global optionGlobalRadix
    optionGlobalDisplay = optionDisplay
    optionGlobalBackup = optionBackup
    optionGlobalRadix = optionRadix


# Ensure no output directory is set for gstlearn
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
    Turns gray (as if locked if the 'condition' is met)
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


def _saveAndDisplay(contents=None, filename="File.NF", flagForceDisplay=False):
    if contents is None:
        return

    if optionGlobalBackup:
        local_filename = f"{optionGlobalRadix}{filename}"
        contents.dumpToNF(local_filename)

    if optionGlobalDisplay or flagForceDisplay:
        contents.display()


# =========================
# Widget to manage messages
# =========================


def WdefineMessage():
    get_message, set_message = mo.state("")
    return get_message, set_message


def WshowMessage(
    WMessages,
    border=False,
    minHeight="100px",
    width="500px",
):
    get_message, _ = WMessages

    style = {
        "padding": "8px",
        "minHeight": minHeight,
        "width": width,
    }

    if border:
        style["border"] = "1px solid #ccc"

    return mo.md(get_message()).style(style)


def WsetMessage(WMessages, *args, sep=" ", newline=True):
    _, set_message = WMessages

    message = sep.join(str(arg) for arg in args)
    if newline:
        message += "<br>"

    set_message(message)


def WaddMessage(WMessages, *args, sep=" ", newline=True):
    get_message, set_message = WMessages

    message = sep.join(str(arg) for arg in args)
    if newline:
        message += "<br>"

    set_message(get_message() + message)


def WclearMessage(WMessages):
    _, set_message = WMessages
    set_message("")


# ================================================================
# Widget to manage a Covariance (prefix = WCov)
# The covariance is classified in a list of 'ncovmax' covariances
# ================================================================


def WdefineOneCovariance(
    ic=0, ncovmax=1, ncovdef=-1, distmax=100, varmax=100, typeRef="Spherical"
):
    """
    Returns the widget for inquiring the parameters for a single Basic structure
    ncovmax: Maximum number of Basic structures
    ncovdef: Number of covariances used by default (if < 0: it is set to 'ncovmax')
    distmax: Maximum distance
    varmax:  Maximum Variance value
    typeRef: Type of covariance used as default
    """
    if ncovdef < 0:
        ncovdef = ncovmax
    distRef = distmax * (ic + 1) / (ncovdef + 1)
    distAux = distRef
    varRef = varmax / ncovdef
    angRef = 0

    WCovUsed = mo.ui.switch(True, label="Basic Structure Used")
    WCovType = mo.ui.dropdown(
        options=_getCovarianceDict(), value=typeRef, label="Structure"
    )
    WCovRange = mo.ui.number(start=None, stop=None, value=distRef, label="Range")
    WCovSill = mo.ui.number(start=0, stop=None, value=varRef, label="Sill")
    WCovAniso = mo.ui.switch(value=False, label="Anisotropy")
    WCovRange2 = mo.ui.number(start=0, stop=None, value=distAux, label="Range Aux.")
    WCovAngle = mo.ui.number(start=0, stop=None, value=angRef, label="Angle")
    return mo.ui.array(
        [WCovUsed, WCovType, WCovRange, WCovSill, WCovAniso, WCovRange2, WCovAngle]
    )


def WshowOneCovariance(WAll, flagTitle=True):
    [WCovUsed, WCovType, WCovRange, WCovSill, WCovAniso, WCovRange2, WCovAngle] = WAll

    WCovTitle = _WgetTitle("Covariance Definition", flagTitle)

    WCovTypeupd = _WLock(WCovType, not WCovUsed.value)
    WCovRangeupd = _WLock(WCovRange, not WCovUsed.value)
    WCovSillupd = _WLock(WCovSill, not WCovUsed.value)
    WCovAnisoupd = _WLock(WCovAniso, not WCovUsed.value)
    WCovRange2upd = _WLock(WCovRange2, not WCovUsed.value or not WCovAniso.value)
    WCovAngleupd = _WLock(WCovAngle, not WCovUsed.value or not WCovAniso.value)
    return mo.ui.array(
        [
            WCovTitle,
            WCovUsed,
            WCovTypeupd,
            WCovRangeupd,
            WCovSillupd,
            WCovAnisoupd,
            WCovRange2upd,
            WCovAngleupd,
        ]
    )


def WgetOneCovariance(WAll):
    [WCovUsed, WCovType, WCovRange, WCovSill, WCovAniso, WCovRange2, WCovAngle] = WAll

    if WCovUsed.value:
        if not WCovAniso.value:
            # Isotropic covariance
            return gl.CovAniso.createIsotropic(
                ctxt=gl.CovContext(1, 2),
                type=gl.ECov.fromKey(WCovType.value),
                range=WCovRange.value,
                sill=WCovSill.value,
                param=1.0,
                flagRange=True,
            )
        else:
            # Anisotropic covariance
            return gl.CovAniso.createAnisotropic(
                ctxt=gl.CovContext(1, 2),
                type=gl.ECov.fromKey(WCovType.value),
                ranges=[WCovRange.value, WCovRange2.value],
                sill=WCovSill.value,
                param=1.0,
                angles=[WCovAngle.value, 0.0],
                flagRange=True,
            )
    else:
        return None


# ========================================
# Widget to manage the list of Covariances
# ========================================


def WdefineCovariances(ncovmax=1, ncovdef=-1, distmax=100, varmax=100):
    """
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    ncovmax: Maximum number of Basic structures
    ncovdef: Number of covariances used by default (if < 0: it is set to 'ncovmax')
    distmax: Maximum distance
    varmax:  Maximum Variance value
    """
    return mo.ui.array(
        [
            WdefineOneCovariance(ic, ncovmax, ncovdef, distmax, varmax)
            for ic in range(ncovmax)
        ]
    )


def WshowCovariances(WAll, flagTitle=True):
    WMTitle = _WgetTitle("Model Definition", flagTitle)
    UI = mo.accordion(
        {
            f"Covariance {ic + 1}": WshowOneCovariance(cov, False)
            for ic, cov in enumerate(WAll)
        }
    )
    return mo.ui.array([WMTitle, UI])


def WgetCovariances(WAll):
    model = gl.Model()
    for cov in WAll:
        cova = WgetOneCovariance(cov)
        if cova is not None:
            model.addCov(cova)
    return model


# ==============================================================
# Widget to manage the list of Basic Structures used for Fitting
# ==============================================================


def WshowBasicList(basic_list, flagTitle=True):
    WTitle = _WgetTitle("Basic Structures for Fitting", flagTitle)
    WList = basic_list["types"]
    return mo.ui.array([WTitle, WList])


# ========================
# Widget to manage a Model
# ========================


def WdefineModel(
    ncovmax=1,
    ncovdef=-1,
    distmax=100,
    varmax=100,
    vario=None,
    deftypes=["Spherical"],
    valdef="Fit",
):
    """
    Returns the array of widgets for inquiring a series of 'ncovmax' basic structures
    ncovmax: Maximum number of Basic structures (used for defaulting range)
    ncovdef: Number of covariances used by default (if < 0: it is set to 'ncovmax')
    distmax: Maximum distance
    varmax:  Maximum Variance value
    vario: Vario used for providing default values (if provided)
    valdef: Defaulted option for Model definition
    """
    if vario is not None:
        distmax = vario.getMaximumDistance()
        varmax = vario.getVar()

    WMChoice = mo.ui.radio(
        options={"Interactive": 1, "Fit": 2, "From NF": 3}, value=valdef
    )
    WInter = WdefineCovariances(
        ncovmax=ncovmax, ncovdef=ncovdef, distmax=distmax, varmax=varmax
    )
    WFitVario = WdefineModelFitVario(deftypes=deftypes)
    WMFromNF = WdefineModelFromNF()

    return mo.ui.array([WMChoice, WInter, WFitVario, WMFromNF])


def WshowModel(WAll, flagTitle=True, gapv=2):
    [WMChoice, WInter, WFitVario, WMFromNF] = WAll
    WMTitle = _WgetTitle("Model Definition", flagTitle)
    option = WMChoice.value

    # Contenu à afficher selon le choix
    if option == 1:
        UI = mo.accordion(
            {f"Covariance {ic + 1}": mo.vstack(WInter[ic]) for ic in range(len(WInter))}
        )
        return mo.vstack([WMTitle, WMChoice, UI], gap=gapv)
    elif option == 2:
        return mo.vstack([WMTitle, WMChoice, *WFitVario], gap=gapv)
    elif option == 3:
        return mo.vstack([WMTitle, WMChoice, *WMFromNF], gap=gapv)
    else:
        return mo.md("Invalid option selected")


def WgetModel(WAll, vario=None):
    [WMChoice, WInter, WFitVario, WMFromNF] = WAll

    option = WMChoice.value

    model = None
    if option == 1:
        model = WgetCovariances(WInter)
    elif option == 2:
        model = WgetModelFitVario(WFitVario, vario)
    elif option == 3:
        model = WgetModelFromNF(WMFromNF)
    else:
        return None

    if model is not None:
        # Add the Universality Condition (always)
        model.setDriftIRF(order=0, nfex=0)

    _saveAndDisplay(model, "Model.NF")

    return model


def WdefineModelFromNF():
    WMFile = mo.ui.file_browser(
        label="Select a 'Model' Neutral File",
        multiple=False,
    )
    # Add filetypes=[".NF", ".ascii"]
    # if you want to filter only NF or ascii files (extension)
    return mo.ui.array([WMFile])


def WgetModelFromNF(WAll):
    [WMFile] = WAll
    filename = WMFile.name()
    if filename is None:
        return None
    return gl.Model.createFromNF(str(WMFile.path(index=0)))


def WdefineModelFitVario(deftypes=["Spherical"]):
    WMTypes = mo.ui.multiselect(options=_getCovarianceDict(), value=deftypes)
    return mo.ui.array([WMTypes])


def WgetModelFitVario(WAll, vario):
    [WMTypes] = WAll
    if vario is None:
        return None

    types = WMTypes.value
    if not types:
        return None

    model = gl.Model.createFromVario(vario, gl.ECov.fromKeys(types))
    return model


# =======================
# Widget to manage a Grid
# =======================


def WdefineGrid(nxdef=100, dxdef=1):
    """
    Returns parameters for a regular 2-D grid
    nxdef: Number of grid meshes (same along X and Y)
    """
    WGridNX = mo.ui.number(start=1, stop=200, value=nxdef)
    WGridNY = mo.ui.number(start=1, stop=200, value=nxdef)
    WGridDX = mo.ui.number(start=1, stop=None, value=dxdef)
    WGridDY = mo.ui.number(start=1, stop=None, value=dxdef)
    WGridX0 = mo.ui.number(start=0, stop=None, value=0)
    WGridY0 = mo.ui.number(start=0, stop=None, value=0)
    return mo.ui.array([WGridNX, WGridNY, WGridDX, WGridDY, WGridX0, WGridY0])


def WshowGrid(WAll, flagTitle=True, gapv=0, gaph=1):
    [WGridNX, WGridNY, WGridDX, WGridDY, WGridX0, WGridY0] = WAll
    WTitle = _WgetTitle("Grid Definition", flagTitle)
    Wgrid = mo.hstack(
        [
            mo.vstack(
                [mo.md("Parameters"), mo.md("Nodes"), mo.md("Mesh"), mo.md("Origin")],
                gap=gapv,
            ),
            mo.vstack(
                [mo.md("along X"), WGridNX, WGridDX, WGridX0], align="end", gap=gapv
            ),
            mo.vstack(
                [mo.md("along Y"), WGridNY, WGridDY, WGridY0], align="end", gap=gapv
            ),
        ],
        gap=gaph,
    )
    return mo.vstack([WTitle, Wgrid], gap=gapv)


def WgetGrid(WAll):
    [WGridNX, WGridNY, WGridDX, WGridDY, WGridX0, WGridY0] = WAll
    grid = gl.DbGrid.create(
        nx=[WGridNX.value, WGridNY.value],
        dx=[WGridDX.value, WGridDY.value],
        x0=[WGridX0.value, WGridY0.value],
    )
    return grid


# ============================
# Widget to manage Simulations
# ============================


def WdefineSimtub(nbtuba=100, nbsimu=1, seed=13134):
    """
    Returns parameters for performing Turning Bands simulations
    nbtuba: Number of Turning Bands
    nbsimu: Number of Simulations
    seed: Seed for random number generator
    """
    WSimNbtuba = mo.ui.number(
        start=1, stop=None, value=nbtuba, label="Number of Turning Bands"
    )
    WSimNbsimu = mo.ui.number(
        start=1, stop=None, value=nbsimu, label="Number of Simulations"
    )
    WSimSeed = mo.ui.number(start=0, stop=None, value=seed, label="Seed")

    WDisplayBinary = mo.ui.switch(False, label="Display in Binary Mode")

    return mo.ui.array([WSimNbtuba, WSimNbsimu, WSimSeed, WDisplayBinary])


def WshowSimtub(WAll, flagTitle=True, gapv=2):
    [WSimNbtuba, WSimNbsimu, WSimSeed, WDisplayBinary] = WAll
    WSimTitle = _WgetTitle("Turning Bands Method", flagTitle)
    widgets = [
        WSimTitle,
        WSimNbtuba,
        WSimNbsimu,
        WSimSeed,
        WDisplayBinary,
    ]

    return mo.vstack(widgets, gap=gapv)


def WgetSimtub(WAll):
    [WSimNbtuba, WSimNbsimu, WSimSeed, WDisplayBinary] = WAll
    return (
        WSimNbtuba.value,
        WSimNbsimu.value,
        WSimSeed.value,
        WDisplayBinary.value,
    )


# ========================
# Widget to manage a Vario
# ========================


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

    WVChoice = mo.ui.radio(options={"Omni": 1, "Multi": 2, "From NF": 3}, value=valdef)
    WVOmni = WdefineVarioParamOmni(nlag=nlag, dlag=dlag)
    WVMulti = WdefineVarioParamMulti(ndir=ndir, nlag=nlag, dlag=dlag)
    WVFromNF = WdefineVarioFromNF()
    return mo.ui.array([WVChoice, WVOmni, WVMulti, WVFromNF])


def WshowVario(WAll, flagTitle=True):
    [WVChoice, WVOmni, WVMulti, WVFromNF] = WAll

    WVTitle = _WgetTitle("Variogram Parameters", flagTitle)
    option = WVChoice.value

    # Sélection du contenu selon le choix
    if option == 1:
        return mo.vstack([WVTitle, WVChoice, *WVOmni])
    elif option == 2:
        return mo.vstack([WVTitle, WVChoice, *WVMulti])
    elif option == 3:
        return mo.vstack([WVTitle, WVChoice, *WVFromNF])
    else:
        return mo.md("Invalid selection")


def WgetVario(WAll, db=None):
    [WVChoice, WVOmni, WVMulti, WVFromNF] = WAll
    option = WVChoice.value

    varioparam = None
    if option == 1:
        varioparam = WgetVarioParamOmni(WVOmni)
    elif option == 2:
        varioparam = WgetVarioParamMulti(WVMulti)
    elif option == 3:
        return WgetVarioFromNF(WVFromNF)
    else:
        return None

    vario = None
    if varioparam is not None and db is not None:
        vario = gl.Vario.computeFromDb(
            varioparam, db, calculType=gl.ECalcVario.VARIOGRAM, verbose=True
        )

    _saveAndDisplay(vario, "Vario.NF")

    return vario


def WdefineVarioParamOmni(nlag=10, dlag=1):
    WVOmniNlag = mo.ui.number(start=1, stop=100, value=nlag, label="Number of Lags")
    WVOmniDlag = mo.ui.number(start=0, stop=100, value=dlag, label="Lag Value")
    WVOmniToldis = mo.ui.number(
        start=0, stop=1, value=0.5, label="Tolerance on Distance"
    )
    WVOmniCylrad = mo.ui.number(start=0, stop=None, value=0, label="Cylinder Radius")
    return mo.ui.array([WVOmniNlag, WVOmniDlag, WVOmniToldis, WVOmniCylrad])


def WshowVarioParamOmni(WAll, flagTitle=True):
    [WVOmniNlag, WVOmniDlag, WVOmniToldis, WVOmniCylrad] = WAll

    WVOmniTitle = _WgetTitle("Variogram Parameters", flagTitle)
    return mo.ui.array(
        [WVOmniTitle, WVOmniNlag, WVOmniDlag, WVOmniToldis, WVOmniCylrad]
    )


def WgetVarioParamOmni(WAll):
    [WVOmniNlag, WVOmniDlag, WVOmniToldis, WVOmniCylrad] = WAll
    if WVOmniCylrad.value > 0:
        varioparam = gl.VarioParam.createOmniDirection(
            nlag=WVOmniNlag.value,
            dlag=WVOmniDlag.value,
            toldis=WVOmniToldis.value,
            cylrad=WVOmniCylrad.value,
        )
    else:
        varioparam = gl.VarioParam.createOmniDirection(
            nlag=WVOmniNlag.value, dlag=WVOmniDlag.value, toldis=WVOmniToldis.value
        )

    return varioparam


def WdefineVarioParamMulti(ndir=4, nlag=10, dlag=1):
    WVMultiNdir = mo.ui.number(
        start=1, stop=10, value=ndir, label="Number of Directions"
    )
    WVMultiNlag = mo.ui.number(start=1, stop=100, value=nlag, label="Number of Lags")
    WVMultiDlag = mo.ui.number(start=0, stop=100, value=dlag, label="Lag Value")
    WVMultiAngref = mo.ui.number(
        start=0, stop=180, value=0.0, label="Reference angle (degree)"
    )
    WVMultiToldis = mo.ui.number(
        start=0, stop=1, value=0.5, label="Tolerance on Distance"
    )
    return mo.ui.array(
        [WVMultiNdir, WVMultiNlag, WVMultiDlag, WVMultiAngref, WVMultiToldis]
    )


def WshowVarioParamMulti(WAll, flagTitle=True):
    [WVMultiNdir, WVMultiNlag, WVMultiDlag, WVMultiAngref, WVMultiToldis] = WAll

    WVMultiTitle = _WgetTitle("Variogram Definition", flagTitle)

    return mo.ui.array(
        [
            WVMultiTitle,
            WVMultiNdir,
            WVMultiNlag,
            WVMultiDlag,
            WVMultiAngref,
            WVMultiToldis,
        ]
    )


def WgetVarioParamMulti(WAll):
    [WVMultiNdir, WVMultiNlag, WVMultiDlag, WVMultiAngref, WVMultiToldis] = WAll
    return gl.VarioParam.createMultiple(
        ndir=WVMultiNdir.value,
        nlag=WVMultiNlag.value,
        dlag=WVMultiDlag.value,
        toldis=WVMultiToldis.value,
        angref=WVMultiAngref.value,
    )


def WdefineVarioFromNF():
    WVFile = mo.ui.file_browser(
        label="Select a 'Vario' Neutral File",
        multiple=False,
    )
    # Add filetypes=[".NF", ".ascii"]
    # if you want to filter only NF or ascii files (extension)
    return mo.ui.array([WVFile])


def WgetVarioFromNF(WAll):
    [WVFile] = WAll
    filename = WVFile.name()
    if filename is None:
        return None
    return gl.Vario.createFromNF(str(WVFile.path(index=0)))


# =====================
# Widget to manage a Db
# =====================


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

    WDChoice = mo.ui.radio(
        options={"From Box": 1, "From Grid": 2, "From CSV": 3, "From NF": 4},
        value=valdef,
    )
    WDFromBox = WdefineDbFromBox(nech, nvarmax, xmin, ymin, xmax, ymax, seed)
    WDFromGrid = WdefineDbFromGrid(nvarmax, nxdef, seed)
    WDFromCSV = WdefineDbFromCSV()
    WDFromNF = WdefineDbFromNF()

    return mo.ui.array([WDChoice, WDFromBox, WDFromGrid, WDFromCSV, WDFromNF])


def WshowDb(WAll, flagTitle=True):
    [WDChoice, WDFromBox, WDFromGrid, WDFromCSV, WDFromNF] = WAll

    WDTitle = _WgetTitle("Data Base Parameters", flagTitle)
    option = WDChoice.value

    if option == 1:
        return mo.vstack([WDTitle, WDChoice, *WDFromBox])
    elif option == 2:
        return mo.vstack([WDTitle, WDChoice, *WDFromGrid])
    elif option == 3:
        return mo.vstack([WDTitle, WDChoice, *WDFromCSV])
    elif option == 4:
        return mo.vstack([WDTitle, WDChoice, *WDFromNF])
    else:
        return None


def WgetDb(WAll):
    [WDChoice, WDFromBox, WDFromGrid, WDFromCSV, WDFromNF] = WAll
    option = WDChoice.value

    if option == 1:
        db = WgetDbFromBox(WDFromBox)
    elif option == 2:
        db = WgetDbFromGrid(WDFromGrid)
    elif option == 3:
        db = WgetDbFromCSV(WDFromCSV)
    elif option == 4:
        db = WgetDbFromNF(WDFromNF)
    else:
        db = None

    _saveAndDisplay(db, "Db.NF")

    return db


# ===========================
# Widget to manage a DbIsoPot
# ===========================


def createDefaultIsoPot(nech=3, seed=13424):
    """
    Create the default database used to define the IsoPotential.
    """
    dbIso = gl.Db.createFillRandom(
        ndat=nech,
        ndim=2,
        nvar=0,
        coormin=[0.0, 0.0],
        coormax=[100.0, 100.0],
        seed=seed,
        flagAddSampleRank=False,
    )
    # Adding the Layer information
    dbIso.addColumnsByConstant(1, 1.0, "Layer", gl.ELoc.LAYER)
    dbIso.setName("x-1", "x")
    dbIso.setName("x-2", "y")
    return dbIso


def WdefineIsoPot(dbIsoPot):
    return WdefineEdit(dbIsoPot)


def WshowIsoPot(WAll, flagTitle=True, gapv=0):
    if WAll is None:
        return None

    comments = mo.md(
        """
        Define the iso-potential values by editing the table below.
        - You can add / remove / modify iso-potential points as needed.
        - Samples of the same IsoPotential must have the same **Layer**.
        """
    )

    return mo.vstack(
        [
            comments,
            WshowEdit(WAll, flagTitle=flagTitle, gapv=gapv),
        ],
        gap=gapv,
    )


def WgetIsoPot(WAll, dbIsoPot):
    return WgetEdit(WAll, dbIsoPot)


def createDefaultGradPot(nech=1, seed=13424):
    """
    Create the default database used to define the Gradient for Potential.
    """
    dbGrad = gl.Db.createFillRandom(
        ndat=nech,
        ndim=2,
        nvar=0,
        coormin=[0.0, 0.0],
        coormax=[100.0, 100.0],
        seed=seed,
        flagAddSampleRank=False,
    )
    # Add the gradient coordinates
    dbGrad.addColumnsByConstant(1, 1.0, "gx", gl.ELoc.G, 0)
    dbGrad.addColumnsByConstant(1, 0.0, "gy", gl.ELoc.G, 1)

    dbGrad.setName("x-1", "x")
    dbGrad.setName("x-2", "y")
    return dbGrad


def WdefineGradPot(dbGradPot):
    return WdefineEdit(dbGradPot)


def WshowGradPot(WAll, flagTitle=True, gapv=0):
    if WAll is None:
        return None

    comments = mo.md(
        """
        Define the gradient values by editing the table below.
        - You can add / remove / modify Gradient Points as needed.
        """
    )

    return mo.vstack(
        [
            comments,
            WshowEdit(WAll, flagTitle=flagTitle, gapv=gapv),
        ],
        gap=gapv,
    )


def WgetGradPot(WAll, dbGradPot):
    return WgetEdit(WAll, dbGradPot)


def createDefaultTgtePot(nech=1, seed=42683):
    """
    Create the default database used to define the Tangent for Potential.
    """
    dbTgte = gl.Db.createFillRandom(
        ndat=nech,
        ndim=2,
        nvar=0,
        coormin=[0.0, 0.0],
        coormax=[100.0, 100.0],
        seed=seed,
        flagAddSampleRank=False,
    )
    # Add the tangent coordinates
    dbTgte.addColumnsByConstant(1, 1.0, "tx", gl.ELoc.TGT, 0)
    dbTgte.addColumnsByConstant(1, 0.0, "ty", gl.ELoc.TGT, 1)

    dbTgte.setName("x-1", "x")
    dbTgte.setName("x-2", "y")
    return dbTgte


def WdefineTgtePot(dbTgtePot, flagUseTangent=False):
    flagTgtePot = mo.ui.checkbox(
        label="Use tangent points",
        value=flagUseTangent,
    )

    WAll = WdefineEdit(dbTgtePot)

    return flagTgtePot, WAll


def WshowTgtePot(flagTgtePot, WAll, flagTitle=True, gapv=0):
    comments = mo.md(
        """
        Define the Tangent values by editing the table below.
        - You can add / remove / modify Tangent Points as needed.
        """
    )

    return mo.vstack(
        [
            flagTgtePot,
            comments,
            WshowEdit(WAll, flagTitle=flagTitle, gapv=gapv),
        ],
        gap=gapv,
    )


def WgetTgtePot(flagTgtePot, WAll, dbTgtePot):
    if not flagTgtePot.value:
        return None

    return WgetEdit(WAll, dbTgtePot)


def WdefineModelPot(
    type="Cubic",
    isotropic=True,
    range=60,
    rangeU=60,
    rangeV=60,
    angle=0.0,
    driftOrder=0,
):
    """
    Returns parameters for defining a Potential Model

    Parameters
    ----------
    type : str
        Covariance type ("Cubic", "Gaussian", "Spline")
    isotropic : bool
        True for isotropic model, False for anisotropic model
    range : float
        Main range value
    rangeU : float
        Anisotropic range along U direction
    rangeV : float
        Anisotropic range along V direction
    angle : float
        Anisotropy angle
    driftOrder : int
        Drift order (0, 1 or 2)
    """

    WType = mo.ui.radio(
        options={
            "Cubic": "Cubic",
            "Gaussian": "Gaussian",
            "Spline": "Spline",
        },
        value=type,
    )

    WIso = mo.ui.radio(
        options={
            "Isotropic": True,
            "Anisotropic": False,
        },
        value="Isotropic" if isotropic else "Anisotropic",
    )

    WRange = mo.ui.number(
        start=0.0,
        value=range,
    )

    WRangeU = mo.ui.number(
        start=0.0,
        value=rangeU,
    )

    WRangeV = mo.ui.number(
        start=0.0,
        value=rangeV,
    )

    WAngle = mo.ui.number(
        value=angle,
    )

    WDriftOrder = mo.ui.dropdown(
        options={
            "Stationary": 0,
            "Order 1": 1,
            "Order 2": 2,
        },
        value={
            0: "Stationary",
            1: "Order 1",
            2: "Order 2",
        }[driftOrder],
    )

    return mo.ui.array(
        [
            WType,
            WIso,
            WRange,
            WRangeU,
            WRangeV,
            WAngle,
            WDriftOrder,
        ]
    )


def WshowModelPot(WAll, flagTitle=True, gapv=0, gaph=1):
    [
        WType,
        WIso,
        WRange,
        WRangeU,
        WRangeV,
        WAngle,
        WDriftOrder,
    ] = WAll

    WModelTitle = _WgetTitle("Model Definition", flagTitle)

    def _line(label, widget):
        return mo.hstack(
            [
                mo.md(f"<div style='width:120px'>{label}</div>"),
                widget,
            ],
            justify="start",
            gap=gaph,
        )

    WModelGrid = [
        _line("Type", WType),
        _line("Geometry", WIso),
        _line("Range", WRange),
    ]

    if not WIso.value:
        WModelGrid.extend(
            [
                _line("Range U", WRangeU),
                _line("Range V", WRangeV),
                _line("Angle", WAngle),
            ]
        )

    WModelGrid.append(_line("Drift Order", WDriftOrder))

    return mo.vstack(
        [
            WModelTitle,
            mo.vstack(
                WModelGrid,
                gap=gapv,
            ),
        ],
        gap=gapv,
    )


def WgetModelPot(WAll):
    [
        WType,
        WIso,
        WRange,
        WRangeU,
        WRangeV,
        WAngle,
        WDriftOrder,
    ] = WAll

    type = WType.value
    if type == "Cubic":
        type = gl.ECov.CUBIC
    elif type == "Gaussian":
        type = gl.ECov.GAUSSIAN
    elif type == "Spline":
        type = gl.ECov.SPLINE2_GC
    else:
        raise ValueError(f"Unknown covariance type: {type}")

    isotropic = WIso.value

    modelRef = None
    if isotropic:
        range = WRange.value
        modelRef = gl.Model.createFromParam(type, range=range)
    else:
        rangeU = WRangeU.value
        rangeV = WRangeV.value
        angle = WAngle.value
        modelRef = gl.Model.createFromParam(
            type, ranges=[rangeU, rangeV], angles=[angle, 0.0]
        )

    covp = gl.CovGradientAnalytic(modelRef.getCovAniso(0))

    model = gl.Model()
    model.setCov(covp)
    model.setDriftIRF(order=int(WDriftOrder.value), nfex=0)

    _saveAndDisplay(model, "Model.NF")

    return model


# ===============================
# Widget to manage a Regular Grid
# ===============================


def WdefineDbFromGrid(nvarmax=1, nxdef=10, seed=14543):
    WDbGridNX = mo.ui.number(start=1, stop=100, value=nxdef, label="NX")
    WDbGridNY = mo.ui.number(start=1, stop=100, value=nxdef, label="NY")
    WDbGridDX = mo.ui.number(start=1, stop=None, value=1, label="DX")
    WDbGridDY = mo.ui.number(start=1, stop=None, value=1, label="DY")
    WDbGridX0 = mo.ui.number(start=0, stop=None, value=0, label="X0")
    WDbGridY0 = mo.ui.number(start=0, stop=None, value=0, label="Y0")
    WDbGridNvar = mo.ui.number(
        start=1, stop=nvarmax, value=1, label="Number of Variables"
    )
    WDbGridPerc = mo.ui.number(start=0, stop=100, value=10, label="Random Displacement")
    WDbGridSeed = mo.ui.number(start=None, stop=None, value=seed, label="Seed")

    return mo.ui.array(
        [
            WDbGridNX,
            WDbGridNY,
            WDbGridDX,
            WDbGridDY,
            WDbGridX0,
            WDbGridY0,
            WDbGridNvar,
            WDbGridPerc,
            WDbGridSeed,
        ]
    )


def WgetDbFromGrid(WAll):
    [
        WDbGridNX,
        WDbGridNY,
        WDbGridDX,
        WDbGridDY,
        WDbGridX0,
        WDbGridY0,
        WDbGridNvar,
        WDbGridPerc,
        WDbGridSeed,
    ] = WAll

    grid = gl.DbGrid.create(
        nx=[WDbGridNX.value, WDbGridNY.value],
        dx=[WDbGridDX.value, WDbGridDY.value],
        x0=[WDbGridX0.value, WDbGridY0.value],
    )
    db = gl.Db.createFromGridRandomized(grid, randperc=WDbGridPerc.value)
    db.addColumnsRandom(WDbGridNvar.value, "z", seed=WDbGridSeed.value)
    return db


def WdefineDbFromBox(
    nech=100, nvarmax=1, xmin=0, ymin=0, xmax=100, ymax=100, seed=14543
):
    WDBoxNech = mo.ui.number(start=1, stop=None, value=nech, label="Number of Samples")
    WDBoxNvar = mo.ui.number(
        start=1, stop=nvarmax, value=1, label="Number of Variables"
    )
    WDBoxXmin = mo.ui.number(start=None, stop=None, value=xmin, label="Minimum along X")
    WDBoxYmin = mo.ui.number(start=None, stop=None, value=ymin, label="Minimum along Y")
    WDBoxXmax = mo.ui.number(start=None, stop=None, value=xmax, label="Maximum along X")
    WDBoxYmax = mo.ui.number(start=None, stop=None, value=ymax, label="Maximum along Y")
    WDBoxSeed = mo.ui.number(start=None, stop=None, value=seed, label="Seed")

    return mo.ui.array(
        [WDBoxNech, WDBoxNvar, WDBoxXmin, WDBoxYmin, WDBoxXmax, WDBoxYmax, WDBoxSeed]
    )


def WgetDbFromBox(WAll):
    [WDBoxNech, WDBoxNvar, WDBoxXmin, WDBoxYmin, WDBoxXmax, WDBoxYmax, WDBoxSeed] = WAll
    return gl.Db.createFillRandom(
        ndat=WDBoxNech.value,
        ndim=2,
        nvar=WDBoxNvar.value,
        coormin=[WDBoxXmin.value, WDBoxYmin.value],
        coormax=[WDBoxXmax.value, WDBoxYmax.value],
        seed=WDBoxSeed.value,
    )


def WdefineDbFromNF():
    WDFile = mo.ui.file_browser(label="Select a 'Db' Neutral File", multiple=False)
    # Add filetypes=[".NF", ".ascii"]
    # if you want to filter only NF or ascii files (extension)
    return mo.ui.array([WDFile])


def WgetDbFromNF(WAll):
    [WDFile] = WAll
    filename = WDFile.name()
    if filename is None:
        return None
    return gl.Db.createFromNF(str(WDFile.path(index=0)))


def WdefineDbFromCSV(
    nameX="Longitude", nameY="Latitude", nameVar="pH", flagEnglishStyle=True
):
    WDCSVnameX = mo.ui.text(label="X Coordinate", value=nameX)
    WDCSVnameY = mo.ui.text(label="Y Coordinate", value=nameY)
    WDCSVnameVar = mo.ui.text(label="Variable Name", value=nameVar)
    WDCSVengStyle = mo.ui.checkbox(label="English Style", value=flagEnglishStyle)
    WDCSVFile = mo.ui.file_browser(label="Select a CSV File", multiple=False)
    # Add filetypes=[".csv"] if you want to filter only CSV files (extension)
    return mo.ui.array([WDCSVnameX, WDCSVnameY, WDCSVnameVar, WDCSVengStyle, WDCSVFile])


def WgetDbFromCSV(WAll, flagHeader=True):
    [WDCSVnameX, WDCSVnameY, WDCSVnameVar, WDCSVengStyle, WDCSVFile] = WAll
    filename = WDCSVFile.name()
    if filename is None:
        return None
    path = WDCSVFile.path(index=0)
    if WDCSVengStyle.value:
        charSep = ","
        charDec = "."
    else:
        charSep = ";"
        charDec = ","
    dataframe = pd.read_csv(
        path,
        sep=charSep,
        decimal=charDec,
        header=0 if flagHeader else None,
        on_bad_lines="warn",
    )
    db = gl.Db_fromPandas(dataframe)
    if db.getNSample() <= 0:
        print("Reading of CSV file failed: Check its Style")
        db = None
    else:
        db.setLocators([WDCSVnameX.value, WDCSVnameY.value], gl.ELoc.X)
        db.setLocator(WDCSVnameVar.value, gl.ELoc.Z)
    return db


# ==============================================
# Widget to manage a Box based on an optional Db
# ==============================================


def WdefineBox(db=None):
    """
    Returns parameters for defining a Box (by meshes only)
    db: Database (optional) for providing default values
    """
    if db is not None:
        box = db.getExtremas()
        longmin = box[0][0]
        longmax = box[0][1]
        latmin = box[1][0]
        latmax = box[1][1]
    else:
        longmin = -180
        longmax = 180
        latmin = -90
        latmax = 90

    WBoxLongMin = mo.ui.number(start=None, stop=None, value=longmin)
    WBoxLongMax = mo.ui.number(start=None, stop=None, value=longmax)
    WBoxLatMin = mo.ui.number(start=None, stop=None, value=latmin)
    WBoxLatMax = mo.ui.number(start=None, stop=None, value=latmax)
    WBoxFlagBackground = mo.ui.checkbox(
        label="Background (if coordinates are Long/Lat)", value=False
    )

    return mo.ui.array(
        [WBoxLongMin, WBoxLongMax, WBoxLatMin, WBoxLatMax, WBoxFlagBackground]
    )


def WshowBox(WAll, flagTitle=True, gapv=0, gaph=1):
    [WBoxLongMin, WBoxLongMax, WBoxLatMin, WBoxLatMax, WBoxFlagBackground] = WAll

    WBoxTitle = _WgetTitle("Box Definition", flagTitle)
    WBoxgrid = mo.hstack(
        [
            mo.vstack(
                [mo.md("Parameters"), mo.md("Minimum"), mo.md("Maximum")], gap=gapv
            ),
            mo.vstack(
                [mo.md("Longitude"), WBoxLongMin, WBoxLongMax], align="end", gap=gapv
            ),
            mo.vstack(
                [mo.md("Latitude"), WBoxLatMin, WBoxLatMax], align="end", gap=gapv
            ),
        ],
        gap=gaph,
    )
    return mo.vstack([WBoxTitle, WBoxgrid, WBoxFlagBackground], gap=gapv)


def WgetBox(WAll):
    [WBoxLongMin, WBoxLongMax, WBoxLatMin, WBoxLatMax, WBoxFlagBackground] = WAll
    box = np.ndarray(shape=(2, 2))
    box[0, 0] = WBoxLongMin.value
    box[0, 1] = WBoxLongMax.value
    box[1, 0] = WBoxLatMin.value
    box[1, 1] = WBoxLatMax.value
    return box, WBoxFlagBackground.value


# ======================================
# Widget to manage a discretization Grid
# ======================================


def WdefineGridN(nxdef=50):
    """
    Returns parameters for defining a discretization Grid
    nxdef: Number of grid meshes (same along X and Y)
    """
    WGridNNX = mo.ui.number(start=1, stop=None, value=nxdef)
    WGridNNY = mo.ui.number(start=1, stop=None, value=nxdef)
    return mo.ui.array([WGridNNX, WGridNNY])


def WshowGridN(WAll, flagTitle=True, gapv=0, gaph=1):
    [WGridNNX, WGridNNY] = WAll
    WGridNTitle = _WgetTitle("Grid Discretization", flagTitle=flagTitle)

    WGridNgrid = mo.hstack(
        [
            mo.vstack([mo.md("Parameters"), mo.md("Nodes")], gap=gapv),
            mo.vstack([mo.md("Along X"), WGridNNX], align="end", gap=gapv),
            mo.vstack([mo.md("Along Y"), WGridNNY], align="end", gap=gapv),
        ],
        gap=gaph,
    )
    return mo.vstack([WGridNTitle, WGridNgrid], gap=gapv)


def WgetGridN(WAll, box):
    [WGridNNX, WGridNNY] = WAll

    nx = WGridNNX.value
    ny = WGridNNY.value
    deltax = box[0, 1] - box[0, 0]
    deltay = box[1, 1] - box[1, 0]
    dx = deltax / (nx - 1)
    dy = deltay / (ny - 1)
    x0 = box[0, 0]
    y0 = box[1, 0]
    return gl.DbGrid.create(nx=[nx, ny], dx=[dx, dy], x0=[x0, y0])


# ==================================
# Widget to manage a Editing of a Db
# ==================================


def WdefineEdit(db):
    if db is None:
        return
    names = db.getAllNames()
    cols = db.getColumnsAsVVD(names)
    df = pd.DataFrame(cols.T, columns=names).round(3)
    Wedit = mo.ui.data_editor(df).form(bordered=False)
    return mo.ui.array([Wedit])


def WshowEdit(WAll, flagTitle=True, gapv=1, gaph=1):
    if WAll is None:
        return None

    [Wedit] = WAll

    WeditTitle = _WgetTitle("Edit Database", flagTitle=flagTitle)

    return mo.vstack(
        [WeditTitle, Wedit],
        gap=gapv,
    )


def WgetEdit(WAll, db):
    if WAll is None:
        return None
    [Wedit] = WAll
    if Wedit.value is None:
        return db

    df = Wedit.value

    if df.shape[0] != db.getNSample():
        db.resizeSamples(df.shape[0])

    names = db.getAllNames()

    for name in names:
        if name in df.columns:
            db.setColumn(df[name], name)

    _saveAndDisplay(db, "Db.NF")

    return db


# =======================
# Widget to manage a Rule
# =======================


def WdefineRule(ruleDef=None, maxFacies=10):
    if ruleDef is None:
        ruleDef = ["S", "T", "F1", "F2", "F3"]

    WRule = mo.ui.text(
        value=", ".join(ruleDef),
        label="Rule",
    )

    WProps = mo.ui.array(
        [
            mo.ui.number(
                start=0.0, stop=None, value=1.0, label=f"Proportion of Facies {i + 1}"
            )
            for i in range(maxFacies)
        ]
    )

    return mo.ui.array([WRule, WProps])


def WshowRule(WAll, flagTitle=True, gapv=0):
    [WRule, WProps] = WAll
    WRuleTitle = _WgetTitle("Rule and Proportions", flagTitle=flagTitle)

    return mo.vstack(
        [
            WRuleTitle,
            WRule,
            mo.md(
                """
                **Proportions**

                You must define the Proportions of each facies.
                Only those included in the Rule will be considered.
                They are normalized internally so that they sum to 1.
                """
            ),
            *WProps,
        ],
        gap=gapv,
    )


def WgetRule(WAll):
    [WRule, WProps] = WAll

    ruleDef = [s.strip() for s in WRule.value.split(",") if s.strip()]
    rule = gl.Rule.createFromNames(ruleDef)

    nfacies = rule.getNFacies() if rule is not None else 1

    raw_props = np.array([float(WProps[i].value) for i in range(nfacies)])

    s = raw_props.sum()

    if s == 0:
        props = np.ones(nfacies) / nfacies
    else:
        props = raw_props / s

    _saveAndDisplay(rule, "Rule.NF")

    ruleprop = gl.RuleProp.createFromRule(rule, props)
    return ruleprop


# =======================
# Widget to manage Layout
# =======================

_LAYOUT_OPTIONS = {
    "data": {"label": "Display Data", "default": True},
    "rule": {"label": "Display Lithotype Rule", "default": True},
    "model": {"label": "Display Model(s)", "default": True},
    "estimation": {"label": "Display Estimation Map(s)", "default": True},
    "stdev": {"label": "Display Standard Deviation Map(s)", "default": True},
    "simulation": {"label": "Display Simulation Map(s)", "default": True},
    "average": {
        "label": "Display Simulation Average and Dispersion Map(s)",
        "default": True,
    },
    "KWeights": {"label": "Display Kriging Weights", "default": True},
    "XWeights": {"label": "Display Cross-Validation", "default": True},
}


def WdefineLayout(nrow=3, ncol=3, width=5, height=5, defaults=None):
    widgets = [
        mo.ui.number(start=1, value=nrow),
        mo.ui.number(start=1, value=ncol),
        mo.ui.number(start=1, value=width),
        mo.ui.number(start=1, value=height),
    ]
    defaults = defaults or {}

    for key, opt in _LAYOUT_OPTIONS.items():
        widgets.append(
            mo.ui.switch(
                defaults.get(key, opt["default"]),
                label=opt["label"],
            )
        )

    return mo.ui.array(widgets)


def WshowLayout(WAll, flagTitle=True, gapv=0, gaph=1):
    WTitle = _WgetTitle("Graphic Layout", flagTitle)

    WLayoutNrow, WLayoutNcol, WLayoutWidth, WLayoutHeight = WAll[:4]
    switches = WAll[4:]

    # --- 2x2 grid for geometry parameters
    geometry = mo.hstack(
        [
            mo.vstack(
                [mo.md("Layout"), mo.md("Along X"), mo.md("Along Y")],
                gap=gapv,
            ),
            mo.vstack(
                [mo.md("Number of cells"), WLayoutNcol, WLayoutNrow],
                align="end",
                gap=gapv,
            ),
            mo.vstack(
                [mo.md("Figure size"), WLayoutWidth, WLayoutHeight],
                align="end",
                gap=gapv,
            ),
        ],
        gap=gaph,
    )

    # --- switches block
    options = mo.vstack(
        switches,
        gap=0.3,
    )

    return mo.vstack(
        [
            WTitle,
            mo.md("### Geometry"),
            geometry,
            mo.md("### Options"),
            options,
        ],
        gap=gapv,
    )


def WgetLayout(WAll, nvar=1, nbsimu=1, ngrf=1, valid=()):
    WLayoutNrow, WLayoutNcol, WLayoutWidth, WLayoutHeight = WAll[:4]
    switches = WAll[4:]

    layout = {
        "nx": WLayoutNrow.value,
        "ny": WLayoutNcol.value,
        "dimx": WLayoutWidth.value,
        "dimy": WLayoutHeight.value,
    }

    # ✔️ map full UI state
    for key, widget in zip(_LAYOUT_OPTIONS.keys(), switches):
        layout[key] = widget.value

    # ✔️ CRITICAL: valid defines API capability
    render_plan = []

    def add_if_allowed(key, count):
        if key in valid and layout.get(key):
            render_plan.append((key, count))

    add_if_allowed("data", nvar)
    add_if_allowed("rule", 1)
    add_if_allowed("model", ngrf)
    add_if_allowed("estimation", nvar)
    add_if_allowed("stdev", nvar)
    add_if_allowed("simulation", nbsimu * nvar)
    add_if_allowed("average", 2 * nvar)
    add_if_allowed("KWeights", 1)
    add_if_allowed("XWeights", 1)

    layout["render_plan"] = render_plan

    return layout


# =========================
# Widget to manage AutoSave
# =========================


def WdefineAutoSave(
    autosave=False, autoDisplay=False, rootname="My", directory=str(pathlib.Path.cwd())
):
    """
    Returns parameters for the AutoSave
    autosave: Whether to save automatically
    autoDisplay: Whether to display automatically
    rootname: Root name of the file
    directory: Directory where to save
    """
    return {
        "autosave": mo.ui.checkbox(
            value=autosave,
            label="Automatic Backup of item as Neutral File",
        ),
        "autoDisplay": mo.ui.checkbox(
            value=autoDisplay,
            label="Automatic Display of item",
        ),
        "rootname": mo.ui.text(
            value=rootname,
            label="Radix of the Neutral File name",
        ),
        "directory": mo.ui.file_browser(
            initial_path=directory,
            multiple=False,
            label="Directory of backup",
            filetypes=[],
        ),
    }


def WshowAutoSave(panel):
    return mo.vstack(
        [
            mo.md(
                """
                **Automatic Backup**

                When this option is enabled, a Neutral File will be created
                automatically after any Item is modified in the interface.
                The file will be saved in the selected directory with the specified root name.
                """,
            ),
            panel["autosave"],
            panel["autoDisplay"],
            panel["rootname"],
            panel["directory"],
        ]
    )


def WgetAutoSave(panel):
    selected = panel["directory"].value

    directory = selected[0].path if selected else str(pathlib.Path.cwd())

    global optionGlobalBackup
    global optionGlobalDisplay
    global optionGlobalRadix
    optionGlobalBackup = panel["autosave"].value
    optionGlobalDisplay = panel["autoDisplay"].value
    optionGlobalRadix = panel["rootname"].value

    return {
        "autosave": panel["autosave"].value,
        "autoDisplay": panel["autoDisplay"].value,
        "rootname": panel["rootname"].value,
        "directory": directory,
    }


# =============================
# Widget to manage Neighborhood
# =============================


def WdefineNeigh(nmaxi=10, radius=10.0, flagUnique=False):
    """
    Returns parameters for the Neighborhood
    """
    Wunique = mo.ui.checkbox(label="Unique Neighborhood", value=flagUnique)
    Wnmaxi = mo.ui.number(start=1, value=nmaxi, label="Maximum number of samples")
    Wradius = mo.ui.number(value=radius, step=0.1, label="Search radius")
    Wnmini = mo.ui.number(start=1, value=1, label="Minimum number of samples")
    Wnsect = mo.ui.number(start=1, value=1, label="Number of sectors")
    Wnsmax = mo.ui.number(value=10, label="Maximum samples per sector")
    Wangle = mo.ui.number(value=0, label="Rotation Angle (degree)")

    return mo.ui.array(
        [
            Wunique,
            Wnmaxi,
            Wradius,
            Wnmini,
            Wnsect,
            Wnsmax,
            Wangle,
        ]
    )


def WshowNeigh(WAll, flagTitle=True, gapv=1):
    [Wunique, Wnmaxi, Wradius, Wnmini, Wnsect, Wnsmax, Wangle] = WAll
    WNeighTitle = _WgetTitle("Neighborhood", flagTitle)
    # mode UNIQUE → on masque tout sauf titre + switch
    if Wunique.value:
        widgets = [
            WNeighTitle,
            Wunique,
        ]
    else:
        widgets = [
            WNeighTitle,
            Wunique,
            Wnmaxi,
            Wradius,
            Wnmini,
            Wnsect,
            Wnsmax,
            Wangle,
        ]
    return mo.vstack(widgets, gap=gapv)


def WgetNeigh(WAll):
    [Wunique, Wnmaxi, Wradius, Wnmini, Wnsect, Wnsmax, Wangle] = WAll
    if Wunique.value:
        return gl.NeighUnique()

    return gl.NeighMoving.create(
        flag_xvalid=False,
        nmaxi=Wnmaxi.value,
        radius=Wradius.value,
        nmini=Wnmini.value,
        nsect=Wnsect.value,
        nsmax=Wnsmax.value,
        angles=[Wangle.value, 0.0],
    )


# =====================================
# Some display functions used in Marimo
# =====================================


def plotData(
    ax,
    db,
    name=None,
    box=None,
    title=None,
    flagLiteral=True,
    flagCst=False,
    flagProj=False,
    flagBackground=False,
    flagTitle=True,
    s=20,
    c="r",
):
    if db is None:
        return None

    if box is not None:
        ax.baseMap(db=db, box=box, flagProj=flagProj)

    if flagLiteral:
        ax.literal(db=db, name=name, fontsize=6, c=c)
    else:
        ax.symbol(db=db, nameSize=name, flagCst=flagCst, s=s, c=c)

    if flagBackground:
        ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:4326")
    if title is None:
        title = name
    if flagTitle and title is not None:
        ax.decoration(title=title)
    ax.geometry(aspect=1)


def plotVario(ax, vario=None, model=None, title=None, showPairs=True):
    ax.varmod(vario, model, showPairs=showPairs)
    if title is None:
        title = "Variogram & Model"
    if title is not None:
        ax.decoration(title=title)


# Plot a grid with optional isolines
# * ax: Marimo axis
# * grid: DbGrid to plot
# * name: Name of the variable to plot (must be in the DbGrid)
# * title: Title of the plot (optional)
# * flagLegend: Whether to display the legend (default: False)
# * flagBinary: Whether to use a binary color map (default: False)
# * nlevel: Number of isoline levels (default: 10, set to 0 for no isolines or flagBinary=True)
# * levels: Optional list of specific levels for isolines (overrides nlevel if provided)
# * flagTitle: Whether to display the title (default: True)
# * levelColor: Color of the isoline levels (default: "black")
# * levelWidth: Width of the isoline levels (default: 0.5)
def plotGrid(
    ax,
    grid,
    name=None,
    rule=None,
    title=None,
    flagLegend=False,
    flagBinary=False,
    nlevel=10,
    levels=None,
    flagTitle=True,
    levelColor="black",
    levelWidth=0.5,
):
    if grid is None:
        return
    if name is None or grid.getColIdx(name) <= 0:
        return
    if title is None:
        title = f"Grid: {name}"
    if flagBinary:
        flagLegend = False

    ax.raster(
        dbgrid=grid,
        name=name,
        rule=rule,
        flagLegend=flagLegend,
        flagBinary=flagBinary,
    )
    if nlevel > 0 and not flagBinary:
        ax.isoline(
            dbgrid=grid,
            name=name,
            nlevel=nlevel,
            levels=levels,
            colors=levelColor,
            linewidths=levelWidth,
            linestyles="solid",
        )

    if flagTitle and title is not None:
        ax.decoration(title=title)
    ax.geometry(aspect=1)


# Plot a Lithotype Rule
# * ax: Marimo axis
# * rule: Rule to plot
# * title: Title of the plot (optional)
# * flagLegend: Whether to display the legend (default: False)
def plotRule(
    ax,
    rule,
    title=None,
    flagLegend=False,
):
    if rule is None:
        return
    ax.rule(
        ruleobj=rule,
        flagLegend=flagLegend,
    )

    if title is not None:
        ax.decoration(title=title)
    ax.geometry(aspect=1)


# Plot a Neighborhood
# * ax: Marimo axis
# * neigh: Neighborhood to plot
# * grid: DbGrid to plot the neighborhood on
# * node: Rank of the node where centering the neighborhood (default: 0)
# * title: Title of the plot (optional)
# * flagCell: Whether to display the cell (default: False)
def plotNeigh(
    ax,
    neigh,
    grid,
    node=0,
    flagCell=False,
    title=None,
):
    if neigh is None:
        return
    ax.neigh(neighobj=neigh, grid=grid, node=node, flagCell=flagCell)

    if title is not None:
        ax.decoration(title=title)
    ax.geometry(aspect=1)


# Plot the Kriging weights
# * ax: Marimo axis
# * target: Data Base containing the target points (used for picking)
# * db: Db which is displayed as a background and also contains the weights
# * node: Rank of the node where centering the neighborhood (default: -1)
# * title: Title of the plot (optional)
def plotWeights(
    ax,
    target,
    db,
    nameW,
    neigh,
    node=-1,
    flagCell=False,
    title=None,
):
    # Display the target points without any decoration (used for picking)
    plotData(ax, target, flagLiteral=False, flagCst=True, s=2, c="red")

    # Display the data points without any decoration (used for picking)
    plotData(ax, db, flagLiteral=False, flagCst=True, s=10, c="blue")
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    if node >= 0 and node < target.getNSample():
        # Display the Weights on the data points involved in the Neighborhood
        plotData(ax, db, name=nameW, flagLiteral=True, c="green")
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Display the Neighborhood characteristics
        plotNeigh(ax, neigh, target, node, flagCell)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if title is not None:
        ax.decoration(title=title)
    ax.geometry(aspect=1)


def getSelectedIndex(plot_widget, db):
    """
    Returns the index of the first selected sample.

    Parameters
    ----------
    plot_widget : marimo plot widget
        Interactive widget containing the selection.
    db : gl.Db
        Database containing coordinates.

    Returns
    -------
        - -1 if no axis is provided;
        - -1 if no selection is active;
        - the first selected index otherwise.
    """

    if plot_widget is None:
        return -1

    sel = plot_widget.value

    if sel is None:
        return -1

    tabx, taby = gp.getCoordinates(db)

    mask = sel.get_mask(tabx, taby)

    if mask is None:
        return -1

    selected = np.flatnonzero(mask)
    return -1 if selected.size == 0 else int(selected[0])
