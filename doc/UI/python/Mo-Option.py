import marimo

__generated_with = "0.19.2"
app = marimo.App(css_file="custom.css")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import matplotlib.pyplot as plt

    import numpy as np
    import pandas as pd

    return gl, gmo, gp, mo


@app.cell(hide_code=True)
def _(gmo):
    gmo.setEnvironment(optionSaveNF=True, optionPrint=False)
    return


@app.cell(hide_code=True)
def _(gmo):
    # Définition du super-widget
    WidgetDb = gmo.WdefineDb()
    return (WidgetDb,)


@app.cell(hide_code=True)
def _(WidgetDb, gmo):
    dbinit = gmo.WgetDb(WidgetDb)
    return (dbinit,)


@app.cell(hide_code=True)
def _(dbinit, gmo):
    WidgetEdit = gmo.WdefineEdit(dbinit)
    return (WidgetEdit,)


@app.cell(hide_code=True)
def _(WidgetEdit, dbinit, gmo):
    db = gmo.WgetEdit(WidgetEdit, dbinit)
    return (db,)


@app.cell(hide_code=True)
def _(db, gmo):
    WidgetView = gmo.WdefineBox(db)
    return (WidgetView,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetGrid = gmo.WdefineGridN()
    return (WidgetGrid,)


@app.cell(hide_code=True)
def _(db, gmo):
    WidgetVario = gmo.WdefineVario(nlag=10, db=db)
    return (WidgetVario,)


@app.cell(hide_code=True)
def _(WidgetVario, db, gmo):
    vario = gmo.WgetVario(WidgetVario, db=db)
    return (vario,)


@app.cell(hide_code=True)
def _(gmo, vario):
    WidgetModel = gmo.WdefineModel(ncovmax=3, vario=vario)
    return (WidgetModel,)


@app.cell(hide_code=True)
def _(
    WidgetDb,
    WidgetGrid,
    WidgetModel,
    WidgetVario,
    WidgetView,
    gl,
    gmo,
    gp,
    mo,
):
    def myaction():
        # Define the Input Db
        db = gmo.WgetDb(WidgetDb)
        if db is None:
            return None

        targetName = db.getNameByLocator(gl.ELoc.Z)
        ndim = db.getNLoc(gl.ELoc.X)
        if ndim != 2:
            print("The 'db' should be 2-D (", ndim, ")")
            return None
        if db.getColIdx(targetName) < 0:
            print("The 'db' should contain the target variable")
            return None

        # Define the output Grid
        box, flagProj = gmo.WgetBox(WidgetView)
        grid = gmo.WgetGridN(WidgetGrid, box)

        # Define the Variogram parameters
        vario = gmo.WgetVario(WidgetVario, db)

        # Define the Model
        model = gmo.WgetModel(WidgetModel, vario)

        # Define Neighborhood (Unique)
        neigh = gl.NeighUnique.create()

        # Perform the Estimation
        if (
            db is not None
            and grid is not None
            and model is not None
            and neigh is not None
        ):
            err = gl.kriging(db, grid, model, neigh)

        fig, ax = gp.init(2, 2, figsize=(8, 8))
        gmo.plotData(ax[0, 0], db, name=targetName, box=box, flagProj=flagProj)
        gmo.plotVario(ax[0, 1], vario, model, showPairs=True)

        gmo.plotGrid(ax[1, 0], grid, name="Kriging.*.estim", flagLegend=True)
        gmo.plotData(ax[1, 0], db, name=targetName, box=box, flagProj=flagProj)
        gmo.plotGrid(ax[1, 1], grid, name="Kriging.*.stdev", flagLegend=True)
        gmo.plotData(ax[1, 1], db, name=targetName, box=box, flagProj=flagProj)
        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def _(
    WidgetDb,
    WidgetEdit,
    WidgetGrid,
    WidgetModel,
    WidgetVario,
    WidgetView,
    gmo,
    mo,
    myaction,
):
    param = mo.ui.tabs(
        {
            "Data": gmo.WshowDb(WidgetDb),
            "Edit": gmo.WshowEdit(WidgetEdit),
            "View": gmo.WshowBox(WidgetView, gapv=0, gaph=1),
            "Grid": gmo.WshowGridN(WidgetGrid),
            "Variogram": gmo.WshowVario(WidgetVario),
            "Model": gmo.WshowModel(WidgetModel),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack(
        [mo.md(""), mo.md(f"Data and its Estimation {mo.as_html(myaction())}")], gap=4
    )

    mo.hstack([param, simu], gap=4)
    return


if __name__ == "__main__":
    app.run()
