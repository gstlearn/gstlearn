import marimo

__generated_with = "0.19.2"
app = marimo.App()


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import matplotlib.pyplot as plt
    import contextily as ctx

    import numpy as np
    import pandas as pd

    return ctx, gl, gmo, gp, mo


@app.cell(hide_code=True)
def _():
    # Global parameters
    nxdef = 100
    return (nxdef,)


@app.cell(hide_code=True)
def _(gmo):
    # Définition du super-widget
    WidgetDb = gmo.WdefineDb()
    return (WidgetDb,)


@app.cell(hide_code=True)
def _(WidgetDb, gmo):
    db = gmo.WgetDb(WidgetDb)
    return (db,)


@app.cell(hide_code=True)
def _(db, gmo):
    WidgetView = gmo.WdefineBox(db)
    return (WidgetView,)


@app.cell(hide_code=True)
def _(gmo, nxdef):
    WidgetGrid = gmo.WdefineGridN(nxdef)
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
    ctx,
    gl,
    gmo,
    gp,
    mo,
):
    def plotVario(ax, vario, model, showPairs=True):
        ax.varmod(vario, model, showPairs=showPairs)
        ax.decoration(title="Variogram & Model")

    def plotData(ax, db, box, targetName, flagProj=False):
        ax.baseMap(db=db, box=box, flagProj=flagProj)
        ax.literal(db=db, name=targetName, fontsize=6)
        if flagProj:
            ctx.add_basemap(
                ax, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:4326"
            )
        ax.decoration(title=targetName)
        ax.geometry(aspect=1)

    def plotEstim(ax, db, grid, targetName, flagProj=False):
        ax.raster(dbgrid=grid, name="Kriging.*.estim", alpha=0.5)
        ax.literal(db=db, name=targetName, fontsize=6)
        if flagProj:
            ctx.add_basemap(
                ax, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:4326"
            )
        ax.decoration(title="Estimation")
        ax.geometry(aspect=1)

    def plotStdev(ax, db, grid, targetName, flagProj=False):
        ax.raster(dbgrid=grid, name="Kriging.*.stdev", alpha=0.5)
        ax.literal(db=db, name=targetName, fontsize=6)
        if flagProj:
            ctx.add_basemap(
                ax, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:4326"
            )
        ax.decoration(title="St. Dev. of Estimation Error")
        ax.geometry(aspect=1)

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
        err = gl.kriging(db, grid, model, neigh)

        fig, ax = gp.init(2, 2, figsize=(10, 10))
        plotData(ax[0, 0], db, box, targetName, flagProj=flagProj)
        plotVario(ax[0, 1], vario, model, showPairs=True)

        plotEstim(ax[1, 0], db, grid, targetName)
        plotStdev(ax[1, 1], db, grid, targetName)
        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def _(
    WidgetDb,
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
