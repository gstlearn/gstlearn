import marimo

__generated_with = "0.11.25"
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
    return ctx, gl, gmo, gp, mo, np, pd, plt


@app.cell(hide_code=True)
def _():
    # Global parameters
    nxdef = 100
    return (nxdef,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetDb = gmo.WdefineDb()
    return (WidgetDb,)


@app.cell(hide_code=True)
def _(WidgetDb, gmo):
    db = gmo.WgetDb(WidgetDb)
    return (db,)


@app.cell(hide_code=True)
def _(db, gmo):
    WidgetZoom = gmo.WdefineBox(db)
    return (WidgetZoom,)


@app.cell(hide_code=True)
def _(gmo, nxdef):
    WidgetGrid = gmo.WdefineGridN(nxdef)
    return (WidgetGrid,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetVario = gmo.WdefineVario(nlag=10, dlag=0.005)
    return (WidgetVario,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetCovList = gmo.WdefineCovList()
    return (WidgetCovList,)


@app.cell(hide_code=True)
def _(WidgetCovList, WidgetDb, WidgetGrid, WidgetVario, ctx, gl, gmo, gp, mo):
    def plotVario(ax, vario, model, showPairs=True):
        ax.varmod(vario, model, showPairs=showPairs)

    def plotData(ax, db, box, targetName, flagProj=False):
        ax.baseMap(db=db, box=box, flagProj=flagProj)
        ax.literal(db=db, name=targetName, fontsize=6)
        if flagProj:
            ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326")
        ax.decoration(title=targetName + " (long/lat)")

    def plotEstim(ax, db, grid, targetName, flagProj=False):
        ax.raster(dbgrid=grid, name="Kriging.*.estim", alpha=0.5)
        ax.literal(db=db, name=targetName, fontsize=6)
        if flagProj:
            ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326")
        ax.decoration(title="Estimation")

    def plotStdev(ax, db, grid, targetName, flagProj=False):
        ax.raster(dbgrid=grid, name="Kriging.*.stdev", alpha=0.5)
        ax.literal(db=db, name=targetName, fontsize=6)
        if flagProj:
            ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik,
                            crs="EPSG:4326")
        ax.decoration(title="St. Dev. of Estimation Error")

    def myaction():

        # Define the Input Db
        db = gmo.WgetDb(WidgetDb)
        if db is None:
            return None

        # Check that the Db is 2-D and contains the variable of interest
        targetName = "pH"
        db.setLocators(["longitude", "latitude"], gl.ELoc.X)
        db.setLocator("sel", gl.ELoc.SEL) # Set the selection (if present in CSV)
        db.setLocator(targetName, gl.ELoc.Z)
        ndim = db.getNLoc(gl.ELoc.X)
        if ndim!= 2 or db.getColIdx(targetName) < 0:
            print("The 'db' should be 2-D (",ndim,") and contain the target variable")
            return None

        # Define the output Grid
        box = db.getExtremas()
        grid = gmo.WgetGridN(WidgetGrid, box)

        # Define the Variogram parameters
        vario = gmo.WgetVario(WidgetVario, db)

        # Define the Model
        model = gmo.WgetCovList(WidgetCovList, vario)

        # Define Neighborhood (Unique)
        neigh = gl.NeighUnique.create()

        # Perform the Estimation
        err = gl.kriging(db, grid, model, neigh)

        fig, ax = gp.init(2, 2, figsize=(10,10))
        plotData(ax[0,0], db, box, targetName)
        plotVario(ax[0,1], vario, model, showPairs=True)

        plotEstim(ax[1,0], db, grid, targetName)
        plotStdev(ax[1,1], db, grid, targetName)
        mo.mpl.interactive(fig)

        return fig
    return myaction, plotData, plotEstim, plotStdev, plotVario


@app.cell(hide_code=True)
def _(
    WidgetCovList,
    WidgetDb,
    WidgetGrid,
    WidgetVario,
    WidgetZoom,
    gmo,
    mo,
    myaction,
):
    param = mo.ui.tabs(
        {
            "Data":       gmo.WshowDb(WidgetDb),
            "Zoom":       gmo.WshowBox(WidgetZoom),
            "Grid":       gmo.WshowGridN(WidgetGrid),
            "Variogram":  gmo.WshowVario(WidgetVario),
            "Model":      gmo.WshowCovList(WidgetCovList),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack(
        [
             mo.md(""),
             mo.md(f"Data and its Estimation{mo.as_html(myaction())}")
        ],
        gap = 4
    )

    mo.hstack([param, simu], gap=4)
    return param, simu


if __name__ == "__main__":
    app.run()
