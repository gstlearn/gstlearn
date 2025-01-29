import marimo

__generated_with = "0.10.9"
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


@app.cell
def _():
    # Global parameters
    nxdef = 100
    return (nxdef,)


@app.cell(hide_code=True)
def _(mo):
    b = mo.ui.file_browser()
    mo.md( f"""
        Choose your CSV file {b}
        """)
    return (b,)


@app.cell(hide_code=True)
def _(b, gl, mo):
    mo.stop(len(b.value) == 0, "You must define a CSV file first")
    filename = b.value[0].path
    csvformat = gl.CSVformat.create(flagHeader=True, charSep=';', charDec=',')
    db = gl.Db.createFromCSV(filename, csvformat)
    db.setLocators(["longitude", "latitude"], gl.ELoc.X)
    db.setLocator("sel", gl.ELoc.SEL) # Set the selection (if present in CSV)

    box = db.getExtremas()

    dbfmt = gl.DbStringFormat.createFromFlags(flag_stats=True)
    db.display()
    return box, csvformat, db, dbfmt, filename


@app.cell(hide_code=True)
def _(box, mo):
    WLongMin = mo.ui.number(start=1, stop=200, value=box[0,0])
    WLongMax = mo.ui.number(start=1, stop=200, value=box[0,1])
    WLatMin  = mo.ui.number(start=1, stop=200, value=box[1,0])
    WLatMax  = mo.ui.number(start=1, stop=200, value=box[1,1])

    mo.vstack([mo.md("Define the Area to be displayed"), 
               mo.hstack([mo.md("Longitude"), WLongMin, WLongMax]),
               mo.hstack([mo.md("Latitude"),  WLatMin,  WLatMax])
              ], align='start')
    return WLatMax, WLatMin, WLongMax, WLongMin


@app.cell(hide_code=True)
def _(WLatMax, WLatMin, WLongMax, WLongMin, box):
    box[0,0] = WLongMin.value
    box[0,1] = WLongMax.value
    box[1,0] = WLatMin.value
    box[1,1] = WLatMax.value
    return


@app.cell(hide_code=True)
def _(box, ctx, db, mo, plt):
    fig1 = plt.figure(figsize=(5,4))

    axref = fig1.add_subplot(1,1,1)
    axref.baseMap(db, box=box, flagProj=True)
    ctx.add_basemap(axref, source=ctx.providers.OpenStreetMap.Mapnik)
    axref.decoration(title="In projected coordinates")
    axref.axis("equal")

    mo.mpl.interactive(fig1)
    return axref, fig1


@app.cell(hide_code=True)
def _(db, gl):
    # Define the target variable
    targetName = "pH"
    db.setLocator(targetName, gl.ELoc.Z)
    return (targetName,)


@app.cell(hide_code=True)
def _(ctx, db, mo, plt, targetName):
    fig2 = plt.figure(figsize=(5,4))

    ax2 = fig2.add_subplot(1,1,1)
    ax2.literal(db, targetName, fontsize=6)
    ctx.add_basemap(ax2, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:4326")
    ax2.decoration(title=targetName + " (long/lat)")
    ax2.axis("equal")

    mo.mpl.interactive(fig2)
    return ax2, fig2


@app.cell(hide_code=True)
def _(mo, nxdef):
    WNX = mo.ui.number(start=1, stop=200, value=nxdef)
    WNY = mo.ui.number(start=1, stop=200, value=nxdef)
    mo.vstack([mo.md("Grid Definition"), WNX, WNY])
    return WNX, WNY


@app.cell(hide_code=True)
def _(WNX, WNY, box, gl):
    deltax = box[0,1] - box[0,0]
    deltay = box[1,1] - box[1,0]
    nx = WNX.value
    ny = WNY.value
    dx = deltax / (nx-1)
    dy = deltay / (ny-1)
    x0 = box[0,0]
    y0 = box[1,0]
    grid = gl.DbGrid.create(nx = [nx,ny], dx = [dx,dy], x0 = [x0, y0])
    return deltax, deltay, dx, dy, grid, nx, ny, x0, y0


@app.cell(hide_code=True)
def _(gmo):
    WidgetVario = gmo.WdefineVario(nlag=10, dlag=0.005)
    return (WidgetVario,)


@app.cell(hide_code=True)
def _(WidgetVario, gmo):
    WVarioLayout = gmo.WshowVario(WidgetVario)
    WVarioLayout
    return (WVarioLayout,)


@app.cell(hide_code=True)
def _(WVarioLayout, WidgetVario, db, gmo):
    vario = gmo.WgetVario(WidgetVario, WVarioLayout.value, db)
    vario
    return (vario,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetCovList = gmo.WdefineCovList()
    return (WidgetCovList,)


@app.cell(hide_code=True)
def _(WidgetCovList, gmo):
    gmo.WshowCovList(WidgetCovList)
    return


@app.cell(hide_code=True)
def _(WidgetCovList, gmo, vario):
    model = gmo.WgetCovList(WidgetCovList, vario)
    return (model,)


@app.cell(hide_code=True)
def _(model, plt, vario):
    fig3 = plt.figure(figsize=(4,3))
    ax3 = fig3.add_subplot(1,1,1)
    ax3.varmod(vario, model, showPairs=True)
    return ax3, fig3


@app.cell(hide_code=True)
def _(db, gl, grid, model):
    # Perform Kriging
    neigh = gl.NeighUnique.create()
    err = gl.kriging(db, grid, model, neigh)
    return err, neigh


@app.cell(hide_code=True)
def _(ctx, db, grid, mo, plt, targetName):
    fig4 = plt.figure(figsize=(7,3))

    ax4a = fig4.add_subplot(1,2,1)
    ax4a.raster(grid, name="Kriging.*.estim", alpha=0.5)
    ax4a.literal(db, targetName, fontsize=6)
    ctx.add_basemap(ax4a, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:4326")
    ax4a.decoration(title="Estimation")
    #ax4a.axis("equal")

    ax4b = fig4.add_subplot(1,2,2)
    ax4b.raster(grid, name="Kriging.*.stdev", alpha=0.5)
    ax4b.literal(db, targetName, fontsize=6)
    ctx.add_basemap(ax4b, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:4326")
    ax4b.decoration(title="St. Dev.")
    #ax4b.axis("equal")

    mo.mpl.interactive(fig4)
    return ax4a, ax4b, fig4


if __name__ == "__main__":
    app.run()
