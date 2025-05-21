import marimo

__generated_with = "0.11.25"
app = marimo.App()


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import altair as alt

    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import gstlearn.document as gdoc

    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    from IPython.display import Markdown
    return Markdown, alt, copy, gdoc, gl, gmo, gp, mo, np, plt


@app.cell(hide_code=True)
def _():
    # Parametrization for the Model
    ncovmax = 2
    distmax = 100
    varmax  = 50
    return distmax, ncovmax, varmax


@app.cell(hide_code=True)
def _():
    # Version gstlearn a faire marcher
    #Markdown(gdoc.loadDoc("Statistics_mean.md"))
    # version decortiquee qui fonctionne
    #filename = "/home/drenard/project_gstlearn/gstlearn/doc/references/Cvv.md"
    #Markdown(filename)
    return


@app.cell(hide_code=True)
def _(distmax, gmo, ncovmax, varmax):
    WidgetModel = gmo.WdefineModel(ncovmax, distmax, varmax)
    return (WidgetModel,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetGrid = gmo.WdefineGrid(100)
    return (WidgetGrid,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetSimtub = gmo.WdefineSimtub()
    return (WidgetSimtub,)


@app.cell(hide_code=True)
def _(WidgetGrid, WidgetModel, WidgetSimtub, gl, gmo, gp):
    def myaction():

        model = gmo.WgetModel(WidgetModel)
        grid = gmo.WgetGrid(WidgetGrid)
        nbtuba, seed = gmo.WgetSimtub(WidgetSimtub)
        err = gl.simtub(None, dbout=grid, model=model, nbtuba=nbtuba, seed=int(seed))

        fig, ax = gp.init(2,1,figsize=(10,14))
        ax[0,0].model(model, hmax=100)
        ax[1,0].raster(grid)
        ax[1,0].axis("equal")
        return fig
    return (myaction,)


@app.cell(hide_code=True)
def _(WidgetGrid, WidgetModel, WidgetSimtub, gmo, mo, myaction):
    param = mo.ui.tabs(
        {
            "Grid":       gmo.WshowGrid(WidgetGrid),
            "Model":      gmo.WshowModel(WidgetModel),
            "Simulation": gmo.WshowSimtub(WidgetSimtub),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack(
        [
             mo.md(""),
             mo.md(f"Model and Simulation{mo.as_html(myaction())}")
        ],
        gap = 4
    )

    mo.hstack([param, simu], gap=4)
    return param, simu


if __name__ == "__main__":
    app.run()
