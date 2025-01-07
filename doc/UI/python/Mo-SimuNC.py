import marimo

__generated_with = "0.10.9"
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
def _(WidgetModel, gmo):
    gmo.WshowModel(WidgetModel)
    return


@app.cell(hide_code=True)
def _(gmo):
    WidgetGrid = gmo.WdefineGrid(100)
    return (WidgetGrid,)


@app.cell(hide_code=True)
def _(WidgetGrid, gmo):
    gmo.WshowGrid(WidgetGrid)
    return


@app.cell(hide_code=True)
def _(gmo):
    WidgetSimtub = gmo.WdefineSimtub()
    return (WidgetSimtub,)


@app.cell(hide_code=True)
def _(WidgetSimtub, gmo):
    gmo.WshowSimtub(WidgetSimtub)
    return


@app.cell(hide_code=True)
def _(WidgetGrid, WidgetModel, WidgetSimtub, gl, gmo, plt):
    def mareaction():

        model = gmo.WgetModel(WidgetModel)
        grid = gmo.WgetGrid(WidgetGrid)
        nbtuba, seed = gmo.WgetSimtub(WidgetSimtub)
        err = gl.simtub(None, dbout=grid, model=model, nbtuba=nbtuba, seed=int(seed))

        fig = plt.figure(figsize=(20,6))
        ax1 = fig.add_subplot(1,2,1)
        ax1.model(model, hmax=100)
        ax2 = fig.add_subplot(1,2,2)
        ax2.raster(grid)
        ax2.axis("equal")
        return fig
    return (mareaction,)


@app.cell(hide_code=True)
def _(mareaction, mo):
    mo.md(f"A non-conditional simulation corresponding to a Model: {mo.as_html(mareaction())}")
    return


if __name__ == "__main__":
    app.run()
