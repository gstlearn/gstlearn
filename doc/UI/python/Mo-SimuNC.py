import marimo

__generated_with = "0.19.2"
app = marimo.App(width="full")


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
    WidgetModel = gmo.WdefineModel(
        ncovmax=2, distmax=100, varmax=50, valdef="Interactive"
    )
    return (WidgetModel,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetGrid = gmo.WdefineGrid()
    return (WidgetGrid,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetSimtub = gmo.WdefineSimtub()
    return (WidgetSimtub,)


@app.cell(hide_code=True)
def _(WidgetGrid, WidgetModel, WidgetSimtub, gl, gmo, gp, mo):
    def myaction():
        grid = gmo.WgetGrid(WidgetGrid)

        model = gmo.WgetModel(WidgetModel)

        nbtuba, nbsimu, seed = gmo.WgetSimtub(WidgetSimtub)

        if model is not None and grid is not None:
            gl.simtub(
                None, dbout=grid, model=model, nbtuba=nbtuba, nbsimu=nbsimu, seed=seed
            )

        fig, ax = gp.init(2, 2, figsize=[8, 8])

        gmo.plotVario(ax[0, 0], model=model)
        name = "Simu"
        if nbsimu > 1:
            name = "Simu.1"
        gmo.plotGrid(ax[1, 0], grid, name=name, title="Simulation #1", flagLegend=True)
        if nbsimu >= 2:
            gmo.plotGrid(
                ax[0, 1], grid, name="Simu.2", title="Simulation #2", flagLegend=True
            )
        if nbsimu >= 3:
            gmo.plotGrid(
                ax[1, 1], grid, name="Simu.3", title="Simulation #3", flagLegend=True
            )
        mo.mpl.interactive(fig)

        print(nbsimu)
        return fig

    return (myaction,)


@app.cell(hide_code=True)
def _(WidgetGrid, WidgetModel, WidgetSimtub, gmo, mo, myaction):
    param = mo.ui.tabs(
        {
            "Grid": gmo.WshowGrid(WidgetGrid),
            "Model": gmo.WshowModel(WidgetModel),
            "Simulation": gmo.WshowSimtub(WidgetSimtub),
        }
    ).style({"minWidth": "400px", "width": "350px"})

    simu = mo.vstack(
        [mo.md(""), mo.md(f"Model and Simulation{mo.as_html(myaction())}")], gap=4
    )

    mo.hstack([param, simu], gap=4)
    return


if __name__ == "__main__":
    app.run()
