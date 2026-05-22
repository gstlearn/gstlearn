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
    WidgetSimtub = gmo.WdefineSimtub(nbsimu=4)
    return (WidgetSimtub,)


@app.cell(hide_code=True)
def _(WidgetGrid, WidgetModel, WidgetSimtub, gl, gmo, gp, mo):
    def myaction():
        grid = gmo.WgetGrid(WidgetGrid)

        model = gmo.WgetModel(WidgetModel)

        nbtuba, nbsimu, seed, flagDisplaySimu, flagDisplayBinary = gmo.WgetSimtub(
            WidgetSimtub
        )

        if model is not None and grid is not None:
            gl.simtub(
                None, dbout=grid, model=model, nbtuba=nbtuba, nbsimu=nbsimu, seed=seed
            )
            grid.statisticsBySample(
                names=["Simu.*"], opers=[gl.EStatOption.MEAN, gl.EStatOption.STDV]
            )

        fig, ax = gp.init(2, 2, figsize=[8, 8])
        fig.suptitle("Model and Simulations", fontsize=15, fontweight="bold")

        # First figure (Model)
        gmo.plotVario(ax[0, 0], model=model)

        # Second Figure (Always first simulation)
        if nbsimu >= 1:
            title = f"Simulation #1/{nbsimu}"
            name = "Simu"
            if nbsimu > 1:
                name = "Simu.S1"
            gmo.plotGrid(
                ax[1, 0],
                grid,
                name=name,
                title=title,
                flagLegend=True,
                flagBinary=flagDisplayBinary,
            )

        # Third figure: Second simulation or Average
        if flagDisplaySimu:
            if nbsimu >= 2:
                title = f"Simulation #2/{nbsimu}"
                gmo.plotGrid(
                    ax[0, 1],
                    grid,
                    name="Simu.S2",
                    title=title,
                    flagLegend=True,
                    flagBinary=flagDisplayBinary,
                )
        else:
            title = "Average of Simulations"
            gmo.plotGrid(
                ax[0, 1],
                grid,
                name="Stats.MEAN",
                title=title,
                flagLegend=True,
                flagBinary=flagDisplayBinary,
            )

        # Fourth figure: Third simulation or Dispersion
        if flagDisplaySimu:
            if nbsimu >= 3:
                title = f"Simulation #3/{nbsimu}"
                gmo.plotGrid(
                    ax[1, 1],
                    grid,
                    name="Simu.S3",
                    title=title,
                    flagLegend=True,
                    flagBinary=flagDisplayBinary,
                )
        else:
            if nbsimu > 1:
                title = "St. Deviation of Simulations"
                gmo.plotGrid(
                    ax[1, 1],
                    grid,
                    name="Stats.STDV",
                    title=title,
                    flagLegend=True,
                    flagBinary=False,
                )
        mo.mpl.interactive(fig)

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

    simu = mo.as_html(myaction())

    mo.hstack([param, simu], gap=4)
    return


if __name__ == "__main__":
    app.run()
