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
    WidgetGrid = gmo.WdefineGrid()
    WidgetSimtub = gmo.WdefineSimtub(nbsimu=4)
    WidgetLayout = gmo.WdefineLayout(3, 3, 3, 3)
    WidgetAutoSave = gmo.WdefineAutoSave()
    return WidgetAutoSave, WidgetGrid, WidgetLayout, WidgetModel, WidgetSimtub


@app.cell(hide_code=True)
def _(
    WidgetAutoSave,
    WidgetGrid,
    WidgetLayout,
    WidgetModel,
    WidgetSimtub,
    gl,
    gmo,
    gp,
    mo,
):
    def myaction():
        autosave = gmo.WgetAutoSave(WidgetAutoSave)
        grid = gmo.WgetGrid(WidgetGrid)
        model = gmo.WgetModel(WidgetModel)
        nbtuba, nbsimu, seed, flagDisplayBinary = gmo.WgetSimtub(WidgetSimtub)

        if model is not None and grid is not None:
            gl.simtub(
                None, dbout=grid, model=model, nbtuba=nbtuba, nbsimu=nbsimu, seed=seed
            )
            grid.statisticsBySample(
                names=["Simu.*"], opers=[gl.EStatOption.MEAN, gl.EStatOption.STDV]
            )

        layout = gmo.WgetLayout(WidgetLayout, 1, nbsimu, 1)
        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        contents = layout["contents"]
        valid = {"model", "simu", "average"}  # Constraint to available outputs
        contents_local = (c for c in contents if c in valid)

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy])
        axes = ax.ravel()

        isimu = 0
        iaverage = 0
        for axi in axes:
            content = next(contents_local, None)

            if content == "model":
                gmo.plotVario(axi, model=model)

            elif content == "simu":
                isimu += 1
                gmo.plotGrid(
                    axi,
                    grid,
                    name="Simu" if nbsimu == 1 else f"Simu.S{isimu}",
                    title=f"Simulation #{isimu}/{nbsimu}",
                    flagLegend=False,
                    flagBinary=flagDisplayBinary,
                )

            elif content == "average":
                iaverage += 1
                gmo.plotGrid(
                    axi,
                    grid,
                    name="Stats.MEAN" if iaverage == 1 else "Stats.STDV",
                    title="Average" if iaverage == 1 else "Dispersion",
                    flagLegend=False,
                    flagBinary=flagDisplayBinary,
                )

            else:
                axi.axis("off")  # ou skip visuel propre

        fig.tight_layout(pad=0)

        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def _(
    WidgetAutoSave,
    WidgetGrid,
    WidgetLayout,
    WidgetModel,
    WidgetSimtub,
    gmo,
    mo,
    myaction,
):
    param = mo.ui.tabs(
        {
            "Grid": gmo.WshowGrid(WidgetGrid),
            "Model": gmo.WshowModel(WidgetModel),
            "Simulation": gmo.WshowSimtub(WidgetSimtub),
            "Layout": gmo.WshowLayout(WidgetLayout),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style({"minWidth": "400px", "width": "350px"})

    simu = mo.as_html(myaction())

    mo.hstack([param, simu], gap=4)
    return


if __name__ == "__main__":
    app.run()
