import marimo

__generated_with = "0.23.14"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def my_imports():
    import marimo as mo

    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import matplotlib.pyplot as plt

    import numpy as np
    import pandas as pd

    gmo.setEnvironment(optionBackup=True, optionDisplay=False)
    return gl, gmo, gp, mo, plt


@app.cell(hide_code=True)
def define_widgets(gmo):
    WidgetModel = gmo.WdefineModel(
        ncovmax=2, distmax=100, varmax=50, valdef="Interactive"
    )
    WidgetGrid = gmo.WdefineGrid()
    WidgetSimtub = gmo.WdefineSimtub(nbsimu=4)
    WidgetLayout = gmo.WdefineLayout(3, 3, 3, 3)
    WidgetAutoSave = gmo.WdefineAutoSave()
    return WidgetAutoSave, WidgetGrid, WidgetLayout, WidgetModel, WidgetSimtub


@app.cell(hide_code=True)
def define_action(
    WidgetAutoSave,
    WidgetGrid,
    WidgetLayout,
    WidgetModel,
    WidgetSimtub,
    gl,
    gmo,
    gp,
    mo,
    plt,
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

        layout = gmo.WgetLayout(
            WidgetLayout,
            nvar=1,
            nbsimu=nbsimu,
            ngrf=1,
            valid=["model", "simulation", "average"],
        )
        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy], squeeze=False)
        axes = ax.ravel()

        i = 0
        isimu = 0
        iaverage = 0

        for name, count in layout.get("render_plan", []):
            for k in range(count):
                if i >= len(axes):
                    break

                axi = axes[i]

                if name == "model":
                    gmo.plotVario(axi, model=model)

                elif name == "simulation":
                    isimu += 1
                    gmo.plotGrid(
                        axi,
                        grid,
                        name="Simu" if nbsimu == 1 else f"Simu.S{isimu}",
                        title=f"Simulation #{isimu}/{nbsimu}",
                        flagLegend=False,
                        flagBinary=flagDisplayBinary,
                    )

                elif name == "average":
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
                    axi.axis("off")

                i += 1

        for axi in axes[i:]:
            axi.axis("off")

        plt.tight_layout(pad=0.2)
        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def render_ui(
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
    ).style({"minWidth": "400px", "width": "400px"})

    simu = mo.vstack([mo.md(""), mo.md(f"{mo.as_html(myaction())}")], gap=0)

    layout = mo.hstack([param, simu], gap=4)
    layout
    return


if __name__ == "__main__":
    app.run()
