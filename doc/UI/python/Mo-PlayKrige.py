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
    WidgetAutoSave = gmo.WdefineAutoSave()
    WidgetDb = gmo.WdefineDb()
    WidgetGrid = gmo.WdefineGrid(nxdef=20, dxdef=5)
    WidgetVario = gmo.WdefineVario(nlag=10)
    WidgetNeigh = gmo.WdefineNeigh(radius=50)
    WidgetWeights = gmo.WdefineWeights()
    WidgetLayout = gmo.WdefineLayout(nrow=2, ncol=3, width=3, height=3)
    return (
        WidgetAutoSave,
        WidgetDb,
        WidgetGrid,
        WidgetLayout,
        WidgetNeigh,
        WidgetVario,
        WidgetWeights,
    )


@app.cell(hide_code=True)
def define_widget_model(WidgetDb, WidgetVario, gmo):
    db0 = gmo.WgetDb(WidgetDb)
    vario = gmo.WgetVario(WidgetVario, db=db0)
    WidgetModel = gmo.WdefineModel(ncovmax=3, vario=vario)
    return (WidgetModel,)


@app.cell(hide_code=True)
def define_action(
    WidgetAutoSave,
    WidgetDb,
    WidgetGrid,
    WidgetLayout,
    WidgetModel,
    WidgetNeigh,
    WidgetVario,
    WidgetWeights,
    gl,
    gmo,
    gp,
    mo,
    plt,
):
    def myaction():
        autosave = gmo.WgetAutoSave(WidgetAutoSave)

        data = gmo.WgetDb(WidgetDb)
        grid = gmo.WgetGrid(WidgetGrid)
        vario = gmo.WgetVario(WidgetVario, data)
        model = gmo.WgetModel(WidgetModel, vario)
        neigh = gmo.WgetNeigh(WidgetNeigh)

        Kindex, Xindex = gmo.WgetWeights(WidgetWeights)

        if (
            data is not None
            and grid is not None
            and model is not None
            and neigh is not None
        ):
            gl.kriging(data, grid, model, neigh)

            if Kindex >= 0 and Kindex < grid.getNSample():
                gl.krigWeights(data, grid, model, neigh, Kindex)

            if Xindex >= 0 and Xindex < data.getNSample():
                gl.xvalidWeights(data, model, neigh, Xindex)

        layout = gmo.WgetLayout(
            WidgetLayout,
            nvar=1,
            nbsimu=0,
            ngrf=1,
            valid=["data", "model", "estimation", "stdev", "KWeights", "XWeights"],
        )
        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy], squeeze=False)
        axes = ax.ravel()

        i = 0
        for name, count in layout.get("render_plan", []):
            for k in range(count):
                if i >= len(axes):
                    break

                axi = axes[i]
                targetName = data.getNameByLocator(gl.ELoc.Z, k)

                if name == "data":
                    gmo.plotData(
                        axi, data, name=targetName, title="Data:" + targetName, c="blue"
                    )

                elif name == "model":
                    gmo.plotVario(axi, vario, model)

                elif name == "estimation":
                    gridName = "Kriging." + targetName + ".estim"
                    gmo.plotGrid(axi, grid, name=gridName, title="Estimation")
                    gmo.plotData(axi, data, name=targetName, flagTitle=False, c="blue")

                elif name == "stdev":
                    gridName = "Kriging." + targetName + ".stdev"
                    gmo.plotGrid(
                        axi, grid, name=gridName, title="St. dev. of Estimation Error"
                    )
                    gmo.plotData(axi, data, name=targetName, flagTitle=False, c="blue")

                elif name == "KWeights":
                    name_kw = "KWeights." + targetName
                    if Kindex >= 0 and Kindex < grid.getNSample():
                        gmo.plotWeights(
                            axi,
                            grid,
                            data,
                            name_kw,
                            neigh,
                            Kindex,
                            title="Kriging Weights",
                        )
                    else:
                        axi.set_title("Kriging Weights")

                elif name == "XWeights":
                    name_xw = "XWeights." + targetName
                    if Xindex >= 0 and Xindex < data.getNSample():
                        gmo.plotWeights(
                            axi,
                            data,
                            data,
                            name_xw,
                            neigh,
                            Xindex,
                            title="Cross-Validation Weights",
                        )
                    else:
                        axi.set_title("Cross-Validation Weights")

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
    WidgetDb,
    WidgetGrid,
    WidgetLayout,
    WidgetModel,
    WidgetNeigh,
    WidgetVario,
    WidgetWeights,
    gmo,
    mo,
    myaction,
):
    param = mo.ui.tabs(
        {
            "Data": gmo.WshowDb(WidgetDb),
            "Grid": gmo.WshowGrid(WidgetGrid),
            "Variogram": gmo.WshowVario(WidgetVario),
            "Model": gmo.WshowModel(WidgetModel),
            "Neigh": gmo.WshowNeigh(WidgetNeigh),
            "Weights": gmo.WshowWeights(WidgetWeights),
            "Layout": gmo.WshowLayout(WidgetLayout),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack([mo.md(""), mo.md(f"{mo.as_html(myaction())}")], gap=0)

    layout = mo.hstack([param, simu], gap=4)
    layout
    return


if __name__ == "__main__":
    app.run()
