import marimo

__generated_with = "0.25.0"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def my_imports():
    import marimo as mo

    import gstlearn as gl
    import gstlearn.gstmarimo as gmo
    import gstlearn.plot as gp
    import matplotlib.pyplot as plt

    import numpy as np
    import pandas as pd

    gmo.setEnvironment(optionBackup=True, optionDisplay=False)
    return gl, gmo, gp, mo, plt


@app.cell(hide_code=True)
def widget_definition(gmo):
    WidgetGrid = gmo.WdefineGrid()
    WidgetSimtub = gmo.WdefineSimtub(nbsimu=1)
    WidgetModel = gmo.WdefineModel(
        ncovmax=2, distmax=100, varmax=50, valdef="Interactive", withFit=False
    )

    options_list = [
        ("simulation", "Display Simulations", True),
        ("vmap", "Display Variogram Map", True),
        ("model", "Display Model(s)", True),
        ("vario", "Display Experimental Variograms", True),
    ]

    WidgetLayout = gmo.WdefineLayout(
        options_list,
        nrow=2,
        ncol=2,
        width=3,
        height=3,
    )
    WidgetAutoSave = gmo.WdefineAutoSave()
    return (
        WidgetAutoSave,
        WidgetGrid,
        WidgetLayout,
        WidgetModel,
        WidgetSimtub,
        options_list,
    )


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
    options_list,
    plt,
):
    def myaction():
        autosave = gmo.WgetAutoSave(WidgetAutoSave)
        grid = gmo.WgetGrid(WidgetGrid)
        model = gmo.WgetModel(WidgetModel)
        nbtuba, nbsimu, seed, flagDisplayBinary = gmo.WgetSimtub(WidgetSimtub)

        n = 50
        ymax = 1.4
        nlag = grid.getNX(0) / 2
        hmax = grid.getDX(0) * grid.getNX(0) / 2

        dbp = gl.DbGrid.create(nx=[1, 1], x0=[n, n])
        db = gl.DbGrid.create(nx=[2 * n, 2 * n])
        gl.simtub(None, db, model, nbtuba=nbtuba, seed=seed)
        db["model"] = 1 - model.evalCovMat(db, dbp).toTL()
        sill = model.getTotalSill(0, 0)
        gmax = sill * ymax if sill > 0 else 1.0

        varioparam = gl.VarioParam.createMultipleFromGrid(db, nlag)
        vario = gl.Vario.computeFromDb(varioparam, db)

        layout = gmo.WgetLayout(
            WidgetLayout,
            options_list,
        )
        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy], squeeze=False)
        axes = ax.ravel()

        i = 0
        isimu = 0

        for content in layout.get("selected", []):
            if content == "vmap":
                if i < len(axes):
                    axes[i].raster(db, "model")
                    axes[i].decoration(title="Model")
                    i += 1

            elif content == "simulation":
                if i < len(axes):
                    axes[i].raster(db, "*Simu")
                    axes[i].decoration(title="Simulation")
                    i += 1

            elif content == "model":
                if i < len(axes):
                    axes[i].model(model, hmax=70, codir=[1, 0], color="red")
                    axes[i].model(model, hmax=70, codir=[0, 1], color="black")
                    axes[i].geometry(xlim=[0, hmax], ylim=[0, gmax])
                    axes[i].decoration(title="Variogram Model")
                    i += 1

            elif content == "vario":
                if i < len(axes):
                    axes[i].variogram(vario, idir=0, color="red")
                    axes[i].variogram(vario, idir=1, color="black")
                    axes[i].geometry(xlim=[0, hmax], ylim=[0, gmax])
                    axes[i].decoration(title="Experimental Variogram")
                    i += 1

        for axi in axes[i:]:
            axi.axis("off")

        plt.tight_layout(pad=0.2)
        mo.mpl.interactive(fig)

        db.deleteColumn("model")
        db.deleteColumn("*Simu")

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
    options_list,
):
    param = mo.ui.tabs(
        {
            "Grid": gmo.WshowGrid(WidgetGrid),
            "Model": gmo.WshowModel(WidgetModel),
            "Simulation": gmo.WshowSimtub(WidgetSimtub),
            "Layout": gmo.WshowLayout(WidgetLayout, options_list),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style({"minWidth": "400px", "width": "400px"})

    simu = mo.vstack([mo.md(""), mo.md(f"{mo.as_html(myaction())}")], gap=0)

    layout = mo.hstack([param, simu], gap=4)
    layout
    return


if __name__ == "__main__":
    app.run()
