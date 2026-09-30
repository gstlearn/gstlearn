import marimo

__generated_with = "0.25.0"
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

    # Définition des options sous forme de tuples (short_name, long_name, default_bool)
    options_list = [
        ("model", "Display Model(s)", True),
        ("simulation", "Display Simulations", True),
        ("average", "Display Simulation Average", False),
        ("dispersion", "Display Simulation Dispersion", False),
    ]

    WidgetLayout = gmo.WdefineLayout(
        options_list,
        nrow=3,
        ncol=3,
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

        if model is not None and grid is not None:
            gl.simtub(
                None, dbout=grid, model=model, nbtuba=nbtuba, nbsimu=nbsimu, seed=seed
            )
            grid.statisticsBySample(
                names=["Simu.*"], opers=[gl.EStatOption.MEAN, gl.EStatOption.STDV]
            )

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
            if content == "model":
                if model is not None and i < len(axes):
                    gmo.plotVario(axes[i], model=model)
                    i += 1

            elif content == "simulation":
                for sim_idx in range(nbsimu):
                    if i < len(axes):
                        isimu += 1
                        gmo.plotGrid(
                            axes[i],
                            grid,
                            name="Simu" if nbsimu == 1 else f"Simu.S{isimu}",
                            title=f"Simulation #{isimu}/{nbsimu}",
                            flagLegend=False,
                            flagBinary=flagDisplayBinary,
                        )
                        i += 1

            elif content == "average":
                if i < len(axes):
                    gmo.plotGrid(
                        axes[i],
                        grid,
                        name="Stats.MEAN",
                        title="Simulation Average",
                        flagLegend=False,
                        flagBinary=flagDisplayBinary,
                    )
                    i += 1

            elif content == "dispersion":
                if i < len(axes):
                    gmo.plotGrid(
                        axes[i],
                        grid,
                        name="Stats.STDV",
                        title="Simulation Dispersion",
                        flagLegend=False,
                        flagBinary=flagDisplayBinary,
                    )
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
