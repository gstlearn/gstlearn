import marimo

__generated_with = "0.24.0"
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
    WidgetModel1 = gmo.WdefineModel(
        ncovmax=2, ncovdef=1, distmax=30, varmax=50, valdef="Interactive"
    )
    WidgetModel2 = gmo.WdefineModel(
        ncovmax=2, ncovdef=1, distmax=30, varmax=50, valdef="Interactive"
    )
    WidgetGrid = gmo.WdefineGrid()
    WidgetSimtub = gmo.WdefineSimtub(nbsimu=7)
    WidgetRule = gmo.WdefineRule()

    # Définition de la liste d'options sous forme de tuples (short_name, long_name, default_bool)
    options_list = [
        ("rule", "Display Lithotype Rule", True),
        ("model", "Display Model(s)", True),
        ("simulation", "Display Simulations", True),
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
        WidgetModel1,
        WidgetModel2,
        WidgetRule,
        WidgetSimtub,
        options_list,
    )


@app.cell(hide_code=True)
def define_action(
    WidgetAutoSave,
    WidgetGrid,
    WidgetLayout,
    WidgetModel1,
    WidgetModel2,
    WidgetRule,
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
        model1 = gmo.WgetModel(WidgetModel1)
        model2 = gmo.WgetModel(WidgetModel2)
        nbtuba, nbsimu, seed, flagDisplayBinary = gmo.WgetSimtub(WidgetSimtub)

        ruleprop = gmo.WgetRule(WidgetRule)

        if model1 is not None and grid is not None and ruleprop is not None:
            gl.simpgs(
                None,
                dbout=grid,
                model1=model1,
                model2=model2,
                ruleprop=ruleprop,
                nbtuba=nbtuba,
                nbsimu=nbsimu,
                seed=seed,
            )

        ngrf = ruleprop.getNGRF() if ruleprop is not None else 1
        rule = ruleprop.getRule() if ruleprop is not None else None

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
        imodel = 0
        isimu = 0

        for content in layout.get("selected", []):
            if content == "rule":
                if i < len(axes):
                    gmo.plotRule(axes[i], rule, flagLegend=True)
                    i += 1

            elif content == "model":
                for m in [model1, model2]:
                    if m is not None and i < len(axes):
                        imodel += 1
                        title = f"Model #{imodel}"
                        gmo.plotVario(axes[i], model=m, title=title)
                        i += 1

            elif content in ("simu", "simulation"):
                for sim_idx in range(nbsimu):
                    if i < len(axes):
                        isimu += 1
                        title = f"Simulation #{isimu}/{nbsimu}"
                        gmo.plotGrid(
                            axes[i],
                            grid,
                            name="Facies" if nbsimu == 1 else f"Facies.S{isimu}",
                            rule=rule,
                            title=title,
                            nlevel=0,
                            flagLegend=False,
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
    WidgetModel1,
    WidgetModel2,
    WidgetRule,
    WidgetSimtub,
    gmo,
    mo,
    myaction,
    options_list,
):
    param = mo.ui.tabs(
        {
            "Grid": gmo.WshowGrid(WidgetGrid),
            "Model1": gmo.WshowModel(WidgetModel1),
            "Model2": gmo.WshowModel(WidgetModel2),
            "Rule": gmo.WshowRule(WidgetRule),
            "Simulation": gmo.WshowSimtub(WidgetSimtub),
            "Layout": gmo.WshowLayout(WidgetLayout, options_list, gapv=1),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack([mo.md(""), mo.md(f"{mo.as_html(myaction())}")], gap=0)

    layout = mo.hstack([param, simu], gap=4)
    layout
    return


if __name__ == "__main__":
    app.run()
