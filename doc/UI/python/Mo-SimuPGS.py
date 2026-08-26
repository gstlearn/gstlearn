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
    WidgetModel1 = gmo.WdefineModel(
        ncovmax=2, ncovdef=1, distmax=30, varmax=50, valdef="Interactive"
    )
    WidgetModel2 = gmo.WdefineModel(
        ncovmax=2, ncovdef=1, distmax=30, varmax=50, valdef="Interactive"
    )
    WidgetGrid = gmo.WdefineGrid()
    WidgetSimtub = gmo.WdefineSimtub(nbsimu=7)
    WidgetRule = gmo.WdefineRule()
    WidgetLayout = gmo.WdefineLayout(3, 3, 3, 3)
    WidgetAutoSave = gmo.WdefineAutoSave()
    return (
        WidgetAutoSave,
        WidgetGrid,
        WidgetLayout,
        WidgetModel1,
        WidgetModel2,
        WidgetRule,
        WidgetSimtub,
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
            nvar=1,
            nbsimu=nbsimu,
            ngrf=ngrf,
            valid=["model", "rule", "simu", "simulation", "average"],
        )
        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy], squeeze=False)
        axes = ax.ravel()

        render_plan = layout.get("render_plan", [])
        plan_items = []
        for name, count in render_plan:
            plan_items.extend([name] * count)

        i = 0
        imodel = 0
        isimu = 0

        for axi in axes:
            content = plan_items[i] if i < len(plan_items) else None

            if content == "rule":
                gmo.plotRule(axi, rule, flagLegend=True)

            elif content == "model":
                imodel += 1
                model = model1 if imodel == 1 else model2
                title = f"Model #{imodel}"
                gmo.plotVario(axi, model=model, title=title)

            elif content in ("simu", "simulation"):
                isimu += 1
                title = f"Simulation #{isimu}/{nbsimu}"
                gmo.plotGrid(
                    axi,
                    grid,
                    name="Facies" if nbsimu == 1 else f"Facies.S{isimu}",
                    rule=rule,
                    title=title,
                    nlevel=0,
                    flagLegend=False,
                )

            else:
                axi.axis("off")

            i += 1

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
):
    param = mo.ui.tabs(
        {
            "Grid": gmo.WshowGrid(WidgetGrid),
            "Model1": gmo.WshowModel(WidgetModel1),
            "Model2": gmo.WshowModel(WidgetModel2),
            "Rule": gmo.WshowRule(WidgetRule),
            "Simulation": gmo.WshowSimtub(WidgetSimtub),
            "Layout": gmo.WshowLayout(WidgetLayout, gapv=1),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack([mo.md(""), mo.md(f"{mo.as_html(myaction())}")], gap=0)

    layout = mo.hstack([param, simu], gap=4)
    layout
    return


if __name__ == "__main__":
    app.run()
