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
def _(
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

        layout = gmo.WgetLayout(WidgetLayout, 1, nbsimu, ngrf)
        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        contents = layout["contents"]
        valid = {"rule", "model", "simu"}
        contents_local = (c for c in contents if c in valid)

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy])
        axes = ax.ravel()

        isimu = 0
        imodel = 0
        for axi in axes:
            content = next(contents_local, None)

            if content == "rule":
                gmo.plotRule(axi, rule, flagLegend=True)

            elif content == "model":
                imodel += 1
                model = model1 if imodel == 1 else model2
                title = f"Model #{imodel}"
                gmo.plotVario(axi, model=model, title=title)

            elif content == "simu":
                isimu += 1

                title = f"Simulation #{isimu}/{nbsimu}"
                gmo.plotGrid(
                    axi,
                    grid,
                    name="Facies" if nbsimu == 1 else f"Facies.S{isimu}",
                    rule=rule,
                    title=f"Simulation #{isimu}/{nbsimu}",
                    nlevel=0,
                    flagLegend=False,
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
    ).style({"minWidth": "470px", "width": "400px"})

    simu = mo.as_html(myaction())

    mo.hstack([param, simu], gap=4)
    return


if __name__ == "__main__":
    app.run()
