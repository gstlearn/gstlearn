import marimo

__generated_with = "0.23.14"
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

    gmo.setEnvironment(optionBackup=False, optionDisplay=False)
    return gl, gmo, gp, mo, np, plt


@app.cell(hide_code=True)
def _(gmo):
    dbIsoPot = gmo.createDefaultIsoPot()
    dbGradPot = gmo.createDefaultGradPot()
    dbTgtePot = gmo.createDefaultTgtePot()

    WidgetMessage = gmo.WdefineMessage()
    WidgetIsoPot = gmo.WdefineIsoPot(dbIsoPot)
    WidgetGradPot = gmo.WdefineGradPot(dbGradPot)
    flagTgtePot, WidgetTgtePot = gmo.WdefineTgtePot(dbTgtePot, True)

    WidgetGrid = gmo.WdefineGrid()
    WidgetModel = gmo.WdefineModelPot()
    return (
        WidgetGradPot,
        WidgetGrid,
        WidgetIsoPot,
        WidgetMessage,
        WidgetModel,
        WidgetTgtePot,
        dbGradPot,
        dbIsoPot,
        dbTgtePot,
        flagTgtePot,
    )


@app.cell(hide_code=True)
def _(
    WidgetGradPot,
    WidgetGrid,
    WidgetIsoPot,
    WidgetMessage,
    WidgetModel,
    WidgetTgtePot,
    dbGradPot,
    dbIsoPot,
    dbTgtePot,
    flagTgtePot,
    gl,
    gmo,
    gp,
    mo,
    np,
    plt,
):
    def myaction():
        gmo.WclearMessage(WidgetMessage)
        gmo.WaddMessage(WidgetMessage, "Potential Interpolation")

        # Define the IsoPot DataBase
        dbIso = gmo.WgetIsoPot(WidgetIsoPot, dbIsoPot.clone())
        gmo.WaddMessage(
            WidgetMessage,
            "Iso-Potential ",
            "- Number of Samples =",
            dbIso.getNSample(),
            "- Number of Layers =",
            dbIso.getNOccurence(gl.ELoc.LAYER),
        )

        # Define the Gradient DataBase
        dbGrad = gmo.WgetGradPot(WidgetGradPot, dbGradPot.clone())
        gmo.WaddMessage(
            WidgetMessage, "Gradients ", "- Number of Samples =", dbIso.getNSample()
        )

        # Define the Tangent DataBase
        dbTgte = gmo.WgetTgtePot(flagTgtePot, WidgetTgtePot, dbTgtePot.clone())
        if dbTgte is not None:
            gmo.WaddMessage(
                WidgetMessage, "Tangents ", "- Number of Samples =", dbIso.getNSample()
            )

        # Define the output Grid
        grid = gmo.WgetGrid(WidgetGrid)

        # Define the Model
        model = gmo.WgetModelPot(WidgetModel)

        # Perform the Estimation
        if dbIso is not None and grid is not None and model is not None:
            err = gl.krigingPotential(
                dbiso=dbIso,
                dbgrd=dbGrad,
                dbtgt=dbTgte,
                dbout=grid,
                model=model,
                flag_pot=True,
                flag_save_data=True,
            )

        # Graphic Representation
        fig, ax = gp.init(1, 1, figsize=(5, 5))
        gmo.WaddMessage(WidgetMessage, " ")

        # Draw the Potential field (represented as raster) together with several contour lines
        gmo.plotGrid(ax, grid, name="Potential", nlevel=10)

        # Add the Layer intercept information
        gmo.plotData(ax, dbIso, name="Layer")
        levels = np.unique(np.round(dbIso.getColumn("Potential"), 4))
        gp.isoline(
            grid,
            name="Potential",
            levels=levels,
            colors="red",
            linewidths=2,
            linestyles="solid",
        )
        gmo.WaddMessage(WidgetMessage, "Potential values at Layer intercepts =", levels)

        # Add the Gradient information
        gp.gradient(dbGrad, color="black", scale=20)
        levels = np.unique(np.round(dbGrad.getColumn("Potential"), 4))
        gp.isoline(
            grid,
            name="Potential",
            levels=levels,
            colors="black",
            linewidths=2,
            linestyles="solid",
        )

        # Add the Tangent (optional) as vectors
        if dbTgte is not None:
            gp.tangent(dbTgte, color="blue", scale=20)
            levels = np.unique(np.round(dbTgte.getColumn("Potential"), 4))
            gp.isoline(
                grid,
                name="Potential",
                levels=levels,
                colors="blue",
                linewidths=2,
                linestyles="solid",
            )

        plt.tight_layout()
        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def _(
    WidgetGradPot,
    WidgetGrid,
    WidgetIsoPot,
    WidgetMessage,
    WidgetModel,
    WidgetTgtePot,
    flagTgtePot,
    gmo,
    mo,
    myaction,
):
    param = mo.ui.tabs(
        {
            "IsoPotential": gmo.WshowIsoPot(WidgetIsoPot, flagTitle=False),
            "Gradient": gmo.WshowGradPot(WidgetGradPot, flagTitle=False),
            "Tangent": gmo.WshowTgtePot(flagTgtePot, WidgetTgtePot, flagTitle=False),
            "Model": gmo.WshowModelPot(WidgetModel),
            "Grid": gmo.WshowGrid(WidgetGrid),
        }
    ).style({"minWidth": "500px", "width": "500px"})

    LeftPanel = mo.vstack(
        [
            param,
            gmo.WshowMessage(WidgetMessage),
        ],
        gap=8,
    )

    Potential = mo.vstack(
        [
            mo.as_html(myaction()),
        ],
        gap=4,
    )

    mo.hstack([LeftPanel, Potential], gap=4)
    return


if __name__ == "__main__":
    app.run()
