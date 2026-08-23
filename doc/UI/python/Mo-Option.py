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
def define_db_widget(gmo):
    WidgetDb = gmo.WdefineDb()
    return (WidgetDb,)


@app.cell(hide_code=True)
def get_initial_db(WidgetDb, gmo):
    dbinit = gmo.WgetDb(WidgetDb)
    return (dbinit,)


@app.cell(hide_code=True)
def define_edit_widget(dbinit, gmo):
    WidgetEdit = gmo.WdefineEdit(dbinit)
    return (WidgetEdit,)


@app.cell(hide_code=True)
def get_edited_db(WidgetEdit, dbinit, gmo):
    db = gmo.WgetEdit(WidgetEdit, dbinit)
    return (db,)


@app.cell(hide_code=True)
def define_view_widget(db, gmo):
    WidgetView = gmo.WdefineBox(db)
    return (WidgetView,)


@app.cell(hide_code=True)
def define_grid_widget(gmo):
    WidgetGrid = gmo.WdefineGridN()
    return (WidgetGrid,)


@app.cell(hide_code=True)
def define_vario_widget(db, gmo):
    WidgetVario = gmo.WdefineVario(nlag=10, db=db)
    return (WidgetVario,)


@app.cell(hide_code=True)
def get_variogram(WidgetVario, db, gmo):
    vario = gmo.WgetVario(WidgetVario, db=db)
    return (vario,)


@app.cell(hide_code=True)
def define_model_widget(gmo, vario):
    WidgetModel = gmo.WdefineModel(ncovmax=3, vario=vario)
    return (WidgetModel,)


@app.cell(hide_code=True)
def define_action(
    WidgetDb,
    WidgetGrid,
    WidgetModel,
    WidgetVario,
    WidgetView,
    gl,
    gmo,
    gp,
    mo,
    plt,
):
    def myaction():
        # Define the Input Db
        db = gmo.WgetDb(WidgetDb)
        if db is None:
            return None

        targetName = db.getNameByLocator(gl.ELoc.Z)
        ndim = db.getNLoc(gl.ELoc.X)
        if ndim != 2:
            print("The 'db' should be 2-D (", ndim, ")")
            return None
        if db.getColIdx(targetName) < 0:
            print("The 'db' should contain the target variable")
            return None

        # Define the output Grid
        box, flagBackground = gmo.WgetBox(WidgetView)
        grid = gmo.WgetGridN(WidgetGrid, box)

        # Define the Variogram parameters
        vario = gmo.WgetVario(WidgetVario, db)

        # Define the Model
        model = gmo.WgetModel(WidgetModel, vario)

        # Define Neighborhood (Unique)
        neigh = gl.NeighUnique.create()

        # Perform the Estimation
        if (
            db is not None
            and grid is not None
            and model is not None
            and neigh is not None
        ):
            err = gl.kriging(db, grid, model, neigh)

        fig, ax = gp.init(2, 2, figsize=(6, 6))
        gmo.plotData(
            ax[0, 0], db, name=targetName, box=box, flagBackground=flagBackground
        )
        gmo.plotVario(ax[0, 1], vario, model, showPairs=True)

        gmo.plotGrid(ax[1, 0], grid, name="Kriging.*.estim", flagLegend=True, nlevel=20)
        gmo.plotData(
            ax[1, 0], db, name=targetName, box=box, flagBackground=flagBackground
        )
        gmo.plotGrid(ax[1, 1], grid, name="Kriging.*.stdev", flagLegend=True)
        gmo.plotData(
            ax[1, 1], db, name=targetName, box=box, flagBackground=flagBackground
        )

        plt.tight_layout(pad=0.2)
        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def render_ui(
    WidgetDb,
    WidgetEdit,
    WidgetGrid,
    WidgetModel,
    WidgetVario,
    WidgetView,
    gmo,
    mo,
    myaction,
):
    param = mo.ui.tabs(
        {
            "Data": gmo.WshowDb(WidgetDb),
            "Edit": gmo.WshowEdit(WidgetEdit),
            "View": gmo.WshowBox(WidgetView, gapv=0, gaph=1),
            "Grid": gmo.WshowGridN(WidgetGrid),
            "Variogram": gmo.WshowVario(WidgetVario),
            "Model": gmo.WshowModel(WidgetModel),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack([mo.md(""), mo.md(f"{mo.as_html(myaction())}")], gap=0)

    layout = mo.hstack([param, simu], gap=4)
    layout
    return


if __name__ == "__main__":
    app.run()
