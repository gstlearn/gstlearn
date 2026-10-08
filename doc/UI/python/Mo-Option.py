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
def define_independent_widgets(gmo):
    WidgetDb = gmo.WdefineDb()
    WidgetGrid = gmo.WdefineGridN()

    options_list = [
        ("data", "Display of Data", True),
        ("model", "Display of Variogram and Model", True),
        ("estimation", "Display Estimation Map", True),
        ("stdev", "Display St. Dev. map", True),
    ]
    WidgetLayout = gmo.WdefineLayout(
        options_list,
        nrow=2,
        ncol=2,
        width=3,
        height=3,
    )
    return WidgetDb, WidgetGrid, WidgetLayout, options_list


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
def define_dependent_widgets(db, gmo):
    WidgetView = gmo.WdefineBox(db)
    WidgetVario = gmo.WdefineVario(nlag=10, db=db)
    return WidgetVario, WidgetView


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
    WidgetLayout,
    WidgetModel,
    WidgetVario,
    WidgetView,
    gl,
    gmo,
    gp,
    mo,
    options_list,
    plt,
):
    def myaction():
        # Define the Input Db
        db = gmo.WgetDb(WidgetDb)
        if db is None:
            return None

        targetName = db.getName({gl.toERole(gl.ELoc.Z), 0})
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

        # Récupération du layout sélectionné
        layout = gmo.WgetLayout(WidgetLayout, options_list)
        selected_options = layout.get("selected", [])

        # Define Neighborhood (Unique)
        neigh = gl.NeighUnique.create()

        # Perform the Estimation si nécessaire
        if (
            ("estimation" in selected_options or "stdev" in selected_options)
            and db is not None
            and grid is not None
            and model is not None
            and neigh is not None
        ):
            err = gl.kriging(db, grid, model, neigh)

        nx = layout.get("nx", 1)
        ny = layout.get("ny", 1)
        dimx = layout.get("dimx", 4)
        dimy = layout.get("dimy", 4)

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy], squeeze=False)
        axes = ax.ravel()

        i = 0
        for content in selected_options:
            if content == "data":
                gmo.plotData(
                    axes[i],
                    db,
                    name=targetName,
                    box=box,
                    flagBackground=flagBackground,
                )
                axes[i].set_title("Data Display")
                i += 1

            elif content == "model":
                gmo.plotVario(axes[i], vario, model, showPairs=True)
                axes[i].set_title("Variogram and Model")
                i += 1

            elif content == "estimation":
                gmo.plotGrid(
                    axes[i],
                    grid,
                    name="Kriging.*.estim",
                    flagLegend=True,
                    nlevel=20,
                )
                gmo.plotData(
                    axes[i],
                    db,
                    name=targetName,
                    box=box,
                    flagBackground=flagBackground,
                )
                axes[i].set_title("Estimation Map")
                i += 1

            elif content == "stdev":
                gmo.plotGrid(axes[i], grid, name="Kriging.*.stdev", flagLegend=True)
                gmo.plotData(
                    axes[i],
                    db,
                    name=targetName,
                    box=box,
                    flagBackground=flagBackground,
                )
                axes[i].set_title("St. Dev. Map")
                i += 1

        for axi in axes[i:]:
            axi.axis("off")

        plt.tight_layout(pad=0.2)
        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def render_ui(
    WidgetDb,
    WidgetEdit,
    WidgetGrid,
    WidgetLayout,
    WidgetModel,
    WidgetVario,
    WidgetView,
    gmo,
    mo,
    myaction,
    options_list,
):
    param = mo.ui.tabs(
        {
            "Data": gmo.WshowDb(WidgetDb),
            "Edit": gmo.WshowEdit(WidgetEdit),
            "View": gmo.WshowBox(WidgetView, gapv=0, gaph=1),
            "Grid": gmo.WshowGridN(WidgetGrid),
            "Variogram": gmo.WshowVario(WidgetVario),
            "Model": gmo.WshowModel(WidgetModel),
            "Layout": gmo.WshowLayout(WidgetLayout, options_list),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack([mo.as_html(myaction())], gap=2)

    layout = mo.hstack([param, simu], gap=4, justify="start", align="start")
    layout
    return


if __name__ == "__main__":
    app.run()
