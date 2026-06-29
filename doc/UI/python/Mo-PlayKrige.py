import marimo

__generated_with = "0.23.11"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import gstlearn.document as gdoc

    import numpy as np
    import matplotlib.pyplot as plt

    gmo.setEnvironment(optionBackup=True, optionDisplay=False)
    return gl, gmo, gp, mo


@app.cell(hide_code=True)
def _(gmo):
    WidgetAutoSave = gmo.WdefineAutoSave()
    WidgetDb = gmo.WdefineDb()
    WidgetGrid = gmo.WdefineGrid(nxdef=20, dxdef=5)
    WidgetVario = gmo.WdefineVario(nlag=10)
    WidgetNeigh = gmo.WdefineNeigh(radius=50)
    WidgetWeights = gmo.WdefineWeights()
    WidgetLayout = gmo.WdefineLayout(nrow=2, ncol=3, width=3, height=3)
    # WidgetLayout = gmo.WdefineLayout(nrow=2, ncol=2, width=3, height=3, defaults={
    #     "estimation": False,
    #     "stdev": False, })
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
def _(WidgetDb, gmo):
    db0 = gmo.WgetDb(WidgetDb)
    return (db0,)


@app.cell(hide_code=True)
def _(WidgetVario, db0, gmo):
    vario = gmo.WgetVario(WidgetVario, db=db0)
    return (vario,)


@app.cell(hide_code=True)
def _(gmo, vario):
    WidgetModel = gmo.WdefineModel(ncovmax=3, vario=vario)
    return (WidgetModel,)


@app.cell(hide_code=True)
def _(
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
):
    def myaction(Kindex_pick=-1, Xindex_pick=-1):
        # Define the autosave option
        autosave = gmo.WgetAutoSave(WidgetAutoSave)

        # Read the elements
        data = gmo.WgetDb(WidgetDb)
        grid = gmo.WgetGrid(WidgetGrid)
        vario = gmo.WgetVario(WidgetVario, data)
        model = gmo.WgetModel(WidgetModel, vario)
        neigh = gmo.WgetNeigh(WidgetNeigh)

        # Read the target nodes (from interface or from Picking)
        Kindex, Xindex = gmo.WgetWeights(WidgetWeights)
        if Kindex_pick >= 0:
            Kindex = Kindex_pick

        if Xindex_pick >= 0:
            Xindex = Xindex_pick

        print("Kindex final = ", Kindex)
        print("Xindex final = ", Xindex)

        # Perform the Estimation
        if (
            data is not None
            and grid is not None
            and model is not None
            and neigh is not None
        ):
            err = gl.kriging(data, grid, model, neigh)
            if Kindex >= 0:
                err = gl.krigWeights(data, grid, model, neigh, Kindex)
            if Xindex >= 0:
                err = gl.xvalidWeights(data, model, neigh, Xindex)

        nvar = data.getNLoc(gl.ELoc.Z)
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
        ax_KWeight = None
        ax_XWeight = None
        for name, count in layout["render_plan"]:
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
                    ax_Picking = axi
                    name = "KWeights." + targetName
                    gmo.plotWeights(
                        axi,
                        grid,
                        data,
                        name,
                        neigh,
                        Kindex,
                        title="Kriging Weights",
                    )

                elif name == "XWeights":
                    name = "XWeights." + targetName
                    gmo.plotWeights(
                        axi,
                        data,
                        data,
                        name,
                        neigh,
                        Xindex,
                        title="Cross-Validation Weights",
                    )

                else:
                    axi.axis("off")

                i += 1

        # hide remaining axes
        for axi in axes[i:]:
            axi.axis("off")

        fig.tight_layout(pad=0)

        return fig, axes, data, grid, ax_Picking

    return (myaction,)


@app.cell
def _(mo):
    Kindex_state, set_Kindex = mo.state(-1)
    Xindex_state, set_Xindex = mo.state(-1)
    return Kindex_state, Xindex_state, set_Kindex, set_Xindex


@app.cell(hide_code=True)
def _(Kindex_state, Xindex_state, myaction):
    fig, axs, data, grid, ax_Picking = myaction(
        Kindex_pick=Kindex_state(),
        Xindex_pick=Xindex_state(),
    )
    return ax_Picking, data, grid


@app.cell(hide_code=True)
def _(ax_Picking, mo):
    if ax_Picking is not None:
        widget_Picking = mo.ui.matplotlib(ax_Picking, debounce=True)
    return (widget_Picking,)


@app.cell(hide_code=True)
def _(
    WidgetDb,
    WidgetGrid,
    WidgetLayout,
    WidgetModel,
    WidgetNeigh,
    WidgetVario,
    WidgetWeights,
    gmo,
    mo,
    widget_Picking,
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
        }
    ).style({"minWidth": "450px", "width": "450px"})

    layout = mo.hstack(
        [
            param,
            widget_Picking,
        ],
        gap=4,
    )

    layout
    return


@app.cell(hide_code=True)
def _(WidgetWeights, data, gmo, grid, set_Kindex, set_Xindex, widget_Picking):
    index_data = gmo.getSelectedIndex(widget_Picking, data)
    index_grid = gmo.getSelectedIndex(widget_Picking, grid)

    Kindex, Xindex = gmo.WgetWeights(WidgetWeights)

    if index_grid >= 0:
        Kindex = index_grid
        set_Kindex(Kindex)

    if index_data >= 0:
        Xindex = index_data
        set_Xindex(Xindex)

    print("Effective Kindex =", Kindex)
    print("Effective Xindex =", Xindex)
    return


if __name__ == "__main__":
    app.run()
