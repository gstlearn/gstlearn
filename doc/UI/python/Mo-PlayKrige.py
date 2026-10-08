import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def cell_imports():
    import marimo as mo
    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import matplotlib.pyplot as plt

    import numpy as np
    import pandas as pd

    gmo.setEnvironment(
        optionBackup=True,
        optionDisplay=False,
    )
    return gl, gmo, gp, mo, plt


@app.cell(hide_code=True)
def cell_define_widgets(gmo):
    WidgetAutoSave = gmo.WdefineAutoSave()
    WidgetDb = gmo.WdefineDb()
    WidgetGrid = gmo.WdefineGrid(
        nxdef=20,
        dxdef=5,
    )
    WidgetVario = gmo.WdefineVario(
        nlag=10,
    )
    WidgetNeigh = gmo.WdefineNeigh(
        radius=50,
    )

    # Définition des options propres à l'application
    options_list = [
        ("data", "Display Data", True),
        ("model", "Display Model", True),
        ("estimation", "Display Estimation Map(s)", True),
        ("stdev", "Display Standard Deviation Map(s)", True),
        ("KWeights", "Display Kriging Weights", True),
        ("XWeights", "Display Cross-Validation", True),
    ]

    WidgetLayout = gmo.WdefineLayout(
        options_list,
        nrow=2,
        ncol=3,
        width=3,
        height=3,
    )
    return (
        WidgetAutoSave,
        WidgetDb,
        WidgetGrid,
        WidgetLayout,
        WidgetNeigh,
        WidgetVario,
        options_list,
    )


@app.cell(hide_code=True)
def cell_define_widget_model(WidgetDb, WidgetVario, gmo):
    model_vario = gmo.WgetVario(
        WidgetVario,
        db=gmo.WgetDb(WidgetDb),
    )

    WidgetModel = gmo.WdefineModel(
        ncovmax=3,
        vario=model_vario,
    )
    return (WidgetModel,)


@app.cell(hide_code=True)
def cell_define_picking_state(mo):
    # État atomique (Kindex, Xindex)
    picking_indices_state, set_picking_indices = mo.state((-1, -1))
    return picking_indices_state, set_picking_indices


@app.cell(hide_code=True)
def cell_get_render_layout(WidgetLayout, gmo, options_list):
    render_layout = gmo.WgetLayout(
        WidgetLayout,
        options_list,
    )
    return (render_layout,)


@app.cell(hide_code=True)
def cell_compute_results(
    WidgetAutoSave,
    WidgetDb,
    WidgetGrid,
    WidgetModel,
    WidgetNeigh,
    WidgetVario,
    gl,
    gmo,
    picking_indices_state,
):
    result_autosave = gmo.WgetAutoSave(WidgetAutoSave)
    result_data = gmo.WgetDb(WidgetDb)
    result_grid = gmo.WgetGrid(WidgetGrid)
    result_vario = gmo.WgetVario(WidgetVario, result_data)
    result_model = gmo.WgetModel(WidgetModel, result_vario)
    result_neigh = gmo.WgetNeigh(WidgetNeigh)

    result_Kindex, result_Xindex = picking_indices_state()

    if (
        result_data is not None
        and result_grid is not None
        and result_model is not None
        and result_neigh is not None
    ):
        gl.kriging(
            result_data,
            result_grid,
            result_model,
            result_neigh,
        )

        if result_Kindex >= 0 and result_Kindex < result_grid.getNSample():
            gl.krigWeights(
                result_data,
                result_grid,
                result_model,
                result_neigh,
                result_Kindex,
            )

        if result_Xindex >= 0 and result_Xindex < result_data.getNSample():
            gl.xvalidWeights(
                result_data,
                result_model,
                result_neigh,
                result_Xindex,
            )
    return (
        result_Kindex,
        result_Xindex,
        result_data,
        result_grid,
        result_model,
        result_neigh,
        result_vario,
    )


@app.cell(hide_code=True)
def cell_make_figure(
    gl,
    gmo,
    gp,
    plt,
    render_layout,
    result_Kindex,
    result_Xindex,
    result_data,
    result_grid,
    result_model,
    result_neigh,
    result_vario,
):
    result_nx = render_layout["nx"]
    result_ny = render_layout["ny"]
    result_dimx = render_layout["dimx"]
    result_dimy = render_layout["dimy"]

    result_fig, result_ax = gp.init(
        result_nx,
        result_ny,
        figsize=[result_ny * result_dimy, result_nx * result_dimx],
        squeeze=False,
    )

    result_axes = result_ax.ravel()
    result_Kaxis = None
    result_i = 0

    n_targets = result_data.getNLoc(gl.ELoc.Z) if result_data is not None else 1

    for result_name in render_layout.get("selected", []):
        for result_k in range(n_targets):
            if result_i >= len(result_axes):
                break

            result_axi = result_axes[result_i]
            result_targetName = (
                result_data.getName(gl.RoleID(gl.toERole(gl.ELoc.Z), result_k))
                if result_data
                else ""
            )

            if result_name == "data":
                gmo.plotData(
                    result_axi,
                    result_data,
                    name=result_targetName,
                    title="Data:" + result_targetName,
                    c="blue",
                )

            elif result_name == "model":
                gmo.plotVario(result_axi, result_vario, result_model)

            elif result_name == "estimation":
                result_gridName = "Kriging." + result_targetName + ".estim"
                gmo.plotGrid(
                    result_axi,
                    result_grid,
                    name=result_gridName,
                    title="Estimation",
                )
                gmo.plotData(
                    result_axi,
                    result_data,
                    name=result_targetName,
                    flagTitle=False,
                    c="blue",
                )

            elif result_name == "stdev":
                result_gridName = "Kriging." + result_targetName + ".stdev"
                gmo.plotGrid(
                    result_axi,
                    result_grid,
                    name=result_gridName,
                    title="St. dev. of Estimation Error",
                )
                gmo.plotData(
                    result_axi,
                    result_data,
                    name=result_targetName,
                    flagTitle=False,
                    c="blue",
                )

            elif result_name == "KWeights":
                result_Kaxis = result_axi
                result_name_kw = "KWeights." + result_targetName
                result_title_kw = (
                    f"Kriging (Index: {result_Kindex})"
                    if result_Kindex >= 0
                    else "Kriging (Pick target)"
                )
                gmo.plotWeights(
                    result_axi,
                    result_grid,
                    result_data,
                    result_name_kw,
                    result_neigh,
                    result_Kindex,
                    title=result_title_kw,
                )

            elif result_name == "XWeights":
                result_name_xw = "XWeights." + result_targetName
                result_title_xw = (
                    f"XValidation (Index: {result_Xindex})"
                    if result_Xindex >= 0
                    else "XValidation (Pick target)"
                )
                gmo.plotWeights(
                    result_axi,
                    result_data,
                    result_data,
                    result_name_xw,
                    result_neigh,
                    result_Xindex,
                    title=result_title_xw,
                )

            else:
                result_axi.axis("off")

            result_i += 1

    for result_axi in result_axes[result_i:]:
        result_axi.axis("off")

    plt.tight_layout(pad=0.2)
    return result_Kaxis, result_fig


@app.cell(hide_code=True)
def cell_make_KWidget(mo, result_Kaxis):
    KWidget = None

    if result_Kaxis is not None:
        KWidget = mo.ui.matplotlib(
            result_Kaxis,
            debounce=True,
        )
    return (KWidget,)


@app.cell(hide_code=True)
def cell_update_indices_from_KWidget(
    KWidget,
    gmo,
    picking_indices_state,
    result_data,
    result_grid,
    set_picking_indices,
):
    if KWidget is not None and bool(KWidget.value):
        selected_grid_index = gmo.getSelectedIndex(KWidget, result_grid)
        selected_data_index = gmo.getSelectedIndex(KWidget, result_data)

        grid_idx = selected_grid_index if selected_grid_index is not None else -1
        data_idx = selected_data_index if selected_data_index is not None else -1

        new_indices = (grid_idx, data_idx)
        current_indices = picking_indices_state()

        # Met à jour l'état uniquement s'il y a un réel changement d'indice sélectionné
        if new_indices != current_indices and new_indices != (-1, -1):
            set_picking_indices(new_indices)
    return


@app.cell(hide_code=True)
def cell_render_ui(
    KWidget,
    WidgetAutoSave,
    WidgetDb,
    WidgetGrid,
    WidgetLayout,
    WidgetModel,
    WidgetNeigh,
    WidgetVario,
    gmo,
    mo,
    options_list,
    result_fig,
):
    ui_param = mo.ui.tabs(
        {
            "Data": gmo.WshowDb(WidgetDb),
            "Grid": gmo.WshowGrid(WidgetGrid),
            "Variogram": gmo.WshowVario(WidgetVario),
            "Model": gmo.WshowModel(WidgetModel),
            "Neigh": gmo.WshowNeigh(WidgetNeigh),
            "Layout": gmo.WshowLayout(WidgetLayout, options_list),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style(
        {
            "minWidth": "380px",
            "width": "380px",
        }
    )

    if KWidget is not None:
        ui_graph = mo.vstack(
            [
                mo.md(""),
                KWidget,
            ],
            gap=0,
        )
    else:
        ui_graph = mo.vstack(
            [
                mo.md(""),
                mo.mpl.interactive(result_fig),
            ],
            gap=0,
        )

    mo.hstack(
        [
            ui_param,
            ui_graph,
        ],
        gap=4,
    )
    return


if __name__ == "__main__":
    app.run()
