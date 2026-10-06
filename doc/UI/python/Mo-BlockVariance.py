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
    return gl, gmo, gp, mo, np, plt


@app.cell(hide_code=True)
def define_independent_widgets(gmo):
    # Model
    WidgetModel = gmo.WdefineModel(ncovmax=2, withFit=False)

    # Grille fine
    WidgetFineGrid = gmo.WdefineGrid()

    # Layout
    options_list = [
        ("fine", "Display the Fine Grid", True),
        ("coarse", "Display the Coarse Grid", True),
    ]
    WidgetLayout = gmo.WdefineLayout(
        options_list,
        nrow=1,
        ncol=2,
        width=3,
        height=3,
    )

    # Messages pour les variances
    WMessages = gmo.WdefineMessage()
    return WMessages, WidgetFineGrid, WidgetLayout, WidgetModel, options_list


@app.cell(hide_code=True)
def explain_state(mo):
    get_show_explain, set_show_explain = mo.state(False)
    return get_show_explain, set_show_explain


@app.cell(hide_code=True)
def explain_buttons(mo, set_show_explain):
    open_explain_btn = mo.ui.button(
        label="Explanations...",
        on_click=lambda _: set_show_explain(True),
    )

    close_explain_btn = mo.ui.button(
        label="❌ Close",
        on_click=lambda _: set_show_explain(False),
    )
    return close_explain_btn, open_explain_btn


@app.cell(hide_code=True)
def explain_modal_builder(
    close_explain_btn,
    get_show_explain,
    gmo,
    mo,
    open_explain_btn,
):
    # Appel simplifié du composant hébergé dans gstmarimo
    modal_overlay = gmo.WexplainMarkdown(
        get_show_modal=get_show_explain,
        close_btn=close_explain_btn,
        filepaths=["Cvv.md", "Variance_Reduction.md"],
        modal_title="Block Variance Theory",
    )

    WidgetExplain = mo.vstack([open_explain_btn, modal_overlay])
    return (WidgetExplain,)


@app.cell(hide_code=True)
def get_fine_grid(WidgetFineGrid, gmo):
    grid_fine = gmo.WgetGrid(WidgetFineGrid)
    return (grid_fine,)


@app.cell(hide_code=True)
def define_coarsening_widget(gmo, grid_fine):
    WidgetCoarseGrid = gmo.WdefineCoarseGrid(grid_fine)
    return (WidgetCoarseGrid,)


@app.cell(hide_code=True)
def define_action(
    WMessages,
    WidgetCoarseGrid,
    WidgetLayout,
    WidgetModel,
    gl,
    gmo,
    gp,
    grid_fine,
    mo,
    np,
    options_list,
    plt,
):
    def myaction():
        model = gmo.WgetModel(WidgetModel)
        if model is None:
            gmo.WclearMessage(WMessages)
            return None

        if grid_fine is None:
            gmo.WclearMessage(WMessages)
            return None

        grid_coarse = gmo.WgetCoarseGrid(WidgetCoarseGrid, grid_fine)
        if grid_coarse is None:
            gmo.WclearMessage(WMessages)
            return None

        grid_local = grid_fine.clone()
        err = gl.simtub(None, grid_local, model, nbtuba=1000)

        # Calcul des statistiques (moyenne par bloc)
        err = gl.dbStatisticsOnGrid(grid_local, grid_coarse, gl.EStatOption.MEAN)

        sim_fine = np.array(grid_local["Simu*"])
        sim_coarse = np.array(grid_coarse["Stats*"])

        # Variances empiriques
        var_fine = float(np.nanvar(sim_fine)) if sim_fine.size > 0 else 0.0
        var_coarse = float(np.nanvar(sim_coarse)) if sim_coarse.size > 0 else 0.0

        # Reduction de variance empirique: Var(Z_fine) - Var(Z_coarse)
        deltaVarEmp = var_fine - var_coarse

        # Reduction de variance theorique: C(0) - Cbar(v, v)
        ndisc = int(WidgetCoarseGrid["ndisc"].value)
        c0 = model.getTotalSill()
        dx_coarse = grid_coarse.getDXs()
        cbar = model.evalCvv(dx_coarse, ndisc=[ndisc, ndisc])
        deltaVarTheo = c0 - cbar

        # Mise a jour du message reactif
        gmo.WsetMessage(
            WMessages,
            f"Theoretical Variance Reduction: **{deltaVarTheo:.3f}**<br>"
            f"Empirical Variance Reduction: **{deltaVarEmp:.3f}**",
        )

        layout = gmo.WgetLayout(WidgetLayout, options_list)

        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy], squeeze=False)
        axes = ax.ravel()

        i = 0
        for content in layout.get("selected", []):
            if content == "fine":
                gmo.plotGrid(
                    axes[i],
                    grid_local,
                    "*Simu",
                    nlevel=0,
                )
                axes[i].axes.axis("off")
                axes[i].set_title("Z(x)")
                axes[i].cell(grid_coarse, color="black", linewidth=0.5)
                i += 1

            elif content == "coarse":
                gmo.plotGrid(
                    axes[i],
                    grid_coarse,
                    "Stats*",
                    nlevel=0,
                )
                axes[i].axes.axis("off")
                axes[i].set_title("Z(v)")
                axes[i].cell(grid_coarse, color="black", linewidth=0.5, shift=0)
                i += 1

        for axi in axes[i:]:
            axi.axis("off")

        plt.tight_layout(pad=0.2)
        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def render_ui(
    WMessages,
    WidgetCoarseGrid,
    WidgetExplain,
    WidgetFineGrid,
    WidgetLayout,
    WidgetModel,
    gmo,
    mo,
    myaction,
    options_list,
):
    tabs = mo.ui.tabs(
        {
            "Fine Grid": gmo.WshowGrid(WidgetFineGrid),
            "Coarse Grid": gmo.WshowCoarseGrid(WidgetCoarseGrid),
            "Model": gmo.WshowModel(WidgetModel),
            "Layout": gmo.WshowLayout(WidgetLayout, options_list),
        },
    )

    fig = myaction()

    msg_rendered = gmo.WshowMessage(
        WMessages, border=False, minHeight="auto", width="100%"
    )

    param = mo.vstack([tabs, WidgetExplain, msg_rendered], gap=2).style(
        {"minWidth": "350px", "width": "350px"}
    )
    simu = mo.vstack([mo.as_html(fig) if fig is not None else mo.md("")], gap=2)

    layout = mo.hstack([param, simu], gap=4, justify="start", align="start")
    layout
    return


if __name__ == "__main__":
    app.run()
