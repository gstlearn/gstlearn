import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def my_imports():
    import marimo as mo

    import gstlearn as gl
    import gstlearn.document as gdoc
    import gstlearn.gstmarimo as gmo
    import gstlearn.plot as gp

    import matplotlib.pyplot as plt
    import numpy as np

    gmo.setEnvironment(optionBackup=True, optionDisplay=False)
    return gl, gmo, gp, mo, plt


@app.cell(hide_code=True)
def define_base_widgets(gmo):
    WidgetAutoSave = gmo.WdefineAutoSave()
    WidgetDb = gmo.WdefineDb()

    # Définition des options propres à l'application Db
    options_list = [
        ("data", "Display Data", True),
    ]

    WidgetLayout = gmo.WdefineLayout(
        options_list,
        nrow=1,
        ncol=1,
        width=4,
        height=4,
    )
    return WidgetAutoSave, WidgetDb, WidgetLayout, options_list


@app.cell(hide_code=True)
def define_edit_widget(WidgetDb, gmo):
    dbinit = gmo.WgetDb(WidgetDb)
    WidgetEdit = gmo.WdefineEdit(dbinit)
    return WidgetEdit, dbinit


@app.cell(hide_code=True)
def update_db(WidgetAutoSave, WidgetEdit, dbinit, gmo):
    autosave = gmo.WgetAutoSave(WidgetAutoSave)
    db = gmo.WgetEdit(WidgetEdit, dbinit)
    return (db,)


@app.cell(hide_code=True)
def render_ui(
    WidgetAutoSave,
    WidgetDb,
    WidgetEdit,
    WidgetLayout,
    db,
    gl,
    gmo,
    gp,
    mo,
    options_list,
    plt,
):
    def myplot():
        if db is None:
            return None

        # Récupération du layout sélectionné par l'utilisateur
        layout = gmo.WgetLayout(WidgetLayout, options_list)
        selected = layout.get("selected", [])

        if not selected:
            return None

        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy], squeeze=False)
        axes = ax.ravel()
        axi = axes[0]

        n_targets = db.getNLoc(gl.ELoc.Z) if db is not None else 1

        if "data" in selected:
            target_name = db.getNameByLocator(gl.ELoc.Z, 0) if db else "z"
            gmo.plotData(axi, db, name=target_name, title=f"Data: {target_name}")
        else:
            axi.axis("off")

        plt.tight_layout(pad=0.2)
        return fig

    param = mo.ui.tabs(
        {
            "Data": gmo.WshowDb(WidgetDb),
            "Edit": gmo.WshowEdit(WidgetEdit),
            "Layout": gmo.WshowLayout(WidgetLayout, options_list),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    fig = myplot()
    simu = mo.mpl.interactive(fig) if fig is not None else mo.md("")
    simu_centered = mo.vstack([simu]).center()
    layout = mo.hstack([param, simu_centered], gap=4)
    layout
    return


if __name__ == "__main__":
    app.run()
