import marimo

__generated_with = "0.23.14"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def my_imports():
    import marimo as mo

    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import gstlearn.document as gdoc

    import numpy as np
    import matplotlib.pyplot as plt

    gmo.setEnvironment(optionBackup=True, optionDisplay=False)
    return gmo, gp, mo


@app.cell(hide_code=True)
def define_base_widgets(gmo):
    WidgetAutoSave = gmo.WdefineAutoSave()
    WidgetDb = gmo.WdefineDb()
    return WidgetAutoSave, WidgetDb


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
def render_ui(WidgetAutoSave, WidgetDb, WidgetEdit, db, gmo, gp, mo):
    def myplot():
        fig, ax = gp.init(figsize=(4, 4))
        gmo.plotData(ax, db, name="z")
        return fig

    # L'accolade du dictionnaire englobe désormais bien les 3 onglets :
    param = mo.ui.tabs(
        {
            "Data": gmo.WshowDb(WidgetDb),
            "Edit": gmo.WshowEdit(WidgetEdit),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style({"minWidth": "500px", "width": "350px"})

    fig = myplot()
    simu = mo.as_html(fig) if fig is not None else mo.md("")
    simu_centered = mo.vstack([simu]).center()
    layout = mo.hstack([param, simu_centered], gap=4)
    layout
    return


if __name__ == "__main__":
    app.run()
