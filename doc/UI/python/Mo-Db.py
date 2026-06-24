import marimo

__generated_with = "0.19.2"
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
    return (gmo,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetAutoSave = gmo.WdefineAutoSave()
    WidgetDb = gmo.WdefineDb()
    return WidgetAutoSave, WidgetDb


@app.cell(hide_code=True)
def _(WidgetDb, gmo):
    dbinit = gmo.WgetDb(WidgetDb)
    WidgetEdit = gmo.WdefineEdit(dbinit)
    return WidgetEdit, dbinit


@app.cell(hide_code=True)
def _(WidgetAutoSave, WidgetEdit, dbinit, gmo):
    autosave = gmo.WgetAutoSave(WidgetAutoSave)
    db = gmo.WgetEdit(WidgetEdit, dbinit)
    return


app._unparsable_cell(
    r"""
    def myplot():
        fig, ax = gp.init(figsize=(4, 4))
        gmo.plotData(ax, db, name="z")
        return fig

    param = mo.ui.tabs(
        {"Data": gmo.WshowDb(WidgetDb),
         "Edit": gmo.WshowEdit(WidgetEdit)},
         "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
    ).style({"minWidth": "700px", "width": "350px"})

    simu = mo.vstack(
        [mo.md(""), mo.md(f"Plotting the Data Base:{mo.as_html(myplot())}")], gap=4
    )

    mo.hstack([param, simu], gap=4)
    """,
    column=None,
    disabled=False,
    hide_code=True,
    name="_",
)


if __name__ == "__main__":
    app.run()
