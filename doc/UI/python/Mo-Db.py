import marimo

__generated_with = "0.19.2"
app = marimo.App(css_file="custom.css")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import gstlearn.document as gdoc

    import numpy as np
    import matplotlib.pyplot as plt

    return gmo, gp, mo


@app.cell(hide_code=True)
def _(gmo):
    gmo.setEnvironment(optionSaveNF=True, optionPrint=False)
    return


@app.cell(hide_code=True)
def _(gmo):
    WidgetDb = gmo.WdefineDb()
    return (WidgetDb,)


@app.cell(hide_code=True)
def _(WidgetDb, gmo):
    dbinit = gmo.WgetDb(WidgetDb)
    return (dbinit,)


@app.cell(hide_code=True)
def _(dbinit, gmo):
    WidgetEdit = gmo.WdefineEdit(dbinit)
    return (WidgetEdit,)


@app.cell(hide_code=True)
def _(WidgetEdit, dbinit, gmo):
    db = gmo.WgetEdit(WidgetEdit, dbinit)
    return (db,)


@app.cell(hide_code=True)
def _(db, gmo, gp):
    def myplot():
        fig, ax = gp.init(figsize=(4, 4))
        gmo.plotData(ax, db, name="z")
        return fig

    return (myplot,)


@app.cell(hide_code=True)
def _(WidgetDb, WidgetEdit, gmo, mo, myplot):
    param = mo.ui.tabs(
        {"Data": gmo.WshowDb(WidgetDb), "Edit": gmo.WshowEdit(WidgetEdit)}
    ).style({"minWidth": "700px", "width": "350px"})

    simu = mo.vstack(
        [mo.md(""), mo.md(f"Plotting the Data Base:{mo.as_html(myplot())}")], gap=4
    )

    mo.hstack([param, simu], gap=4)
    return


if __name__ == "__main__":
    app.run()
