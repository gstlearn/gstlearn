import marimo

__generated_with = "0.11.25"
app = marimo.App()


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import altair as alt

    import gstlearn as gl
    import gstlearn.plot as gp
    import gstlearn.gstmarimo as gmo
    import gstlearn.document as gdoc

    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    from IPython.display import Markdown
    return Markdown, alt, copy, gdoc, gl, gmo, gp, mo, np, plt


@app.cell(hide_code=True)
def _(gmo):
    WidgetDb = gmo.WdefineDb()
    return (WidgetDb,)


@app.cell(hide_code=True)
def _(WidgetDb, gmo):
    db = gmo.WgetDb(WidgetDb)
    return (db,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetVario = gmo.WdefineVario()
    return (WidgetVario,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetCovList = gmo.WdefineCovList()
    return (WidgetCovList,)


@app.cell
def _(WidgetCovList, WidgetDb, WidgetVario, gmo, gp):
    def myaction():

        # Define the data base
        db = gmo.WgetDb(WidgetDb)
        if db is None:
            return

        # Define the Variogram
        vario = gmo.WgetVario(WidgetVario, db)
        if vario is None:
            return

        # Define the Model
        model = gmo.WgetCovList(WidgetCovList, vario)
        if model is None:
            return

        fig, ax = gp.init(2,1,figsize=(10,16))
        ax[0,0].symbol(db)
        ax[0,0].decoration(title="Data Location")

        if vario is not None:
            if model is None:
                ax[1,0].variogram(vario, idir=-1)
            else:
                ax[1,0].varmod(vario, model)
            ax[1,0].decoration(title="Variogram and Model")

        return fig
    return (myaction,)


@app.cell(hide_code=True)
def _(WidgetCovList, WidgetDb, WidgetVario, gmo, mo, myaction):
    param = mo.ui.tabs(
        {
            "Data":       gmo.WshowDb(WidgetDb),
            "Variogram":  gmo.WshowVario(WidgetVario),
            "Model":      gmo.WshowCovList(WidgetCovList),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack(
        [
             mo.md(""),
             mo.md(f"Variogram and fitted Model{mo.as_html(myaction())}")
        ],
        gap = 0
    )

    mo.hstack([param, simu], gap=2)
    return param, simu


if __name__ == "__main__":
    app.run()
