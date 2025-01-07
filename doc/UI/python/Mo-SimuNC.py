import marimo

__generated_with = "0.10.9"
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
    WDbLayout = gmo.WshowDb(WidgetDb)
    WDbLayout
    return (WDbLayout,)


@app.cell(hide_code=True)
def _(WDbLayout, WidgetDb, gmo):
    db = gmo.WgetDb(WidgetDb, WDbLayout.value)
    return (db,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetVario = gmo.WdefineVario()
    return (WidgetVario,)


@app.cell(hide_code=True)
def _(WidgetVario, gmo):
    WVarioLayout = gmo.WshowVario(WidgetVario)
    WVarioLayout
    return (WVarioLayout,)


@app.cell(hide_code=True)
def _(WVarioLayout, WidgetVario, db, gmo):
    vario = gmo.WgetVario(WidgetVario, WVarioLayout.value, db)
    return (vario,)


@app.cell(hide_code=True)
def _(gmo):
    WidgetCovList = gmo.WdefineCovList()
    return (WidgetCovList,)


@app.cell(hide_code=True)
def _(WidgetCovList, gmo):
    gmo.WshowCovList(WidgetCovList)
    return


@app.cell(hide_code=True)
def _(WidgetCovList, gmo, vario):
    model = gmo.WgetCovList(WidgetCovList, vario)
    return (model,)


@app.cell(hide_code=True)
def _(db, model, plt, vario):
    def myplot():

        fig = plt.figure(figsize=(20,6))
        if db is not None:
            ax1 = fig.add_subplot(1,2,1)
            ax1.symbol(db)

        if vario is not None and model is None:
            ax2 = fig.add_subplot(1,2,2)
            ax2.variogram(vario, idir=-1)

        if vario is not None and model is not None:
            ax2 = fig.add_subplot(1,2,2)
            ax2.varmod(vario, model)

        return fig
    return (myplot,)


@app.cell(hide_code=True)
def _(mo, myplot):
    mo.md(f"Data Base and Experimental Variogram {mo.as_html(myplot())}")
    return


if __name__ == "__main__":
    app.run()
