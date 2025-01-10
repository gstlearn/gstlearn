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
def _(WDbLayout, WidgetDb, gmo, gp):
    def myplot():
        db = gmo.WgetDb(WidgetDb, WDbLayout.value)

        fig = None
        if db is not None:
            fig, ax = gp.initGeographic(figsize=[4,4])
            ax.symbol(db)
        return fig
    return (myplot,)


@app.cell(hide_code=True)
def _(mo, myplot):
    mo.md(f"Plot of the Db {mo.as_html(myplot())}")
    return


if __name__ == "__main__":
    app.run()
