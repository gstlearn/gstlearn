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
def _(WidgetDb, gmo, gp):
    def myplot():
        db = gmo.WgetDb(WidgetDb)

        fig = None
        if db is not None:
            fig, ax = gp.init(figsize=[4,4])
            ax.symbol(db)
        return fig
    return (myplot,)


@app.cell(hide_code=True)
def _(WidgetDb, gmo, mo, myplot):
    param = mo.ui.tabs(
        {
            "Data":       gmo.WshowDb(WidgetDb),
        }
    ).style({"minWidth": "350px", "width": "350px"})

    simu = mo.vstack(
        [
             mo.md(""),
             mo.md(f"Plotting the Data Base:{mo.as_html(myplot())}")
        ],
        gap = 4
    )

    mo.hstack([param, simu], gap=4)
    return param, simu


if __name__ == "__main__":
    app.run()
