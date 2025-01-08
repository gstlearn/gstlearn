import marimo

__generated_with = "0.10.9"
app = marimo.App()


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import gstlearn as gl
    import gstlearn.plot as gp
    import matplotlib.pyplot as plt
    import contextily as ctx

    import numpy as np
    import pandas as pd
    return ctx, gl, gp, mo, np, pd, plt


@app.cell(hide_code=True)
def _(mo):
    b = mo.ui.file_browser()
    mo.md( f"""
        Choose your file {b}
        """)
    return (b,)


@app.cell(hide_code=True)
def _(b, gl):
    db = gl.Db.createFromNF(b.value[0].path)
    box = db.getExtremas()
    return box, db


@app.cell(hide_code=True)
def _(box, mo):
    value = [box[0,0], box[0,1]]
    Wlongitude = mo.ui.range_slider(start=-20, stop=20, step=0.1, value=value, 
                                      label="Longitude")
    value = [box[1,0], box[1,1]]
    Wlatitude = mo.ui.range_slider(start=30, stop=50, step=0.1, value=value, 
                                      label="Latitude")

    mo.vstack([mo.md("Define the Area to be displayed"), 
               Wlongitude, Wlatitude])
    return Wlatitude, Wlongitude, value


@app.cell(hide_code=True)
def _(Wlatitude, Wlongitude, box, ctx, db, plt):
    box[0,0] = Wlongitude.value[0]
    box[0,1] = Wlongitude.value[1]
    box[1,0] = Wlatitude.value[0]
    box[1,1] = Wlatitude.value[1]

    fig = plt.figure(figsize=(20,6))

    ax1 = fig.add_subplot(1,2,1)
    ax1.baseMap(db)
    ctx.add_basemap(ax1, source=ctx.providers.OpenStreetMap.Mapnik)
    ax1.decoration(title="Based on Db contents")

    ax2 = fig.add_subplot(1,2,2)
    ax2.baseMap(db, box=box)
    ctx.add_basemap(ax2, source=ctx.providers.OpenStreetMap.Mapnik)
    ax2.decoration(title="Based on your Area Definition")

    plt.show()
    return ax1, ax2, fig


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
