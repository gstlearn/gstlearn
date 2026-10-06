import marimo

__generated_with = "0.25.0"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def my_imports():
    import marimo as mo

    import gstlearn as gl
    import gstlearn.gstmarimo as gmo
    import gstlearn.plot as gp
    import gstlearn.document as gdoc
    import matplotlib.pyplot as plt

    import numpy as np
    import pandas as pd
    from scipy.spatial.distance import cdist

    gmo.setEnvironment(optionBackup=True, optionDisplay=False)
    return cdist, gdoc, gl, gmo, gp, mo, np, plt


@app.cell(hide_code=True)
def widget_definition(gmo):
    # Cloud widget definition
    WidgetCloud = gmo.WdefineCloud(
        nlag=10,
        lag=70.0,
        tol=50.0,
        angle=125.0,
        tolangle=20.0,
        cylrad=0.0,
        as_grid=False,
        lagindex=0,
    )

    # Plot choices list with default statuses
    options_list = [
        ("data", "Data Base", True),
        ("mapcloud", "Cloud as a Map", True),
        ("vcloud", "Variogram Cloud", True),
        ("vario", "Experimental Omnidirectional Variogram", True),
        ("vmap", "Variogram Map", True),
        ("vmapNb", "Variogram Map (Nb of pairs)", True),
    ]

    WidgetLayout = gmo.WdefineLayout(
        options_list,
        nrow=2,
        ncol=3,
        width=3,
        height=3,
    )

    WidgetAutoSave = gmo.WdefineAutoSave()
    return WidgetAutoSave, WidgetCloud, WidgetLayout, options_list


@app.cell(hide_code=True)
def read_data(cdist, gdoc, gl, np):
    temp_nf = gdoc.loadData("Scotland", "Scotland_Temperatures.NF")
    db = gl.Db.createFromNF(temp_nf)
    db.setLocator("Janvier", gl.ELoc.Z)

    # Retrieve Coordinates and Variable from db
    zp = db["*temp"]

    # Reduce to defined samples
    ind = np.argwhere(np.invert(np.isnan(zp)))
    ztab = np.atleast_2d(
        zp[ind].reshape(
            -1,
        )
    )

    # Matrix of coordinates
    xtab = db["Longitude"][ind].reshape(
        -1,
    )
    ytab = db["Latitude"][ind].reshape(
        -1,
    )
    coords = np.array([xtab, ytab]).T

    # Matrix of distance flattened in a vector
    dist = cdist(coords, coords, metric="euclidean").reshape(
        -1,
    )

    # Matrix of gamma_{ij} flattened in a vector
    gamma = (
        0.5
        * (ztab - ztab.T).reshape(
            -1,
        )
        ** 2
    )
    variance = np.var(ztab)

    # Compute vectors
    deltaX = np.atleast_2d(xtab) - np.atleast_2d(xtab).T
    deltaY = np.atleast_2d(ytab) - np.atleast_2d(ytab).T
    deltax = deltaX.reshape(
        -1,
    )
    deltaY = deltaY.reshape(
        -1,
    )
    return db, deltaX, deltaY, dist, gamma, variance


@app.cell(hide_code=True)
def _(gl, np):
    def getSelectedPairs(varioparam, db, lagindex):
        selected = gl.buildDbFromVarioParam(db, varioparam)
        ind = np.where(selected["Lag"] == lagindex)
        s1 = selected["Sample-1"][ind]
        s2 = selected["Sample-2"][ind]
        vect = db.getIncrements(s1, s2)
        return vect[0], vect[1]

    return (getSelectedPairs,)


@app.cell(hide_code=True)
def define_action(
    WidgetAutoSave,
    WidgetCloud,
    WidgetLayout,
    db,
    deltaX,
    deltaY,
    dist,
    gamma,
    getSelectedPairs,
    gl,
    gmo,
    gp,
    mo,
    np,
    options_list,
    plt,
    variance,
):
    def myaction():

        # Retrieve parameter values from widgets
        autosave = gmo.WgetAutoSave(WidgetAutoSave)
        cloud_params = gmo.WgetCloud(WidgetCloud)
        layout = gmo.WgetLayout(WidgetLayout, options_list)

        # Retrieve grid layout configuration
        nx = layout["nx"]
        ny = layout["ny"]
        dimx = layout["dimx"]
        dimy = layout["dimy"]

        # Variogram calculation parameters
        nlag = cloud_params["nlag"]
        lag = cloud_params["lag"]
        tol = cloud_params["tol"]
        omni = cloud_params["omni"]
        angle = cloud_params["angle"]
        cylrad = cloud_params["cylrad"]
        tolangle = cloud_params["tolangle"]
        as_grid = cloud_params["as_grid"]
        lagindex = cloud_params["lagindex"]

        # Compute the directional variogram
        if omni:
            varioParam = gl.VarioParam.createOmniDirection(nlag, lag, tol)
        else:
            direction1 = gl.DirParam.create(
                nlag, lag, tol, tolangle, angle2D=angle, cylrad=cylrad
            )
            direction2 = gl.DirParam.create(
                nlag, lag, tol, tolangle, angle2D=angle + 90, cylrad=cylrad
            )
            varioParam = gl.VarioParam()
            varioParam.addDir(direction1)
            varioParam.addDir(direction2)
        vario = gl.Vario.computeFromDb(varioParam, db)

        # Compute the variogram Map
        nxy = 20
        vmap = gl.vmapFromDb(db, nxx=[nxy, nxy])

        # Initialize figure layout
        fig, ax = gp.init(nx, ny, figsize=[ny * dimx, nx * dimy], squeeze=False)
        axes = ax.ravel()
        bound = 500

        flagOneLag = lagindex > 0 and lagindex < nlag
        flagCylinder = cylrad > 0

        i = 0
        for content in layout.get("selected", []):
            if content == "data":
                if i < len(axes):
                    gmo.plotData(axes[i], db)
                    axes[i].decoration(title="Data Base")
                    i += 1

            elif content == "mapcloud":
                if i < len(axes):
                    plt.sca(axes[i])
                    axes[i].scatter(deltaX, deltaY, s=0.01)
                    axes[i].scatter(0, 0)

                    # Highlight the pairs corresponding to lagindex
                    if flagOneLag:
                        xs, ys = getSelectedPairs(varioParam, db, lagindex)
                        axes[i].scatter(
                            xs, ys, s=0.08, label="Selected pairs", c="yellow"
                        )

                    # Draw the circle representation
                    for j in range(nlag):
                        mini, center, maxi = gp.lagDefine(j, lag, tol)
                        gp.drawCircles(mini, maxi)

                    # Draw the Cone representation
                    gp.drawDir(angle, "black")
                    gp.drawDir(angle + tolangle, "red")
                    gp.drawDir(angle - tolangle, "red")

                    # Draw the Cylinder representation
                    if flagCylinder:
                        gp.drawCylrad(angle, cylrad, col="purple")

                    axes[i].axhline(y=0, color="black")
                    axes[i].axvline(x=0, color="black")
                    axes[i].set_adjustable("box")
                    axes[i].set_xlim(left=-bound, right=bound)
                    axes[i].set_ylim(bottom=-bound, top=bound)
                    axes[i].decoration(title="Data Cloud Map")
                    i += 1

            elif content == "vcloud":
                if i < len(axes):
                    axes[i].clear()
                    if as_grid:
                        axes[i].hist2d(dist, gamma, bins=50)
                    else:
                        axes[i].scatter(dist, gamma, s=0.01)
                        for j in range(nlag):
                            mini, center, maxi = gp.lagDefine(j, lag, 0.0)
                            axes[i].axvline(x=center, ymin=0.0, c="black")
                    axes[i].axhline(y=variance, color="black", linestyle="dashed")
                    axes[i].set_adjustable("box")
                    axes[i].set_xlim(left=0)
                    axes[i].set_ylim(bottom=0)
                    axes[i].decoration(title="Variogram Cloud")
                    i += 1

            elif content == "vario":
                if i < len(axes):
                    axes[i].variogram(vario, idir=-1, showPairs=True)
                    axes[i].decoration(title="Exp. Variogram (both directions)")
                    i += 1

            elif content == "vmap":
                if i < len(axes):
                    plt.sca(axes[i])
                    axes[i].raster(vmap, name="*.Var", flagLegend=True)
                    mx = np.min(vmap["x1"]) - vmap.getDX(0) / 2
                    Mx = np.max(vmap["x1"]) + vmap.getDX(0) / 2
                    my = np.min(vmap["x2"]) - vmap.getDX(1) / 2
                    My = np.max(vmap["x2"]) + vmap.getDX(1) / 2

                    for j in range(vmap.getNX(0) + 1):
                        x = vmap.getX0(0) + (j - 0.5) * vmap.getDX(0)
                        axes[i].plot([x, x], [my, My], color="black", linewidth=0.5)
                    for j in range(vmap.getNX(1) + 1):
                        y = vmap.getX0(1) + (j - 0.5) * vmap.getDX(1)
                        axes[i].plot([mx, Mx], [y, y], color="black", linewidth=0.5)
                    axes[i].set_adjustable("box")
                    axes[i].set_xlim(left=-500, right=500)
                    axes[i].set_ylim(bottom=-500, top=500)
                    axes[i].decoration(title="Variogram Map")
                    i += 1

            elif content == "vmapNb":
                if i < len(axes):
                    plt.sca(axes[i])
                    axes[i].raster(vmap, name="*.Nb", flagLegend=True)
                    mx = np.min(vmap["x1"]) - vmap.getDX(0) / 2
                    Mx = np.max(vmap["x1"]) + vmap.getDX(0) / 2
                    my = np.min(vmap["x2"]) - vmap.getDX(1) / 2
                    My = np.max(vmap["x2"]) + vmap.getDX(1) / 2

                    for j in range(vmap.getNX(0) + 1):
                        x = vmap.getX0(0) + (j - 0.5) * vmap.getDX(0)
                        axes[i].plot([x, x], [my, My], color="black", linewidth=0.5)
                    for j in range(vmap.getNX(1) + 1):
                        y = vmap.getX0(1) + (j - 0.5) * vmap.getDX(1)
                        axes[i].plot([mx, Mx], [y, y], color="black", linewidth=0.5)
                    axes[i].set_adjustable("box")
                    axes[i].set_xlim(left=-500, right=500)
                    axes[i].set_ylim(bottom=-500, top=500)
                    axes[i].decoration(title="Variogram Map (Nb of pairs)")
                    i += 1

        # Hide unused subplot axes
        for axi in axes[i:]:
            axi.axis("off")

        plt.tight_layout(pad=0.2)
        mo.mpl.interactive(fig)

        return fig

    return (myaction,)


@app.cell(hide_code=True)
def render_ui(
    WidgetAutoSave,
    WidgetCloud,
    WidgetLayout,
    gmo,
    mo,
    myaction,
    options_list,
):
    param = mo.ui.tabs(
        {
            "Cloud": gmo.WshowCloud(WidgetCloud),
            "Layout": gmo.WshowLayout(WidgetLayout, options_list),
            "AutoSave": gmo.WshowAutoSave(WidgetAutoSave),
        }
    ).style({"minWidth": "400px", "width": "400px"})

    simu = mo.vstack([mo.md(""), mo.md(f"{mo.as_html(myaction())}")], gap=0)

    global_layout = mo.hstack([param, simu], gap=4)
    global_layout
    return


if __name__ == "__main__":
    app.run()
