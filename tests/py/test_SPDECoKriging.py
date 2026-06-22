import gstlearn as gl
import gstlearn.plot as gp
import numpy as np

# %% General parameters
flag_plot = False
ndim = 2
nvar = 2
order = 1
gl.OptCst.defineByKey("NTROW", -1)
gl.ASerializable.setPrefixName("test_SPDECoKriging-")
exts = ["estim", "stdev"]
params = gl.SPDEParam()
params.setNMC(100)

# %% Creating the Model
model = gl.Model(nvar, ndim)
model.addCovFromParam(
    gl.ECov.MATERN, param=1, range=10, sills=np.array([[50, 10], [10, 40]])
)
model.addCovFromParam(gl.ECov.MATERN, param=1, range=20, sills=([[30, -10], [-10, 10]]))
model.addCovFromParam(gl.ECov.NUGGET, sills=np.array([[4, -2], [-2, 6]]))
totalSill = np.max(model.getTotalSills().toTL())

means = np.zeros(nvar)
if order < 0:
    means = [12, -24]

model.setDriftIRF(order)
model.setMeans(means)

# %% Data Bases
dat = gl.Db.createFillRandom(
    ndat=100, ndim=2, nvar=0, coormin=[0, 0], coormax=[100, 100], seed=234
)
err = gl.simtub(None, dat, model, namconv=gl.NamingConvention("Data"))
if order == 1:
    coeffX = +0.5
    coeffY = -0.5
    dat["Data.V1"] = dat["Data.V1"] + dat["x-1"] * coeffX + means[0]
    dat["Data.V1"] = dat["Data.V1"] + dat["x-2"] * coeffY + means[1]

# %% Output Grid
grid = gl.DbGrid.create(nx=[50, 50], dx=[2, 2])

# %% Meshing
nx1 = [70, 70]
mesh1 = gl.MeshETurbo(nx1, [2.0, 2.0], [-20, -20], [], False)
meshes = gl.VectorMeshes([mesh1, mesh1])

############################################
# %% Traditional Kriging (used as Reference)
############################################

gl.mestitle(1, "Co-Kriging (traditional)")
err = gl.kriging(dat, grid, model, gl.NeighUnique(), flag_std=True)
gl.dbStatisticsMono(grid, ["Kriging.*"]).display()

##################################
# %% SPDE Kriging (matrix version)
##################################

gl.mestitle(1, "Co-Kriging using SPDE (Matrix)")
err = gl.krigingSPDE(
    dat,
    grid,
    model,
    True,
    True,
    1,
    meshes,
    params=params,
    namconv=gl.NamingConvention("KM"),
)
gl.dbStatisticsMono(grid, ["KM.*"]).display()
for iext in exts:
    for ivar in range(nvar):
        name1 = gl.NC.getNameEncoded("Kriging", dat, ivar, nvar, 1, 1, iext)
        name2 = gl.NC.getNameEncoded("KM", dat, ivar, nvar, 1, 1, iext)
        print(
            "Difference with",
            name1,
            "(Matrix) = "
            + str(np.round(np.max(np.abs(grid[name1] - grid[name2])) / totalSill, 5)),
        )

#######################################
# %% SPDE Kriging (matrix-free version)
#######################################

gl.mestitle(1, "Co-Kriging using SPDE (Matrix-Free)")
err = gl.krigingSPDE(
    dat,
    grid,
    model,
    True,
    True,
    0,
    meshes,
    params=params,
    namconv=gl.NamingConvention("KF"),
)
gl.dbStatisticsMono(grid, ["KF.*"]).display()
for iext in exts:
    for ivar in range(nvar):
        name1 = gl.NC.getNameEncoded("Kriging", dat, ivar, nvar, 1, 1, iext)
        name2 = gl.NC.getNameEncoded("KF", dat, ivar, nvar, 1, 1, iext)
        print(
            "Difference with",
            name1,
            " = "
            + str(np.round(np.max(np.abs(grid[name1] - grid[name2])) / totalSill, 5)),
        )

for iext in exts:
    for ivar in range(nvar):
        name1 = gl.NC.getNameEncoded("KF", dat, ivar, nvar, 1, 1, iext)
        name2 = gl.NC.getNameEncoded("KM", dat, ivar, nvar, 1, 1, iext)
        print(
            "Difference between Matrix-Free and Matrix ("
            + iext.replace(".", "")
            + ") for variable "
            + iext,
            ivar + 1,
            " = "
            + str(np.round(np.max(np.abs(grid[name1] - grid[name2])) / totalSill, 5)),
        )
grid.dumpToNF("grid.NF")

###############
# Various plots
###############
if flag_plot:
    # Display the result per variable for Traditional Kriging
    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            name = gl.NC.getNameEncoded("Kriging", dat, ivar, nvar, 1, 1, iext)
            gp.raster(grid, name)
            gp.decoration(title=name + " (Traditional)")
            gp.close()

    # Display the result per variable for SPDE Kriging (matrix)
    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            name = gl.NC.getNameEncoded("KM", dat, ivar, nvar, 1, 1, iext)
            gp.raster(grid, name)
            gp.decoration(title=name + " (SPDE Matrix)")
            gp.close()

    # Display the result per variable for SPDE Kriging (matrix-free)
    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init(flagEqual=True)
            name = gl.NC.getNameEncoded("KF", dat, ivar, nvar, 1, 1, iext)
            gp.raster(grid, name)
            gp.decoration(title=name + " (SPDE Matrix-Free)")
            gp.close()

    # Comparing the Krigings
    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init()
            gp.correlation(
                grid,
                gl.NC.getNameEncoded("KM", dat, ivar, nvar, 1, 1, iext),
                gl.NC.getNameEncoded("Kriging", dat, ivar, nvar, 1, 1, iext),
                bissLine=True,
                bissColor="blue",
                bins=100,
                cmin=1,
            )
            gp.decoration(
                title="Comparing Kriging" + iext + " - Variable " + str(ivar + 1),
                xlabel="SPDE (Matrix)",
                ylabel="Traditional",
            )
            gp.close()

    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init()
            gp.correlation(
                grid,
                gl.NC.getNameEncoded("KF", dat, ivar, nvar, 1, 1, iext),
                gl.NC.getNameEncoded("Kriging", dat, ivar, nvar, 1, 1, iext),
                bissLine=True,
                bissColor="blue",
                bins=100,
                cmin=1,
            )
            gp.decoration(
                title="Comparing Kriging" + iext + " - Variable " + str(ivar + 1),
                xlabel="SPDE (Matrix-Free)",
                ylabel="Traditional",
            )
            gp.close()

    for iext in exts:
        for ivar in range(nvar):
            fig, ax = gp.init()
            gp.correlation(
                grid,
                gl.NC.getNameEncoded("KF", dat, ivar, nvar, 1, 1, iext),
                gl.NC.getNameEncoded("KM", dat, ivar, nvar, 1, 1, iext),
                bissLine=True,
                bissColor="blue",
                bins=100,
                cmin=1,
            )
            gp.decoration(
                title="Comparing Matrix-Free/Matrix"
                + iext
                + " - Variable "
                + str(ivar + 1),
                xlabel="SPDE (Matrix-Free)",
                ylabel="SPDE (Matrix)",
            )
            gp.close()
