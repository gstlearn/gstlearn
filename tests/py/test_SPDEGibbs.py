import gstlearn as gl
import numpy as np

gl.law_set_random_seed(67)

# Define data
db = gl.Db()
db["x1"] = [3.25, 1.4, 5, 5.0]
db["x2"] = [2.6, 4.5, 3.0, 8.0]
db.setLocators(["x*"], gl.ELoc.X)
ndat = db.getNSample()

db["value"] = [0, 1, 2, -5.0]
db.setLocators(["value"], gl.ELoc.Z)

db["l_bound"] = [np.nan, -np.inf, -5, 10.0]
db["u_bound"] = [np.nan, -5, 10, np.inf]

db["err"] = [0.001, 0.002, 0.003, 0.004]
db.setLocators(["err"], gl.ELoc.Z)

# Define SPDE parameters
model = gl.Model.createFromParam(gl.ECov.MATERN, range=10, param=1)
modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET, sills=1e-2)  # default 0.01
modelNugg.getCovAniso(0).makeSillNoStatDb("err", 0, 0, db)
invnuggetop = gl.InvNuggetOp(db, modelNugg)

# SPDE mesh
x_max = 10
nx = 20
dx = x_max / (nx - 1)
dbgrid = gl.DbGrid.create(nx=[nx, nx], dx=[dx, dx])
mesh = gl.MeshETurbo(dbgrid)
meshes = gl.VectorMeshes([mesh])
proj_matrix = gl.ProjMatrix(db, mesh)

# SPDE operator
spdeChol = gl.SPDE(db, model, True, False, 1)
spdeChol.setMeshes(True, meshes)
spdeChol.setInvNoise(invnuggetop)
err = spdeChol.makeReady(True)
print("makeReady:", err)

SpdeOp = spdeChol.getSPDEOp()

# Gibbs simulations only on data
out_gibbs_data = SpdeOp.simCondGibbs(
    db["value"],
    db["l_bound"],
    db["u_bound"],
    projK=proj_matrix,
    projS=proj_matrix,
    nIter=5,
)
print("--- out_gibbs_data ---")
print(np.round(out_gibbs_data, 3))

# Gibbs simulations only on all mesh
out_gibbs_mesh = SpdeOp.simCondGibbs(db["value"], db["l_bound"], db["u_bound"], nIter=5)
print("--- out_gibbs_mesh ---")
print(np.round(out_gibbs_mesh, 3))
