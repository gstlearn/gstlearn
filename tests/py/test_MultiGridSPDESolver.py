import gstlearn as gl
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

np.random.seed(4567)
dim = 2
model = gl.Model.createFromParam(gl.ECov.MATERN, param=0.5, range=1, flagRange=False)
gridcoarse = gl.DbGrid.create(nx=[25, 25, 25][0:dim], dx=[0.2, 0.2, 0.2][0:dim])
gridfine = gl.DbGrid.create(nx=[50, 50, 50][0:dim], dx=[0.1, 0.1, 0.1][0:dim])
mg = gl.MultiGridSPDE(model.getCovAniso(0))
gl.OptCustom.define("ompthreads", 10)
import time

time0 = time.time()
mat = mg.buildProlongator(gridfine, gridcoarse, 1)
mesh_c = gl.MeshETurbo(gridcoarse)


gl.OptCustom.define("ompthreads", 10)

alpha = 2
nx = [200, 200, 50][0:dim]
rangev = 90

gridfine = gl.DbGrid.create(nx=nx, dx=[1.0, 1.0, 1.0][0:dim])
gl.defineDefaultSpace(gl.ESpaceType.RN, dim)
model = gl.Model.createFromParam(
    gl.ECov.MATERN, param=alpha - dim / 2, range=rangev, flagRange=True
)
mesh = gl.MeshETurbo(gridfine)
Qtot = gl.PrecisionOp(mesh, model.getCovAniso(0), False)
truesol = np.random.normal(0, 1, size=gridfine.getNSample())
rhs = Qtot.evalDirect(truesol)

nlevels = 5
nrings = 2
niterpower = 10
ratiochebmin = 0.01
ratiochebmax = 1.5
smoothIter = 1

print("Starting multigrid solver...")
solver = gl.createMultiGridSolverSPDE(
    model.getCovAniso(0),
    gridfine,
    nlevels_max=nlevels,
    n_rings=nrings,
    niterPower=niterpower,
    ratiochebmin=ratiochebmin,
    ratiochebmax=ratiochebmax,
    smoothIter=smoothIter,
)

solver.getMaxEigenValues()
solver.display()


def mg(x):
    return solver.vCycle(x, 0 * x, 0)


history = []


def callback(xk):
    res = rhs - Qtot.evalDirect(xk)
    # history.append(np.linalg.norm(res) / np.linalg.norm(rhs))
    print("Residual norm:", np.round(np.linalg.norm(res) / np.linalg.norm(rhs), 4))
    history.append([1])
    print(len(history))
    return


M_op = LinearOperator(
    shape=[Qtot.getSize(), Qtot.getSize()], matvec=mg, dtype=np.float64
)
direct = LinearOperator(
    shape=[Qtot.getSize(), Qtot.getSize()],
    matvec=lambda x: Qtot.evalDirect(x),
    dtype=np.float64,
)

x, info = cg(
    direct, rhs, M=M_op, maxiter=2000, callback=callback, atol=1e-12, rtol=1e-6
)

n_iters = len(history)
if info == 0:
    print("-----------------------------------------------")
    print(f"Convergence reached in {n_iters} iterations.")
    print("-----------------------------------------------")
else:
    print(f"Failed after {n_iters} iterations.")
