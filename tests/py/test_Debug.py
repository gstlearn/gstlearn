# This test is meant to check the Global estimation
# either using the Arithmetic Mean algorithm or the Kriging one
# The global estimate returns a whole structure which is partly printed here
import gstlearn as gl

np = 100
seed = 14313
data = gl.Db.createFillRandom(np, 2, 0)
model = gl.Model.createFromParam(gl.ECov.EXPONENTIAL, 0.5, 1.0)
model.setMeans(0)
model.setDriftIRF(-1)
gl.simtub(None, data, model, None, 1, seed)
nburn = 10
niter = 30
nbsimu = 5
gl.gibbs_sampler(
    data,
    model,
    nbsimu,
    seed,
    nburn,
    niter,
    True,
    True,
    False,
    False,
    False,
    0,
    5.0,
    False,
    False,
    True,
)
data.getContentsAsTable().display()
