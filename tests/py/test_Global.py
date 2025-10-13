# This test is meant to check the Global estimation 
# either using the Arithmetic Mean algorithm or the Kriging one
# The global estimate returns a whole structure which is partly printed here
import gstlearn as gl

# Define a 2-D Data Set in [0,1)x[0,1]
ndata = 10
db = gl.Db.createFillRandom(ndata)
db.display()

# Define a (coarse) grid covering the square
nx = 5
grid = gl.DbGrid.create([nx,nx], [1/nx, 1/nx])
grid.display()

# Define a Model
range = 0.3
model = gl.Model.createFromParam(gl.ECov.SPHERICAL, range)
model.display()

verbose = True

# Arithmetic method: we also print
Gres = gl.global_arithmetic(db, grid, model, 0, verbose)
print("Printing some results from the Output structure")
print("Number of weights = ", Gres.ntot)
print("Weights = ", Gres.weights)

# Kriging solution
Gres = gl.global_kriging(db, grid, model, 0, verbose)
print("Printing some results from the Output structure")
print("Number of weights = ", Gres.ntot)
print("Weights = ", Gres.weights)

