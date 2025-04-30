# %%
import gstlearn as gl
import matplotlib.pyplot as plt
import numpy as np

# %% General parameters
flag_plot = False
ndim = 2
nvar = 2
nbsimu = 2

# %% Model (multivariate) for the field 
if nvar == 1:
    sills = np.array([1.])
    epsNugget = 0.0001

if nvar == 2:
    sills = np.array([[80,-50],[-50,40]])
    epsNugget = 0.001

model = gl.Model.createFromParam(gl.ECov.MATERN,param = 1, range=20,
                                sills = sills)

# %% Data creation (2 variables)
dat = gl.Db.createFillRandom(ndat=100, ndim = 2,nvar = 0,
                             coormin=[20,20],coormax=[80,80],seed=234)

## Coordinates are rounded to check if simulations are exact
dat["x-1"] = np.round(dat["x-1"])
dat["x-2"] = np.round(dat["x-2"])
gl.simtub(None,dat,model,nbtuba = 1000)
Z = dat["Simu*"].T.reshape(-1) #Data vector (two variables in a single vector)

# %% Resolution Grid
grid = gl.DbGrid.create([141,141] ,dx = [1,1], x0 = [-20,-20])
nt = grid.getNSample()

# %% Meshing 
mesh1 = gl.MeshETurbo(grid, False)
meshes = gl.VectorMeshes([mesh1])

# %% Projections

# Projection operators
# (2 to mimic the possibility to have several convolution operators)
AM1 = gl.ProjMatrix(dat,mesh1)
AM2 = gl.ProjMatrix(dat,mesh1)

#Total projection operator
if nvar == 1:
    vectproj = gl.VVectorConstIProj([[AM1]])

if nvar == 2:
    vectproj = gl.VVectorConstIProj([[AM1,None],[None,AM2]])

AM = gl.ProjMulti(vectproj)

# %% Simulation (Matrix) performed with gstlearn

gl.law_set_random_seed(1242)
params = gl.SPDEParam.create(epsNugget = epsNugget)
resultMatGM = gl.simulateSPDENew(dat, grid, model, meshes, nbsimu, 0, params)

SimuGM = {}
for i in range(nbsimu):
    SimuGM[i+1] = {}
    SimuGM[i+1][1] = resultMatGM[i][0:nt]
    if nvar == 2:
        SimuGM[i+1][2] = resultMatGM[i][(nt):(2*nt)]

# Printing Statistics
print("Mean of the gstlearn simulations")
for i in range(nbsimu):
    for j in range(nvar):
        print(np.round(SimuGM[i+1][j+1].mean(),4))

# %% Simulation (with Matrix) is performed by hand

#Noise Operator (here from a nugget but could be from another SPDE)
if nvar == 1:
    sills = np.array([0.0001])

if nvar == 2:
    sills = np.array([[.001,0.],[.0,.001]])

gl.law_set_random_seed(1242)
modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET, sills=sills)
Qop       = gl.PrecisionOpMulti(model,meshes,True)
invnoise  = gl.buildInvNugget(dat,modelNugg)
invnoisep = gl.MatrixSymmetricSim(invnoise)
spdeop    = gl.SPDEOp(Qop, AM, invnoisep)
resultMatH = spdeop.simCond(Z)
resultMatH1 = spdeop.simCond(Z)

SimuH = {}
SimuH[1] = {}
SimuH[1][1] = resultMatH[0:nt]
if nvar == 2:
    SimuH[1][2] = resultMatH[(nt):(2*nt)]
if nbsimu == 2:
    SimuH[2] = {}
    SimuH[2][1] = resultMatH1[0:nt]
    if nvar == 2:
        SimuH[2][2] = resultMatH1[(nt):(2*nt)]

# Printing Statistics
print("Mean of the manual simulations")
for i in range(nbsimu):
    for j in range(nvar):
        print(np.round(SimuH[i+1][j+1].mean(),4))

# %% Various plots

if flag_plot:

    # Display the simulations for all variables and all simulations
    for i in range(nbsimu):
        for j in range(nvar):
            plt.imshow(SimuGM[i+1][j+1].reshape(grid.getNXs()))
            plt.title("S#"+str(i+1)+" V#"+str(j+1)+" (gstlearn)")
            plt.show()
    
    # Checking the simulations are different per variable
    if nvar == 2:
        plt.scatter(SimuGM[1][1],SimuGM[1][2],s=1)
        plt.title("Comparing two Variables for S#1 (gstlearn)")
        plt.show()

    # Checking the simulations are different per simulation
    if nbsimu == 2:
        plt.scatter(SimuGM[1][1], SimuGM[2][1], s=1)
        plt.title("Comparing two Simulations for V#1 (gstlearn)")
        plt.show()
    
    # Comparing Manual and Gstlearn simulations for all variables / simulations
    for i in range(nbsimu):
        for j in range(nvar):
            plt.scatter(SimuH[i+1][j+1], SimuGM[i+1][j+1], s=1)
            plt.title("Comparing Simulations S#"+str(i+1)+" V#"+str(j+1))
            plt.xlabel("(Manual)")
            plt.ylabel("(gstlearn)")
            plt.show()

# %% Checking the exactness of conditional simulations
grid["v1"] = SimuGM[1][1]
if nvar == 2:
    grid["v2"] = SimuGM[1][2]

# %%
gl.migrate(grid,dat,"v1")
dat.setName("Migrate","m1")
if nvar == 2 :
    gl.migrate(grid,dat,"v2")
    dat.setName("Migrate","m2")

# %%
if flag_plot:
    if nvar == 1:
        ax = plt.scatter(dat["Simu"], dat["m1"], s=1)
        plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
        plt.title("Checking Data vs. Simu. Cond. V#1 (gstlearn)")
        plt.show()
    
    if nvar == 2:
        ax = plt.scatter(dat["Simu.1"],dat["m1"],s=1)
        plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
        plt.title("Checking Data vs. Simu. Cond. V#1 (gstlearn)")
        plt.show()
        
        plt.scatter(dat["Simu.2"],dat["m2"],s=1)
        plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
        plt.title("Checking Data vs. Simu. Cond. V#2 (gstlearn)")
        plt.show()
