# %%
import gstlearn as gl
import matplotlib.pyplot as plt
import numpy as np

# %% General parameters
flag_plot = True
ndim = 2
nvar = 2
nbsimu = 2

# %% Model (multivariate) for the field 
if nvar == 1:
    sills = np.array([1.])
    epsNugget = 0.0001
    sillsNugg = np.array([0.0001])

if nvar == 2:
    sills = np.array([[80,-50],[-50,40]])
    epsNugget = 0.001
    sillsNugg = np.array([[.001,0.],[.0,.001]])

model = gl.Model.createFromParam(gl.ECov.MATERN,param = 1, range=20,
                                 sills = sills)
modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET, sills = sillsNugg)

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

print("Simulation using gstlearn (with Matrix) -> SGM")
gl.law_set_random_seed(1242)
params = gl.SPDEParam.create(epsNugget = epsNugget)
resultMat = gl.simulateSPDENew(dat, grid, model, meshes, nbsimu, 1, params)

SimuGM = {}
for i in range(nbsimu):
    SimuGM[i+1] = {}
    SimuGM[i+1][1] = resultMat[i][0:nt]
    if nvar == 2:
        SimuGM[i+1][2] = resultMat[i][(nt):(2*nt)]

# Printing Statistics
print("Mean of the gstlearn simulations")
for i in range(nbsimu):
    for j in range(nvar):
        print(np.round(SimuGM[i+1][j+1].mean(),4))
        
# %% Simulation (Matrix-free) performed with gstlearn

print("Simulation using gstlearn (Matrix Free) -> SGF")
gl.law_set_random_seed(1242)
params = gl.SPDEParam.create(epsNugget = epsNugget)
resultMat = gl.simulateSPDENew(dat, grid, model, meshes, nbsimu, 0, params)

SimuGF = {}
for i in range(nbsimu):
    SimuGF[i+1] = {}
    SimuGF[i+1][1] = resultMat[i][0:nt]
    if nvar == 2:
        SimuGF[i+1][2] = resultMat[i][(nt):(2*nt)]

# Printing Statistics
print("Mean of the gstlearn simulations")
for i in range(nbsimu):
    for j in range(nvar):
        print(np.round(SimuGF[i+1][j+1].mean(),4))

# %% Simulation (with Matrix) is performed by hand

print("Simulation using gstlearperformed by Hand (with Matrix) -> SHF")
gl.law_set_random_seed(1242)
Qop       = gl.PrecisionOpMulti(model,meshes,True)
invnoise  = gl.buildInvNugget(dat,modelNugg)
invnoisep = gl.MatrixSymmetricSim(invnoise)
spdeop    = gl.SPDEOp(Qop, AM, invnoisep)
resultMatH = spdeop.simCond(Z)
resultMatH1 = spdeop.simCond(Z)

SimuHF = {}
SimuHF[1] = {}
SimuHF[1][1] = resultMatH[0:nt]
if nvar == 2:
    SimuHF[1][2] = resultMatH[(nt):(2*nt)]
if nbsimu == 2:
    SimuHF[2] = {}
    SimuHF[2][1] = resultMatH1[0:nt]
    if nvar == 2:
        SimuHF[2][2] = resultMatH1[(nt):(2*nt)]

# Printing Statistics
print("Mean of the manual simulations")
for i in range(nbsimu):
    for j in range(nvar):
        print(np.round(SimuHF[i+1][j+1].mean(),4))

# %% Various plots

if flag_plot:

    # Display the GF simulations for all variables and all simulations
    for i in range(nbsimu):
        for j in range(nvar):
            plt.imshow(SimuGF[i+1][j+1].reshape(grid.getNXs()))
            plt.title("SGF#"+str(i+1)+" V#"+str(j+1)+" (gstlearn)")
            plt.show()

    # Display the GM simulations for all variables and all simulations
    for i in range(nbsimu):
        for j in range(nvar):
            plt.imshow(SimuGM[i+1][j+1].reshape(grid.getNXs()))
            plt.title("SGM#"+str(i+1)+" V#"+str(j+1)+" (gstlearn)")
            plt.show()
    
    # Checking the simulations are different per variable
    if nvar == 2:
        plt.scatter(SimuGF[1][1],SimuGF[1][2],s=1)
        plt.title("Comparing two Variables for SGF#1 (gstlearn)")
        plt.xlabel("Variable #1")
        plt.ylabel("Variable #2")
        plt.show()

    # Checking the simulations are different per simulation
    if nbsimu == 2:
        plt.scatter(SimuGF[1][1], SimuGF[2][1], s=1)
        plt.title("Comparing two Simulations SGF for V#1 (gstlearn)")
        plt.xlabel("Simulation #1")
        plt.ylabel("Simulation #2")
        plt.show()
    
    # Comparing Manual to SGF simulations for all variables / simulations
    for i in range(nbsimu):
        for j in range(nvar):
            plt.scatter(SimuHF[i+1][j+1], SimuGF[i+1][j+1], s=1)
            plt.title("Comparing Simulation SGF#"+str(i+1)+" V#"+str(j+1)+" to SHF")
            plt.xlabel("(Manual SHF)")
            plt.ylabel("(gstlearn GF)")
            plt.show()

# %% Checking the exactness of conditional simulations
grid["v1"] = SimuGF[1][1]
if nvar == 2:
    grid["v2"] = SimuGF[1][2]

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
