
# %%
import gstlearn as gl
import matplotlib.pyplot as plt
import numpy as np
manual = False
nvar = 2
ndat = 100
rangeval = 20
param = 1
nx = [141,141]
dx = [1,1]
x0 = [-20,-20]
#Model (multivariate) for the field 
if nvar == 1:
    sills = np.array([1.])
if nvar == 2:
    sills = np.array([[80,-50],[-50,40]])


model = gl.Model.createFromParam(gl.ECov.MATERN,param = param, range=rangeval,
                                sills = sills)
#Data creation (2 variables)
dat = gl.Db.createFillRandom(ndat,ndim = 2,nvar = 0,coormin=[20,20],coormax=[80,80],seed=234)

## Coordinates are rounded to check if simulations are exact
dat["x-1"] = np.round(dat["x-1"])
dat["x-2"] = np.round(dat["x-2"])
gl.simtub(None,dat,model,nbtuba = 1000)
Z = dat["Simu*"].T.reshape(-1) #Data vector (two variables in a single vector)

#Resolution Grid 
grid = gl.DbGrid.create(nx = nx, dx = dx, x0 = x0)

#Meshes 
mesh1 = gl.MeshETurbo(grid, False)
print(mesh1)

# %%

#Projection operators (2 to mimic the possibility to have several convolution operators)
A1 = gl.ProjMatrix(dat,mesh1)
A2 = gl.ProjMatrix(dat,mesh1)

#Total projection operator
if nvar == 1:
    vectproj = gl.VVectorConstIProj([[A1]])
if nvar == 2:
    vectproj = gl.VVectorConstIProj([[A1,None],[None,A2]])

Atot = gl.ProjMulti(vectproj)

#Precision operator (multivariate)
vmesh = gl.VectorMeshes([mesh1])
stencil = True
Qtot = gl.PrecisionOpMulti(model,vmesh,stencil)


#Noise Operator (here from a nugget but could be from another SPDE)
#Note that it has to be multivariate if nvar > 1
if nvar == 1:
    modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET,
                                        sills = np.array([0.0001]))

if nvar == 2:
    modelNugg = gl.Model.createFromParam(gl.ECov.NUGGET,
                                        sills = np.array([[.001,0.],[.0,.001]]))


Qnoise = gl.buildInvNugget(dat,modelNugg)
NoiseOp = gl.MatrixSquareSymmetricSim(Qnoise)


# %%
spdeop = gl.SPDEOp(Qtot,Atot,NoiseOp)


# %%
result = spdeop.simCond(Z)

# %%
nt = grid.getNSample()
v1 = result[0:nt]
if nvar == 2:
    v2 = result[(nt):(2*nt)]
plt.imshow(v1.reshape(grid.getNXs()))
if manual:
    plt.show()
if nvar == 2:
    plt.imshow(v2.reshape(grid.getNXs()))

# %%
if nvar == 2:
    plt.scatter(v1,v2,s=1)

# %%
grid["v1"] = v1
if nvar == 2:
    grid["v2"] = v2


# %%
gl.migrate(grid,dat,"v1")
dat.setName("Migrate","m1")
if nvar == 2 :
    gl.migrate(grid,dat,"v2")
    dat.setName("Migrate","m2")

# %%
if nvar == 1:
    ax = plt.scatter(dat["Simu"], dat["m1"], s=1)
    plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")

if nvar == 2:
    ax = plt.scatter(dat["Simu.1"],dat["m1"],s=1)
    plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")
    if manual:
        plt.show()
    plt.scatter(dat["Simu.2"],dat["m2"],s=1)
    plt.plot(ax.axes.get_xbound(),ax.axes.get_xbound(),c="r")



