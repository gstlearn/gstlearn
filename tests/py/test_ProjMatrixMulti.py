# %%
import gstlearn as gl
db1 = gl.Db()
db1["x1"] = [0.1, 0.2]
db1["x2"] = [0.1, 0.2]
db1.setLocators(["x*"],gl.ELoc.X)

db2 = gl.Db()
db2["x1"] = [0.1, 0.2,0.3]
db2["x2"] = [0.1, 0.2,0.3]
db2.setLocators(["x*"],gl.ELoc.X)

mesh1 = gl.MeshETurbo([3,3],[1/2,1/2])
mesh2 = gl.MeshETurbo([4,4],[1/3,1/3])
a11 = gl.ProjMatrix(db1,mesh1)
a12 = gl.ProjMatrix(db1,mesh2)

a21 = gl.ProjMatrix(db2,mesh1)
a22 = gl.ProjMatrix(db2,mesh2)

# %%
print("Tests message when arguments are inconsistent")
print("---------")
u = gl.VectorConstProjMatrix([a11,a22])
a=gl.ProjMultiMatrix.create(u,2)
print("---------")
u = gl.VectorConstProjMatrix([a11,None])
a=gl.ProjMultiMatrix.create(u,2)

# %%
u = gl.VVectorConstProjMatrix([[a11,a12],[None,None],[a11,a12]])
a = gl.ProjMultiMatrix(u)
print("-----------")

# %%
u = gl.VVectorConstProjMatrix([[a11,a12],[None,None,a12],[a11,a12]])
a = gl.ProjMultiMatrix(u)
print("-----------")

# %%
u = gl.VVectorConstProjMatrix([[a11,a12],[a21,a12],[a11,a12]])
a = gl.ProjMultiMatrix(u)
print("-----------")

# %%
u = gl.VVectorConstProjMatrix([[a11,a12],[a21,a22],[a11,a11]])
a = gl.ProjMultiMatrix(u)
print("-----------")

# %%
u = gl.VVectorConstProjMatrix([[a11,a12],[a21,a22],[a21,None]])
a = gl.ProjMultiMatrix(u)

# %%
print("Test polymorphism")
u = gl.VVectorConstProjMatrix([[a11,a12],[a21,a22],[a21,None]])
a = gl.ProjMultiMatrix(u)
print("-------------")
u = gl.VVectorConstIProjMatrix([[a11,a12],[a21,a22],[a21,None]])
a = gl.ProjMulti(u)

# %%
print("Test Dimensions")
napex = a.getApexNumber()
napexmanual = a11.getApexNumber()+a12.getApexNumber()
print(napex - napexmanual)
npoint = a.getPointNumber()
npointmanual = a11.getPointNumber()+2 * a21.getPointNumber()
print(npoint - npointmanual)

# %%
import numpy as np
aref = [[a11,a12],[a21,a22],[a21,None]]
uI = gl.VVectorConstIProjMatrix(aref)
uM = gl.VVectorConstProjMatrix(aref)
aI = gl.ProjMulti(uI)
aM = gl.ProjMultiMatrix(uM)
np.random.seed(134)
gaussA = np.random.normal(size = napex)
vect = gl.VectorEigen(gaussA)
result = gl.VectorEigen(npoint)
aI.mesh2point(vect,result)
res = np.array([result.getValue(i) for i in range(result.size())])

# %%
resnumpyM2P = []
for i in range(len(u)):
    temp = 0
    ind = 0
    for j in range(len(u[i])):
        mat = uM[i][j]
        if mat != None:
            mat = mat.toTL()
            temp += mat @ gaussA[ind:ind+mat.shape[1]]
            ind+=mat.shape[1]
    resnumpyM2P += [temp]
resnumpyM2P = np.concatenate(resnumpyM2P)

# %%
print("------------")
print("Test mesh2Point (matrix free)")
print(np.sum(np.abs(resnumpyM2P-res)))

# %%
np.random.seed(134)
gaussP = np.random.normal(size = npoint)
vect = gl.VectorEigen(gaussP)
result = gl.VectorEigen(napex)
aI.point2mesh(vect,result)
res = np.array([result.getValue(i) for i in range(result.size())])

# %%
resnumpyP2M = []
for i in range(len(u[0])):
    temp = 0
    ind = 0
    for j in range(len(u)):
        mat = uM[j][i]
        if mat != None:
            mat = mat.toTL()
            temp += mat.T @ gaussP[ind:ind+mat.shape[0]]
            ind+=mat.shape[0]
    resnumpyP2M += [temp]
resnumpyP2M = np.concatenate(resnumpyP2M)

# %%
print("------------")
print("Test point2mesh (matrix-free)")
print(np.sum(np.abs(resnumpyP2M-res)))

# %%
print("------------")
print("With sparse matrices")


# %%
print("mesh2point from the matrix in python (with toTL())")
np.round(np.sum(np.abs(aM.toTL() @ gaussA - resnumpyM2P)),15)

# %%
print("point2mesh from the matrix in python (with toTL())")
np.sum(np.abs(aM.toTL().T @ gaussP - resnumpyP2M))

# %%
vect = gl.VectorEigen(gaussA)
result = gl.VectorEigen(npoint)
aM.mesh2point(vect,result)
res = np.array([result.getValue(i) for i in range(result.size())])
print("mesh2point from the matrix computed in C++")
print(np.round(np.sum(np.abs(res - resnumpyM2P)),15))



# %%
vect = gl.VectorEigen(gaussP)
result = gl.VectorEigen(napex)
aM.point2mesh(vect,result)
res = np.array([result.getValue(i) for i in range(result.size())])
print("point2mesh from the matrix computed in C++")
print(np.sum(np.abs(res - resnumpyP2M)))


# %%



