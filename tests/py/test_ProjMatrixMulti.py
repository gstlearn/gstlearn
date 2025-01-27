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
print("---------------------------------------------------------------------")
print("---------------------------------------------------------------------")
print("A - Tests message when arguments are inconsistent:")
print("---------------------------------------------------------------------")
print("I. Creation of the vector<vector<Projmatrix*>> from a vector<ProjMatrix*>")
print("---------------------------------------------------------------------")
print("Test1: with a nullptr element:")
print("--------")
u = gl.VectorConstProjMatrix([a11,None])
a=gl.ProjMultiMatrix.create(u,2)
print("---------------------")
print("Test2: with different Point Numbers:")
print("--------")
u = gl.VectorConstProjMatrix([a11,a22])
a=gl.ProjMultiMatrix.create(u,2)

# %%
print("---------------------------------------------------------------------")
print("II. Creation of ProjMulti")
print("---------------------------------------------------------------------")
print("Test3 : full line of nullptr:")
print("--------")
u = gl.VVectorConstProjMatrix([[a11,a12],[None,None],[a11,a12]])
a = gl.ProjMultiMatrix(u)
# %%
print("---------------------")
print("Test4: Rows with different number of elements:")
print("--------")
u = gl.VVectorConstProjMatrix([[a11,a12],[None,None,a12],[a11,a12]])
a = gl.ProjMultiMatrix(u)

# %%
print("---------------------")
print("Test5: Point Numbers different on a row:")
print("--------")
u = gl.VVectorConstProjMatrix([[a11,a12],[a21,a12],[a11,a12]])
a = gl.ProjMultiMatrix(u)
# %%
print("---------------------")
print("Test6: Apex Numbers different on a column:")
print("--------")
u = gl.VVectorConstProjMatrix([[a11,a12],[a21,a22],[a11,a11]])
a = gl.ProjMultiMatrix(u)
# %%
print("---------------------")
print("Test7: empty matrix:")
print("--------")
u = gl.VVectorConstProjMatrix([])
a = gl.ProjMultiMatrix(u)
print("---------------------------------------------------------------------")
# %%
print("End of inconsistency tests")
print("---------------------------------------------------------------------")

# %%
print("B - Various tests")
print("---------------------------------------------------------------------")
print("Test 8: polymorphism:")
print("---------------------")
u = gl.VVectorConstProjMatrix([[a11,a12],[a21,a22],[a21,None]])
a = gl.ProjMultiMatrix(u)
print("Test case 1 passed")
u = gl.VVectorConstIProj([[a11,a12],[a21,a22],[a21,None]])
a = gl.ProjMulti(u)
print("Test case 2 passed")
print("---------------------")# %%
print("Test 9: Dimensions")
print("---------------------")
u = gl.VVectorConstIProj([[a11,a12],[a21,a22],[a21,None]])
a = gl.ProjMulti(u)
print("Number of differences for nlatent " + str(len(u[0])-a.getNLatent()))
print("Number of differences for nlatent " + str(u.size()-a.getNVariable()))
print(a.getNVariable())
print(a.getNLatent())
napex = a.getNApex()
napexmanual = a11.getNApex()+a12.getNApex()
print("Number of differences for apex number " + str(napex - napexmanual))
npoint = a.getNPoint()
npointmanual = a11.getNPoint()+2 * a21.getNPoint()
print("Number of differences for point number " + str(npoint - npointmanual))

# %%
print("---------------------------------------------------------------------")
print("C - Test Operators")
print("---------------------------------------------------------------------")

import numpy as np
aref = [[a11,a12],[a21,a22],[a21,None]]
uI = gl.VVectorConstIProj(aref)
uM = gl.VVectorConstProjMatrix(aref)
aI = gl.ProjMulti(uI)
aM = gl.ProjMultiMatrix(uM)
np.random.seed(134)
gaussA = np.random.normal(size = napex)
vect = gl.VectorDouble(gaussA)
result = gl.VectorDouble(npoint)
aI.mesh2point(vect,result)
res = np.array([result[i] for i in range(result.size())])

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
print("---------------------")
print("Test 10: mesh2Point (matrix free):")
print("---------------------")
print("Difference with manual computation " + str(np.sum(np.abs(resnumpyM2P-res))))

# %%
np.random.seed(134)
gaussP = np.random.normal(size = npoint)
vect = gl.VectorDouble(gaussP)
result = gl.VectorDouble(napex)
aI.point2mesh(vect,result)
res = np.array([result[i] for i in range(result.size())])

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
print("---------------------")
print("Test 11: point2mesh (matrix-free):")
print("---------------------")
print("Difference with manual computation " + str(np.sum(np.abs(resnumpyP2M-res))))

# %%
print("---------------------")
print("With sparse matrices:")
print("---------------------")

# %%
print("Test 12: mesh2point from the matrix in python (with toTL()):")
print("---------------------")
print("Difference with manual computation " + \
      str(np.round(np.sum(np.abs(aM.getProj().toTL() @ gaussA - resnumpyM2P)),15)))


# %%
print("---------------------")
print("Test 13: point2mesh from the matrix in python (with toTL()):")
print("---------------------")
print("Difference with manual computation " + \
      str(np.round(np.sum(np.abs(aM.getProj().toTL().T @ gaussP - resnumpyP2M)),15)))

# %%
vect = gl.VectorDouble(gaussA)
result = gl.VectorDouble(npoint)
aM.mesh2point(vect,result)
res = np.array([result[i] for i in range(result.size())])
print("---------------------")
print("Test 14: mesh2point from the matrix computed in C++:")
print("---------------------")
print("Difference with manual computation " + \
      str(np.round(np.sum(np.abs(res - resnumpyM2P)),15)))



# %%
vect = gl.VectorDouble(gaussP)
result = gl.VectorDouble(napex)
aM.point2mesh(vect,result)
res = np.array([result[i] for i in range(result.size())])
print("---------------------")
print("Test 15: point2mesh from the matrix computed in C++:")
print("---------------------")
print("Difference with manual computation " + \
      str(np.round(np.sum(np.abs(res - resnumpyP2M)),15)))
print("---------------------")
print("End of tests")
# %%



