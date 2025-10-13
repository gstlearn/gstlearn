# %%
import gstlearn as gl
import numpy as np

ndig = 13
np.random.seed(12334)

# %%
A = gl.MatrixSparse(4,4)
A.fillRandom()
A.prodMat(A.transpose())

Apython = A.toTL().todense()
AChol = np.linalg.cholesky(Apython)
x = np.random.normal(size=4)
simuinv = np.linalg.solve(AChol.T,x)
simudir = AChol@x
logdet = np.log(np.linalg.det(Apython))
temp = np.ravel(A.toTL().todense())
random = np.random.normal(size = temp.shape)
B = gl.MatrixSquare(4)
B.setValues(temp)
C = gl.MatrixDense(4,4)
C.setValues(temp)
D = gl.MatrixSymmetric(C)
E = gl.MatrixDense(4,5)
E.fillRandom()
F = gl.MatrixSquare(4)
F.fillRandom()
M = gl.MatrixSymmetric(F)
# %%
def test(mat,inverse,name= ""):
    print("------------------")
    print("-- Test --" + name + " with inverse " + str(inverse))
    ua = gl.MatrixSymmetricSim(mat,inverse)
    if ua.isEmpty():
        return
    a = ua.evalSimulate(x)
    if inverse:
        y = simuinv
    else :
        y = simudir
    print(np.round(np.max(np.abs(a-y)),ndig))

# %%
inverse = True
test(A,inverse,"A")
test(B,inverse,"B")
test(C,inverse,"C")
test(D,inverse,"D")
test(E,inverse,"E")
test(F,inverse,"F")
inverse = False
test(A,inverse,"A")
test(B,inverse,"B")
test(C,inverse,"C")
test(D,inverse,"D")
test(E,inverse,"E")
test(F,inverse,"F")