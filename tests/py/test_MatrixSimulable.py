# %%
import gstlearn as gl
import numpy as np

ndig = 13
np.random.seed(12334)

# %%
A = gl.MatrixSparse(4,4)
A.fillRandom()
A.prodMatInPlace(A.transpose())

Apython = A.toTL().todense()
AChol = np.linalg.cholesky(Apython)
x = np.random.normal(size=4)
simuinv = np.linalg.solve(AChol.T,x)
simudir = AChol@x
logdet = np.log(np.linalg.det(Apython))
temp = np.ravel(A.toTL().todense())

B = gl.MatrixSquareGeneral(4)
B.setValues(temp)
C = gl.MatrixRectangular(4,4)
C.setValues(temp)
D = gl.MatrixSquareSymmetric(C)
E = gl.MatrixRectangular(4,5)
E.fillRandom()

def test(A,inverse,name= ""):
    print("------------------")
    print("-- Test --" + name + " with inverse " + str(inverse))
    ua = gl.MatrixSquareSymmetricSim(A,inverse)
    if ua.isEmpty():
        return
    print("Sparse " + str(ua.isSparse()))
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

inverse = False
test(A,inverse,"A")
test(B,inverse,"B")
test(C,inverse,"C")
test(D,inverse,"D")
test(E,inverse,"E")


