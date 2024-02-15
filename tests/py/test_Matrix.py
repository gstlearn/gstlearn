import numpy as np
import sys
import os
import gstlearn as gl


def reset_to_initial_contents(M, MRR, MSG, MSS, MSP):
  MRR.setValues(M.getValues())
  MSG.setValues(M.getValues())
  MSS.setValues(M.getValues())
  MSP.setValues(M.getValues())

print("Cloning Matrix of integers")
mati = gl.MatrixInt(2,3)
mati.setValues([1, 2, 3, 4, 5, 6])
mati.display()
mati2 = mati.clone()
mati2.display()

print("Cloning Matrix of doubles")
matd = gl.MatrixRectangular(2,3)
matd.setValues([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
matd.display()
matd2 = matd.clone()
matd2.display()

gl.law_set_random_seed(32432)
nrow = 7 # For these tests, the matrix MUST be square (ncol = nrow)
ncol = 7
proba = 0.4 # Probability to set values to 0 (making matrix sparse)

# We create a square symmetrical matrix (not necessarily sparse)

MR = gl.MatrixRectangular(nrow, ncol)
for icol in range(ncol):
  for irow in range(nrow):
    value = gl.law_gaussian()
    tirage = gl.law_uniform(0., 1.)
    if tirage < proba:
      value = 0.
    MR.setValue(irow, icol, value)

# The symmetric matrix is obtained as t(MR) %*% MR . M is symmetric

MRt = MR.transpose() # Using cloneable feature

M = gl.MatrixFactory.prodMatMat(MRt, MR)

print("Matrix M")
M.display()

# Create the different matrix formats (by conversion or extraction)

# To a rectangular matrix
MRR = gl.MatrixRectangular(nrow,ncol)
MRR.setValues(M.getValues())
print("Matrix MRR")
MRR.display()

# To a square general matrix
MSG = gl.MatrixSquareGeneral(M)
print("Matrix MSG")
MSG.display()

# To a square symmetric matrix
MSS = gl.MatrixSquareSymmetric(M)
print("Matrix MSS")
MSS.display()

# To a sparse matrix
MSP = gl.createFromAnyMatrix(M)
print("Matrix MSP")
MSP.display()

# Creating a Diagonal matrix (from M)
cst = gl.law_gaussian()

#
# Adding a constant to the diagonal of a matrix
#
addendum = 1.432

gl.mestitle(0,"Adding a constant value to the diagonal of a matrix")
reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

MRR.addScalarDiag(addendum)
MSG.addScalarDiag(addendum)
print("Are results for MRR and MSG similar: ",MRR.isSame(MSG))
MSS.addScalarDiag(addendum)
print("Are results for MRR and MSS similar: ",MRR.isSame(MSS))
MSP.addScalarDiag(addendum)
print("Are results for MRR and MSP similar: ",MRR.isSame(MSP))

#
# Multiplying the matrix by a constant
#
multiply = 3.2

gl.mestitle(0,"Multiplying a Matrix by a constant")
reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

MRR.prodScalar(multiply)
MSG.prodScalar(multiply)
print("Are results for MRR and MSG similar: ",MRR.isSame(MSG))
MSS.prodScalar(multiply)
print("Are results for MRR and MSS similar: ",MRR.isSame(MSS))
MSP.prodScalar(multiply)
print("Are results for MRR and MSP similar: ",MRR.isSame(MSP))

#
# Adding a constant to a matrix
# Note: This does not make sense for sparse or diagonal matrices
#

gl.mestitle(0,"Adding a constant value to the whole matrix")
reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

MRR.addScalar(addendum)
MSG.addScalar(addendum)
print("Are results for MRR and MSG similar: ",MRR.isSame(MSG))
MSS.addScalar(addendum)
print("Are results for MRR and MSS similar: ",MRR.isSame(MSS))

#
 # Linear combination
 #
gl.mestitle(0,"Linear combination of matrices")
reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

cx =  1.3
cy = -0.3

MRR.addMatInPlace(MRR,cx,cy)
MSG.addMatInPlace(MSG,cx,cy)
print("Are results for MRR and MSG similar: ",MRR.isSame(MSG))
MSS.addMatInPlace(MSS,cx,cy)
print("Are results for MRR and MSS similar: ",MRR.isSame(MSS))
MSP.addMatInPlace(MSP,cx,cy)
print("Are results for MRR and MSP similar: ",MRR.isSame(MSP))

#
# Extraction of a Vector
# All the tests are not performed on all the matrix types
#
gl.mestitle(0,"Extracting Vectors from Matrix")
reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

print("MRR and MSP matrices are used as Reference")
MRR.display()
Vref = MRR.getDiagonal()
V1 = MSP.getDiagonal()
gl.print_vector("Main Diagonal",0,len(Vref),Vref)
print("Are results for MRR and MSP similar: ",gl.VH.isSame(Vref,V1))
Vref = MRR.getDiagonal(1)
V1 = MSP.getDiagonal(1)
gl.print_vector("Second Diagonal Below",0,len(Vref),Vref)
print("Are results for MRR and MSP similar: ",gl.VH.isSame(Vref,V1))
Vref = MRR.getDiagonal(-2)
V1 = MSP.getDiagonal(-2)
gl.print_vector("Third Diagonal Above",0,len(Vref),Vref)
print("Are results for MRR and MSP similar: ",gl.VH.isSame(Vref,V1))
Vref = MRR.getRow(2)
V1 = MSP.getRow(2)
gl.print_vector("Third Row",0,len(Vref),Vref)
print("Are results for MRR and MSP similar: ",gl.VH.isSame(Vref,V1))
Vref = MRR.getColumn(3)
V1 = MSP.getColumn(3)
gl.print_vector("Fourth Column",0,len(Vref),Vref)
print("Are results for MRR and MSP similar: ",gl.VH.isSame(Vref,V1))

#
# Product of the matrix by a vector
#

gl.mestitle(0,"Product of the matrix by a vector")
reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

Vref = gl.VectorDouble(np.zeros(nrow))
V2 = gl.VectorDouble(np.zeros(nrow))
MRR.prodMatVecInPlace(V1, Vref)
MSG.prodMatVecInPlace(V1, V2)
print("Are results for MRR and MSG similar: ",gl.VH.isSame(np.array(Vref.getVector()),np.array(V2.getVector())))
MSS.prodMatVecInPlace(V1, V2)
print("Are results for MRR and MSS similar: ",gl.VH.isSame(np.array(Vref.getVector()),np.array(V2.getVector())))
MSP.prodMatVecInPlace(V1, V2)
print("Are results for MRR and MSP similar: ",gl.VH.isSame(np.array(Vref.getVector()),np.array(V2.getVector())))

#
# Linear solver
#

gl.mestitle(0,"Matrix Linear Solver")
reset_to_initial_contents(M, MRR, MSG, MSS, MSP)
V3 = gl.VectorDouble(np.zeros(nrow))
print("Solve X from A*X=B. Compute A*X and compare with B")

MSS.solve(V1, V2)
MSS.prodMatVecInPlace(np.array(V2.getVector()), V3)
print("Are results correct for MSS: ",gl.VH.isSame(V1,np.array(V3.getVector())))
MSP.solve(V1, V2)
MSP.prodMatVecInPlace(np.array(V2.getVector()), V3)
print("Are results correct for MSP: ",gl.VH.isSame(V1,np.array(V3.getVector())))

#
# Inversion
#

gl.mestitle(0,"Matrix Inversion")
reset_to_initial_contents(M, MRR, MSG, MSS, MSP)
print("Calculate B=A^{-1}. Compute A*B and compare to Identity")

MSGref = MSG.clone() # Used to perform A*A-1 and check Identity

# TODO : This doesn't work - no more identity !!!
MSG.invert()
Res = gl.MatrixFactory.prodMatMat(MSG, MSGref)
#print(Res)
print("Are results correct for MSG: ",Res.isIdentity())

MSS.invert()
Res = gl.MatrixFactory.prodMatMat(MSS, MSGref)
#print(Res)
print("Are results correct for MSS: ",Res.isIdentity())

MSP.invert()
Res = gl.MatrixFactory.prodMatMat(MSP, MSGref)
#print(Res)
print("Are results correct for MSP: ",Res.isIdentity())

print("Test successfully performed")
