# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))
set.seed(32421)

#
# This part is to demonstrate the use of assessors when manipulating matrix
#

print("Testing Matrix")

reset_to_initial_contents <- function(M, MRR, MSG, MSS, MSP){
  MRR$setValues(M$getValues())
  MSG$setValues(M$getValues())
  MSS$setValues(M$getValues())
  MSP$setValues(M$getValues())
  NULL
}

print("Cloning Matrix of integers")
mati  = MatrixInt(2,3)
err   = mati$setValues(c(1, 2, 3, 4, 5, 6))
err   = mati$display()
mati2 = mati$clone()
err   = mati2$display()

print("Cloning Matrix of doubles")
matd  = MatrixRectangular(2,3)
err   = matd$setValues(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0))
err   = matd$display()
matd2 = matd$clone()
err   = matd2$display()

set.seed(32432)
nrow = 7 # For these tests, the matrix MUST be square (ncol = nrow)
ncol = 7
proba = 0.4 # Probability to set values to 0 (making matrix sparse)

# We create a square symmetrical matrix (not necessarily sparse)

MR = MatrixRectangular(nrow, ncol)
for (icol in 1:ncol) {
  for (irow in 1:nrow) {
   value = rnorm(n = 1, mean = 0, sd = 1.0)
   tirage = runif(n = 1, min = 0., max = 1.)
   if (tirage < proba) {
     value = 0.
    }
   err = MR$setValue(irow-1, icol-1, value)
  }
}

# The symmetric matrix is obtained as t(MR) %*% MR . M is symmetric

MRt = MR$transpose() # Using cloneable feature

M = prodMatrix(MRt, MR)

msg = paste0("Matrix M is symmetric? ", M$isSymmetric())
print(msg)
err = M$display()

# Create the different matrix formats (by conversion or extraction)

# To a rectangular matrix
MRR = MatrixRectangular(nrow, ncol)
err = MRR$setValues(M$getValues())
print("Matrix MRR")
err = MRR$display()

# To a square general matrix
MSG = MatrixSquareGeneral(M)
print("Matrix MSG")
err = MSG$display()

# To a square symmetric matrix
MSS = MatrixSquareSymmetric(M)
print("Matrix MSS")
err = MSS$display()

# To a sparse matrix
MSP = toSparse(M)
print("Matrix MSP")
err = MSP$display()

# Creating a Diagonal matrix (from M)
cst = rnorm(n = 1, mean = 0, sd = 1.0)

#
# Adding a constant to the diagonal of a matrix
#
addendum = 1.432

err = mestitle(0,"Adding a constant value to the diagonal of a matrix")
err = reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

err = MRR$addScalarDiag(addendum)
err = MSG$addScalarDiag(addendum)
print(paste0("Are results for MRR and MSG similar: ", MRR$isSame(MSG)))
err = MSS$addScalarDiag(addendum)
print(paste0("Are results for MRR and MSS similar: ", MRR$isSame(MSS)))
err = MSP$addScalarDiag(addendum)
print(paste0("Are results for MRR and MSP similar: ", MRR$isSame(MSP)))

#
# Multiplying the matrix by a constant
#
multiply = 3.2

err = mestitle(0,"Multiplying a Matrix by a constant")
err = reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

err = MRR$prodScalar(multiply)
err = MSG$prodScalar(multiply)
print(paste0("Are results for MRR and MSG similar: ", MRR$isSame(MSG)))
err = MSS$prodScalar(multiply)
print(paste0("Are results for MRR and MSS similar: ", MRR$isSame(MSS)))
err = MSP$prodScalar(multiply)
print(paste0("Are results for MRR and MSP similar: ", MRR$isSame(MSP)))

#
# Adding a constant to a matrix
# Note: This does not make sense for sparse or diagonal matrices
#
err = mestitle(0,"Adding a constant value to the whole matrix")
err = reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

err = MRR$addScalar(addendum)
err = MSG$addScalar(addendum)
print(paste0("Are results for MRR and MSG similar: ", MRR$isSame(MSG)))
err = MSS$addScalar(addendum)
print(paste0("Are results for MRR and MSS similar: ", MRR$isSame(MSS)))

#
# Linear combination
#
err = mestitle(0,"Linear combination of matrices")
err = reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

cx =  1.3
cy = -0.3

err = MRR$linearCombination(cx,cy,MRR)
err = MSG$linearCombination(cx,cy,MSG)
print(paste0("Are results for MRR and MSG similar: ", MRR$isSame(MSG)))
err = MSS$linearCombination(cx,cy,MSS)
print(paste0("Are results for MRR and MSS similar: ", MRR$isSame(MSS)))
err = MSP$linearCombination(cx,cy,MSP)
print(paste0("Are results for MRR and MSP similar: ", MRR$isSame(MSP)))

#
# Extraction of a Vector
# All the tests are not performed on all the matrix types
#
err = mestitle(0,"Extracting Vectors from Matrix")
err = reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

print("MRR and MSP matrices are used as Reference")
err  = MRR$display()
Vref = MRR$getDiagonal()
V1   = MSP$getDiagonal()

print(paste0("Main Diagonal:"))
print(Vref)
print(paste0("Are results for MRR and MSP similar: ", VectorHelper_isSame(Vref, V1)))

Vref = MRR$getDiagonal(1)
V1   = MSP$getDiagonal(1)
print(paste0("Second Diagonal Below:"))
print(Vref)
print(paste0("Are results for MRR and MSP similar: ", VectorHelper_isSame(Vref, V1)))

Vref = MRR$getDiagonal(-2)
V1   = MSP$getDiagonal(-2)
print(paste0("Third Diagonal Above:"))
print(Vref)
print(paste0("Are results for MRR and MSP similar: ", VectorHelper_isSame(Vref, V1)))

Vref = MRR$getRow(2)
V1   = MSP$getRow(2)
print(paste0("Third Row:"))
print(Vref)
print(paste0("Are results for MRR and MSP similar: ",  VectorHelper_isSame(Vref, V1)))

Vref = MRR$getColumn(3)
V1   = MSP$getColumn(3)
print(paste0("Fourth Column:"))
print(Vref)
print(paste0("Are results for MRR and MSP similar: ",  VectorHelper_isSame(Vref, V1)))

#
# Product of the matrix by a vector
#
err = mestitle(0,"Product of the matrix by a vector")
err = reset_to_initial_contents(M, MRR, MSG, MSS, MSP)

# initialisation
print(paste0("nrow = ", nrow))
Vref = VectorDouble(nrow)
V2   = VectorDouble(nrow)
err = MRR$prodVector(V1, Vref)
err = MSG$prodVector(V1, V2)
print(paste0("Are results for MRR and MSG similar: ",  VectorHelper_isSame(Vref, V2)))
err = MSS$prodVector(V1, V2)
print(paste0("Are results for MRR and MSS similar: ",  VectorHelper_isSame(Vref, V2)))
err = MSP$prodVector(V1, V2)
print(paste0("Are results for MRR and MSP similar: ",  VectorHelper_isSame(Vref, V2)))

#
# Linear solver
#

err = mestitle(0,"Matrix Linear Solver")
err = reset_to_initial_contents(M, MRR, MSG, MSS, MSP)
V3 = VectorDouble(rep(0.0, nrow))
print(paste0("Solve X from A*X=B. Compute A*X and compare with B"))

err = MSS$solve(b = V1, x = V2)
err = MSS$prodVector(V2, V3)
print(paste0("Are results correct for MSS: ",  VectorHelper_isSame(V1, V3)))

err = MSP$solve(b = V1, x = V2)
err = MSP$prodVector(V2, V3)
print(paste0("Are results correct for MSP: ", VectorHelper_isSame(V1, V3)))

#
# Inversion
#

err = mestitle(0,"Matrix Inversion")
err = reset_to_initial_contents(M, MRR, MSG, MSS, MSP)
print(paste0("Calculate B=A^{-1}. Compute A*B and compare to Identity"))

MSGref = MSG$clone() # Used to perform A*A-1 and check Identity

err = MSG$invert()
Res = prodMatrix(MSG, MSGref)
print(paste0("Are results correct for MSG: ", Res$isIdentity()))

err = MSS$invert()
Res = prodMatrix(MSS, MSGref)
print(paste0("Are results correct for MSS: ", Res$isIdentity()))

err = MSP$invert()
Res = prodMatrix(MSP, MSGref)
print(paste0("Are results correct for MSP: ", Res$isIdentity()))

print(paste0("Test successfully performed"))
