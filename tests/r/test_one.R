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

M = MatrixFactory_prodMatMat(MRt, MR)

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
MSP = createFromAnyMatrix(M)
print("Matrix MSP")
err = MSP$display()

