# This test is meant to test the function extractDiag() applied to a PrecisionOp
# whether the matrix is stored explicitely (PrecisionOpMatrix) or not (PrecisionOp)

import gstlearn as gl
import time

gl.OptCst.define(gl.ECst.NTCOL, -1)
gl.OptCst.define(gl.ECst.NTROW, -1)

# The first part is used to check the method extractDiag() on a small Data Set
# depending on whether we use the Matrix-Free PrecisionOp or not
# Due to the small size of the data set, we can affort displaying the results

gl.mestitle(0,"Small Data Set")
nx = 3
print("(Based on the Turbo Mahing derived from a", nx ,"x", nx, "grid)")

# Creating the environment
mesh = gl.MeshETurbo([nx, nx],[],[],[],False)
model = gl.Model.createFromParam(gl.ECov.MATERN, param=1, range=nx/2)
print("Number of apices = ", mesh.getNApices())

# Creating the PrecisionOp with Sparse matrix implemented
QOpCs = gl.PrecisionOpMatrix(mesh, model.getCova(0))

# Printing the complete Precision Matrix to visualize the Diagonal
Qmat = QOpCs.getQ()
#Qmat.display()

ref1 = QOpCs.extractDiag()
print("Using the version with explicit Matrix")
print(ref1)

# Creating the PrecisionOp without Sparse matrix implemented
Qop = gl.PrecisionOp(mesh, model.getCova(0))
ref2 = Qop.extractDiag()
print("Using the Matrix-free version")
print(ref2)

# Comparing the two resulting vectors
if gl.VH.isEqualExtended(ref1, ref2, string="Vector are not equal"):
    print("--> Vectors are identical")

# The second part is used to check the method extractDiag() on a large Data Set
# depending on whether we use the Matrix-Free PrecisionOp or not.
# It is meant to evaluate the elapsed time

print_Time = False
if not print_Time:
    print("\nNote: Time is not printed in order not to generate conflicts in non-regression files)")

gl.mestitle(0,"Large Data Set")
nx = 350
print("(Based on the Turbo Meshing derived from a", nx ,"x", nx, "grid)")

# Creating the environment
mesh = gl.MeshETurbo([nx, nx],[],[],[],False)
model = gl.Model.createFromParam(gl.ECov.MATERN, param=1, range=nx/2)
print("Number of apices = ", mesh.getNApices())

# Creating the PrecisionOpMatrix with Sparse matrix implemented
start_time = time.time()
QOpCs = gl.PrecisionOpMatrix(mesh, model.getCova(0))
ref1 = QOpCs.extractDiag()
checkpointOpCs = str((time.time() - start_time))
if print_Time:
    print("Time with explicit Matrix =", checkpointOpCs)

# Creating the Matrix-free PrecisionOp
start_time = time.time()
Qop = gl.PrecisionOp(mesh, model.getCova(0))
ref2 = Qop.extractDiag()
checkpointOp = str((time.time() - start_time))
if print_Time:
    print("Time Matrix-free version =",checkpointOp)

# Comparing the two resulting vectors
if gl.VH.isEqualExtended(ref1, ref2, string="Vector are not equal"):
    print("--> Vectors are identical")

# The third part is used to check the method extractDiag() on a 3-D Data Set
# We do not compare with the Matrix solution which rapidly becomes untractable
# It is meant to evaluate the elapsed time

gl.mestitle(0,"3-D Data Set")
nx = 50
print("(Based on the Turbo Meshing derived from a", nx ,"x", nx, "x", nx, "grid)")

# Creating the environment
gl.defineDefaultSpace(gl.ESpaceType.RN, 3)
mesh = gl.MeshETurbo([nx, nx, nx],[],[],[],False)
model = gl.Model.createFromParam(gl.ECov.MATERN, param=0.5, range=nx/2)
print("Number of apices = ", mesh.getNApices())

# Creating the Matrix-free PrecisionOp
start_time = time.time()
Qop = gl.PrecisionOp(mesh, model.getCova(0))
ref = Qop.extractDiag()
checkpointOp3D = str((time.time() - start_time))
if print_Time:
    print("Time Matrix-free version =",checkpointOp3D)

