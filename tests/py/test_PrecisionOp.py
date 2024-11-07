# This test is meant to test the function "extractDiag()" applied to a PrecisionOp
# whether the matrix is stored explicitely (PrecisionOpCs) or not (PrecisionOp)

import gstlearn as gl
gl.OptCst.define(gl.ECst.NTCOL, -1)
gl.OptCst.define(gl.ECst.NTROW, -1)

# Creating the environment
nx = 3
mesh = gl.MeshETurbo([3,3])
model = gl.Model.createFromParam(gl.ECov.MATERN, param=1, range=2)

# Creating the PrecisionOp with Sparse matrix implemented
QOpCs = gl.PrecisionOpCs(mesh, model.getCova(0))

# Extracting the ShiftOp Matrix 'S'
S = QOpCs.getShiftOp().getS()
print("S Matrix")
print(S)

# Extracting the Precision Matrix 'Q'
Qmat = QOpCs.getQ()
print("Precision Matrix")
print(Qmat)

# Extracting the Diagonal from Q with SCIPY
ref1 = Qmat.toTL().diagonal()
print("Diagonal extracted from Q using Scipy")
print(ref1)

# Extracting the Diagonal from Q with gstlearn
ref2 = Qmat.extractDiag()
print("Diagonal extracted from Q using gstlearn")
print(ref2)

# Extracting the Diagonal with the direct function in gstlearn
ref3 = QOpCs.extractDiag()
print("Diagonal extract from PrecisionOpCs using gstlearn")
print(ref3)

# Creating the PrecisionOp without Sparse matrix implemented
Qop = gl.PrecisionOp(mesh, model.getCova(0))

# Extracting the Diagonal with the direct function in gstlearn
ref4 = Qop.extractDiag()
print("Diagonal extracted from PrecisionOp using gstlearn")
print(ref4)
