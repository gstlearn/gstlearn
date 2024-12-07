# %%
import gstlearn as gl
import numpy as np

np.random.seed(123)
dbg = gl.DbGrid.create(nx = [11,11])
mesh = gl.MeshETurbo(dbg)
model = gl.Model.createFromParam(gl.ECov.MATERN,ranges = [3.,4.],param=1)
cova = model.getCova(0)
Smat  = gl.ShiftOpMatrix(mesh,cova)    
Ssten = gl.ShiftOpStencil(mesh,cova)

# %%
newvar = np.random.normal(size=dbg.getSampleNumber())
newvar[dbg["x1"]==0] = 0
newvar[dbg["x1"]==10] = 0
newvar[dbg["x2"]==0] = 0
newvar[dbg["x2"]==10] = 0

result = Smat.evalDirect(newvar)
#On doit patcher les valeurs aux bords sur le résultats
result[dbg["x1"]==0] = 0
result[dbg["x1"]==10] = 0
result[dbg["x2"]==0] = 0
result[dbg["x2"]==10] = 0

# %%
resultNew = Ssten.evalDirect(newvar)

# %%
print("Difference =", np.max(np.abs(resultNew-result)))

