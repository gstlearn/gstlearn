# %%
import gstlearn as gl
import numpy as np

np.random.seed(123)
dbg = gl.DbGrid.create(nx = [11,11])
ncell = dbg.getSampleNumber()
seltab = np.ones(ncell)
middle = int(ncell / 2)
seltab[middle] = 0
flagSel = True
meshref = gl.MeshETurbo(dbg)
if flagSel:
    dbg.addColumns(seltab, "sel", gl.ELoc.SEL)
    print("Pixel", middle, "has been masked off")

meshsel = gl.MeshETurbo(dbg)
model = gl.Model.createFromParam(gl.ECov.MATERN,ranges = [3.,4.],param=1)
cova = model.getCova(0)

# %%
newvar = np.random.normal(size=dbg.getSampleNumber())
newvar[dbg["x1"]==00] = 0
newvar[dbg["x1"]==10] = 0
newvar[dbg["x2"]==00] = 0
newvar[dbg["x2"]==10] = 0
if flagSel:
    newvar[middle]    = 0

# %%
Smat  = gl.ShiftOpMatrix(meshref,cova)
result = Smat.evalDirect(newvar)
result[dbg["x1"]==00] = 0
result[dbg["x1"]==10] = 0
result[dbg["x2"]==00] = 0
result[dbg["x2"]==10] = 0
if flagSel:
    result[middle]    = 0

# %%
Ssten = gl.ShiftOpStencil(meshsel,cova)
resultNew = Ssten.evalDirect(newvar)

# %%
print("Difference =", np.max(np.abs(resultNew-result)))


# %%
