# %% [markdown]
#         

# %%
import gstlearn as gl
import numpy as np
dig = 13
scales = [2.,3.]
scaleT = 5.3
coords1 = [12,3,1]
coords2 = [4,5,2]
space1D = gl.SpaceRN.create(1)
sep = 1
covtemp = gl.Model.createFromParam(gl.ECov.EXPONENTIAL,ranges =scaleT,flagRange = False, space = space1D).getCova(0)
covspat = gl.Model.createFromParam(gl.ECov.EXPONENTIAL,ranges =scales,flagRange = False).getCova(0)
gneiting = gl.CorGneiting(covspat.getCorAniso(),covtemp.getCorAniso(),sep)
p1 = gl.SpacePoint(gneiting.getSpaceSh())
p2 = gl.SpacePoint(gneiting.getSpaceSh())
p1.setCoords(coords1)
p2.setCoords(coords2)
cres = gneiting.eval(p1,p2)
print("Gneiting eval " + str(cres))

# %%
p1_0 = p1.spacePointOnSubspace(0) 
print("Displaying the first two coordinates (spatial) of the first space point")
p1_0.display() #Problem, it should display only the first two coordinates


# %%
p1_1 = p1.spacePointOnSubspace(1)
print("Displaying the last coordinate (temporal) of the first space point")

p1_1.display()

# %%
p2_0 = p2.spacePointOnSubspace(0)
print("Displaying the first two coordinates (spatial) of the second space point")
p2_0.display()  #Problem, it should display only the first two coordinates

# %%
p2_1 = p2.spacePointOnSubspace(1)
print("Displaying the last coordinate (temporal) of the second space point")

p2_1.display()

# %%
ct = covtemp.eval(p1_1,p2_1) 
print("Difference for temporal covariance " + str(ct - np.exp(-np.abs(coords1[2]-coords2[2])/scaleT)))

# %%
al = sep/2

covspatCopy = gl.CovAniso(covspat)
covspatCopy.setScale(0,covspat.getScale(0)/ct**al)
covspatCopy.setScale(1,covspat.getScale(1)/ct**al)

cs = covspatCopy.eval(p1_0,p2_0)

delta = [np.abs(coords1[0]-coords2[0]),np.abs(coords1[1]-coords2[1])]
diff = cs - np.exp(-np.sqrt((ct**al * delta[0]/scales[0])**2 
                         +  (ct**al * delta[1]/scales[1])**2))
print("Difference for spatial covariance (corrected) " 
      + str(np.round(diff,dig)))

# %%
print("Difference for Gneiting cov " +str(cres - cs * ct))



