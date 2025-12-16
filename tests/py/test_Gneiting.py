import gstlearn as gl
import numpy as np

dig = 13
scales = [2.0, 3.0]
scaleT = 5.3
coords1 = [12, 3, 1]
coords2 = [4, 5, 2]
space1D = gl.SpaceRN.create(1)
sep = 1
covtemp = gl.Model.createFromParam(
    gl.ECov.EXPONENTIAL, ranges=scaleT, flagRange=False, space=space1D
)
ct = covtemp.getCovAniso(0)
ct0 = ct.getCorAniso()
covspat = gl.Model.createFromParam(gl.ECov.EXPONENTIAL, ranges=scales, flagRange=False)
cs = covspat.getCovAniso(0)
cs0 = cs.getCorAniso()
gneiting = gl.CorGneiting(cs0, ct0, sep)
p1 = gl.SpacePoint(gneiting.getSpace())
p2 = gl.SpacePoint(gneiting.getSpace())
p1.setCoords(coords1)
p2.setCoords(coords2)
cres = gneiting.evalCov(p1, p2)
print("Gneiting eval " + str(cres))

p1_0 = p1.spacePointOnSubspace(0)
print("Displaying the first two coordinates (spatial) of the first space point")
p1_0.display()  # Problem, it should display only the first two coordinates


p1_1 = p1.spacePointOnSubspace(1)
print("Displaying the last coordinate (temporal) of the first space point")

p1_1.display()

p2_0 = p2.spacePointOnSubspace(0)
print("Displaying the first two coordinates (spatial) of the second space point")
p2_0.display()  # Problem, it should display only the first two coordinates

p2_1 = p2.spacePointOnSubspace(1)
print("Displaying the last coordinate (temporal) of the second space point")

p2_1.display()

covt = covtemp.evalCov(p1_1, p2_1)
print(
    "Difference for temporal covariance "
    + str(covt - np.exp(-np.abs(coords1[2] - coords2[2]) / scaleT))
)

al = sep / 2

covspatCopy = gl.CovAniso(cs)
covspatCopy.setScaleDim(0, cs.getScale(0) / covt**al)
covspatCopy.setScaleDim(1, cs.getScale(1) / covt**al)

covs = covspatCopy.evalCov(p1_0, p2_0)


delta = [np.abs(coords1[0] - coords2[0]), np.abs(coords1[1] - coords2[1])]
diff = covs - np.exp(
    -np.sqrt(
        (covt**al * delta[0] / scales[0]) ** 2 + (covt**al * delta[1] / scales[1]) ** 2
    )
)
print("Difference for spatial covariance (corrected) " + str(np.round(diff, dig)))

print("Difference for Gneiting cov " + str(cres - covs * covt))
