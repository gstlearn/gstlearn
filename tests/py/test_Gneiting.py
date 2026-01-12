import gstlearn as gl
import numpy as np

# creating the space-time
space1D = gl.SpaceRN.create(1)
space2D = gl.SpaceRN.create(2)

# creating the time trace covariance
ctxt1D = gl.CovContext(1, space1D)
alpha = 1.0
beta = 1.0
params = [alpha, beta * 2 / 2]
scaleT = 5.3
angleT = 0.0
corT = gl.CorAniso.create(ctxt1D, gl.ECov.CAUCHY_GEN, params, scaleT, angleT, False)

# creating the space trace covariance
ctxt2D = gl.CovContext(1, space2D)
nu = [0.5]
kappa = [1.0]
scales = [2.0, 3.0]
angles = [0, 0]
corS = gl.CorGaussianMixture.create(
    gl.ECov.MATERN, ctxt2D, nu, kappa, scales, angles, False
)

sep = 1
gneiting = gl.CorGneiting(corS, corT, sep)
print("Space dimension of Gneiting Covariance = " + str(gneiting.getNDim()))

coords1 = [12, 3, 1]
coords2 = [4, 5, 2]
p1 = gl.SpacePoint(gneiting.getSpace())
p2 = gl.SpacePoint(gneiting.getSpace())
p1.setCoords(coords1)
p2.setCoords(coords2)
cres = gneiting.evalCov(p1, p2)
print("Gneiting eval " + str(cres))

p1_0 = p1.spacePointOnSubspace(0)
print("Displaying the first two coordinates (spatial) of the first space point")
p1_0.display()

p1_1 = p1.spacePointOnSubspace(1)
print("Displaying the last coordinate (temporal) of the first space point")

p1_1.display()

p2_0 = p2.spacePointOnSubspace(0)
print("Displaying the first two coordinates (spatial) of the second space point")
p2_0.display()

p2_1 = p2.spacePointOnSubspace(1)
print("Displaying the last coordinate (temporal) of the second space point")

p2_1.display()

# temporal trace covariance
covt_gstl = corT.evalCov(p1_1, p2_1)
covt_calc = 1 / (
    (1 + (np.abs(coords1[2] - coords2[2]) / scaleT) ** params[0]) ** params[1]
)
print(
    "Difference for temporal covariance = "
    + str(np.round(np.abs(covt_gstl - covt_calc), 6))
)

# spatial trace covariance
delta = [np.abs(coords1[0] - coords2[0]), np.abs(coords1[1] - coords2[1])]
covs_gstl = corS.evalCov(p1_0, p2_0)
covs_calc = np.exp(-np.sqrt((delta[0] / scales[0]) ** 2 + (delta[1] / scales[1]) ** 2))
print(
    "Difference for spatial covariance = "
    + str(np.round(np.abs(covs_gstl - covs_calc), 6))
)

# gneiting covariance
al = sep / 2
scale = covt_gstl**al
covspatCopy = gl.CorGaussianMixture(corS)
covspatCopy.setScaleGneiting(scale)
covs = covspatCopy.evalCov(p1_0, p2_0)

print(
    "Difference for Gneiting covariance = "
    + str(np.round(np.abs(cres - covs * covt_gstl), 6))
)
