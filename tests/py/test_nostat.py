import gstlearn as gl
import gstlearn.test as gt
import numpy as np

# Creating a data base with non-stationary parameters
db = gl.DbGrid.create([2,2])
db["scales1"] = np.arange(1,5)
db["scales2"] = 2*np.arange(1,5)
db["angles"]  = (20+10*np.arange(1,5))*1

db.setLocators(["scales*","angles"],gl.ELoc.NOSTAT)

# Creating the non-stationary Model
m = gl.Model.createFromParam(gl.ECov.EXPONENTIAL)
nostat = gl.NoStatArrayCov(["S1","S2","A1"],db)
m.getCova(0).addNoStat(nostat)

# Calculating the data covariance matrix (using gstlearn)
covmat_gstlearn = m.evalCovMatrix(db).toTL()

# Calculating the data covariance matrix (by hand)
x = db["x1"]
y = db["x2"]
covmat_hand = np.empty(shape=(4,4))

def rotmat(theta):
    thp = theta/180 * np.pi
    return np.array([[np.cos(thp),-np.sin(thp)],[np.sin(thp),np.cos(thp)]])

for i in range(4):
    R = rotmat(db["angles"][i])
    Sigmai = R@np.diag([db["scales1"][i]**2,db["scales2"][i]**2])@R.T
    deti = np.linalg.det(Sigmai)**(1/4)
    
    for j in range(4):
            R = rotmat(db["angles"][j])
            Sigmaj = R@np.diag([db["scales1"][j]**2,db["scales2"][j]**2])@R.T
            detj = np.linalg.det(Sigmaj)**(1/4)
            
            SigmaM = (Sigmai + Sigmaj)/2.
            detM = np.linalg.det(SigmaM)**(1/2)
            
            delta = [x[j]-x[i],y[j]-y[i]]
            xt = np.linalg.solve(SigmaM,delta)
            norm = np.sqrt(np.sum(xt * delta))
            ratio = deti * detj / detM
            covmat_hand[i,j] = ratio * np.exp(-norm)

# Comparing the results
gt.checkEqualityVector(covmat_gstlearn, covmat_hand)
