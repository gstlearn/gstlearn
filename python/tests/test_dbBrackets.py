import gstlearn as gl
import numpy as np

np.random.seed(124)

a = gl.Db.createFromGrid([2,2],[1.,1.])
a.display()

a["var1"] = np.random.normal(size=4)
a.display()

print(a["var1"])

a["var1"] = np.random.normal(size=4)

print(a["var1"])

a["var2"] = np.random.normal(size =4)
a.display()

print(a["var*"])

a["var*"]=a["var*"]>0

print(a["var*"])

a["newvar"] = np.random.normal(size = (4,3))

a.display()

print(a["newvar*"])

v = a["newvar*"]
v[0,0]=None

a["newvar*"] = v

print(a["newvar*"])

