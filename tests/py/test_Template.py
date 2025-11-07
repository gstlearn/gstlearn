# Script in python
# It is meant to check the Template exported in Python via Swig

# Next we load the packages used by the subsequent script, such as *gstlearn* (note that any function belonging to this library will have a prefix *gl*. Similarly, we load the graphic extension (called *gslearn.plot*).

import gstlearn as gl

# - created using VH functions
tab = gl.VectorHelper.sequenceVD(0, 9., 1.)
tab[3]           = gl.TEST
print("Type of VD created using VH function")
print(type(tab))

# - created using the VD constructor only
tab = gl.VectorDouble([2,3,4])
print("Type of VD created using VD constructor")
print(type(tab))
print("Identifying the vector type (double)")
tab.identify()

# - created using VH function
itab = gl.VectorHelper.sequence(10)
itab[3]        = gl.ITEST
print("Type of VD created using VH function")
print(type(itab))

# - created using the VD constructor only
itab = gl.VectorInt({1, 2, 4, 6})
print("Type of VD created using VD constructor")
print(type(itab))
print("Identifying a vector of Id")
itab.identify()

# - manipulation of VD
tab.dump("Initial Vector of Double")
tab.add(tab)
tab.dump("Add this to itself")
tab.addCst(12)
tab.dump("Add a constant to this")
tabaux = tab.addVec(tab)
tabaux.dump("Add 'this' and return a new VD")

# - manipulation of VI
itab.dump("Initial Vector of Id")
itab.add(itab)
itab.dump("Add this to itself")
itab.addCst(12)
itab.dump("Add a constant to this")
itabaux = itab.addVec(itab)
itabaux.dump("Add 'this' and return a new VI")
