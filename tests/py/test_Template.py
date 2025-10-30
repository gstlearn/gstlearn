# Script in python
# It is meant to check the Template exported in Python via Swig

# Next we load the packages used by the subsequent script, such as *gstlearn* (note that any function belonging to this library will have a prefix *gl*. Similarly, we load the graphic extension (called *gslearn.plot*).

import gstlearn as gl

# Establish a sequence of real values
tab = gl.VectorHelper.sequenceVD(0, 9., 1.)
tab[3]           = gl.TEST
gl.VectorHelper.dump("Creating Vector of Double", tab)

# Establish a sequence of integer values
itab = gl.VectorHelper.sequence(10)
itab[3]        = gl.ITEST
gl.VectorHelper.dump("Creating Vector of Integer", itab)

# Use the traditional function cumul() in its real version (polymorphism)
rvalue = gl.VectorHelper.cumul(tab)
print("Cumul Double = ", rvalue)

# Use the traditional function cumul() in its integer version (polymorphism)
ivalue = gl.VectorHelper.cumul(itab)
print("Cumul Integer = ", ivalue)

# Use the template function cumul2() in its real version
#rvalueNew = gl.VectorHelper.cumul2(tab)
#print("Cumul Neutral (tab) = ", rvalueNew)

# Use the template function cumul2() in its integer version
#ivalueNew = gl.VectorHelper.cumul2(itab)
#print("Cumul Neutral (itab) = ", ivalueNew)
