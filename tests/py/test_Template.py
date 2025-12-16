# Script in python
# It is meant to check the Template exported in Python via Swig

import gstlearn as gl
import numpy as np

# - created using the VD constructor only
rtab1 = gl.VectorDouble([2.1, 3.2, 4.1, np.nan, 8.3])
rtab1.dump("rtab1")
rtab1.identify()
print("Moyenne = ", np.round(rtab1.mean(), 4))
print("Mediane = ", rtab1.median())
print("Minimum = ", rtab1.minimum())
print("Maximum = ", rtab1.maximum())
print("Variance = ", np.round(rtab1.variance(), 4))
print("St. dev. = ", np.round(rtab1.stdv(), 4))

rtab2 = gl.VectorDouble([-1.1, np.nan, 6.1, 4.1, -0.2])
rtab2.dump("rtab2")
rtab1.add(rtab2)
rtab1.dump("'rtab1 = 'rtab1 + 'rtab2'")
rtab1.addCst(12.0)
rtab1.dump("'rtab1' = 'rtab1' + 12")
tabaux = rtab1.addVec(rtab2)
tabaux.dump("'out' = 'rtab1' + 'rtab2'")

# - created using the VD constructor only
itab1 = gl.VectorInt({1, 2, 4, 6, 3, 5})
itab1.dump("itab1")
itab1.identify()
print("Moyenne = ", itab1.mean())
print("Mediane = ", itab1.median())
print("Variance = ", np.round(itab1.variance(), 4))
print("St. dev. = ", np.round(itab1.stdv(), 4))

itab2 = gl.VectorInt({-5, 4.0, 2, 0, 1, 8})
itab2.dump("itab2")
itab1.add(itab2)
itab1.dump("'itab1 = 'itab1 + 'itab2'")
itab1.addCst(12)
itab1.dump("'itab1' = 'itab1' + 12")
itabaux = itab1.addVec(itab2)
itabaux.dump("'out' = 'itab1' + 'itab2'")
