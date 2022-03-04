import gstlearn as gl
import os
import sys

# Redirection

filename = os.path.splitext(os.path.basename(__file__))[0] + '.out'
sys.stdout = open(filename,'w')

# Enums are global unique object instances
d1 = gl.ENeigh.UNIQUE
print ("d1=", "Enum#", d1.getValue(), ":", d1.getKey(), "|", d1.getDescr())
# Test defaut value
d2 = gl.ENeigh()
print ("d2=", "Enum#", d2.getValue(), ":", d2.getKey(), "|", d2.getDescr())
# Create from an integer
d3 = gl.ENeigh.fromValue(2)
print ("d3=", "Enum#", d3.getValue(), ":", d3.getKey(), "|", d3.getDescr())

# Use the iterator (getIterator is a static function, the returned iterator is static also)
it = gl.ENeigh.getIterator()
while it.hasNext() :
    print ("Enum#", it.getValue(), ": ", it.getKey(), "|", it.getDescr())
    it.toNext()

# Test existence (static function)
nei = "SUPERMOVING"
if gl.ENeigh.existsKey(nei) :
    print(nei, "exists")
else :
    print(nei, "doesn't exists")

# No switch in python
if d1 == gl.ENeigh.UNIQUE or d1 == gl.ENeigh.MOVING :
    print ("d1=", d1.getDescr(), "is often used!")
else:
    print ("d1=", d1.getDescr(), "is rarely used!")

print("Test successfully performed")
