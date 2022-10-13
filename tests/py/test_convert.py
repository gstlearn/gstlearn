#
# This file is meant to test the writing of the contents of a Data Base
# into a specific file (such as ArcGis)
#

import gstlearn as gl

# Create a grid: 60 by 40 meshes, with a square mesh of size 10

grid = gl.DbGrid.create([60,40],[10,10])

# Write the contents of the file into an external file using ArcGis format

aof = gl.GridArcGis("ArcGis.dat",grid)
aof.setCol(0)
if aof.isAuthorized():
    aof.writeInFile()

print("Test successfully performed")
