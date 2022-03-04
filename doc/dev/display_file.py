# This script is meant to visualize the contents of any serialized object 
# from the current Directory

import sys
import os
import gstlearn as gl
import gstlearn.plot as gp
import matplotlib.pyplot as plt
from attr._make import NOTHING

args = sys.argv
if len(args) < 2:
    print("This script should be called according to the following syntax:")
    print("  display_file.py filename")
    print("- filename: Name of the Serialized file")
    exit()
filename = args[1]

# Get the Type of the File and stop if returned empty (file not found)
filetype = gl.ASerializable.getFileIdentity(filename)
if filetype == "":
    exit()

if filetype == "Db":
    db = gl.Db.createFromNF(filename,False)
    dbfmt = gl.DbStringFormat()
    dbfmt.setParams(gl.FLAG_VARS)
    db.display(dbfmt)
           
elif filetype == "DbGrid":
    dbgrid = gl.DbGrid.createFromNF(filename,False)
    dbfmt = gl.DbStringFormat()
    dbfmt.setParams(gl.FLAG_VARS)
    dbgrid.display(dbfmt)
            
elif filetype == "Vario":
    vario = gl.Vario.createFromNF(filename,False)
    vario.display()
    
elif filetype == "Model":
    model = gl.Model.createFromNF(filename,False)
    model.display()
    
elif filetype == "Rule":
    rule = gl.Rule.createFromNF(filename,False)
    rule.display()
    
elif filetype == "Table":
    table = gl.Table.createFromNF(filename,False)
    table.display()

elif filetype == "Polygon":
    poly = gl.Polygons.createFromNF(filename,False)
    poly.display()
        
else:
    print("Unknown type")
