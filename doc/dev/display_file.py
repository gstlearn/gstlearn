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

def checkValidPointer(pointer):
    if not pointer:
        print(" ")
        print("The object has not been read correctly")
        print("Procedure is stopped")
        exit()
        
# Get the Type of the File and stop if returned empty (file not found)
filetype = gl.ASerializable.getFileIdentity(filename)
if filetype == "":
    exit()

if filetype == "Db":
    db = gl.Db.createFromNF(filename,False)
    checkValidPointer(db)
    dbfmt = gl.DbStringFormat()
    dbfmt.setFlags(flag_vars=True, flag_stats=True)
    db.display(dbfmt)
           
elif filetype == "DbGrid":
    dbgrid = gl.DbGrid.createFromNF(filename,False)
    checkValidPointer(dbgrid)
    dbfmt = gl.DbStringFormat()
    dbfmt.setFlags(flag_vars=True, flag_stats=True)
    dbgrid.display(dbfmt)
            
elif filetype == "Vario":
    vario = gl.Vario.createFromNF(filename,False)
    checkValidPointer(vario)
    vario.display()
    
elif filetype == "Model":
    model = gl.Model.createFromNF(filename,False)
    checkValidPointer(model)
    model.display()
    
elif filetype == "Rule":
    rule = gl.Rule.createFromNF(filename,False)
    checkValidPointer(rule)
    rule.display()
    
elif filetype == "Table":
    table = gl.Table.createFromNF(filename,False)
    checkValidPointer(table)
    table.display()

elif filetype == "Polygon":
    poly = gl.Polygons.createFromNF(filename,False)
    checkValidPointer(poly)
    poly.display()
        
elif filetype == "MeshETurbo":
    mesh = gl.MeshETurbo.createFromNF(filename, False)
    checkValidPointer(mesh)
    mesh.display()
    
else:
    print("This type of file is UNKNOWN in display_file:", filetype)
