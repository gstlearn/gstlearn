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
filetype = gl.ASerializable.getFileIdentify(filename)
if filetype == "":
    exit()

if filetype == "Db":
    
    # Case of Db file
    db = gl.Db.createFromNF(filename,False)
    dbfmt = gl.DbStringFormat()
    dbfmt.setParams(gl.FLAG_VARS)
    db.display(dbfmt)
            
elif filetype == "Vario":
    
    # Case of Vario file
    vario = gl.Vario.createFromNF(filename,False)
    vario.display()
    
elif filetype == "Model":
    
    # Case of Model file
    model = gl.Model.createFromNF(filename,False)
    model.display()
    
elif filetype == "Rule":
    
    # Case of Rule file
    rule = gl.Rule(filename,False)
    rule.display()
    
elif filetype == "Table":

    # Case of a Table
    table = gl.Table(filename,False)
    table.display()

elif filetype == "Polygon":
    
    # Case of a polygon
    poly = gl.Polygons.createFromNF(filename,False)
    poly.display()
        
else:
    print("Unknown type")
