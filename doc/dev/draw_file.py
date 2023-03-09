# This script is meant to visualize any serialized object 
# from the current Directory

import sys
import os
import gstlearn as gl
import gstlearn.plot as gp
import matplotlib.pyplot as plt
from attr._make import NOTHING
from pandas.core.sorting import nargsort
from pandas.core.indexing import check_deprecated_indexers

def invalid(ranks, number):
    last = number - 1
    for rank in ranks:
        if int(rank) > last:
            print("Incorrect 'rank' argument (",int(rank),"): it should be smaller than",number)
            return True
    return False

def checkValidPointer(pointer):
    if not pointer:
        print(" ")
        print("The object has not been read correctly")
        print("Procedure is stopped")
        exit()
        
args = sys.argv
nargs = len(args)
if nargs < 2:
    print("This script should be called according to the following syntax:")
    print("  draw_file.py filename [ranks]")
    print("- filename: Name of the Serialized file")
    print("- ranks: list of ranks")
    exit()
filename = args[1]
fileaux = ""
if nargs > 2:
    fileaux = args[2]
ranks = args[2:nargs]

# Get the Type of the File and stop if returned empty (file not found)
filetype = gl.ASerializable.getFileIdentity(filename)
if filetype == "":
    exit()
filetaux = gl.ASerializable.getFileIdentity(fileaux)

if filetype == "Db":
    db = gl.Db.createFromNF(filename,False)
    checkValidPointer(db)
    if invalid(ranks, db.getColumnNumber()): 
        exit()
    if len(ranks) == 0:
        name = db.getLastName()
        flagDb = True
    elif len(ranks) == 1:
        name = db.getNameByColIdx(int(ranks[0]))
        flagDb = True
    elif len(ranks) == 2:
        nameX = db.getNameByColIdx(int(ranks[0]))
        nameY = db.getNameByColIdx(int(ranks[1]))
        flagDb = False
    else:
        print("Number of Variable ranks should be 0, 1 or 2")
        exit()
    if flagDb:
        ax = gp.point(db, name)
        ax.decoration(title=name)
        plt.show()
    else:
        gp.correlation(db, nameX, nameY)
        plt.show()
            
elif filetype == "DbGrid":
    dbgrid = gl.DbGrid.createFromNF(filename,False)
    checkValidPointer(dbgrid)
    if invalid(ranks, dbgrid.getColumnNumber()): exit()
    if len(ranks) == 0:
        name = dbgrid.getLastName()
    elif len(ranks) == 1:
        name = dbgrid.getNameByColIdx(int(ranks[0]))
    else:
        print("Number of Variable rank should be 0 or 1")
        exit()
        
    if dbgrid.getNDim() > 1:
        ax = gp.grid(dbgrid, name)
        ax.decoration(title=name)
        plt.show()
    else:
        ax = gp.grid1D(dbgrid, name)
        ax.decoration(title=name)
        plt.show()
    
elif filetype == "Vario":
    vario = gl.Vario.createFromNF(filename,False)
    checkValidPointer(vario)
    
    if filetaux == "Model":
        model = gl.Model.createFromNF(fileaux,False)
        checkValidPointer(model)
        gp.varmod(vario, model)
        plt.show()
    else:
        gp.vario(vario)
        plt.show()
    
elif filetype == "Model":
    model = gl.Model.createFromNF(filename,False)
    checkValidPointer(model)
    
    if filetaux == "Vario":
        vario = gl.Vario.createFromNF(fileaux,False)
        checkValidPointer(vario)
        gp.varmod(vario, model)
        plt.show()
    else:
        gp.model(model)
        plt.show()
    
elif filetype == "Rule":
    rule = gl.Rule.createFromNF(filename,False)
    checkValidPointer(rule)
    gp.rule(rule)
    plt.show()
    
elif filetype == "Table":
    table = gl.Table.createFromNF(filename,False)
    checkValidPointer(table)
    if invalid(ranks, table.getNCols()): exit()
    ax = gp.table(table,ranks)
    ax.decoration(title=filename)
    plt.show()

elif filetype == "Polygon":
    poly = gl.Polygons.createFromNF(filename,False)
    checkValidPointer(poly)
    gp.polygon(poly,colorPerSet=True,flagFace=False)
    plt.show()
      
elif filetype == "MeshETurbo":
    mesh = gl.MeshETurbo.createFromNF(filename, False)
    checkValidPointer(mesh)
    gp.mesh(mesh)
    plt.show()
 
else:
    print("This type of file is UNKNOWN in draw_file:", filetype)
