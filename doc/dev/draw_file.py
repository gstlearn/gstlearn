# This script is meant to visualize any serialized object 
# from the current Directory

import sys
import os
import gstlearn as gl
import gstlearn.plot as gp
import matplotlib.pyplot as plt
from attr._make import NOTHING
from pandas.core.sorting import nargsort

def invalid(ranks, number):
    last = number - 1
    for rank in ranks:
        if int(rank) > last:
            print("Incorrect 'rank' argument (",int(rank),"): it should be smaller than",number)
            return True
    return False

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
        gp.point(db, name, title=name, end_plot=True)
    else:
        gp.correlation(db, nameX, nameY, end_plot=True)
            
elif filetype == "DbGrid":
    dbgrid = gl.DbGrid.createFromNF(filename,False)
    if invalid(ranks, dbgrid.getColumnNumber()): exit()
    if len(ranks) == 0:
        name = dbgrid.getLastName()
    elif len(ranks) == 1:
        name = dbgrid.getNameByColIdx(int(ranks[0]))
    else:
        print("Number of Variable rank should be 0 or 1")
        exit()
        
    if dbgrid.getNDim() > 1:
        gp.grid(dbgrid, name, title=name, end_plot=True)
    else:
        gp.grid1D(dbgrid, name, title=name, end_plot=True)
    
elif filetype == "Vario":
    vario = gl.Vario.createFromNF(filename,False)
    
    if filetaux == "Model":
        model = gl.Model.createFromNF(fileaux,False)
        gp.varmod(vario, model, end_plot=True)
    else:
        gp.vario(vario,end_plot=True)
    
elif filetype == "Model":
    model = gl.Model.createFromNF(filename,False)
    
    if filetaux == "Vario":
            vario = gl.Vario.createFromNF(fileaux,False)
            gp.varmod(vario, model, end_plot=True)
    else:
        gp.model(model,end_plot=True)
    
elif filetype == "Rule":
    rule = gl.Rule.createFromNF(filename,False)
    gp.rule(rule,end_plot=True)
    
elif filetype == "Table":
    table = gl.Table.createFromNF(filename,False)
    if invalid(ranks, table.getColNumber()): exit()
    gp.table(table,ranks,end_plot=True,title=filename)

elif filetype == "Polygon":
    poly = gl.Polygons.createFromNF(filename,False)
    gp.polygon(poly,colorPerSet=True,flagFace=True,end_plot=True)
       
else:
    print("Unknown type")
