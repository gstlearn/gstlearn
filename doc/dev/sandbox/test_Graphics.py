# Regular Packages
import numpy as np
import sys
import os
import gstlearn as gl
import gstlearn.plot as gp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Parameters
# This test does not make sense as a non-regression test.
# Nevertheless it tells how to produce several figures simultaneously
# Turn the relevant flags ON when used interactively

draw_grid  = False
draw_point = False

# Create representation grid

grid = gl.DbGrid.create([50,50])

# Create a Data Set

data = gl.Db.createFromBox(100,grid.getCoorMinimum(),grid.getCoorMaximum())

# Create the model

model = gl.Model.createFromParam(gl.ECov.CUBIC,15)

# Non conditionnal simulations (4)

err = gl.simtub(None,grid,model,None,4)

# Representation of the 4 simulations of the Grid

if draw_grid:
	norm = mpl.colors.Normalize(vmin=0, vmax=1)
	fig, axs = plt.subplots(2,2,figsize=(30,15))
	ecr = 1
	for i in range(2):
		for j in range(2):
			ax = gp.grid(grid,"Simu."+str(ecr),ax=axs[j,i],flagColorBar=False)
			ecr = ecr + 1
#	cbar = fig.colorbar(im,ax=axs[0:,:],location='right',shrink=0.5)
	plt.show()

# Representation of two figures on the Data

if draw_point:
	fig, axs = plt.subplots(1,2,figsize=(15,15))
	for i in range(2):
		ax = gp.point(data,ax=axs[i])
	plt.show()

print("Test successfully performed")
