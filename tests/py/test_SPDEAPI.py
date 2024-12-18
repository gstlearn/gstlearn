#!/usr/bin/env python
# coding: utf-8

# # Checking Robustness of SPDE API

# In[1]:


import gstlearn as gl
import gstlearn.plot as gp
import numpy as np
import matplotlib.pyplot as plt

import gstlearn.document as gdoc

gdoc.setNoScroll()


# This script is meant to demonstrate the robustness of the SPDE interface against:
# - undefined data
# - selection on input data
# - combination of both
# - mask on the output grid
# - data located outside the grid
# 
# We also compare the result with the results of a traditional Kriging through a scatter plot between the two estimations performed on the Grid.

# In[2]:


range = 0.3
sill = 3.
model = gl.Model.createFromParam(gl.ECov.MATERN, range, sill)
model.display()


# In[3]:


def visualize(name):
    return # Do nothing in a PY file
    fig, ax = gp.initGeographic()
    ax.raster(grid, name=name, flagLegend=True)
    ax.symbol(data, nameSize="z", c='blue', s=100)
    plt.show()


# In[4]:


def estimate():
    err = gl.krigingSPDE(data, grid, model, namconv="SPDE")
    visualize("SPDE.*")
    
    neighU = gl.NeighUnique()
    err = gl.kriging(data, grid, model, neighU, flag_std=False, namconv="Kriging")
    gp.correlation(grid, "SPDE.*", "Kriging.*", bins=100, cmin=0.001)
    
    corr = grid.getCorrelation("SPDE.*", "Kriging.*")
    print("Correlation = ", round(corr,3))


# ## All values are defined

# In[5]:


grid = gl.DbGrid.create([100,100], [0.01, 0.01])

data = gl.Db.createFillRandom(ndat=10)
print(data[:])


# In[6]:


estimate()


# ## Some undefined data

# In[7]:


grid = gl.DbGrid.create([100,100], [0.01, 0.01])

data = gl.Db.createFillRandom(ndat=10)
data[2, "z"] = np.nan
data[5, "z"] = np.nan
data[7, "z"] = np.nan
data.display()
print(data[:])


# In[8]:


estimate()


# ## Some undefined data and a Mask on Data

# In[9]:


grid = gl.DbGrid.create([100,100], [0.01, 0.01])

data = gl.Db.createFillRandom(ndat=10, selRatio=0.3)
data[2, "z"] = np.nan
data[5, "z"] = np.nan
data[7, "z"] = np.nan
data.display()
print(data[:])


# In[10]:


estimate()


# ## Some undefined data and a Mask on Data and Grid

# In[11]:


grid = gl.DbGrid.create([100,100], [0.01, 0.01])
grid.addSelectionRandom(0.9)

data = gl.Db.createFillRandom(ndat=10, selRatio=0.3)
data[2, "z"] = np.nan
data[5, "z"] = np.nan
data[7, "z"] = np.nan
data.display()
print(data[:])


# In[12]:


estimate()


# ## Some Data located outside the Grid

# In[13]:


grid = gl.DbGrid.create([100,100], [0.01, 0.01])

data = gl.Db.createFillRandom(ndat=10)
data[2, "x-1"] = 1.5
data[5, "x-2"] = -0.8
data.display()
print(data[:])


# In[14]:


estimate()


# ## Estimating at (missing) data location

# In[15]:


data = gl.Db.createFillRandom(ndat=10)
data[2, "z"] = np.nan
data[5, "z"] = np.nan
data[7, "z"] = np.nan
data.display()
print(data[:])


# In[16]:


err = gl.krigingSPDE(data, data, model, namconv="Missing")
print(data[:])


# We can check that kriging using SPDE produces an almost exact interpolator. It allows filling the missing Data.
