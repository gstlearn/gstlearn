#!/usr/bin/env python
# coding: utf-8

#python3 -m pip install mlxtend

import numpy as np
import matplotlib.pyplot as plt
from mlxtend.plotting import scatterplotmatrix
from mlxtend.data import iris_data

import gstlearn as gl
import gstlearn.plot as gp
import gstlearn.test as gt


# We use the data coming from Fisher (1936) and called *iris*

X, y = iris_data()

ndat = 120
Xtrain = X[range(ndat),:]
Xtest = X[range(ndat+1,X.shape[0]),:]


# ## Load the values inside a Data Base

# We load the data into the Data Base of gstlearn. This operation creates a Data base with 5 columns (the first one, named*rank* is generated automatically; the other four columns come from the external file).
# Note that we also construct the flag (*dbfmt*) used to produce the complete printout of the Db contents.

nech = Xtrain.shape[0]
ncol = Xtrain.shape[1]
names = gl.generateMultipleNames("x",ncol, delim="")
locatornames = gl.generateMultipleNames("z",ncol, delim="")
db = gl.Db.createFromSamples(nech=Xtrain.shape[0], tab=Xtrain.flatten(), names=names, locatorNames=locatornames)
dbfmt = gl.DbStringFormat.createFromFlags(flag_array=True, flag_resume=False, flag_vars=False)
gl.OptCst.defineByKey("NTCOL",-1)
db

# PCA in gstlearn

# Construct the PCA and print its contents.

pcag = gl.PCA()
pcag.pca_compute(db, True)
pcag.display()

# Transform the initial data into factors.

err = pcag.dbZ2F(db)
db

# Transform the vector back to the data.

err = pcag.dbF2Z(db)
db

# Compare the initial variables with the ones obtained by back-transforming the factors

gt.checkEqualityVector(db["x1"],db["F2Z.Z2F.x1"])
gt.checkEqualityVector(db["x2"],db["F2Z.Z2F.x2"])
gt.checkEqualityVector(db["x3"],db["F2Z.Z2F.x3"])
gt.checkEqualityVector(db["x4"],db["F2Z.Z2F.x4"])

# Get the vector of variance ratio.

varratiog = pcag.getVarianceRatio()

# PCA directly in Python

## On the original data set

# Calculate the Covariance matrix (dividing by n-1)

means = np.mean(Xtrain,axis=0)
Xc = Xtrain-means
cova = 1/(ndat-1) * Xc.T@Xc

# Computing the Variance Ratio

eig = np.linalg.eig(cova)
varratio = eig[0]/np.sum(eig[0])

# Comparing with gstlearn

gt.checkEqualityVector(varratio, varratiog, message="Error in Variance Ratio")

# Transform the original data set

Xt = Xc@eig[1]

# Comparing with gstlearn

gt.checkEqualityVector(Xt, db["Z2F*"], flagAbsolute=True, message="Transformed original data set")

# Transform another data set

# Creating the Db for the other data set

nech = Xtest.shape[0]
ncol = Xtest.shape[1]
names = gl.generateMultipleNames("x",ncol, delim="")
locatornames = gl.generateMultipleNames("z",ncol, delim="")
dbtest = gl.Db.createFromSamples(nech=Xtest.shape[0], tab=Xtest.flatten(), names=names, locatorNames=locatornames)
dbtest

# Applying the PCA using gstlearn

err = pcag.dbZ2F(dbtest)

# Applying the Python PCA on the secondary data set

Xtestnm =(Xtest-means)@eig[1]

# Comparing with gstlearn

gt.checkEqualityVector(Xtestnm, dbtest["Z2F*"], flagAbsolute=True, message="Transformed another data set")
