#!/usr/bin/env python
# coding: utf-8

import gstlearn as gl
import numpy as np

def compute(i,proj):
    mv = gl.VectorDouble(proj.getNApex())
    pv = gl.VectorDouble(proj.getNPoint())
    for k in range(mv.size()):
        mv[k]=0
    mv[i]=1.
    proj.mesh2point(mv,pv)
    return pv

def computeT(i,proj):
    mv = gl.VectorDouble(proj.getNApex())
    pv = gl.VectorDouble(proj.getNPoint())
    for k in range(pv.size()):
        pv[k]=0
    pv[i]=1.
    proj.point2mesh(pv,mv)
    return mv

wavelet = [-0.000000,-0.010629,0.079015,0.209873,
           0.278554,0.257160,0.116605,-0.171370,-0.599146,
           -1.070201,-1.412208,-1.424379,-0.922758,0.193134,
           1.757338,3.207395,3.810365,3.207395,1.757338,
           0.193134,-0.922758,-1.424379,-1.412208,-1.070201,
           -0.599146,-0.171370,0.116605,0.257160,0.278554,
           0.209873,0.079015,-0.010629,-0.000000]

test  = gl.DbGrid.create([5,5,40])
projc = gl.ProjConvolution(wavelet,test,nodeRes2D=[3,3])
nrow  = test.getNSample()
ncol  = projc.getResolutionGrid().getNSample()
A = np.zeros(shape=(nrow,ncol))
B = np.zeros(shape=(ncol,nrow))

for k in range(ncol):
    v = compute(k,projc)
    A[:,k]=[i for i in v]

for k in range(nrow):
    v = computeT(k,projc)
    B[:,k]=[i for i in v]

# The result of the following chunk has to be zero :

value = np.sum(np.abs(A-B.T))
print("Difference between two calculations (should be close to zero) =", value)
assert value==0

