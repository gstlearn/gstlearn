################################################################################
#                                                                              #
#                         gstlearn Python package                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# License: BSD 3-clause                                                        #
#                                                                              #
################################################################################
import numpy as np
import gstlearn as gl

import segyio
from segyio import BinField
from segyio import TraceField
from mpl_toolkits.axes_grid1.inset_locator import AnchoredZoomLocator
from pickle import NONE

def open(filename):
    f = segyio.open(filename, ignore_geometry=False)
    return f

def summary(f, verbose=False):
    ntot = f.tracecount
    nxline = f.xlines[-1] - f.xlines[0] + 1
    niline = f.ilines[-1] - f.ilines[0] + 1
    
    if verbose:
        print("Characteristics of the SEGY file")
        print("- Number of traces = ", ntot)
        print("- Number of XLines = ", nxline)
        print("- Number of ILines = ", niline)

    return ntot, nxline, niline


def getScale(head):
    
    scale = head[0].__getitem__(TraceField.SourceGroupScalar)
    if scale < 0:
        scale = 1. / (-scale)
    return scale

def getCoordinates(f, rank):
    
    head = f.header
    
    scale = getScale(head)
    x = head[rank].__getitem__(TraceField.SourceX) * scale 
    y = head[rank].__getitem__(TraceField.SourceY) * scale
    
    return x, y
    
def getDistance(f, r1, r2):
    x1, y1 = getCoordinates(f, r1)
    x2, y2 = getCoordinates(f, r2)
    
    d1 = x2 - x1
    d2 = y2 - y1
    print("Distance between ", r1, "and", r2, " = ", 
          np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2))


def getCornerPoints(f, verbose=False):
    
    # Get global characteristics
    ntot, nxline, niline = summary(f, verbose=False)
    
    xcorner = np.ones(4)
    ycorner = np.ones(4)

    # Coordinates of the origin of the grid
    xcorner[0], ycorner[0] = getCoordinates(f, 0)

    # Coordinates of the end of column
    xcorner[1], ycorner[1] = getCoordinates(f, ntot-nxline)

    # Coordinates of the end of line
    xcorner[2], ycorner[2] = getCoordinates(f, nxline-1)

    # Coordinates of the opposite corner
    xcorner[3], ycorner[3] = getCoordinates(f, ntot-1)

    # Optional printout
    if verbose:
        print("Coordinates of the corners (according to XLine, ILine indices)")
        print("- XL=",0,      " IL=", 0,      " : ", xcorner[0], ycorner[0])
        print("- XL=",0,      " IL=", niline, " : ", xcorner[1], ycorner[1])
        print("- XL=",nxline, " IL=", 0,      " : ", xcorner[2], ycorner[2])
        print("- XL=",nxline, " IL=", niline, " : ", xcorner[3], ycorner[3])
        
    return xcorner, ycorner

def getGridCharacteristics(f, verbose=False):
    
    # Get global characteristics
    ntot, nxline, niline = summary(f, verbose=False)
    
    # Get vertical information
    bin = f.bin
    nz  = bin.__getitem__(BinField.Samples)
    dz  = bin.__getitem__(BinField.Samples) / 1000
    z00 = -dz * nz
    
    # Coordinates of the several key nodes
    x00, y00 = getCoordinates(f, 0)
    x01, y01 = getCoordinates(f, 1)
    xn1, yn1 = getCoordinates(f, nxline)
    xny, yny = getCoordinates(f, ntot-nxline)
    
    # Get horizontal information
    dy = np.sqrt((x01 - x00) ** 2 + (y01 - y00) ** 2)
    dx = np.sqrt((xn1 - x00) ** 2 + (yn1 - y00) ** 2)
    
    # Get the rotation angle
    deltax = xny - x00
    deltay = yny - y00
    theta  = np.arctan(deltay/deltax)
    thetaD = theta * 180 / np.pi
    
    nx = niline
    ny = nxline
    
    # Printout (optional)
    if verbose:
        print("Grid characteristics")
        print("- Origin:", x00, y00, z00)
        print("- Mesh  :", dx, dy, dz)
        print("- Count :", nx, ny, nz)
        print("- AngleD:", thetaD)
    
    return x00, y00, z00, dx, dy, dz, nx, ny, nz, theta, thetaD

def getCornersFromGrid(nx, ny, x0, y0, dx, dy, theta):

    deltax = dx * (nx - 1)
    deltay = dy * (ny - 1)

    # Reconstruction for check
    xrec = np.ones(4)
    yrec = np.ones(4)

    # The origin coincides (by construction)
    xrec[0] = x0
    yrec[0] = y0

    # The end of the column
    xrec[1] = x0 + deltax * np.cos(theta)
    yrec[1] = y0 + deltax * np.sin(theta)

    # The end of the line
    theta2 = theta + np.pi / 2
    xrec[2] = xrec[0] + deltay * np.cos(theta2)
    yrec[2] = yrec[0] + deltay * np.sin(theta2)

    # The opposite corner
    xrec[3] = xrec[1] + deltay * np.cos(theta2)
    yrec[3] = yrec[1] + deltay * np.sin(theta2)
    
    return xrec, yrec

def defineGrid(f, verbose=False):
    # Get the grid characteristics
    x0, y0, z0, dx, dy, dz, nx, ny, nz, theta, thetaD = getGridCharacteristics(f)
    
    # Create the grid
    nxvec = [nx,ny]
    grid = gl.DbGrid.create(nx = nxvec, dx = [dx,dy], x0 = [x0,y0], angles = [thetaD,0])
    
    # Add a variable which contains the order of the traces
    rankTrace = gl.Grid.generateGridIndices(nxvec, "+x2+x1")
    grid["rankTrace"] = rankTrace
    
    # Optional display
    if verbose:
        grid.display()
        
    return grid

def readAllTraceHeaders(f):
    
    head = f.header
    scale = getScale(head)
    
    ind = range(len(head))
    x=[head[i].__getitem__(TraceField.SourceX) * scale for i in ind]
    y=[head[i].__getitem__(TraceField.SourceY) * scale for i in ind]
    
    return x, y

def checkCompatible(fileSEGYs):
    nfileSEGY = len(fileSEGYs)
    if nfileSEGY <= 1:
        return True
    
    # Getting the file characteristics of the first file (used as reference)
    
    f0 = open(fileSEGYs[0])
    x0r, y0r, z0r, dxr, dyr, dzr, nxr, nyr, nzr, thetar, thetaDr = getGridCharacteristics(f0)
    
    for iseg in range(1,nfileSEGY):
        fi = open(fileSEGYs[iseg])
        x0, y0, z0, dx, dy, dz, nx, ny, nz, theta, thetaD = getGridCharacteristics(fi)
        
        if x0 != x0r or y0 != y0r or z0 != z0r or dx != dxr or dy != dyr or dz != dzr or nx != nxr or ny != nyr or nz != nzr or thetaD != thetaDr:
            print("SEGY File #", iseg+1, " ...")
            _, _, _, _, _, _, _, _, _, _, _ = getGridCharacteristics(fi, verbose= True)
            print("... is not compatible with SEGY File # 1")
            _, _, _, _, _, _, _, _, _, _, _ = getGridCharacteristics(f0, verbose= True)           
            return False
    
    return True

def create3DGrid(fileSEGYs, dbsegy, limits2D=None, topName = None, botName = None, 
                 limitZ = None, verbose=False):
    
    # Get the information on the vertical extension of SEGY file
    f = open(fileSEGYs[0])
    _, _, z0, _, _, dz, _, _, nz, _, _ = getGridCharacteristics(f)
    f.close()
    z00 = z0
    zmax = z0 + (nz - 1) * dz
    
    # Characteristics of the initial SEGY (based on the SEGY grid 'dbsegy')
    nx = dbsegy.getNXs()
    dx = dbsegy.getDXs()
    x0 = dbsegy.getX0s()
    angles = dbsegy.getAngles()

    # Calculating the vertical extension (if not specified in input)
    if limitZ is None:
        limitZ = np.empty(2)
        limitZ[0] = 0
        limitZ[1] = nz
    
        if botName is not None:
            botArray = dbsegy.getColumn(botName)
            zbot = gl.VectorHelper.minimum(botArray)
            if zbot > z0:
                limitZ[0] = (zbot - z0) // dz

    
        if topName is not None:
            topArray = dbsegy.getColumn(topName)
            ztop = gl.VectorHelper.maximum(topArray)
            if ztop < zmax:
                limitZ[1] = (ztop - z0 + dz/2) // dz + 1
    
    nx.resize(3)
    if limits2D is not None:
        limitX = limits2D[0]
        limitY = limits2D[1]
    else:
        limitX = np.empty(2)
        limitX[0] = 0
        limitX[1] = nx[0]
        limitY = np.empty(2)
        limitY[0] = 0
        limitY[1] = nx[1]

    nx[0] = limitX[1] - limitX[0]
    nx[1] = limitY[1] - limitY[0]
    nx[2] = limitZ[1] - limitZ[0]
    
    dx.resize(3)
    dx[2] = dz 
    
    indg = np.empty(3)
    indg[0] = limitX[0]
    indg[1] = limitY[0]
    indg[2] = limitZ[0]
    x0 = dbsegy.indicesToCoordinate(indg)
    x0.resize(3)
    x0[2] = z0 + dz * limitZ[0]
    
    if verbose:
        print("Creating a 3-D:")
        print("- Origin :", x0)
        print("- Mesh   :", dx)
        print("- Count  :", nx)
        print("- Angle  :", angles[0])
        
        print("- Limits along X :", limitX)
        print("- Limits along Y :", limitY)
        print("- Limits along Z :", limitZ)
    
    # Creating the 3-D Grid
    grid = gl.DbGrid.create(nx=nx, x0=x0, dx=dx, angles=angles)
    
    # Copy the traces
    nfileSEGY = len(fileSEGYs)
    for iseg in range(nfileSEGY):
        mat = segyio.tools.cube(fileSEGYs[iseg])
        matred = mat[int(limitX[0]):int(limitX[1]),:,:][:,int(limitY[0]):int(limitY[1]),:][:,:,int(limitZ[0]):int(limitZ[1])]
        matred = matred.transpose(2,1,0)  
        matred = matred.reshape(-1)
    
        name = "SEGY." + str(iseg+1)
        grid[name] = matred
    
    return grid
