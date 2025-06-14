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

import sys
import segyio
from segyio import BinField
from segyio import TraceField

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
    dz  = bin.__getitem__(BinField.Interval) / 1000
    z00 = -dz * (nz-1) - f.header[0][segyio.tracefield.TraceField.DelayRecordingTime]
    
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

def create2DGrid(fileSEGY, verbose=False):
    
    # Open the SEGY file
    f = open(fileSEGY)

    # Retrieve the grid characteristics
    x0, y0, z0, dx, dy, dz, nx, ny, nz, theta, thetaD = getGridCharacteristics(f)

    # Close the file
    f.close()
    
    # Create the grid
    nxvec = [nx,ny]
    grid = gl.DbGrid.create(nx = nxvec, dx = [dx,dy], x0 = [x0,y0], angles = [thetaD,0])
    
    # Add a variable which contains the order of the traces
    rankTrace = gl.Grid.gridIndices(nxvec, "+x2+x1")
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

def checkCompatible(fileSEGYs, verbose=False):
    nfileSEGY = len(fileSEGYs)
    if nfileSEGY <= 1:
        return True
    
    # Optional printout
    if verbose:
        print("Printout of the", nfileSEGY, "SEGY files")
        for iseg in range(0,nfileSEGY):
            print("SEG File #", iseg+1)
            fi = open(fileSEGYs[iseg])
            _, _, _, _, _, _, _, _, _, _, _ =  getGridCharacteristics(fi, verbose=True)
  
    # Getting the file characteristics of the first file (used as reference)
    f0 = open(fileSEGYs[0])
    x0r, y0r, z0r, dxr, dyr, dzr, nxr, nyr, nzr, thetar, thetaDr = getGridCharacteristics(f0)

    # Compare with the other SEGY files
    for iseg in range(1,nfileSEGY):
        fi = open(fileSEGYs[iseg])
        x0, y0, z0, dx, dy, dz, nx, ny, nz, theta, thetaD = getGridCharacteristics(fi)
        
        if x0 != x0r or y0 != y0r or z0 != z0r or dx != dxr or dy != dyr or dz != dzr or nx != nxr or ny != nyr or nz != nzr or thetaD != thetaDr:
            print("\nWarning: SEGY File #", iseg+1, " ...")
            _, _, _, _, _, _, _, _, _, _, _ = getGridCharacteristics(fi, verbose= True)
            print("... is not 'strictly' compatible with SEGY File # 1.")
            _, _, _, _, _, _, _, _, _, _, _ = getGridCharacteristics(f0, verbose= True) 

            if nx == nxr and ny == nyr and nz == nzr:
                print("However, the count of grid nodes being similar, the difference is considered as minor.")
                print("The SEGYS are considered as compatible: the information of SEGY#1 is used")
                return True
                  
            return False
    
    return True


def create3DGrid(fileSEGYs, dblabel, topName = None, botName = None, limitZ = None, 
                 verbose=False):
    '''
    Create a 3D grid from the set of SEGY files and the 2D surface information ('dbsegy')

    fileSEGYs: Vector of SEGY files
    dblabel: 2D grid file containing the surface information
    topName: Name of the Top surface (in 'dbsegy')
    botName: Name of the Bottom surface (in 'dbsegy')
    limitZ: Vector giving the vertical indices (optional)
    restrictHorizontal: True if the horizontal extension must be restricted to the
       only pixels where top and bottom are defined and correctly ordered
    
    Remarks: 
    - 'topName' and 'botName' are used (when both defined and correctly ordered) 
      for the vertical extension of 3D grid
    - If 'limitZ' is not defined, they are used to derive the vertical extension.
    '''

    # Get the information on the vertical extension of SEGY file
    f = open(fileSEGYs[0])
    _, _, z0, _, _, dz, _, _, nz, _, _ = getGridCharacteristics(f)
    f.close()
    z00 = z0
    zmax = z0 + (nz - 1) * dz
    
    # Extract the limits for horizontal extension
    Limits2D = dblabel.getLimitsFromVariableExtend(topName, botName)
    limitX = Limits2D[0]
    limitY = Limits2D[1]

    # Build the 2-D grid restricted to the area of interest
    dbsegy2D = gl.DbGrid.createSubGrid(dblabel, Limits2D, False)

    # Characteristics of the initial SEGY (based on the SEGY grid 'dbsegy')
    nx = dbsegy2D.getNXs()
    dx = dbsegy2D.getDXs()
    x0 = dbsegy2D.getX0s()
    angles = dbsegy2D.getAngles()
    
    nx[0] = limitX[1] - limitX[0]
    nx[1] = limitY[1] - limitY[0]

    # Get the grid of maximum extension (even if variables are not defined)
    iuidSel = dbsegy2D.setSelectionFromVariableExtend(topName, botName)
    
    # Calculating the vertical extension (if not specified in input)
    if limitZ is None:
        limitZ = [0, nz]
            
        if botName is not None:
            botArray = dbsegy2D.getColumn(botName, useSel=True)
            zbot = gl.VectorHelper.minimum(botArray)
            if zbot > z0:
                limitZ[0] = int((zbot - z0) // dz)
    
        if topName is not None:
            topArray = dbsegy2D.getColumn(topName, useSel=True)
            ztop = gl.VectorHelper.maximum(topArray)
            if ztop < zmax:
                limitZ[1] = int((ztop - z0 + dz/2) // dz + 1)

    # Number of nodes
    nx.resize(3)
    nx[2] = limitZ[1] - limitZ[0]
    
    # Size of the grid meshes
    dx.resize(3)
    dx[2] = dz 
    
    # Origin of the 3D grid
    x0.resize(3)
    x0[2] = z0 + dz * limitZ[0]
    
    flagIncorrect = False
    if nx[0] <= 0 or nx[1] <= 0 or nx[2] <= 0:
        flagIncorrect = True

    if verbose or flagIncorrect:
        print("Creating the 3-D between",botName,"and",topName)
        print("- Origin :", x0)
        print("- Mesh   :", dx)
        print("- Count  :", nx)
        print("- Angle  :", angles[0])
        print("- Limits along Z :", limitZ)
        print("(This definition takes the SEGY grid characteristics into account)")

    if flagIncorrect:
        sys.exit("Error in the 3-D Grid characteristics")
    
    # Creating the 3-D Grid
    dbsegy3D = gl.DbGrid.create(nx=nx, x0=x0, dx=dx, angles=angles, flagAddCoordinates=False)
    
    # Extract the traces and copy them to the output 3D file        mat = segyio.tools.cube(fileSEGYs[iseg])

    nfileSEGY = len(fileSEGYs)
    flipforward = True
    flipback = False
    for iseg in range(nfileSEGY):
        mat = segyio.tools.cube(fileSEGYs[iseg])
        if flipforward:
            print("Flip forward is performed")
            mat = np.flip(mat, 2)
       
        matred = mat[int(limitX[0]):int(limitX[1]),:,:][:,int(limitY[0]):int(limitY[1]),:][:,:,int(limitZ[0]):int(limitZ[1])]
        
        if flipback:
           print("Flip backward is perfomed")
           np.flip(matred,2)
        matred = matred.transpose(2,1,0)  
        matred = matred.reshape(-1)
    
        name = "SEGY." + str(iseg+1)
        if verbose:
            print("Processing variable", name)

        dbsegy3D[name] = matred
    
    return dbsegy2D, dbsegy3D

def restrict2DGrid(dblabel, topName, botName):
    '''
    Create a new 2D grid restricted the bounding box of the samples
    where both 'top' and 'bottom' variables are defined and correctly ordered
    '''
    
    # Extract the limits for horizontal extension
    Limits2D = dblabel.getLimitsFromVariableExtend(topName, botName)
    
    return gl.DbGrid.createSubGrid(dblabel, Limits2D, False)

def SaveSegy(source, destination, prob) :

    # Read the extension of the input array 'prob1'
    nX, nY, nSamples = prob.shape
    with segyio.open(source, ignore_geometry=True) as src:
        spec = segyio.tools.metadata(src)
        spec.samples = spec.samples[:nSamples]
        with segyio.create(destination, spec) as dst:
            dst.text[0] = src.text[0]
            dst.bin = src.bin
            dst.bin.update(hns=len(spec.samples))
            dst.header = src.header
            for h in dst.header :
                h[segyio.TraceField.TRACE_SAMPLE_COUNT] = nSamples
            dst.trace = src.trace

        # Here we assign the at each segy trace the probability

            iGlobal = 0
            for iy in range(nY) :
                for ix in range(nX) :
                    dst.trace[iGlobal] = np.array(prob[ix, iy], dtype=np.float32)
                    iGlobal += 1
