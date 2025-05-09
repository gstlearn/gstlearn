################################################################################
#                                                                              #
#                         gstlearn Python package                              #
#                                                                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://gstlearn.org                                                #
# License: BSD 3-clause                                                        #
#                                                                              #
################################################################################

import numpy                as np
import gstlearn             as gl
from numpy import pi
from matplotlib.animation import adjusted_figsize

try:
    import plotly.graph_objects as go
except ModuleNotFoundError as ex:
    msg = ("Python dependency 'plotly' not found.\n"
          "To install it alongside gstlearn, please run `pip install gstlearn[plot]'")
    raise ModuleNotFoundError(msg) from ex

def getCscale():
    cscale = [
        [0.0, '#313695'],
        [0.07692307692307693, '#3a67af'],
        [0.15384615384615385, '#5994c5'],
        [0.23076923076923078, '#84bbd8'],
        [0.3076923076923077, '#afdbea'],
        [0.38461538461538464, '#d8eff5'],
        [0.46153846153846156, '#d6ffe1'],
        [0.5384615384615384, '#fef4ac'],
        [0.6153846153846154, '#fed987'],
        [0.6923076923076923, '#fdb264'],
        [0.7692307692307693, '#f78249'],
        [0.8461538461538461, '#e75435'],
        [0.9230769230769231, '#cc2727'],
        [1.0, '#a50026']
        ]
    return cscale

def __invalidFileDimension(db, ndim):
    
    if db.getNDim() != ndim:
        print("This representation is only designed for Db of dimension", ndim)
        return True
    return False

def __linearInterpolate(values, Min=0.0, Max=1.0, flagNoInterpolate = False):
    if flagNoInterpolate:
        return values
    values = np.asarray(values)
    vmin = np.min(values)
    vmax = np.max(values)
    if vmax == vmin:
        return np.full_like(values, (Min + Max) / 2.0)
    # Linear interpolation between Min and Max
    return Min + (values - vmin) * (Max - Min) / (vmax - vmin)

def SurfaceOnMesh(mesh, intensity=None, cscale=None, color='white', opacity=0.50, 
                  **plot_args):
    '''
    Represent a Function defined on a Mesh
    
    plot_args Arguments passed to Mesh3d()
    '''
    
    tab = np.array(mesh.getEmbeddedCoordinatesPerApex())
    meshes = np.array(mesh.getMeshesAsVVI())
    
    if cscale is None:
        cscale = getCscale()
    
    if intensity is None:
        intensity = np.zeros(mesh.getNApices())
    else:
        if type(intensity) == gl.VectorDouble:
            intensity = np.array(intensity.getVector())

    surface = go.Mesh3d(x=tab[0,:], y=tab[1,:], z=tab[2,:], 
                        color=color, colorbar_title='z', name='y',
                        i=meshes[:,0],j=meshes[:,1],k=meshes[:,2],
                        colorscale=cscale, intensity=intensity,
                        opacity=opacity, 
                        **plot_args
                        )

    return surface
    
def ScatterOnMesh(points, meshes, intensity=None, cscale=None, color='white', 
                  opacity=0.50, **plot_args):
    '''
    Represent a Function defined on a Mesh defined by its apices and meshes
    
    plot_args Arguments passed to Mesh3d()
    '''
    
    napices = meshes.shape[0]
    if cscale is None:
        cscale = getCscale()
    
    if intensity is None:
        intensity = np.zeros(napices)
    else:
        if type(intensity) == gl.VectorDouble:
            intensity = np.array(intensity.getVector())

    surface = go.Mesh3d(x=points[:,0], y=points[:,1], z=points[:,2], 
                        color=color, colorbar_title='z', name='y',
                        i=meshes[:,0],j=meshes[:,1],k=meshes[:,2],
                        colorscale=cscale, intensity=intensity,
                        opacity=opacity, 
                        **plot_args)

    return surface
    
def Meshing(mesh, color='black', width=1, **plot_args):
    '''
    Represent the contents of a Mesh object
    
    plot_args Arguments passed to dict(type='scatter3d')
    '''
    
    xs = list()
    ys = list()
    zs = list()
    for i in range(mesh.getNMeshes()):
        a = np.array(mesh.getEmbeddedCoordinatesPerMesh(i))
        xs.extend(a[:,0].tolist() + [None])
        ys.extend(a[:,1].tolist() + [None])
        zs.extend(a[:,2].tolist() + [None])

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    meshing = dict(type='scatter3d',x=xs, y=ys, z=zs, mode='lines',
                   line=dict(color=color, width=width), **plot_args)
    return meshing
    
def ScatterOnDb(db, mode='lines', color='black', width=1, 
                m_symbol = 'circle', m_color='black', m_line = 'black', 
                m_size=15, m_width=2,
                **plot_args):
    
    if __invalidFileDimension(db, 3):
        return None
                      
    meshing = dict(type='scatter3d',
                   x=db.getOneCoordinate(0), 
                   y=db.getOneCoordinate(1),
                   z=db.getOneCoordinate(2),
                   mode=mode,marker_symbol=m_symbol,
                   marker_line_color=m_line, marker_color=m_color, 
                   marker_line_width=m_width, marker_size=m_size,
                   line=dict(color=color, width=width),
                   **plot_args)
    return meshing

def Scatter(x, y, z, mode='lines', color='black', width=1, 
            m_symbol = 'circle', m_color='black', m_line = 'black', m_size=15, m_width=2,
            **plot_args):
    '''
    Represent a set of points in 3-D
    
    x, y, z Vector of coordinates (same dimension)
    mode can be 'lines' or 'markers'
    color Color assigned to the lines (if mode == 'lines')
    width Width of the line (if mode == 'lines')
    m_symbol Type of symbol (if mode == 'markers')
    m_color  Color assigned to the symbol (if mode == 'markers')
    m_line   Type of line for symbol (if mode == 'markers')
    m_size   Size of the symbol (if mode == 'markers')
    m_width  Width of the line for symbol (if mode == 'markers')
    
    plot_args Arguments passed to dict(type='scatter3d')
    '''
    meshing = dict(type='scatter3d',x=x, y=y, z=z, mode=mode,
                   marker_symbol=m_symbol,
                   marker_line_color=m_line, marker_color=m_color, 
                   marker_line_width=m_width, marker_size=m_size,
                   line=dict(color=color, width=width),
                   **plot_args)
    return meshing
    
def ScatterOnSphere(long, lat, mode='lines', color='black', width=1, 
                    m_symbol = 'circle', m_color='black', m_line = 'black', 
                    m_size=15, m_width=2, dilate=1, 
                    **plot_args):
    '''
    Represent a set of points on the Sphere
    
    plot_args Arguments passed to Scatter()
    '''
    tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
    meshing = Scatter(tab[0,:], tab[1,:], tab[2,:], mode=mode, 
                      color=color, width=width,
                      m_symbol=m_symbol, m_color=m_color, m_line=m_line, 
                      m_size=m_size, m_width=m_width,
                      **plot_args)

    return meshing

def Line(x, y, z, color='black', width=1, 
         **plot_args):
    
    '''
    Represent a set of points koined by a Line 
    
    plot_args Arguments passed to dict(type='scatter3d')
    '''
    line = dict(type='scatter3d',x=x, y=y, z=z, mode='lines',
                line=dict(color=color, width=width),
                **plot_args
                )
    return line
    
def PolygonOnSphere(poly, flagClose=False, 
                    color='black', width=1, dilate=1,
                    **plot_args):
    '''
    Represent a Polygon object on the Sphere
    
    plot_args Arguments passed to dict(type='scatter3d')
    '''

    xs = list()
    ys = list()
    zs = list()

    for i in range(poly.getNPolyElem()):
        a = poly.getX(i)
        b = poly.getY(i)
        tab = np.array(gl.GH.convertLongLatTo3D(a, b, dilate, np.nan))
        xp = tab[0,:]
        yp = tab[1,:]
        zp = tab[2,:]
        xs.extend(list(xp) + [None])
        ys.extend(list(yp) + [None])
        zs.extend(list(zp) + [None])
    
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    boundaries=dict(type='scatter3d', x=xs, y=ys, z=zs, mode='lines', 
                    line=dict(color=color, width=width),
                    **plot_args)
    return boundaries

def SliceOnDbGrid(grid, name, section=0, rank=0, useSel=False, 
                  cmin = None, cmax = None):


    if __invalidFileDimension(grid, 3):
        return None
                      
    shape = list(grid.getNXs())
    shape.pop(section)
    vect = grid.getSlice(name, section, rank, useSel)
    # The next two lines are meant to let NAN values be represented as transparent
    zc = vect[2]
    zc[np.isnan(vect[3])] = np.nan
    x = np.array(vect[0]).reshape(shape)
    y = np.array(vect[1]).reshape(shape)
    z = np.array(zc).reshape(shape)
    values = np.array(vect[3]).reshape(shape)
    
    slice = go.Surface(x=x, y=y, z=z, surfacecolor=values, 
                       coloraxis='coloraxis', cmin = cmin, cmax = cmax)
    return slice

def Slice3DOnDbGrid(grid, name, corner=None, 
                       useSel=False, cmin=None, cmax=None):
    
    '''
    Represent a series of three slices (XoY, YoZ, XoZ)
    performed around the corner (given by its grid indices)

    Returns a go.Figure element
    '''
    if __invalidFileDimension(grid, 3):
        return None
    if corner is None:
        corner = grid.getNXs() / 2
    
    data = [SliceOnDbGrid(grid,name,0,corner[2]),
            SliceOnDbGrid(grid,name,1,corner[1]),
            SliceOnDbGrid(grid,name,2,corner[0])
       ]
    fig = go.Figure(data=data)
    return fig
    
def IsoSurfaceOnDbGrid(grid, name, useSel=False, levels=None, 
                       colorscale='BlueRed', isomin=0, isomax=1, surface_count = 1, 
                       showlegend=False):
    
    if __invalidFileDimension(grid, 3):
        return None
                      
    shape = list(grid.getNXs())

    x = grid.getOneCoordinate(0, useSel).reshape(shape)
    y = grid.getOneCoordinate(1, useSel).reshape(shape)
    z = grid.getOneCoordinate(2, useSel).reshape(shape)
    values = grid.getColumn( name, useSel).reshape(shape)
    
    surfaces = go.Isosurface(x=x.flatten(), y=y.flatten(), z=z.flatten(), 
                             value = values.flatten(), 
                             isomin = isomin, isomax = isomax,
                             surface_count = surface_count, 
                             colorscale=colorscale,
                             showscale = showlegend, 
                             caps = dict(x_show=False, y_show=False)
                            )
    return surfaces
   
def SurfaceOnDbGrid(grid, name, useSel=False, showscale=False, **plot_args):
    
    if __invalidFileDimension(grid, 2):
        return None
    
    shape = list(np.flip(grid.getNXs()))
    values = grid.getColumn(name, useSel).reshape(shape)
    
    surface = go.Surface(z=values, showscale=showscale, opacity=0.9, **plot_args)
    
    return surface
    
def PointDb(db, nameColor=None, nameSize=None, useSel=True, 
            color='black', size=3, sizeMin=1, sizeMax=3, 
            nColors=10, flagNoColorInterpolate=False,
            opacity=1, posX=0, posY=1, posZ=2,
            fromLongLat = False, dilate = 1,
            **plot_args): 
    '''
    Represent a set of Points contained in a Db
    
    plot_args Arguments passed to Scatter3d()
    '''
    if fromLongLat:
        if __invalidFileDimension(db, 2):
            return None
        long = db.getOneCoordinate(0, useSel)
        lat  = db.getOneCoordinate(1, useSel)
        tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
        x = tab[0,:]
        y = tab[1,:]
        z = tab[2,:]
    else:
        if __invalidFileDimension(db, 3):
            return None
        x = db.getOneCoordinate(posX, useSel)
        y = db.getOneCoordinate(posY, useSel)
        z = db.getOneCoordinate(posZ, useSel)

    if nameColor is not None:
        colors = __linearInterpolate(db.getColumn(nameColor, useSel), 
                                     1, nColors, flagNoColorInterpolate)
    else:
        colors = color
    
    if nameSize is not None:
        sizes = __linearInterpolate(db.getColumn(nameSize, useSel), 
                                    sizeMin, sizeMax)
    else:
        sizes = size
    
    object = go.Scatter3d(x=x, y=y, z=z, mode='markers',
                          marker = dict(size = sizes,
                                        color = colors,
                                        colorscale ='Viridis',
                                        opacity = opacity
                                   ),
                          **plot_args)
    return object

def GradientDb(db, useSel=True, colorscale='Blues', sizemode='absolute', size=2, 
               **plot_args): 
    '''
    Represent a set of Gradients contained in a Db
    
    plot_args Arguments passed to Cone()
    '''

    if __invalidFileDimension(db, 3):
        return None
                      
    x = db.getOneCoordinate(0, useSel)
    y = db.getOneCoordinate(1, useSel)
    z = db.getOneCoordinate(2, useSel)
    
    gx = db.getGradient(0, useSel)
    gy = db.getGradient(1, useSel)
    gz = db.getGradient(2, useSel)
    
    if len(gx) <= 0 or len(gy) <= 0 or len(gz) <= 0:
        print("Gradient components must be present")
        return
    
    objects = go.Cone(x=x.flatten(), y=y.flatten(), z=z.flatten(), 
                      u=gx.flatten(), v=gy.flatten(), w=gz.flatten(),
                      colorscale = colorscale, sizemode=sizemode,
                      sizeref = size, 
                      **plot_args)
    
    return objects

def TangentDb(db, useSel=True, colorscale='Blues', sizemode='absolute', size=2, 
               **plot_args): 
    '''
    Represent a set of Tangents contained in a Db
    
    plot_args Arguments passed to Cone()
    '''

    if __invalidFileDimension(db, 3):
        return None
                      
    x = db.getOneCoordinate(0, useSel)
    y = db.getOneCoordinate(1, useSel)
    z = db.getOneCoordinate(2, useSel)
    
    tx = db.getTangent(0, useSel)
    ty = db.getTangent(1, useSel)
    tz = db.getTangent(2, useSel)
    
    x = np.concatenate((x, x)) 
    y = np.concatenate((y, y)) 
    z = np.concatenate((z, z)) 
    
    tx = np.concatenate((tx, -tx)) 
    ty = np.concatenate((ty, -ty)) 
    tz = np.concatenate((tz, -tz)) 
    
    if len(tx) <= 0 or len(ty) <= 0 or len(tz) <= 0:
        print("Tangent components must be present")
        return
    
    objects = go.Cone(x=x.flatten(), y=y.flatten(), z=z.flatten(), 
                      u=tx.flatten(), v=ty.flatten(), w=tz.flatten(),
                      colorscale = colorscale, sizemode=sizemode, anchor='tail',
                      sizeref = size, 
                      **plot_args)
                    
    return objects

def Equator(ndisc = 360, color='black', width=3, dilate=1, 
            **plot_args):
    '''
    Represent the Ecuador on a Sphere
    
    plot_args Arguments passed to Line()
    '''
    
    long = np.arange(0,ndisc+1) * 360. / ndisc
    lat  = np.zeros(ndisc+1)
    
    tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
    line = Line(tab[0,:], tab[1,:], tab[2,:], color=color, width=width,
                **plot_args
                )

    return line

def Meridians(angle=10, ndisc=360, color = 'black', width=1, dilate=1.,
              **plot_args):
    '''
    Represent the Meridians on a Sphere
    
    plot_args Arguments passed to Line()
    '''

    number = (int) (360 / angle)  
    xs = list()
    ys = list()
    zs = list()
    
    for i in range(number):
        lat = (np.arange(0,ndisc+1) - ndisc / 2.) * 180. / ndisc
        long = np.zeros(ndisc+1)
        long.fill(i * angle)
        tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
        xp = tab[0,:]
        yp = tab[1,:]
        zp = tab[2,:]
        xs.extend(list(xp) + [None])
        ys.extend(list(yp) + [None])
        zs.extend(list(zp) + [None])
        
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    line = Line(x=xs, y=ys, z=zs, color=color, width=width,
                **plot_args
                )
    return line

def Parallels(angle = 10, ndisc=360, color='black', width=1, dilate=1.):
    
    number = (int) (180 / angle)  
    xs = list()
    ys = list()
    zs = list()
    
    for i in range(number+1):
        long = np.arange(0,ndisc+1) * 360. / ndisc
        lat  = np.zeros(ndisc+1)
        lat.fill((i - number/2) * angle)
        tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
        xp = tab[0,:]
        yp = tab[1,:]
        zp = tab[2,:]
        xs.extend(list(xp) + [None])
        ys.extend(list(yp) + [None])
        zs.extend(list(zp) + [None])
        
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    line = Line(x=xs, y=ys, z=zs, color=color, width=width)
    
    return line

def Pole(sizeref = 1000, dilate=1.3):
    
    long = np.zeros(1)
    lat = np.ones(1) * 90
    tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
    pole = go.Cone(
        u=[0],v=[0],w=[1],
        x=tab[0,:],y=tab[1,:],z=tab[2,:],
        sizemode="absolute",
        sizeref=sizeref,
        anchor="tip",
        showscale=False)
    return pole

def PolarAxis(color='black', width=3, dilate=1.2):
    
    long = np.zeros(2)
    lat = np.zeros(2)
    lat[0] = -90.
    lat[1] = 90.
    tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
    
    line = Line(tab[0,:], tab[1,:], tab[2,:], color=color, width=width)

    return line
