import plotly.graph_objects as go
import numpy                as np
import gstlearn             as gl
from numpy import pi


def getCscale():
    cscale = [[0.0, '#313695'],
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
 [1.0, '#a50026']]
    return cscale

def SurfaceOnMesh(mesh, intensity=None, cscale=None, color='lightpink', opacity=0.50):
    
    tab = np.array(mesh.getEmbeddedApexCoordinates())
    meshes = np.array(mesh.getMeshes())-1
    
    if cscale is None:
        cscale = getCscale()
    
    if intensity is None:
        intensity = np.zeros(mesh.getNApices())

    surface = go.Mesh3d(x=tab[0,:], y=tab[1,:], z=tab[2,:], 
                        color=color, colorbar_title='z', intensity=intensity,
                        i=meshes[:,0],j=meshes[:,1],k=meshes[:,2],
                        name='y', colorscale=cscale, showscale=False)

    return surface
    
def Meshing(mesh, color='black', width=1):
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

    meshing = dict(type='scatter3d',x=xs, y=ys, z=zs, 
                   mode='lines',
                   line=dict(color=color, width=width)
                   )
    return meshing
    
def Scatter(x, y, z, mode='lines', color='black', width=1, 
            m_symbol = 'circle', m_color='black', m_line = 'black', m_size=15, m_width=2):
    ''' 
    mode can be 'lines' or 'markers'
    '''
    meshing = dict(type='scatter3d',x=x, y=y, z=z, 
                   mode=mode,
                   marker_symbol=m_symbol,
                   marker_line_color=m_line, marker_color=m_color, marker_line_width=m_width, 
                   marker_size=m_size,
                   line=dict(color=color, width=width)
                   )
    return meshing
    
def ScatterOnSphere(long, lat, mode='lines', color='black', width=1, 
                    m_symbol = 'circle', m_color='black', m_line = 'black', m_size=15, m_width=2,
                    dilate=1):
    tab = np.array(gl.util_convert_longlat(long, lat, dilate))
    meshing = scatter(tab[0,:], tab[1,:], tab[2,:], mode=mode, color=color, width=width,
                      m_symbol=m_symbol, m_color=m_color, m_line=m_line, m_size=m_size, m_width=m_width)
    return meshing

def Line(x, y, z, color='black', width=1):
    line = dict(type='scatter3d',x=x, y=y, z=z, 
                   mode='lines',
                   line=dict(color=color, width=width),
                   showlegend=False
                   )
    return line
    
def LineOnSphere(long, lat, color='black', width=1, dilate=1.):
    tab = np.array(gl.util_convert_longlat(long, lat, dilate))
    line = Line(tab[0,:], tab[1,:], tab[2,:], color=color, width=width)
    return line

def PolygonOnSphere(poly, flagClose=False, color='black', width=1):
    xs = list()
    ys = list()
    zs = list()

    for i in range(poly.getPolySetNumber()):
        a = poly.getX(i)
        b = poly.getY(i)
        tab = np.array(gl.util_convert_longlat(a, b))
        xp = tab[0,:]
        yp = tab[1,:]
        zp = tab[2,:]
        xs.extend(list(xp) + [None])
        ys.extend(list(yp) + [None])
        zs.extend(list(zp) + [None])
    
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    boundaries=dict(type='scatter3d', x=xs, y=ys, z=zs,
               mode='lines', line=dict(color=color, width=width)
              )
    return boundaries

def SliceOnDbGrid3D(grid, name, section=0, rank=0, usesel=False):
    shape = list(grid.getNXs())
    shape.pop(section)
    vect = grid.getSlice("Simu", section, rank, usesel)
    x = np.array(vect[0]).reshape(shape)
    y = np.array(vect[1]).reshape(shape)
    z = np.array(vect[2]).reshape(shape)
    values = np.array(vect[3]).reshape(shape)
    
    slice = go.Surface(x=x, y=y, z=z, surfacecolor=values,
                    coloraxis='coloraxis')
    return slice
   
def Equator(ndisc = 360, color='black', width=3, dilate=1.):
    long = np.arange(0,ndisc+1) * 360. / ndisc
    lat  = np.zeros(ndisc+1)
    line = LineOnSphere(long, lat, color=color, width=width, dilate=dilate)
    return line

def Meridians(angle=10, ndisc=360, color = 'black', width=1, dilate=1.):
    number = (int) (360 / angle)  
    xs = list()
    ys = list()
    zs = list()
    
    for i in range(number):
        lat = (np.arange(0,ndisc+1) - ndisc / 2.) * 180. / ndisc
        long = np.zeros(ndisc+1)
        long.fill(i * angle)
        tab = np.array(gl.util_convert_longlat(long, lat, dilate))
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

def Parallels(angle = 10, ndisc=360, color='black', width=1, dilate=1.):
    number = (int) (180 / angle)  
    xs = list()
    ys = list()
    zs = list()
    
    for i in range(number+1):
        long = np.arange(0,ndisc+1) * 360. / ndisc
        lat  = np.zeros(ndisc+1)
        lat.fill((i - number/2) * angle)
        tab = np.array(gl.util_convert_longlat(long, lat, dilate))
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
    tab = np.array(gl.util_convert_longlat(long, lat, dilate))
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
    tab = np.array(gl.util_convert_longlat(long, lat, dilate))
    
    line = Line(tab[0,:], tab[1,:], tab[2,:], color=color, width=width)

    return line
