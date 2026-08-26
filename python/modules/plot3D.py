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

from matplotlib import axis
import numpy as np
import gstlearn as gl

try:
    import plotly.graph_objects as go
except ModuleNotFoundError as ex:
    msg = (
        "Python dependency 'plotly' not found.\n"
        "To install it alongside gstlearn, please run `pip install gstlearn[plot]'"
    )
    raise ModuleNotFoundError(msg) from ex


def getCscale():
    cscale = [
        [0.0, "#313695"],
        [0.07692307692307693, "#3a67af"],
        [0.15384615384615385, "#5994c5"],
        [0.23076923076923078, "#84bbd8"],
        [0.3076923076923077, "#afdbea"],
        [0.38461538461538464, "#d8eff5"],
        [0.46153846153846156, "#d6ffe1"],
        [0.5384615384615384, "#fef4ac"],
        [0.6153846153846154, "#fed987"],
        [0.6923076923076923, "#fdb264"],
        [0.7692307692307693, "#f78249"],
        [0.8461538461538461, "#e75435"],
        [0.9230769230769231, "#cc2727"],
        [1.0, "#a50026"],
    ]
    return cscale


def __invalidFileDimension(db, ndim):
    if db.getNDim() != ndim:
        print("This representation is only designed for Db of dimension", ndim)
        return True
    return False


def __linearInterpolate(values, Min=0.0, Max=1.0, flagNoInterpolate=False):
    if flagNoInterpolate:
        return values
    values = np.asarray(values)
    vmin = np.min(values)
    vmax = np.max(values)
    if vmax == vmin:
        return np.full_like(values, (Min + Max) / 2.0)
    # Linear interpolation between Min and Max
    return Min + (values - vmin) * (Max - Min) / (vmax - vmin)


def SurfaceOnMesh(
    mesh, intensity=None, cscale=None, color="white", opacity=0.50, **plot_args
):
    """
    Represent a Function defined on a Mesh

    plot_args Arguments passed to Mesh3d()
    """

    tab = np.array(mesh.getEmbeddedCoordinatesPerApex())
    meshes = np.array(mesh.getMeshesAsVVI())

    if cscale is None:
        cscale = getCscale()

    if intensity is None:
        intensity = np.zeros(mesh.getNApices())
    else:
        if type(intensity) == gl.VectorDouble:
            intensity = np.array(intensity.getVector())

    surface = go.Mesh3d(
        x=tab[0, :],
        y=tab[1, :],
        z=tab[2, :],
        color=color,
        colorbar_title="z",
        name="y",
        i=meshes[:, 0],
        j=meshes[:, 1],
        k=meshes[:, 2],
        colorscale=cscale,
        intensity=intensity,
        opacity=opacity,
        **plot_args,
    )

    return surface


def ScatterOnMesh(
    points,
    meshes,
    intensity=None,
    cscale=None,
    color="white",
    opacity=0.50,
    **plot_args,
):
    """
    Represent a Function defined on a Mesh defined by its apices and meshes

    plot_args Arguments passed to Mesh3d()
    """

    napices = meshes.shape[0]
    if cscale is None:
        cscale = getCscale()

    if intensity is None:
        intensity = np.zeros(napices)
    else:
        if type(intensity) == gl.VectorDouble:
            intensity = np.array(intensity.getVector())

    surface = go.Mesh3d(
        x=points[:, 0],
        y=points[:, 1],
        z=points[:, 2],
        color=color,
        colorbar_title="z",
        name="y",
        i=meshes[:, 0],
        j=meshes[:, 1],
        k=meshes[:, 2],
        colorscale=cscale,
        intensity=intensity,
        opacity=opacity,
        **plot_args,
    )

    return surface


def Meshing(mesh, color="black", width=1, **plot_args):
    """
    Represent the contents of a Mesh object

    plot_args Arguments passed to dict(type='scatter3d')
    """

    xs = list()
    ys = list()
    zs = list()
    for i in range(mesh.getNMeshes()):
        a = np.array(mesh.getEmbeddedCoordinatesPerMesh(i))
        xs.extend(a[:, 0].tolist() + [None])
        ys.extend(a[:, 1].tolist() + [None])
        zs.extend(a[:, 2].tolist() + [None])

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    meshing = dict(
        type="scatter3d",
        x=xs,
        y=ys,
        z=zs,
        mode="lines",
        line=dict(color=color, width=width),
        **plot_args,
    )
    return meshing


def ScatterOnDb(
    db,
    mode="lines",
    color="black",
    width=1,
    m_symbol="circle",
    m_color="black",
    m_line="black",
    m_size=15,
    m_width=2,
    showscale=False,
    **plot_args,
):
    if __invalidFileDimension(db, 3):
        return None

    meshing = dict(
        type="scatter3d",
        x=db.getVecCoordinate(0),
        y=db.getVecCoordinate(1),
        z=db.getVecCoordinate(2),
        mode=mode,
        marker_symbol=m_symbol,
        marker_line_color=m_line,
        marker_color=m_color,
        marker_line_width=m_width,
        marker_size=m_size,
        line=dict(color=color, width=width),
        showscale=showscale,
        **plot_args,
    )
    return meshing


def Scatter(
    x,
    y,
    z,
    mode="lines",
    color="black",
    width=1,
    m_symbol="circle",
    m_color="black",
    m_line="black",
    m_size=15,
    m_width=2,
    **plot_args,
):
    """
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
    """
    meshing = dict(
        type="scatter3d",
        x=x,
        y=y,
        z=z,
        mode=mode,
        marker_symbol=m_symbol,
        marker_line_color=m_line,
        marker_color=m_color,
        marker_line_width=m_width,
        marker_size=m_size,
        line=dict(color=color, width=width),
        **plot_args,
    )
    return meshing


def ScatterOnSphere(
    long,
    lat,
    mode="lines",
    color="black",
    width=1,
    m_symbol="circle",
    m_color="black",
    m_line="black",
    m_size=15,
    m_width=2,
    dilate=1,
    **plot_args,
):
    """
    Represent a set of points on the Sphere

    plot_args Arguments passed to Scatter()
    """
    tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
    meshing = Scatter(
        tab[0, :],
        tab[1, :],
        tab[2, :],
        mode=mode,
        color=color,
        width=width,
        m_symbol=m_symbol,
        m_color=m_color,
        m_line=m_line,
        m_size=m_size,
        m_width=m_width,
        **plot_args,
    )

    return meshing


def Line(x, y, z, color="black", width=1, **plot_args):
    """
    Represent a set of points koined by a Line

    plot_args Arguments passed to dict(type='scatter3d')
    """
    line = dict(
        type="scatter3d",
        x=x,
        y=y,
        z=z,
        mode="lines",
        line=dict(color=color, width=width),
        **plot_args,
    )
    return line


def PolygonOnSphere(
    poly, flagClose=False, color="black", width=1, dilate=1, **plot_args
):
    """
    Represent a Polygon object on the Sphere

    plot_args Arguments passed to dict(type='scatter3d')
    """

    xs = list()
    ys = list()
    zs = list()

    for i in range(poly.getNPolyElem()):
        a = poly.getX(i)
        b = poly.getY(i)
        tab = np.array(gl.GH.convertLongLatTo3D(a, b, dilate, np.nan))
        xp = tab[0, :]
        yp = tab[1, :]
        zp = tab[2, :]
        xs.extend(list(xp) + [None])
        ys.extend(list(yp) + [None])
        zs.extend(list(zp) + [None])

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    boundaries = dict(
        type="scatter3d",
        x=xs,
        y=ys,
        z=zs,
        mode="lines",
        line=dict(color=color, width=width),
        **plot_args,
    )
    return boundaries


def SliceOnDbGrid(grid, name, section=0, rank=0, useSel=False, cmin=None, cmax=None):
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

    slice = go.Surface(
        x=x, y=y, z=z, surfacecolor=values, coloraxis="coloraxis", cmin=cmin, cmax=cmax
    )
    return slice


def Slice3DOnDbGrid(grid, name, corner=None, useSel=False, cmin=None, cmax=None):
    """
    Represent a series of three slices (XoY, YoZ, XoZ)
    performed around the corner (given by its grid indices)

    Returns a go.Figure element
    """
    if __invalidFileDimension(grid, 3):
        return None
    if corner is None:
        corner = grid.getNXs() / 2

    data = [
        SliceOnDbGrid(grid, name, 0, corner[2]),
        SliceOnDbGrid(grid, name, 1, corner[1]),
        SliceOnDbGrid(grid, name, 2, corner[0]),
    ]
    fig = go.Figure(data=data)
    return fig


def IsoSurfaceOnDbGrid(
    grid,
    name,
    useSel=False,
    levels=None,
    colorscale="BlueRed",
    isomin=0,
    isomax=1,
    surface_count=1,
    showlegend=False,
):
    if __invalidFileDimension(grid, 3):
        return None

    shape = list(grid.getNXs())

    x = grid.getVecCoordinate(0, useSel).reshape(shape)
    y = grid.getVecCoordinate(1, useSel).reshape(shape)
    z = grid.getVecCoordinate(2, useSel).reshape(shape)
    values = grid.getColumn(name, useSel).reshape(shape)

    surfaces = go.Isosurface(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value=values.flatten(),
        isomin=isomin,
        isomax=isomax,
        surface_count=surface_count,
        colorscale=colorscale,
        showscale=showlegend,
        caps=dict(x_show=False, y_show=False),
    )
    return surfaces


def SurfaceOnDbGrid(
    grid, name, useSel=False, showscale=False, opacity=0.9, **plot_args
):
    if __invalidFileDimension(grid, 2):
        return None

    shape = list(np.flip(grid.getNXs()))
    values = grid.getColumn(name, useSel).reshape(shape)

    surface = go.Surface(z=values, showscale=showscale, opacity=opacity, **plot_args)

    return surface


def PointDb(
    db,
    nameColor=None,
    nameSize=None,
    useSel=True,
    color="black",
    size=3,
    sizeMin=1,
    sizeMax=3,
    nColors=10,
    flagNoColorInterpolate=False,
    opacity=1,
    posX=0,
    posY=1,
    posZ=2,
    fromLongLat=False,
    dilate=1,
    **plot_args,
):
    """
    Represent a set of Points contained in a Db

    plot_args Arguments passed to Scatter3d()
    """
    if fromLongLat:
        if __invalidFileDimension(db, 2):
            return None
        long = db.getVecCoordinate(0, useSel)
        lat = db.getVecCoordinate(1, useSel)
        tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
        x = tab[0, :]
        y = tab[1, :]
        z = tab[2, :]
    else:
        if __invalidFileDimension(db, 3):
            return None
        x = db.getVecCoordinate(posX, useSel)
        y = db.getVecCoordinate(posY, useSel)
        z = db.getVecCoordinate(posZ, useSel)

    if nameColor is not None:
        colors = __linearInterpolate(
            db.getColumn(nameColor, useSel), 1, nColors, flagNoColorInterpolate
        )
    else:
        colors = color

    if nameSize is not None:
        sizes = __linearInterpolate(db.getColumn(nameSize, useSel), sizeMin, sizeMax)
    else:
        sizes = size

    object = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode="markers",
        marker=dict(size=sizes, color=colors, colorscale="Viridis", opacity=opacity),
        **plot_args,
    )
    return object


def GradientDb(
    db, useSel=True, colorscale="Blues", sizemode="absolute", size=2, **plot_args
):
    """
    Represent a set of Gradients contained in a Db

    plot_args Arguments passed to Cone()
    """

    if __invalidFileDimension(db, 3):
        return None

    x = db.getVecCoordinate(0, useSel)
    y = db.getVecCoordinate(1, useSel)
    z = db.getVecCoordinate(2, useSel)

    gx = db.getGradient(0, useSel)
    gy = db.getGradient(1, useSel)
    gz = db.getGradient(2, useSel)

    if len(gx) <= 0 or len(gy) <= 0 or len(gz) <= 0:
        print("Gradient components must be present")
        return

    objects = go.Cone(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        u=gx.flatten(),
        v=gy.flatten(),
        w=gz.flatten(),
        colorscale=colorscale,
        sizemode=sizemode,
        sizeref=size,
        **plot_args,
    )
    return objects


def TangentDb(
    db, useSel=True, colorscale="Blues", sizemode="absolute", size=2, **plot_args
):
    """
    Represent a set of Tangents contained in a Db

    plot_args Arguments passed to Cone()
    """

    if __invalidFileDimension(db, 3):
        return None

    x = db.getVecCoordinate(0, useSel)
    y = db.getVecCoordinate(1, useSel)
    z = db.getVecCoordinate(2, useSel)

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

    objects = go.Cone(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        u=tx.flatten(),
        v=ty.flatten(),
        w=tz.flatten(),
        colorscale=colorscale,
        sizemode=sizemode,
        anchor="tail",
        sizeref=size,
        **plot_args,
    )

    return objects


def Equator(ndisc=360, color="black", width=3, dilate=1, **plot_args):
    """
    Represent the Ecuador on a Sphere

    plot_args Arguments passed to Line()
    """

    long = np.arange(0, ndisc + 1) * 360.0 / ndisc
    lat = np.zeros(ndisc + 1)

    tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
    line = Line(tab[0, :], tab[1, :], tab[2, :], color=color, width=width, **plot_args)

    return line


def Meridians(angle=10, ndisc=360, color="black", width=1, dilate=1.0, **plot_args):
    """
    Represent the Meridians on a Sphere

    plot_args Arguments passed to Line()
    """

    number = (int)(360 / angle)
    xs = list()
    ys = list()
    zs = list()

    for i in range(number):
        lat = (np.arange(0, ndisc + 1) - ndisc / 2.0) * 180.0 / ndisc
        long = np.zeros(ndisc + 1)
        long.fill(i * angle)
        tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
        xp = tab[0, :]
        yp = tab[1, :]
        zp = tab[2, :]
        xs.extend(list(xp) + [None])
        ys.extend(list(yp) + [None])
        zs.extend(list(zp) + [None])

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    line = Line(x=xs, y=ys, z=zs, color=color, width=width, **plot_args)
    return line


def Parallels(angle=10, ndisc=360, color="black", width=1, dilate=1.0):
    number = (int)(180 / angle)
    xs = list()
    ys = list()
    zs = list()

    for i in range(number + 1):
        long = np.arange(0, ndisc + 1) * 360.0 / ndisc
        lat = np.zeros(ndisc + 1)
        lat.fill((i - number / 2) * angle)
        tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
        xp = tab[0, :]
        yp = tab[1, :]
        zp = tab[2, :]
        xs.extend(list(xp) + [None])
        ys.extend(list(yp) + [None])
        zs.extend(list(zp) + [None])

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    line = Line(x=xs, y=ys, z=zs, color=color, width=width)

    return line


def Pole(sizeref=1000, dilate=1.3):
    long = np.zeros(1)
    lat = np.ones(1) * 90
    tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))
    pole = go.Cone(
        u=[0],
        v=[0],
        w=[1],
        x=tab[0, :],
        y=tab[1, :],
        z=tab[2, :],
        sizemode="absolute",
        sizeref=sizeref,
        anchor="tip",
        showscale=False,
    )
    return pole


def PolarAxis(color="black", width=3, dilate=1.2):
    long = np.zeros(2)
    lat = np.zeros(2)
    lat[0] = -90.0
    lat[1] = 90.0
    tab = np.array(gl.GH.convertLongLatTo3D(long, lat, dilate, np.nan))

    line = Line(tab[0, :], tab[1, :], tab[2, :], color=color, width=width)

    return line


def Sphere(
    radius, center, color="red", opacity=0.5, ntheta=180, showscale=False, **plot_args
):
    if radius <= 0:
        return None

    # Center
    x0, y0, z0 = center

    # angles
    phi = np.linspace(0, np.pi, ntheta)  # latitude
    theta = np.linspace(0, 2 * np.pi, ntheta)  # longitude

    phi, theta = np.meshgrid(phi, theta)

    # coordonnées cartésiennes
    x = radius * np.sin(phi) * np.cos(theta) + x0
    y = radius * np.sin(phi) * np.sin(theta) + y0
    z = radius * np.cos(phi) + z0

    # créer la surface
    sphere = go.Surface(
        x=x,
        y=y,
        z=z,
        opacity=opacity,
        showscale=showscale,
        surfacecolor=np.ones_like(z),
        colorscale=[[0, color], [1, color]],
        **plot_args,
    )

    return sphere


def Cone(
    apex,
    axis,
    angle_deg,
    orientation=1,
    length=5,
    color="orange",
    opacity=0.5,
    n_theta=180,
    **plot_args,
):
    if angle_deg == gl.TEST or angle_deg >= 90 or angle_deg <= 0:
        return None

    x0, y0, z0 = apex
    axis_vec = np.array(axis)
    if orientation == -1:
        axis_vec = -axis_vec

    # normaliser l’axe
    axis_vec /= np.linalg.norm(axis_vec)

    # vecteurs perpendiculaires
    if np.allclose(axis_vec[:2], 0):
        v = np.array([1, 0, 0])
    else:
        v = np.cross(axis_vec, [0, 0, 1])
        v /= np.linalg.norm(v)
    u = np.cross(axis_vec, v)

    # angle en radians
    alpha = np.deg2rad(angle_deg)

    # cercle de base
    r = length * np.tan(alpha)  # rayon à l’extrémité du cône
    theta = np.linspace(0, 2 * np.pi, n_theta)
    circle = (
        apex
        + axis_vec * length
        + r * (np.outer(np.cos(theta), u) + np.outer(np.sin(theta), v))
    )

    # points du cône (triangles vers le sommet)
    x = np.hstack([circle[:, 0], np.full(n_theta, x0)])
    y = np.hstack([circle[:, 1], np.full(n_theta, y0)])
    z = np.hstack([circle[:, 2], np.full(n_theta, z0)])

    # Mesh3d simple
    cone = go.Mesh3d(
        x=x, y=y, z=z, color=color, opacity=opacity, showscale=False, **plot_args
    )

    return cone


def Cylinder(
    center,
    axis,
    radius=1,
    length=5,
    color="green",
    opacity=0.5,
    n_theta=180,
    n_height=20,
    **plot_args,
):
    if radius == gl.TEST or radius <= 0.0:
        return None

    center = np.array(center, dtype=float)
    axis_vec = np.array(axis, dtype=float)
    axis_vec /= np.linalg.norm(axis_vec)

    # vecteurs orthogonaux
    if abs(axis_vec[0]) < 0.1 and abs(axis_vec[1]) < 0.1:
        arbitrary = np.array([0, 1, 0])
    else:
        arbitrary = np.array([0, 0, 1])

    u = np.cross(axis_vec, arbitrary)
    u /= np.linalg.norm(u)
    v = np.cross(axis_vec, u)  # vecteur orthogonal complet

    # coordonnées paramétriques
    theta = np.linspace(0, 2 * np.pi, n_theta)
    t = np.linspace(-length / 2, length / 2, n_height)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    X = np.zeros((n_height, n_theta))
    Y = np.zeros((n_height, n_theta))
    Z = np.zeros((n_height, n_theta))

    for i, ti in enumerate(t):
        circle = (
            np.array(center)
            + axis_vec * ti
            + radius * (np.outer(cos_theta, u) + np.outer(sin_theta, v))
        )
        X[i, :] = circle[:, 0]
        Y[i, :] = circle[:, 1]
        Z[i, :] = circle[:, 2]

    Cylindre = go.Mesh3d(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        color=color,
        opacity=opacity,
        showscale=False,
        **plot_args,
    )

    return Cylindre


def Bench(
    center,
    bench,
    orientation=1,
    length=5,
    color="lightblue",
    opacity=0.5,
    n1=10,
    n2=10,
    **plot_args,
):
    if bench == gl.TEST or bench <= 0.0:
        return None

    codir = (0, 0, 1)
    center = np.array(center, dtype=float)
    codir = np.array(codir, dtype=float)
    codir /= np.linalg.norm(codir)  # vecteur unitaire perpendiculaire aux plans

    # calcul des centres des deux plans
    plane_center = center + (bench if orientation == 1 else -bench) * codir

    # vecteurs directeurs pour les plans (orthogonaux à codir)
    if np.allclose(codir[:2], 0):
        arbitrary = np.array([1, 0, 0])
    else:
        arbitrary = np.array([0, 0, 1])

    dir1 = np.cross(codir, arbitrary)
    dir1 /= np.linalg.norm(dir1)
    dir2 = np.cross(codir, dir1)

    # fonction interne pour créer un plan
    s = np.linspace(-length / 2, length / 2, n1)
    t = np.linspace(-length / 2, length / 2, n2)
    S, T = np.meshgrid(s, t)
    X = plane_center[0] + S * dir1[0] + T * dir2[0]
    Y = plane_center[1] + S * dir1[1] + T * dir2[1]
    Z = plane_center[2] + S * dir1[2] + T * dir2[2]
    Plane = go.Mesh3d(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        color=color,
        opacity=opacity,
        showscale=False,
        **plot_args,
    )

    return Plane
