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

# TODO : 1. geopandas datasets are soon deprecated
#        2. geopandas dependency should be optionnal (issue #68)
#            => Do not restore without reading this: https://github.com/gstlearn/gstlearn/issues/68


try:
    import geopandas as gpd
    from shapely.geometry import Polygon, Point
except ModuleNotFoundError as ex:
    msg = ("Python dependencies 'geopandas' and 'shapely' not found.\n"
          "To install them alongside gstlearn, please run `pip install gstlearn[plot]'")
    raise ModuleNotFoundError(msg) from ex

''' Generate GeometryArray of shapely Point geometries from x, y(, z) coordinates.
    x, y: Data coordinates
    crsFrom: Coordinate Reference System of the geometry objects. Can be anything accepted 
    by pyproj.CRS.from_user_input(), such as an authority string (eg “EPSG:4326”) 
    It returns a GeometryArray.
    '''
    
def proj(x, y, crsFrom="EPSG:4326", crsTo="EPSG:2154"):
    
    # Creating a new Geometry from coordinates using input projection
    a = gpd.points_from_xy(x, y, crs=crsFrom)

    # Conversion into the output projection system
    a = a.to_crs(crsTo)

    return a

def projGP(geometry, crsFrom="EPSG:4326", crsTo="EPSG:3857"):
    '''
    Create a GeoPanda DataFrame in order to perform a projection between two CRS
    '''
    gdf = gpd.GeoDataFrame(geometry=geometry, crs=crsFrom) 
    gdf = gdf.to_crs(crsTo) 
    return gdf

'''
    Define the world contour line, possibly intersectd by a Bounding Box
    minx, miny, maxx, maxy: A set of 4 tuple defining the bounding box
    Example: minx=(-10,36), miny=(-1,36), maxx=(-1,49), maxy=(-10,49)
    csrTo: Coordinate Reference System in which the country lines must be projected
'''
def world(minx=None, miny=None, maxx=None, maxy=None, crsTo="EPSG:2154"):

    # Reading the world contour lines
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    if minx is not None:
        # Creating the Bounding Box
        polyGpd = Polygon([minx, miny, maxx, maxy])

        # Intersecting The Country map and the Bounding Box
        world['geometry'] = world.geometry.intersection(polyGpd)

    # Projection
    world = world.to_crs(crsTo) 
    
    return world
    