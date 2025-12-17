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

try:
    import pyproj as ppr
except ModuleNotFoundError as ex:
    msg = (
        "Python dependencies 'pyproj' not found.\n"
        "To install it alongside gstlearn, please run `pip install gstlearn[plot]'"
    )
    raise ModuleNotFoundError(msg) from ex

""" Generate projected x', y' coordinates from x, y coordinates.
    x, y: Data coordinates
    crsFrom: Coordinate Reference System of the input coordinates. Can be anything accepted
    by pyproj package, such as an authority string (eg “EPSG:4326”) 
    crsTo: Coordinate Reference System of the output coordinates. Can be anything accepted
    by pyproj package, such as an authority string (eg “EPSG:2154”)
    Returns: Pair of vector of projected coordinates
    """


def proj(x, y, crsFrom="EPSG:4326", crsTo="EPSG:2154"):
    a = ppr.Transformer.from_crs(crsFrom, crsTo, always_xy=True).transform(x, y)
    return a
