#!/usr/bin/env python
# encoding: utf-8
"""
file: tiftools.py

   Provides several useful functions
   for converting geotif format to tt3 format.

Contains:
   get_topo:     downloads topo file from GeoClaw repository on web.
   topo1writer:  create files of topotype1 from synthetic topo function.
   topo2writer:  create files of topotype2 from synthetic topo function.
   gcdist:       computes great circle distance between points on sphere.
   dx_from_gcdist:  inverts gcdist function at a give latitude.

   tif2tt3
   rasterio2tt3

Authors: Katy Barnhart
"""

import os

import numpy as np

def tif2tt3(inpath, outpath, band=1, strformat=r"%7.5e", nanerror=True):
    """
    Convert a tif file to a tt3 format file.

    Depends on rasterio to read in the tif. Raises a value error if rasterio
    package is not installed.

    inpath : string indicating path of input tif
    outpath : string indicating path of output .tt3
    band : band from tif to write out (presumes that band exists)
    strformat : string format to write in tt3 (default "%7.5e")
    nanerror : if True (default) raise an error if there are any nodata values
        in the input tif band.
    """
    try:
        import rasterio
    except:
        raise ImportError("rasterio dependency not met.")
    with rasterio.open(inpath) as src:
        rasterio2tt3(src.read(band), src.transform, src.meta["nodata"], outpath, strformat=strformat, nanerror=nanerror)


def rasterio2tt3(array, transform, nodata, outpath, strformat=r"%7.5e", nanerror=True):
    """
    Convert formats used by rasterio to a tt3 format file.

    Depends on rasterio to read in the tif. Raises a value error if rasterio
    package is not installed.

    array : numpy array representing 2d data read in by rasterio.
    transform : rasterio transfrom object representing spatial extent and grid
        cell size.
    nodata : nodata value to check against.
    outpath : string indicating path of output .tt3
    strformat : string format to write in tt3 (default "%7.5e")
    nanerror : if True (default) raise an error if there are any nodata values
        in the input tif band.
    """
    try:
        import rasterio.transform
    except:
        raise ImportError("rasterio dependency not met.")
    # write out in tt3 format.
    if np.any(np.isnan(array)) and nanerror:
        raise ValueError("array cannot have nan, to supress this message set nanerror=False")
    nrows = array.shape[0]
    ncols = array.shape[1]
    west, south, _, _ = rasterio.transform.array_bounds(nrows, ncols, transform)
    dx = transform[0]
    xllcenter = west + (dx / 2)
    yllcenter = south + (dx / 2)
    with open(outpath, "w") as f:

        f.write("%6i ncols\n" % ncols)
        f.write("%6i nrows\n" % nrows)
        f.write("%22.15e xllcenter\n" % xllcenter)
        f.write("%22.15e yllcenter\n" % yllcenter)
        f.write("%22.15e cellsize\n" % dx,)
        f.write("%22.15e nodata_value\n" % nodata)
        outstr = strformat + "   "
        for i in range(nrows):
            for j in range(ncols):
                f.write(outstr % (array[i, j]))
            f.write("\n")
