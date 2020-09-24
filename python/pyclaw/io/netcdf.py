#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing a NetCDF output file

Routines for reading and writing a NetCDF output file via either
    - netcdf4-python - http://code.google.com/p/netcdf4-python/
    - pupynere - http://pypi.python.org/pypi/pupynere/
    
These interfaces are very similar so if a different module needs to be used,
it can more than likely be inserted with a minimal of effort.

This module will first try to import the netcdf4-python module which is based
on the compiled libraries and failing that will attempt to import the pure
python interface pupynere which requires no libraries.

To install the netCDF 4 library, please see:
    http://www.unidata.ucar.edu/software/netcdf/
    
:Authors:
    Kyle T. Mandli (2009-02-17) Initial version
    Josh Jacobs (2011-04-22) NetCDF 3 Support
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#      Copyright (C) 2011 J. Jacobs
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD)
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================
import os, sys
import logging

import pyclaw.solution

logger = logging.getLogger("io")

# Import appropriate netcdf package
use_netcdf3 = (
    True  # use this implementation as default...even though it does not have groups
)
use_netcdf4 = False
use_pupynere = False
try:
    from Scientific.IO import NetCDF

    use_netcdf3 = True
except:
    try:
        from scipy.io import netcdf as NetCDF
    except:
        print("*** Could not import NetCDF from Scientific.IO or scipy.io")
if not use_netcdf3:
    try:
        import netCDF4

        use_netcdf4 = True
    except:
        pass
if not (use_netcdf3 or use_netcdf4):
    try:
        import pupynere

        use_pupynere = True
    except:
        error_msg = (
            "Could not import netCDF3,netCDF4 or Pupynere, please install "
            + "one of the available modules for netcdf files.  Refer to this "
            + "modules doc_string for more information."
        )
        # raise Exception(error_msg)
        print(error_msg)


def write_netcdf(
    solution, frame, path, file_prefix="fort", write_aux=False, options={}
):
    r"""
    Write out a NetCDF data file representation of solution
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw object to be 
       output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name. ``default = 'claw'``
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out. ``default = False``     
     - *options* - (dict) Optional argument dictionary, see 
       `NetCDF Option Table`_
    
    .. _`NetCDF Option Table`:
    
    +-------------------------+----------------------------------------------+
    | Key                     | Value                                        |
    +=========================+==============================================+
    | description             | Dictionary of key/value pairs that will be   |
    |                         | attached to the root group as attributes,    |
    |                         | i.e. {'time':3}                              |
    +-------------------------+----------------------------------------------+
    | format                  | Can be one of the following netCDF flavors:  |
    |                         | NETCDF3_CLASSIC, NETCDF3_64BIT,              |
    |                         | NETCDF4_CLASSIC, and NETCDF4                 |
    |                         | ``default = NETCDF4``                        |
    +-------------------------+----------------------------------------------+
    | clobber                 | if True (Default), file will be overwritten, |
    |                         | if False an exception will be raised         |
    +-------------------------+----------------------------------------------+
    | zlib                    | if True, data assigned to the Variable       |
    |                         | instance is compressed on disk.              |
    |                         | ``default = False``                          |
    +-------------------------+----------------------------------------------+
    | complevel               | the level of zlib compression to use (1 is   |
    |                         | the fastest, but poorest compression, 9 is   |
    |                         | the slowest but best compression).  Ignored  |
    |                         | if zlib=False.  ``default = 6``              |
    +-------------------------+----------------------------------------------+
    | shuffle                 | if True, the HDF5 shuffle filter is applied  |
    |                         | to improve compression. Ignored if           |
    |                         | zlib=False. ``default = True``               |
    +-------------------------+----------------------------------------------+
    | fletcher32              | if True (default False), the Fletcher32      |
    |                         | checksum algorithm is used for error         |
    |                         | detection.                                   |
    +-------------------------+----------------------------------------------+
    | contiguous              | if True (default False), the variable data   |
    |                         | is stored contiguously on disk.  Setting to  |
    |                         | True for a variable with an unlimited        |
    |                         | dimension will trigger an error.             |
    |                         | ``default = False``                          |
    +-------------------------+----------------------------------------------+
    | chunksizes              | Can be used to specify the HDF5 chunksizes   |
    |                         | for each dimension of the variable. A        |
    |                         | detailed discussion of HDF chunking and I/O  |
    |                         | performance is available here. Basically,    |
    |                         | you want the chunk size for each dimension   |
    |                         | to match as closely as possible the size of  |
    |                         | the data block that users will read from the |
    |                         | file. chunksizes cannot be set if            |
    |                         | contiguous=True.                             |
    +-------------------------+----------------------------------------------+
    | least_significant_digit | If specified, variable data will be          |
    |                         | truncated (quantized). In conjunction with   |
    |                         | zlib=True this produces 'lossy', but         |
    |                         | significantly more efficient compression.    |
    |                         | For example, if least_significant_digit=1,   |
    |                         | data will be quantized using around          |
    |                         | (scale*data)/scale, where scale = 2**bits,   |
    |                         | and bits is determined so that a precision   |
    |                         | of 0.1 is retained (in this case bits=4).    |
    |                         | ``default = None``, or no quantization.      |
    +-------------------------+----------------------------------------------+
    | endian                  | Can be used to control whether the data is   |
    |                         | stored in little or big endian format on     | 
    |                         | disk. Possible values are little, big or     |
    |                         | native (default). The library will           |
    |                         | automatically handle endian conversions when |
    |                         | the data is read, but if the data is always  |
    |                         | going to be read on a computer with the      |
    |                         | opposite format as the one used to create    |
    |                         | the file, there may be some performance      |
    |                         | advantage to be gained by setting the        |
    |                         | endian-ness.                                 |
    +-------------------------+----------------------------------------------+
    | fill_value              | If specified, the default netCDF _FillValue  |
    |                         | (the value that the variable gets filled     |
    |                         | with before any data is written to it) is    |
    |                         | replaced with this value. If fill_value is   |
    |                         | set to False, then the variable is not       |
    |                         | pre-filled.                                  |
    +-------------------------+----------------------------------------------+
    
    .. note:: 
        The zlib, complevel, shuffle, fletcher32, contiguous, chunksizes and
        endian keywords are silently ignored for netCDF 3 files that do not 
        use HDF5.
        
    """

    # Option parsing
    option_defaults = {
        "format": "NETCDF4",
        "zlib": False,
        "complevel": 6,
        "shuffle": True,
        "fletcher32": False,
        "contiguous": False,
        "chunksizes": None,
        "endian": "native",
        "least_significant_digit": None,
        "fill_value": None,
        "clobber": True,
        "description": {},
    }
    for (k, v) in list(option_defaults.items()):
        if k in options:
            exec("%s = options['%s']" % (k, k))
        else:
            exec("%s = v" % k)

    # Filename
    filename = os.path.join(path, "%s.q%s.nc" % (file_prefix, str(frame).zfill(4)))
    if use_netcdf3:
        import numpy

        # Open new file
        f = NetCDF.NetCDFFile(
            filename, "w"
        )  # if I want not to clobber this mode needs to be 'a'

        # Loop through description dictionary and add the attributes to the
        # root group
        for (k, v) in list(description.items()):
            exec("f.%s = %s" % (k, v))

        # Create Global Dimensions
        f.createDimension("timedimension", None)  # NoneType is unlimited
        f.createDimension("meqn", solution.grid.meqn)

        # Create Global variables
        time = f.createVariable("timedimension", "f", ("timedimension",))
        ngrids = f.createVariable("ngrids", "d", ())
        naux = f.createVariable("naux", "d", ())
        ndim = f.createVariable("ndim", "d", ())

        # For each grid, write out attributes
        for grid in solution.grids:
            if solution.grids.index(grid) == 0:
                time.assignValue(getattr(grid, "t"))
                ngrids.assignValue(
                    len(solution.grids)
                )  # I am not quite sure about these
                naux.assignValue(grid.naux)  #
                ndim.assignValue(len(grid.dimensions))  #

            # Create dimensions for q (and aux)
            dim_names = []
            for dim in grid.dimensions:
                f.createDimension(dim.name + str(grid.gridno), dim.n)
                dim_names.append(dim.name)

            # Write q array
            tmp_names = dim_names + ["meqn", "timedimension"]
            index_str = ",".join([":" for name in tmp_names])
            q = f.createVariable("grid_" + str(grid.gridno), "d", tmp_names)
            exec("q[%s] = grid.q" % index_str)

            # General grid properties
            for attr in ["gridno", "level"]:  # ,'mbc']:
                setattr(q, attr, getattr(grid, attr))
            setattr(q, "dim_names", dim_names)
            for dim in grid.dimensions:
                setattr(q, dim.name + ".lower", dim.lower)
                setattr(q, dim.name + ".d", dim.d)

        #            # Write out aux
        #            if grid.maux > 0 and write_aux:
        #                dim_names[-1] = 'maux'
        #                subgroup.createDimension('maux',grid.maux)
        #                aux = subgroup.createVariable('aux','f8',dim_names,
        #                                            zlib,complevel,shuffle,fletcher32,
        #                                            contiguous,chunksizes,endian,
        #                                            least_significant_digit,fill_value)
        #                exec("aux[%s] = grid.aux" % index_str)

        f.close()
    elif use_netcdf4:
        # Open new file
        f = netCDF4.Dataset(filename, "w", clobber=clobber, format=format)

        # Loop through description dictionary and add the attributes to the
        # root group
        for (k, v) in list(description.items()):
            exec("f.%s = %s" % (k, v))

        # For each grid, write out attributes
        for grid in solution.grids:
            # Create group for this grid
            subgroup = f.createGroup("grid%s" % grid.gridno)

            # General grid properties
            for attr in ["t", "meqn", "gridno", "level", "mbc"]:
                setattr(subgroup, attr, getattr(grid, attr))

            # Write out dimension names
            setattr(subgroup, "dim_names", grid.name)

            # Create dimensions for q (and aux)
            for dim in grid.dimensions:
                subgroup.createDimension(dim.name, dim.n)
                # Write other dimension attributes
                for attr in [
                    "n",
                    "lower",
                    "d",
                    "upper",
                    "mthbc_lower",
                    "mthbc_upper",
                    "units",
                ]:
                    if hasattr(dim, attr):
                        if getattr(dim, attr) is not None:
                            attr_name = "%s.%s" % (dim.name, attr)
                            setattr(subgroup, attr_name, getattr(dim, attr))
            subgroup.createDimension("meqn", grid.meqn)

            # Write q array
            dim_names = grid.name
            dim_names.append("meqn")
            index_str = ",".join([":" for name in dim_names])
            q = subgroup.createVariable(
                "q",
                "f8",
                dim_names,
                zlib,
                complevel,
                shuffle,
                fletcher32,
                contiguous,
                chunksizes,
                endian,
                least_significant_digit,
                fill_value,
            )
            exec("q[%s] = grid.q" % index_str)

            # Write out aux
            if grid.maux > 0 and write_aux:
                dim_names[-1] = "maux"
                subgroup.createDimension("maux", grid.maux)
                aux = subgroup.createVariable(
                    "aux",
                    "f8",
                    dim_names,
                    zlib,
                    complevel,
                    shuffle,
                    fletcher32,
                    contiguous,
                    chunksizes,
                    endian,
                    least_significant_digit,
                    fill_value,
                )
                exec("aux[%s] = grid.aux" % index_str)

        f.close()
    elif use_pupynere:
        logging.critical("Pupynere support has not been implemented yet.")
        raise IOError("Pupynere support has not been implemented yet.")
    else:
        err_msg = "No netcdf python modules available."
        logging.critical(err_msg)
        raise Exception(err_msg)


def read_netcdf(
    solution, frame, path="./", file_prefix="fort", read_aux=True, options={}
):
    r"""
    Read in a NetCDF data files into solution
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw object to be 
       output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name.  ``default = 'claw'``
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out.  ``default = False``     
     - *options* - (dict) Optional argument dictionary, unused for reading.
    """

    # Option parsing
    option_defaults = {}
    for (k, v) in list(option_defaults.items()):
        if k in options:
            exec("%s = options['%s']" % (k, k))
        else:
            exec("%s = v" % k)

    # Filename
    filename = os.path.join(path, "%s.q%s.nc" % (file_prefix, str(frame).zfill(4)))
    # jj-2010.02.04-added support for netcdf classic format
    if use_netcdf3:
        from Scientific.IO import NetCDF
        import numpy

        # Open file
        print(filename)
        f = NetCDF.NetCDFFile(filename, "r")
        # Each file written by the fortran code has
        # Dimensions:
        #           timedimension : UNLIMITED
        #           meqn          : The number of equations
        #           dimx_<gridno> : X dimension for grid number <gridno>
        #           dimy_<gridno> : Y dimension for grid number <gridno>
        # Variables:
        #           timedimension : Stores the time of the frame
        #           ngrids        : Number of grids in this frame
        #           naux          : Number of Auxilary Variables
        #           ndim          : Number of Dimensions in the frame
        #           grid_<gridno> : A grid of (dimx,dimy,meqn)
        # Attributes:
        # (grid_<no>) gridno      : The number of this grid <grid_no>
        #           level         : The AMR level
        #           dim_names     : a list of dimensions [dimx,dimy]
        #           dim<>.low     : The lowest dimension value
        #           dim<>.d       : The distance between grid points

        time = numpy.squeeze(f.variables["timedimension"].getValue())
        ngrids = numpy.squeeze(f.variables["ngrids"].getValue())
        naux = numpy.squeeze(f.variables["naux"].getValue())
        ndims = numpy.squeeze(f.variables["ndim"].getValue())
        meqn = f.dimensions["meqn"]
        #        print time,ngrids,naux,ndims,meqn
        for var_name in list(f.variables.keys()):
            if var_name[:4] == "grid":
                var = f.variables[var_name]
                gridno = numpy.squeeze(getattr(var, "gridno"))
                # print gridno
                # Construct each dimension
                dimensions = []

                # Read in dimension attribute to keep dimension order
                dim_names = eval(getattr(var, "dim_names"))
                for dim_name in dim_names:
                    dim_n = f.dimensions[dim_name + "_" + str(gridno).zfill(2)]
                    dim_lower = numpy.squeeze(getattr(var, "%s.lower" % dim_name))
                    dim_d = numpy.squeeze(getattr(var, "%s.d" % dim_name))
                    # print dim_lower,dim_d
                    dim_upper = dim_lower + dim_d * dim_n
                    dim = pyclaw.solution.Dimension(
                        dim_name, dim_lower, dim_upper, dim_n
                    )

                    #                 # Optional attributes
                    #                 # jj-2011/02/04--these are boundary conditions... don't know how to handle yet.
                    #                for attr in ['mthbc_lower','mthbc_upper','units']:
                    #                    attr_name = "%s.%s" % (dim_name,attr)
                    #                    if hasattr(subgroup,attr_name):
                    #                        setattr(dim,attr,getattr(subgroup, "%s.%s" % (dim_name,attr)))

                    dimensions.append(dim)

                # Create grid
                grid = pyclaw.solution.Grid(dimensions)

                # General grid properties
                for attr in ["t", "meqn", "gridno", "level"]:
                    if attr in ["gridno", "level"]:
                        setattr(grid, attr, numpy.squeeze(getattr(var, attr)))
                    elif attr == "t":
                        setattr(grid, attr, time)
                    elif attr == "meqn":
                        setattr(grid, attr, meqn)
                # Read in q
                # index_str = ','.join( [':' for i in xrange(grid.ndim+1)] )
                # grid_str='grid_'+str(gridno).zfill(2)
                # exec("tmp_grid = f.variables['"+grid_str+"'][0,%s]" % index_str)
                tmp_grid = numpy.array(var.getValue())[0]  # .variables[grid_str][0,:]
                # print type(tmp_grid),numpy.shape(tmp_grid)
                # swap axes from column major to row major order
                tmp_len = len(tmp_grid.shape)
                for i in range(0, tmp_len / 2):
                    tmp_grid = tmp_grid.swapaxes(i, tmp_len - (i + 1))

                grid.q = tmp_grid
                # Read in aux if applicable
                # jj--2011.02.04--Not going to do this yet.
                ##            if read_aux and subgroup.dimensions.has_key('maux'):
                ##                exec("grid.aux = subgroup.variables['aux_<grid_no>'][%s]" % index_str)

                solution.grids.append(grid)
        solution.grids.sort(key=lambda grid: grid.level)
        f.close()
    elif use_netcdf4:
        # Open file
        f = netCDF4.Dataset(filename, "r")

        # We only expect subgroups of grids, otherwise we need to put some
        # sort of conditional here
        for subgroup in list(f.groups.values()):
            # Construct each dimension
            dimensions = []

            # Read in dimension attribute to keep dimension order
            dim_names = getattr(subgroup, "dim_names")
            for dim_name in dim_names:
                dim = pyclaw.solution.Dimension(
                    dim_name,
                    getattr(subgroup, "%s.lower" % dim_name),
                    getattr(subgroup, "%s.upper" % dim_name),
                    getattr(subgroup, "%s.n" % dim_name),
                )
                # Optional attributes
                for attr in ["mthbc_lower", "mthbc_upper", "units"]:
                    attr_name = "%s.%s" % (dim_name, attr)
                    if hasattr(subgroup, attr_name):
                        setattr(
                            dim, attr, getattr(subgroup, "%s.%s" % (dim_name, attr))
                        )
                dimensions.append(dim)

            # Create grid
            grid = pyclaw.solution.Grid(dimensions)

            # General grid properties
            for attr in ["t", "meqn", "gridno", "level"]:
                setattr(grid, attr, getattr(subgroup, attr))

            # Read in q
            index_str = ",".join([":" for i in range(grid.ndim + 1)])
            exec("grid.q = subgroup.variables['q'][%s]" % index_str)

            # Read in aux if applicable
            if read_aux and "maux" in subgroup.dimensions:
                exec("grid.aux = subgroup.variables['aux'][%s]" % index_str)

            solution.grids.append(grid)
        solution.grids.sort(key=lambda grid: grid.level)

        f.close()
    elif use_pupynere:
        logging.critical("Pupynere support has not been implemented yet.")
        raise IOError("Pupynere support has not been implemented yet.")
    else:
        err_msg = "No netcdf python modules available."
        logging.critical(err_msg)
        raise Exception(err_msg)


def read_netcdf_t(frame, path="./", file_prefix="fort"):
    # Option parsing
    option_defaults = {}
    for (k, v) in list(option_defaults.items()):
        if k in options:
            exec("%s = options['%s']" % (k, k))
        else:
            exec("%s = v" % k)

    # Filename
    filename = os.path.join(path, "%s.q%s.nc" % (file_prefix, str(frame).zfill(4)))
    #    print filename
    if use_netcdf3:
        from Scientific.IO import NetCDF
        import numpy

        # Open file
        f = NetCDF.NetCDFFile(filename, "r")
        t = numpy.squeeze(f.variables["timedimension"].getValue())
        print(t)
        return t
