#!/usr/bin/python
"""
file: fortconvert
=========
   Provides several functions
   for manipulating output data in AMR files, fort.qXXXX
   into other formats of grids with non-overlapping data (scattered, uniform, topotype)

   Contains:

       convertfortdir (convert entire directory or files from directory)
       fort2topotype    (create file in the form of topotype for a single component of q or b.)
       fort2xyqscatter (finest grid data for entire domain as x,y,q scattered data)
       fort2refined  (covert fort.qXXXX with amr into a single uniform grid with highest level spacing)
       fort2uniform  (covert fort.qXXXX with amr into a single uniform grid with user defined spacing)
       fort2griddata
       fort2griddata_vector
       fort2griddata_framenumbers
       fort2list (used to read fort data...4.x compatible...)
       forttheaderread
       fortqheaderread
       forttheaderwrite
       fortqheaderwrite
       pointfromfort
       intersection
       griddata2fort
       array2fort


    Future plans:
        probably should use the pyclaw solution class rather than solutionlist
        or at least provide options for both.

        provide conversion to clawpack netcdf and possibly other binary formats.

    Note:
        this was done expeditiously and might be inefficient or bug infested
        (KRB 2022.06.14): x,y locations considered here represent CELL CENTERS
        rather than lower left corner of cells.

    David George dgeorge@uw.edu 2014.05.20

"""

import copy
import os
import string

import numpy as np

import dclaw.topotools as gt
#import dclaw.netcdf_tools as gn


# ================================================================================
def convertfortdir(
    outputtype,
    nplots="fort.nplot",
    outputname="fort.q",
    components="all",
    outdir=None,
    fortdir=None,
    parallel=False,
    num_cores=1,
    **kwargs
):

    """
    convert an entire directory of fort.q files to another form

    arguments:
    ---------
        outputtype: 'scattered','topotype' , 'fortrefined' , 'fortuniform'
                    Note: only scattered option retains actual non-interpolated output, generally.
                     --scattered: highest level data as columns x,y,q
                     --topotype: like a DEM with a single component of q (useful for GIS software)
                     --fortrefined: standard clawpack output style, but a single grid at highest resolution
                     --fortuniform: standard clawpack output style, but a single grid with user defined grid parameters.
        nplots: integer for the number of plots, list of plots or fort.nplot filename
            a sinlge integer will be all plots up to that integer
            for individual plots make nplots a list
        outputname: the base name of the output files...filenames will be appended with frameno
        components: for scattered and fort a list of components q, or 'all'
            For topotype it must be an integer to select the single component of q, or 'topo' to get eta-h,
                or 'depth' to get nodata_value where h=0...for plotting purposes.
        outdir: specify a directory for files to go in relative to cwd
        fortdir: location of fort.q files if not ./
        parallel: bool flag whether to use parallel processing (needs joblib, default false)
        num_cores: num cores to use if doing parallel (default 1)
        topotype: topotype to output to. If provided when using "fortrefined" or "fortuniform",
                    output will get written to the specified topotype instead of to the standard
                    fort.qXXX file type. A fort.tXXXX will also get written. The geotiff will have
                    the name fort_qXXXX.tif
        **kwargs: for topotype:  xll,yll,cellsize,ncols,nrows,topotype = 1,2,3,'gdal', 'gtif' (for standard gdal/esri header)
                for fortuniform: xlower,ylower,mx,my
                if using topotype="gtif" can provide epsg=XXXX (EPSG code for CRS)
                if kwargs are omitted the grid parameters are taken from the fort file with finest level spacing


    """

    if not fortdir:
        fortdir = "."

    if not outdir:
        outdir = "."
    elif not os.path.isdir(outdir):
        os.system(("mkdir " + outdir))

    curdir = os.path.abspath(os.path.curdir)
    fortdir = os.path.abspath(fortdir)
    outdir = os.path.abspath(outdir)

    write_level = kwargs.get("write_level", False)
    east = kwargs.get("east", None)
    west = kwargs.get("west", None)
    south = kwargs.get("south", None)
    north = kwargs.get("north", None)

    bilinear = kwargs.get("bilinear", True)

    if east is not None and west is not None:
        assert east > west
    if north is not None and south is not None:
        assert north > south

    epsg = kwargs.get("epsg", None)

    if outputtype == "topotype":
        try:
            xll = kwargs["xll"]
            yll = kwargs["yll"]
            nrows = kwargs["nrows"]
            ncols = kwargs["ncols"]
            cellsize = kwargs["cellsize"]
        except:
            print("using grid parameters from fort.qXXXX files")
            xll = None
            yll = None
            nrows = None
            ncols = None
            cellsize = None

        try:
            topotype = kwargs["topotype"]
        except:
            topotype = 2

    if outputtype == "fortuniform":
        try:
            xlower = kwargs["xlower"]
            ylower = kwargs["ylower"]
            xupper = kwargs["xupper"]
            yupper = kwargs["yupper"]
            mx = kwargs["mx"]
            my = kwargs["my"]

        except:
            print(
                "for outputtype==fortuniform you must provide xlower,ylower,xupper,yupper,mx,my as kwargs"
            )
            raise ValueError()

    if isinstance(nplots, str):
        nplotfile = os.path.join(fortdir, nplots)
        fin = open(nplotfile, "r")
        nplots = fin.readline()
        nplots = int(nplots)
        nplots = np.arange(nplots + 1)
        fin.close()
    elif isinstance(nplots, int):
        nplots = np.arange(nplots + 1)
        print(("converting frames 0 - %s" % (nplots)))
        print("to convert individual frames, call with nplots = a list of frames")
    else:
        nplots = np.array(nplots, dtype=int)

    os.chdir(fortdir)

    arg_list = []

    for frameno in nplots:
        numstring = str(10000 + frameno)
        framenostr = numstring[1:]
        forttname = "fort.t" + framenostr
        fortqname = "fort.q" + framenostr
        # print(('converting '+os.path.join(fortdir,fortqname)))
        outfname = os.path.join(outdir, outputname + framenostr)

        # print(('writing to ' +outfname))

        if outputtype == "scattered":
            _func = fort2xyqscattered
            arg_list.append([frameno, outfname, components])

        elif outputtype == "topotype":
            _func = fort2topotype
            arg_list.append(
                [
                    frameno,
                    outfname,
                    fortdir,
                    xll,
                    yll,
                    cellsize,
                    ncols,
                    nrows,
                    components,
                    topotype,
                    bilinear
                ]
            )

        elif outputtype == "fortrefined":
            outfortt = os.path.join(outdir, "fort.t" + framenostr)
            topotype = kwargs.get("topotype", None)
            _func = fort2refined

            arg_list.append(
                [
                    frameno,
                    outfname,
                    outfortt,
                    components,
                    topotype,
                    write_level,
                    west,
                    east,
                    south,
                    north,
                    epsg,
                    bilinear,
                ]
            )

        elif outputtype == "fortuniform":
            outfortt = os.path.join(outdir, "fort.t" + framenostr)
            topotype = kwargs.get("topotype", None)
            _func = fort2uniform

            arg_list.append(
                [
                    frameno,
                    outfname,
                    outfortt,
                    xlower,
                    xupper,
                    ylower,
                    yupper,
                    mx,
                    my,
                    components,
                    topotype,
                    write_level,
                    epsg,
                    bilinear,
                ]
            )

    # now run in parallel based on func and arg list
    if parallel:
        try:
            from joblib import Parallel, delayed
        except ImportError:
            raise ImportError("joblib needed for parallel functionality")
        Parallel(n_jobs=num_cores)(delayed(_func)(*args) for args in arg_list)
    else:
        # loop through prepared args.
        for args in arg_list:
            print(("Converting {}".format(args[1])))
            _func(*args)

    # return to curdir if changed.
    os.chdir(curdir)


# =============================================================================
def fort2xyqscattered(framenumber, outfile=None, components="all"):
    """
    convert a fort.qXXXX amr file into a scattered data with columns x,y,q...
    q can be 1-meqn columns according to the list 'components'
    data is taken from the finest of grid intersections

    data are returned with x,y coordinates reflecting the cell center.

    arguments
    ----------
        framenumber : of fort.qXXXX
        outfile: name or file handle or None
            if outfile=None routine returns a numpy array
        components: list of q components eg. [1,3,5] or 'all' for all components
    """

    if isinstance(outfile, str):
        fout = open(outfile, "w")
    else:
        fout = outfile

    numstring = str(10000 + framenumber)
    framenostr = numstring[1:]
    forttname = "fort.t" + framenostr
    fortqname = "fort.q" + framenostr

    solutionlist = fort2list(fortqname, forttname)

    if components == "all":
        qlst = np.arange(solutionlist[0]["meqn"])
    else:
        qlst = np.array(components, dtype=int) - 1

    levels = solutionlist[0]["AMR_maxlevel"]

    # note that solutionlist is ordered from highest levels to lowest.
    for grid in solutionlist:
        if grid["AMR_level"] == levels:  # highest level...data assumed nonoverlapping
            x = grid["xlow"] + grid["dx"] * (0.5 + (np.arange(grid["mx"], dtype=float))) # converted from lower left to cell center
            y = grid["ylow"] + grid["dy"] * (0.5 + (np.arange(grid["my"], dtype=float)))
            Q = grid["data"][:, qlst]
            (X, Y) = np.meshgrid(x, y)
            X = np.reshape(X, (grid["mx"] * grid["my"], 1))
            Y = np.reshape(Y, (grid["mx"] * grid["my"], 1))

            try:
                XYQ = np.vstack((XYQ, np.hstack((X, Y, Q))))
            except:
                XYQ = np.hstack((X, Y, Q))
            gridrect = [grid["xlow"], grid["ylow"], grid["xupper"], grid["yupper"]]
        else:  # not highest level
            # check for intersection
            gridrect = [grid["xlow"], grid["ylow"], grid["xupper"], grid["yupper"]]
            gridoverlap = False
            for rectangle in rectangles:
                if intersection(gridrect, rectangle):
                    gridoverlap = True
                    break
            if not gridoverlap:  # take all data
                x = grid["xlow"] + grid["dx"] * (
                    0.5 + (np.arange(grid["mx"], dtype=float))
                )
                y = grid["ylow"] + grid["dy"] * (
                    0.5 + (np.arange(grid["my"], dtype=float))
                )
                Q = grid["data"][:, qlst]
                (X, Y) = np.meshgrid(x, y)
                X = np.reshape(X, (grid["mx"] * grid["my"], 1))
                Y = np.reshape(Y, (grid["mx"] * grid["my"], 1))
                XYQ = np.vstack((XYQ, np.hstack((X, Y, Q))))
            else:  # need to loop to find points
                row = 0
                for j in range(grid["my"]):
                    y = grid["ylow"] + grid["dy"] * (0.5 + float(j))
                    for i in range(grid["mx"]):
                        x = grid["xlow"] + grid["dx"] * (0.5 + float(i))
                        q = grid["data"][row, qlst]
                        # check point
                        gridoverlap = False
                        for rectangle in rectangles:
                            if intersection([x, y, x, y], rectangle):
                                gridoverlap = True
                                break
                        if not gridoverlap:
                            XYQ = np.vstack((XYQ, np.hstack((x, y, q))))
                        row = row + 1
        try:
            rectangles = np.vstack((rectangles, gridrect))
        except:
            rectangles = gridrect

    if not outfile:
        return XYQ
    else:
        np.savetxt(fout, XYQ)
        fout.close()


# ==============================================================================
def fort2uniform(
    framenumber,
    outfortq,
    outfortt,
    xlow,
    xhi,
    ylow,
    yhi,
    mx,
    my,
    components="all",
    topotype=None,
    write_level=False,
    epsg=None,
    bilinear=True,
):
    """
    convert fort.qXXXX with AMR data into fort.qXXXX with data on a uniform single grid.
    Resolution is user defined.
    Format is still standard clawpack fort.qXXXX /fort.tXXXX and can be plotted with clawpack utilities
    Should call with outdir being a new directory to keep original fort.q/fort.t files. Names are the same

    x,y of results are provided as cell centers.

    arguments
    ----------
        framenumber : of fort.qXXXX
        outfortq: name of output fort.qXXXX file
        outfortt: name of output fort.tXXXX file
            if above =None routine returns a numpy array
        xlower,xupper,ylower,yupper,mx,my: output grid parameters
        components: list of q components eg. [1,3,5] or 'all' for all components
        xlow
        xhi
        ylow
        yhi
        mx
        my
        components = "all"
        topotype
        write_level
        epsg
    """

    numstring = str(10000 + framenumber)
    framenostr = numstring[1:]
    forttname = "fort.t" + framenostr
    fortqname = "fort.q" + framenostr

    solutionlist = fort2list(fortqname, forttname)

    if components == "all":
        qlst = np.arange(solutionlist[0]["meqn"])
    else:
        qlst = np.array(components, dtype=int) - 1

    dx = float((xhi - xlow) / mx)
    dy = float((yhi - ylow) / my)

    # print(dx, dy, mx, my, xlow, xhi, ylow, yhi )
    fortheader = {}
    fortheader["grid_number"] = 1
    fortheader["AMR_level"] = 1
    fortheader["mx"] = mx
    fortheader["my"] = my
    fortheader["xlow"] = xlow
    fortheader["ylow"] = ylow
    fortheader["dx"] = dx
    fortheader["dy"] = dy
    fortheader["naux"] = solutionlist[0]["naux"]
    fortheader["time"] = solutionlist[0]["time"]
    fortheader["ndim"] = solutionlist[0]["ndim"]
    fortheader["ngrids"] = 1
    if components == "all":
        fortheader["meqn"] = solutionlist[0]["meqn"]
    else:
        fortheader["meqn"] = len(components)

    if isinstance(outfortt, str):
        foutt = open(outfortt, "w")
    else:
        foutf = outfortt

    if topotype is not None:
        # write t file.
        forttheaderwrite(fortheader, foutt)
        foutt.close()

        # prepare to return an array, OR write to topo file.
        Q = np.empty((mx * my, len(qlst)))
        if write_level:
            source_level = np.empty((my, mx))

        for j in range(my):
            y = ylow + (j + 0.5) * dy # here x and y are cell centers based on 0.5 dx and dy adjustement.
            for i in range(mx):
                x = xlow + (i + 0.5) * dx
                qv, lev = pointfromfort((x, y), solutionlist, bilinear=bilinear)
                qout = qv[qlst]
                Q[j * mx + i] = qout

                if write_level:
                    source_level[j, i] = lev

        # if topotype is specified, write out as topotype instead of
        # array or standard fort.q
        if topotype is not None:

            xv = np.array(xlow + dx * np.arange(mx))
            yv = np.array(ylow + dy * np.arange(my))
            (X, Y) = np.meshgrid(xv, yv)
            Y = np.flipud(Y)

            if topotype == "gtif":
                outfile = outfortq.replace("fortq.", "fortq_") + ".tif"
                # this should only change the file name.

                # manipulate shape.order of Q
                # written row by row, so shape into my, mx, meqn
                # reorder into my, mx, meqn by shift axis so meq is at front,
                # finally flip along axis 1 so that up is up.

                Q_out = np.flip(
                    np.moveaxis(Q.reshape((my, mx, len(qlst))), (0, 1, 2), (1, 2, 0)),
                    axis=1,
                )
                if write_level:
                    source_level = np.flipud(source_level)
                    Q_out = np.concatenate(
                        (np.atleast_3d(Q_out), source_level.reshape((1, my, mx))),
                        axis=0,
                    )

                gt.griddata2gtif(
                    X,
                    Y,
                    Q_out,
                    outfile,
                    epsg=epsg,
                )

            elif topotype == "netcdf":
                import dclaw.netcdf_tools as gn
                outfile = outfortq.replace("fortq.", "fortq_") + ".nc"
                # this should only change the file name.

                # manipulate shape.order of Q
                # written row by row, so shape into my, mx, meqn
                # reorder into my, mx, meqn by shift axis so meq is at front,
                # finally flip along axis 1 so that up is up.

                Q_out = np.flip(
                    np.moveaxis(Q.reshape((my, mx, len(qlst))), (0, 1, 2), (1, 2, 0)),
                    axis=1,
                )
                if write_level:
                    source_level = np.flipud(source_level)
                    Q_out = np.concatenate(
                        (np.atleast_3d(Q_out), source_level.reshape((1, my, mx))),
                        axis=0,
                    )

                gn.griddata2netcdf(
                    fortheader["time"],
                    X,
                    Y,
                    Q_out,
                    outfile,
                    qlst,
                    write_level,
                    epsg=None
                )

            else:
                if fortheader["meqn"] > 1:
                    raise ValueError(
                        "refined/uniform to topo only can take 1 element, for all use gtif"
                    )

                outfile = outfortq
                if topotype == "gdal":
                    gt.griddata2topofile(
                        X, Y, Q, ".tmpfile", nodata_value_out=nodata_value
                    )
                    infile = ".tmpfile"
                    gt.esriheader(infile, outfile)
                    os.system("rm .tmpfile")
                else:
                    gt.griddata2topofile(
                        X, Y, Q, outfile, topotype, nodata_value_out=nodata_value
                    )
        else:
            return fortheader, Q

    else:
        if isinstance(outfortq, str):
            foutq = open(outfortq, "w")
        else:
            foutq = outfortq

        forttheaderwrite(fortheader, foutt)
        foutt.close()
        fortqheaderwrite(fortheader, foutq, closefile=False)

        for j in range(my):
            foutq.write("\n")
            y = ylow + (j + 0.5) * dy # here also we are extracting at cell centers.
            for i in range(mx):
                x = xlow + (i + 0.5) * dx
                qv, lev = pointfromfort((x, y), solutionlist, bilinear=bilinear)
                qout = qv[qlst]
                for q in qout:
                    foutq.write("%s " % float(q))
                foutq.write("\n")
        foutq.close()


# ==============================================================================
def fort2refined(
    framenumber,
    outfortq,
    outfortt,
    components="all",
    topotype=None,
    write_level=False,
    west=None,
    east=None,
    south=None,
    north=None,
    epsg=None,
    bilinear=True,
):
    """
    convert fort.qXXXX with AMR data into fort.qXXXX with data on a uniform single grid.
    Resolution is at that of the highest level grids.
    Format is still standard clawpack fort.qXXXX /fort.tXXXX and can be plotted with clawpack utilities
    Should call with outdir being a new directory to keep original fort.q/fort.t files. Names are the same

    future plans:
                  this routine could be more efficient by directly assigning the highest level data
                  for simplicity, now it just loops point by point.
    arguments
    ----------
        framenumber : of fort.qXXXX
        outfortq: name of output fort.qXXXX file
        outfortt: name of output fort.tXXXX file
            if =None routine returns a numpy array
        components: list of q components eg. [1,3,5] or 'all' for all components
        topotype:
        write_level: bool, to write the level from which values are taken.
        west=None,
        east=None,
        south=None,
        north=None,
        epsg=None,
    """
    numstring = str(10000 + framenumber)
    framenostr = numstring[1:]
    forttname = "fort.t" + framenostr
    fortqname = "fort.q" + framenostr

    solutionlist = fort2list(fortqname, forttname)

    if components == "all":
        qlst = np.arange(solutionlist[0]["meqn"])
    else:
        qlst = np.array(components, dtype=int) - 1

    levels = solutionlist[0]["AMR_maxlevel"]
    xlow = solutionlist[0]["xlowdomain"]
    ylow = solutionlist[0]["ylowdomain"]
    xhi = solutionlist[0]["xhidomain"]
    yhi = solutionlist[0]["yhidomain"]
    dx = solutionlist[0]["dx"]
    dy = solutionlist[0]["dy"]

    mx = int((xhi - xlow) / dx)
    my = int((yhi - ylow) / dy)
    # if any of nsew is not none crop extent.
    if east is not None or west is not None:
        xs = np.linspace(xlow, xhi, mx)
        if east is not None:
            xhi = np.max(xs[xs < east])
        if west is not None:
            xlow = np.min(xs[xs > west])
        mx = int((xhi - xlow) / dx)
    if north is not None or south is not None:
        ys = np.linspace(ylow, yhi, my)
        if north is not None:
            yhi = np.max(ys[ys < north])
        if south is not None:
            ylow = np.min(ys[ys > south])
        my = int((yhi - ylow) / dy)

    # fort2uniform will extract at cell centers based on the bounding rectangle defined.
    return fort2uniform(
        framenumber,
        outfortq,
        outfortt,
        xlow,
        xhi,
        ylow,
        yhi,
        mx,
        my,
        components,
        topotype,
        write_level,
        epsg,
    )


# ==============================================================================
def fort2topotype(
    framenumber, outfile, fortdir, xll, yll, cellsize, ncols, nrows, m=1, topotype=2, bilinear=True,
):
    """
    convert data in a fort file of framenumber = XXXX, ie fort.qXXXX
    to a DEM style file of topotype=1,2 or 3, or gdal (gdal/esri header) with uniform spacing.
    m is the component of q, ie. the column in the fort.qXXXX file
    for all components of q, see fort2xyscattered or amr2single

    specifying xll, yll, cellsize, ncols, nrows will result in interpolation of AMR data

    without xll, yll, cellsize, ncols, nrows...the domain will be the same as input frame
    and coarse resolution.
    """

    numstring = str(10000 + framenumber)
    framenostr = numstring[1:]
    forttname = os.path.join(fortdir, "fort.t" + framenostr)
    fortqname = os.path.join(fortdir, "fort.q" + framenostr)

    nodata_value = -9999.0

    solutionlist = fort2list(fortqname, forttname)

    if not xll:
        dx = solutionlist[0]["dx"]
        dy = solutionlist[0]["dy"]
        xll = solutionlist[0]["xlowdomain"]
        yll = solutionlist[0]["ylowdomain"]
        xhi = solutionlist[0]["xhidomain"]
        yhi = solutionlist[0]["yhidomain"]

        cellsize = min(dx, dy)
        ncols = int(np.floor((xhi - xll) / cellsize))
        nrows = int(np.floor((yhi - yll) / cellsize))

        xll = xll + 0.5 * cellsize # Here we adjust to get cell centers.
        yll = yll + 0.5 * cellsize
        xhi = xhi - 0.5 * cellsize
        yhi = yhi - 0.5 * cellsize

    meqn = solutionlist[0]["meqn"]

    xv = np.array(xll + cellsize * np.arange(ncols))
    yv = np.array(yll + cellsize * np.arange(nrows))

    (X, Y) = np.meshgrid(xv, yv)

    Y = np.flipud(Y)

    Q = np.zeros(np.shape(X))

    if m == "topo":
        for j in range(ncols):
            xp = X[0, j]
            for i in range(nrows):
                yp = Y[i, 0]
                qv, lev = pointfromfort((xp, yp), solutionlist, bilinear=bilinear)
                Q[i, j] = qv[meqn - 1] - qv[0]

    elif m == "depth":
        for j in range(ncols):
            xp = X[0, j]
            for i in range(nrows):
                yp = Y[i, 0]
                qv, lev = pointfromfort((xp, yp), solutionlist, bilinear=bilinear)
                depth = qv[0]
                if depth <= 1.0e-3:
                    depth = nodata_value
                Q[i, j] = depth

    elif m == "eta":
        for j in range(ncols):
            xp = X[0, j]
            for i in range(nrows):
                yp = Y[i, 0]
                qv, lev = pointfromfort((xp, yp), solutionlist, bilinear=bilinear)
                eta = qv[meqn - 1]
                if qv[0] <= 1.0e-3:
                    eta = nodata_value
                Q[i, j] = eta

    elif m == "all":
        Q = np.zeros((ncols * nrows, meqn))
        for i in range(nrows):
            yp = Y[i, 0]
            for j in range(ncols):
                xp = X[0, j]
                k = i * ncols + j
                qv, lev = pointfromfort((xp, yp), solutionlist, bilinear=bilinear)
                Q[k] = qv

    else:
        for j in range(ncols):
            xp = X[0, j]
            for i in range(nrows):
                yp = Y[i, 0]
                qv, lev = pointfromfort((xp, yp), solutionlist, bilinear=bilinear)
                Q[i, j] = qv[m - 1]

    if m == "all":
        headerstring = """


                    """
        np.savetxt(outfile, Q)  # ,header=headerstring) save for new numpy version
    elif topotype == "gdal":
        gt.griddata2topofile(X, Y, Q, ".tmpfile", nodata_value_out=nodata_value)
        infile = ".tmpfile"
        gt.esriheader(infile, outfile)
        os.system("rm .tmpfile")
    else:
        gt.griddata2topofile(X, Y, Q, outfile, topotype, nodata_value_out=nodata_value)


# ==============================================================================
def griddata2fort(
    X,
    Y,
    Q,
    qoutputfile,
    toutputfile,
    time=0.0,
    naux=0,
    ndim=2,
    grid_number=1,
    AMR_level=1,
):
    """
    !!!!WIP -- NOT TESTED!!!
    griddata2fort creates a clawpack fort file from grid datasets
    Q - numpy array, shape(Q)=(mx,my,meqn)
    """

    Qshape = np.shape(Q)
    if len(Qshape) == 2:
        my = len(Q[:, 0])
        mx = len(Q[0, :])
        meqn = 1
        Qfort = np.reshape(Q, (mx * my, 1))
    elif len(Qshape) == 3:
        my = len(Q[:, 0, 0])
        mx = len(Q[0, :, 0])
        meqn = Qshape[2]
        Qfort = np.reshape(Q, (mx * my, 1))
    xlow = X[0, 0] # FLAG DAVE: This probably needs adjusting for cell center vs lower left.
    ylow = Y[-1, 0]
    dx = X[0, 1] - X[0, 0]
    dy = Y[0, 0] - Y[1, 0]

    array2fort(
        Qfort,
        outfortqfile,
        outforttfile,
        xlow,
        ylow,
        mx,
        my,
        dx,
        dy,
        time,
        naux,
        ndim,
        grid_number,
        AMR_level,
    )

    return


# ==============================================================================
def array2fort(
    Q,
    fort_q_outfile,
    fort_t_outfile,
    xlow,
    ylow,
    mx,
    my,
    dx,
    dy,
    time=0.0,
    naux=0,
    ndim=2,
    grid_number=1,
    AMR_level=1,
):
    """
    array2fortfile creates a clawpack fort file (single grid) for a single array, Q.
    Assumed that shape(Q) = (mx*my,meqn); ie: array Q is not reshaped.
    """

    meqn = np.shape(Q)[1] + 1

    # creat header for toutputfile
    forttheader = forttheader = {}
    forttheader["time"] = time
    forttheader["meqn"] = meqn
    forttheader["ngrids"] = 1
    forttheader["naux"] = naux
    forttheader["ndim"] = ndim
    forttheaderwrite(forttheader, fort_t_outfile, closefile=True)

    fortqheader = {}
    fortqheader["xlow"] = xlow
    fortqheader["ylow"] = ylow
    fortqheader["dx"] = dx
    fortqheader["dy"] = dy
    fortqheader["mx"] = mx
    fortqheader["my"] = my
    fortqheader["grid_number"] = 1
    fortqheader["AMR_level"] = 1

    foutq = fortqheaderwrite(fortqheader, fort_q_outfile, closefile=False)
    for j in range(my):
        foutq.write("\n")
        for i in range(mx):
            rowind = j * mx + i
            qout = Q[rowind, :]
            for q in qout:
                foutq.write("%s " % float(q))
            foutq.write("\n")
    foutq.close()

    return


# ==============================================================================
def fort2griddata(fortqname, forttname, m=1, bilinear=True):
    """
    convert data in a fort file ie fort.qXXXX
    to numpy arrays X,Y,Q (single gridded data)
    m is the component of q, ie. the column in the fort.qXXXX file
    """

    solutionlist = fort2list(fortqname, forttname)
    fortqheader = fortqheaderread(fortqname)
    xll = fortqheader["xlow"]
    yll = fortqheader["ylow"]
    dx = fortqheader["dx"]
    dy = fortqheader["dy"]
    ncols = fortqheader["mx"]
    nrows = fortqheader["my"]

    xv = np.array(xll + dx * np.arange(ncols)) # FLAG DAVE
    yv = np.array(yll + dy * np.arange(nrows)) # 0.5 DX, DY adjustement?

    (X, Y) = np.meshgrid(xv, yv)

    Y = np.flipud(Y)

    Q = np.zeros(np.shape(X))

    for j in range(ncols):
        xp = X[0, j]
        for i in range(nrows):
            yp = Y[i, 0]
            qv = pointfromfort((xp, yp), solutionlist, bilinear=bilinear)
            Q[i, j] = qv[m - 1]

    return X, Y, Q


# ==============================================================================

# ==============================================================================
def fort2griddata_vector(fortqname, forttname, meqn=7, bilinear=True):
    """
    convert data in a fort file ie fort.qXXXX
    to numpy arrays X,Y,Q (single gridded data)
    m is the component of q, ie. the column in the fort.qXXXX file
    """

    solutionlist = fort2list(fortqname, forttname)
    fortqheader = fortqheaderread(fortqname)
    xll = fortqheader["xlow"]
    yll = fortqheader["ylow"]
    dx = fortqheader["dx"]
    dy = fortqheader["dy"]
    ncols = fortqheader["mx"]
    nrows = fortqheader["my"]

    xv = np.array(xll + dx * np.arange(ncols)) # FLAG DAVE
    yv = np.array(yll + dy * np.arange(nrows)) # 0.5 DX, DY adjustement?

    (X, Y) = np.meshgrid(xv, yv)

    Y = np.flipud(Y)

    Xshape = np.shape(X)
    # meqn = 7
    mq = meqn + 1
    Qshape = (Xshape[0], Xshape[1], mq)
    Q = np.zeros(Qshape)

    for j in range(ncols):
        xp = X[0, j]
        for i in range(nrows):
            yp = Y[i, 0]
            qv = pointfromfort((xp, yp), solutionlist, bilinear=bilinear)
            Q[i, j, :] = qv
            # import pdb;pdb.set_trace()

    return X, Y, Q


# ==============================================================================
def fort2griddata_framenumbers(framenumber, fortdir, m=1, bilinear=True):
    """
    convert data in a fort file of framenumber = XXXX, ie fort.qXXXX
    to numpy arrays X,Y,Q (single gridded data)
    m is the component of q, ie. the column in the fort.qXXXX file
    """

    try:
        import dclaw.topotools as gt
    except:
        import geoclaw.topotools as gt

    if isinstance(framenumber, str):
        numstring = framenumber
    elif isinstance(framenumber, int):
        numstring = str(10000 + framenumber)

    framenostr = numstring[1:]
    forttname = "fort.t" + framenostr
    fortqname = "fort.q" + framenostr

    solutionlist = fort2list(fortqname, forttname)

    xv = np.array(xll + cellsize * np.arange(ncols)) # FLAG DAVE
    yv = np.array(yll + cellsize * np.arange(nrows)) # 0.5 DX, DY adjustement?

    (X, Y) = np.meshgrid(xv, yv)

    Y = np.flipud(Y)

    Q = np.zeros(np.shape(X))

    for j in range(ncols):
        xp = X[0, j]
        for i in range(nrows):
            yp = Y[i, 0]
            qv, lev = pointfromfort((xp, yp), solutionlist, bilinear=bilinear)
            Q[i, j] = qv[m - 1]

    return X, Y, Q


# =========================================================================================
def forttheaderread(inputfile, closefile=True):

    """
    read the header in inputfile (fort.tXXX) and place in dictionary forttheader to be returned.
    The header is of the following form with columns containing the forttheader value and keyword respectively.


        float time
        int   meqn
        int   ngrids
        int   naux
        int   ndim
    """

    if isinstance(inputfile, str):
        fin = open(inputfile, "r")
    else:
        fin = inputfile

    forttheader = {}

    field = fin.readline().split()
    forttheader["time"] = float(field[0])

    field = fin.readline().split()
    forttheader["meqn"] = int(field[0])

    field = fin.readline().split()
    forttheader["ngrids"] = int(field[0])

    field = fin.readline().split()
    forttheader["naux"] = int(field[0])

    field = fin.readline().split()
    forttheader["ndim"] = int(field[0])

    if closefile:
        fin.close()
        return forttheader
    else:
        return fin, forttheader

    # end fortqheader===============================================================================


# ============================================================================
def forttheaderwrite(forttheader, outputfile, closefile=True):

    """
    forttheaderwrite opens an ascii topography data file and writes the header
    using the dictionary "forttheader"

    The header is of the following form with columns containing value and key respectively.

        float time
        int   meqn
        int   ngrids
        int   naux
        int   ndim


    if closefile==True: the file is closed. Otherwise return the open file object.
    """

    if isinstance(outputfile, str):
        fout = open(outputfile, "w")
    else:
        fout = outputfile

    fout.write("%s %s\n" % (float(forttheader["time"]), "time"))
    fout.write("%s %s\n" % (forttheader["meqn"], "meqn"))
    fout.write("%s %s\n" % (int(forttheader["ngrids"]), "ngrids"))
    fout.write("%s %s\n" % (int(forttheader["naux"]), "naux"))
    fout.write("%s %s\n" % (int(forttheader["ndim"]), "ndim"))

    if closefile:
        fout.close()
    else:
        return fout
    # end headerwriter=========================================================================


# =========================================================================================
def fortqheaderread(inputfile, closefile=True):

    """
    read an individual header in inputfile (fort.qXXXX) and place in dictionary fortqheader to be returned.
    The header is of the following form with columns containing the fortqheader value and keyword respectively.

        int grid_number
        int AMR_level
        int mx
        int my
        float xlow
        float ylow
        float dx
        float dy
    """
    if isinstance(inputfile, str):
        fin = open(inputfile, "r")
    else:
        fin = inputfile

    fortqheader = {}

    field = fin.readline().split()
    fortqheader["grid_number"] = int(field[0])

    field = fin.readline().split()
    fortqheader["AMR_level"] = int(field[0])

    field = fin.readline().split()
    fortqheader["mx"] = int(field[0])

    field = fin.readline().split()
    fortqheader["my"] = int(field[0])

    field = fin.readline().split()
    fortqheader["xlow"] = float(field[0])

    field = fin.readline().split()
    fortqheader["ylow"] = float(field[0])

    field = fin.readline().split()
    fortqheader["dx"] = float(field[0])

    field = fin.readline().split()
    fortqheader["dy"] = float(field[0])

    if closefile:
        fin.close()
        return fortqheader
    else:
        return fin, fortqheader

    # end fortqheader===============================================================================


# ============================================================================
def fortqheaderwrite(fortqheader, outfile, closefile=True):

    """
    fortqheaderwrite writes the header in fort.qXXXX files
    using the dictionary "fortqheader"

    outputfile can be a filename or file handle

    The header is of the following form with columns containing value and key respectively.

        int grid_number
        int AMR_level
        int mx
        int my
        float xlow
        float ylow
        float dx
        float dy


    if closefile==True: the file is closed. Otherwise return the open file object.
    """

    if isinstance(outfile, str):
        fout = open(outfile, "w")
    else:
        fout = outfile

    fout.write("%s %s\n" % (int(fortqheader["grid_number"]), "grid_number"))
    fout.write("%s %s\n" % (int(fortqheader["AMR_level"]), "AMR_level"))
    fout.write("%s %s\n" % (int(fortqheader["mx"]), "mx"))
    fout.write("%s %s\n" % (int(fortqheader["my"]), "my"))
    fout.write("%s %s\n" % (float(fortqheader["xlow"]), "xlow"))
    fout.write("%s %s\n" % (float(fortqheader["ylow"]), "ylow"))
    fout.write("%s %s\n" % (float(fortqheader["dx"]), "dx"))
    fout.write("%s %s\n" % (float(fortqheader["dy"]), "dy"))

    if closefile:
        fout.close()
    else:
        return fout
    # end headerwriter=========================================================================


# ================================================================================
def fort2list(fortqname, forttname):
    """
    read fort.qXXXX file and make a list, solutionlist. Each element of the list is
    a dictionary representation of an individual grid (all data and relevant grid params)
    """

    fortt = forttheaderread(forttname)
    fin = open(fortqname, "r")

    ngrids = fortt["ngrids"]
    solutionlist = []

    for gridinlist in range(ngrids):
        fin, fortq = fortqheaderread(fin, closefile=False)
        griddict = {}
        griddict.update(fortt)
        griddict.update(fortq)
        griddict["xupper"] = griddict["xlow"] + griddict["mx"] * griddict["dx"]
        griddict["yupper"] = griddict["ylow"] + griddict["my"] * griddict["dy"]
        rows = griddict["mx"] * griddict["my"]
        grid = np.zeros((rows, griddict["meqn"]))
        row = 0
        for my in range(griddict["my"]):
            fin.readline()
            for mx in range(griddict["mx"]):
                line1 = fin.readline().split()
                # line2 = fin.readline().split()
                grid[row, :] = np.array((line1), dtype=float)
                row = row + 1

        fin.readline()
        # copy might be inefficient...is a deep copy needed? introduced in debugging
        griddict["data"] = copy.copy(grid)
        solutionlist.append(griddict)
    fin.close()

    xlowdomain = np.inf
    ylowdomain = np.inf
    xhidomain = -np.inf
    yhidomain = -np.inf
    maxlevel = 0
    for gridinlist in range(ngrids):
        xlowdomain = min(xlowdomain, solutionlist[gridinlist]["xlow"])
        ylowdomain = min(ylowdomain, solutionlist[gridinlist]["ylow"])
        xhidomain = max(xhidomain, solutionlist[gridinlist]["xupper"])
        yhidomain = max(yhidomain, solutionlist[gridinlist]["yupper"])
        maxlevel = max(maxlevel, solutionlist[gridinlist]["AMR_level"])
    for gridinlist in range(ngrids):
        solutionlist[gridinlist]["xlowdomain"] = xlowdomain
        solutionlist[gridinlist]["ylowdomain"] = ylowdomain
        solutionlist[gridinlist]["xhidomain"] = xhidomain
        solutionlist[gridinlist]["yhidomain"] = yhidomain
        solutionlist[gridinlist]["AMR_maxlevel"] = maxlevel

        orderedlist = sorted(solutionlist, key=lambda k: k["AMR_level"])
        orderedlist.reverse()

    return orderedlist


# ===============================================================================
def pointfromfort(point, solutionlist, bilinear=True):
    """
    for a point (x,y) return the solution vector q determined from the
    best grid available for that point.
    future plans: array from fort takes numpy grid arrays X,Y
    """

    xp = point[0]
    yp = point[1]

    griddict = solutionlist[0]
    dintersection = (
        (xp >= griddict["xlowdomain"])
        & (xp <= griddict["xhidomain"])
        & (yp >= griddict["ylowdomain"])
        & (yp <= griddict["yhidomain"])
    )

    if not dintersection:
        print("ERROR: point outside of domain:")
        print(("point x= %d y=%d" % (point)))
        print(
            (
                "domain x bounds: %d -- %d"
                % (griddict["xlowdomain"], griddict["xhidomain"])
            )
        )
        print(
            (
                "domain y bounds: %d -- %d"
                % (griddict["ylowdomain"], griddict["yhidomain"])
            )
        )
        raise SystemExit(0)

    intersection = False
    grid = {}
    for trygrid in solutionlist:
        intersection = (
            (xp >= trygrid["xlow"])
            & (xp <= trygrid["xupper"])
            & (yp >= trygrid["ylow"])
            & (yp <= trygrid["yupper"])
        )
        if intersection:
            grid = copy.copy(trygrid)
            break

    try:
        xlow = grid["xlow"]
        ylow = grid["ylow"]
        yhi = grid["yupper"]
        xhi = grid["xupper"]
        dx = grid["dx"]
        dy = grid["dy"]
        mx = grid["mx"]
        my = grid["my"]
        data = grid["data"]
        level = grid["AMR_level"]
    except KeyError:

        print(("point is possibly on amr grid edge: x= %s y=%s" % (point)))
        print(("intersection? %s" % (intersection)))
        print("taking data from adjacent grid")
        eps = 1e-5

        for trygrid in solutionlist:
            intersection = (
                (xp + eps >= trygrid["xlow"])
                & (xp - eps <= trygrid["xupper"])
                & (yp + eps >= trygrid["ylow"])
                & (yp - eps <= trygrid["yupper"])
            )

            if intersection:
                grid = copy.copy(trygrid)
                break
        try:
            xlow = grid["xlow"]
            ylow = grid["ylow"]
            yhi = grid["yupper"]
            xhi = grid["xupper"]
            dx = grid["dx"]
            dy = grid["dy"]
            mx = grid["mx"]
            my = grid["my"]
            data = grid["data"]
            level = grid["AMR_level"]
        except KeyError:

            print(("point is possibly on amr grid edge: x= %s y=%s" % (point)))
            domain = (
                griddict["xlowdomain"],
                griddict["xhidomain"],
                griddict["ylowdomain"],
                griddict["yhidomain"],
            )
            print(("Domain xlow,xhi,ylow,yhi: [%s , %s] , [%s , %s]" % (domain)))
            print(("intersection? %s" % (intersection)))
            print("quitting in protest.")

            raise SystemExit(0)

    # actual i,j values for a given grid (note: this is not the row in data array)
    i1 = int(np.floor((xp - xlow) / dx)) + 1
    i2 = int(np.ceil((xp - xlow) / dx)) + 1
    j1 = int(np.floor((yp - ylow) / dy)) + 1
    j2 = int(np.ceil((yp - ylow) / dy)) + 1

    # catch small rounding errors in case point is very close to grid edge
    i1 = max(i1, 1)
    i1 = min(i1, mx)
    i2 = max(i2, 1)
    i2 = min(i2, mx)
    j1 = max(j1, 1)
    j1 = min(j1, my)
    j2 = max(j2, 1)
    j2 = min(j2, my)

    # x and y values of the four surrounding points in the grid
    # KRB note: I think these are the lower left corners of the cell.
    # all need increasing by 1/2 dx or dy. # does this mean that all the
    # topo writers (non-gtif) would then need adjusting... Since they all
    # assume that x and y are the lower left corners of the grid cell in
    # question.
    xl = xlow + (i1 - 1) * dx + (dx/2)
    xr = xl + dx + (dx/2)
    yl = ylow + (j1 - 1) * dy+ (dy/2)
    yu = yl + dy+ (dy/2)

    # indices into the data array
    ijll = (j1 - 1) * mx + i1
    ijlr = (j1 - 1) * mx + i2
    ijul = (j2 - 1) * mx + i1
    ijur = (j2 - 1) * mx + i2

    # solution vector at the four corners surrounding (xp,yp)
    qll = data[ijll - 1, :]
    qlr = data[ijlr - 1, :]
    qul = data[ijul - 1, :]
    qur = data[ijur - 1, :]

    # bilinear interpolation to (xp,yp)
    if bilinear:
        q = (
            qll * (xr - xp) * (yu - yp)
            + qlr * (xp - xl) * (yu - yp)
            + qul * (xr - xp) * (yp - yl)
            + qur * (xp - xl) * (yp - yl)
        )
        q = q / (dx * dy)
    else:
        # don't interpolate, instead, take the native grid resolution
        # value. Choose whichever of ll, lr, ul, lr is closest based on
        # the relative areas is bigger. https://en.wikipedia.org/wiki/Bilinear_interpolation

        all = (xr - xp) * (yu - yp)
        alr = (xp - xl) * (yu - yp)
        aul = (xr - xp) * (yp - yl)
        aur = (xp - xl) * (yp - yl)

        maxa = max((all, alr, aur, aul))

        if all == maxa:
            q = qll
        elif alr == maxa:
            q = qlr
        elif aul == maxa:
            q = qul
        else:
            q = qur

    return q, level


# ===============================================================================
def intersection(rectangle1, rectangle2):
    """
    return True of False if rectangle1 and rectangle2 intersect
    Note: rectangles may be single points as well, with xupper=xlower etc.
    arguments:
       rectangleX: list [xlower,ylower,xupper,yupper]


    """
    xl1 = rectangle1[0]
    yl1 = rectangle1[1]
    xu1 = rectangle1[2]
    yu1 = rectangle1[3]

    xl2 = rectangle2[0]
    yl2 = rectangle2[1]
    xu2 = rectangle2[2]
    yu2 = rectangle2[3]

    nonintersection = (xl1 > xu2) | (xl2 > xu1) | (yl1 > yu2) | (yl2 > yu1)
    intersection = not nonintersection

    return intersection
