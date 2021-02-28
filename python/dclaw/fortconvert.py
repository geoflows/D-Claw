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
       forttheaderread
       fortqheaderread
       pointfromfort
       fort2list (used to read fort data...4.x compatible...)


    Future plans:
        probably should use the pyclaw solution class rather than solutionlist
        or at least provide options for both.

        provide conversion to clawpack netcdf and possibly other binary formats.

    Note:
        this was done expeditiously and might be inefficient or bug infested

    David George dgeorge@uw.edu 2014.05.20

"""

import copy
import os
import string

import numpy as np

import dclaw.topotools as gt


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
        nplots: integer for the number of plots, list of plots or fort.nplot filename
            a sinlge integer will be all plots up to that integer
            for individual plots make nplots a list
        outputtype: 'scattered','topotype' , 'fortrefined' , 'fortuniform'
                    Note: only scattered option retains actual non-interpolated output, generally.
                     --scattered: highest level data as columns x,y,q
                     --topotype: like a DEM with a single component of q (useful for GIS software)
                     --fortrefined: standard clawpack output style, but a single grid at highest resolution
                     --fortuniform: standard clawpack output style, but a single grid with user defined grid parameters.
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
            raise SystemExit(0)

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
                ]
            )

        elif outputtype == "fortrefined":
            outfortt = os.path.join(outdir, "fort.t" + framenostr)
            topotype = kwargs.get("topotype", None)
            _func = fort2refined

            arg_list.append([frameno, outfname, outfortt, components, topotype])

        elif outputtype == "fortuniform":
            outfortt = os.path.join(outdir, "fort.t" + framenostr)
            topotype = kwargs.get("topotype", None)
            func = fort2uniform

            _arg_list.append(
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
                ]
            )

    # now run in parallel based on func and arg list
    if parallel:
        try:
            from joblib import Parallel, delayed

            Parallel(n_jobs=num_cores)(delayed(_func)(*args) for args in arg_list)
        except ImportError:
            raise ImportError("joblib needed for parallel functionality")

    else:
        # loop through prepared args.
        for args in arg_list:
            print("Converting {}".format(args[1]))
            _func(*args)

    # return to curdir if changed.
    os.chdir(curdir)


# =============================================================================
def fort2xyqscattered(framenumber, outfile=None, components="all"):
    """
    convert a fort.qXXXX amr file into a scattered data with columns x,y,q...
    q can be 1-meqn columns according to the list 'components'
    data is taken from the finest of grid intersections

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
            x = grid["xlow"] + grid["dx"] * (0.5 + (np.arange(grid["mx"], dtype=float)))
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
):
    """
    convert fort.qXXXX with AMR data into fort.qXXXX with data on a uniform single grid.
    Resolution is user defined.
    Format is still standard clawpack fort.qXXXX /fort.tXXXX and can be plotted with clawpack utilities
    Should call with outdir being a new directory to keep original fort.q/fort.t files. Names are the same

    arguments
    ----------
        framenumber : of fort.qXXXX
        outfortq: name of output fort.qXXXX file
        outfortt: name of output fort.tXXXX file
            if above =None routine returns a numpy array
        xlower,xupper,ylower,yupper,mx,my: output grid parameters
        components: list of q components eg. [1,3,5] or 'all' for all components
    """

    if isinstance(outfortq, str):
        foutq = open(outfortq, "w")
    else:
        foutq = outfortq

    if isinstance(outfortt, str):
        foutt = open(outfortt, "w")
    else:
        foutf = outfortt

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

    if (not outfortq) or (topotype is not None):
        Q = np.empty((mx * my, len(qlst)))

        for j in range(my):
            y = ylow + (j + 0.5) * dy
            for i in range(mx):
                x = xlow + (i + 0.5) * dx
                qv = pointfromfort((x, y), solutionlist)
                qout = qv[qlst]
                Q[j * mx + i] = qout

        if topotype is not None:
            raise NotImplementedError()
        else:
            return fortheader, Q
    else:
        forttheaderwrite(fortheader, foutt)
        foutt.close()
        fortqheaderwrite(fortheader, foutq, closefile=False)

        for j in range(my):
            foutq.write("\n")
            y = ylow + (j + 0.5) * dy
            for i in range(mx):
                x = xlow + (i + 0.5) * dx
                qv = pointfromfort((x, y), solutionlist)
                qout = qv[qlst]
                for q in qout:
                    foutq.write("%s " % float(q))
                foutq.write("\n")
        foutq.close()


# ==============================================================================
def fort2refined(framenumber, outfortq, outfortt, components="all", topotype=None):
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
    """

    if isinstance(outfortq, str):
        foutq = open(outfortq, "w")
    else:
        foutq = outfortq

    if isinstance(outfortt, str):
        foutt = open(outfortt, "w")
    else:
        foutf = outfortt

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

    if (not outfortq) or (topotype is not None):

        # prepare to return an array, OR write to topo file.

        Q = np.empty((mx * my, len(qlst)))
        for j in range(my):
            y = ylow + (j + 0.5) * dy
            for i in range(mx):
                x = xlow + (i + 0.5) * dx
                qv = pointfromfort((x, y), solutionlist)
                qout = qv[qlst]
                Q[j * mx + i] = qout

        # if topotype is specified, write out as topotype instead of
        # array or standard fort.q
        if topotype is not None:
            xv = np.array(xlow + dx * np.arange(mx))
            yv = np.array(ylow + dy * np.arange(my))
            (X, Y) = np.meshgrid(xv, yv)
            Y = np.flipud(Y)

            if topotype == "gtif":
                outfile = outfortq.replace(".", "_") + ".tif"
                gt.griddata2gtif(
                    X,
                    Y,
                    np.flip(np.moveaxis(Q.reshape((mx, my, len(qlst))), (0, 1, 2), (1, 2, 0)), axis=1),
                    # reshape into nr, nc, neq, and shift axis so neq is at front,
                    # finally flip along axis 1 so that up is up.
                    outfile,
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
        forttheaderwrite(fortheader, foutt)
        foutt.close()
        fortqheaderwrite(fortheader, foutq, closefile=False)

        for j in range(my):
            foutq.write("\n")
            y = ylow + (j + 0.5) * dy
            for i in range(mx):
                x = xlow + (i + 0.5) * dx
                qv = pointfromfort((x, y), solutionlist)
                qout = qv[qlst]
                for q in qout:
                    foutq.write("%s " % float(q))
                foutq.write("\n")
        foutq.close()


# ==============================================================================
def fort2topotype(
    framenumber, outfile, fortdir, xll, yll, cellsize, ncols, nrows, m=1, topotype=2
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

        xll = xll + 0.5 * cellsize
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
                qv = pointfromfort((xp, yp), solutionlist)
                Q[i, j] = qv[meqn - 1] - qv[0]

    elif m == "depth":
        for j in range(ncols):
            xp = X[0, j]
            for i in range(nrows):
                yp = Y[i, 0]
                qv = pointfromfort((xp, yp), solutionlist)
                depth = qv[0]
                if depth <= 1.0e-3:
                    depth = nodata_value
                Q[i, j] = depth

    elif m == "eta":
        for j in range(ncols):
            xp = X[0, j]
            for i in range(nrows):
                yp = Y[i, 0]
                qv = pointfromfort((xp, yp), solutionlist)
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
                qv = pointfromfort((xp, yp), solutionlist)
                Q[k] = qv

    else:
        for j in range(ncols):
            xp = X[0, j]
            for i in range(nrows):
                yp = Y[i, 0]
                qv = pointfromfort((xp, yp), solutionlist)
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
def fort2griddata(framenumber, xll, yll, cellsize, ncols, nrows, m=1):
    """
    convert data in a fort file of framenumber = XXXX, ie fort.qXXXX
    to numpy arrays X,Y,Q (single gridded data)
    m is the component of q, ie. the column in the fort.qXXXX file
    """

    try:
        import geotools.topotools as gt
    except:
        import clawpack.geoclaw.topotools as gt

    numstring = str(10000 + framenumber)
    framenostr = numstring[1:]
    forttname = "fort.t" + framenostr
    fortqname = "fort.q" + framenostr

    solutionlist = fort2list(fortqname, forttname)

    xv = np.array(xll + cellsize * np.arange(ncols))
    yv = np.array(yll + cellsize * np.arange(nrows))

    (X, Y) = np.meshgrid(xv, yv)

    Y = np.flipud(Y)

    Q = np.zeros(np.shape(X))

    for j in range(ncols):
        xp = X[0, j]
        for i in range(nrows):
            yp = Y[i, 0]
            qv = pointfromfort((xp, yp), solutionlist)
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
def pointfromfort(point, solutionlist):
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
    except:
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
        except:
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
    xl = xlow + (i1 - 1) * dx
    xr = xl + dx
    yl = ylow + (j1 - 1) * dy
    yu = yl + dy

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
    q = (
        qll * (xr - xp) * (yu - yp)
        + qlr * (xp - xl) * (yu - yp)
        + qul * (xr - xp) * (yp - yl)
        + qur * (xp - xl) * (yp - yl)
    )
    q = q / (dx * dy)

    return q


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
