import argparse
import glob
import os

import fiona
import numpy as np
import rasterio
import rasterio.transform
from rasterio import features
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

from dclaw.fortconvert import convertfortdir, fort2list
from dclaw.get_data import (get_amr2ez_data, get_dig_data, get_region_data,
                            get_tsunami_data)


def main():
    parser = argparse.ArgumentParser(
        prog="gridded_maxval",
        description=(
            "Calculate maximum value over a amrclaw region or in the box bounded "
            "by north, west, east, and south. All timesteps are used. First fort "
            "format output is converted to .tif. Then a maximum is calculated. "
            "Cptionally, a shapefile of the extent is generated."
        ),
        usage="%(prog)s [options]",
    )

    parser.add_argument(
        "-wd", "--wdir", nargs="?", help="Working directory", default=".", type=str
    )

    parser.add_argument(
        "-od",
        "--odir",
        nargs="?",
        help="Directory within wd containing fort files",
        default="_output",
        type=str,
    )

    parser.add_argument(
        "-gd",
        "--gdir",
        nargs="?",
        help="Directory within --wdir to place gridded files",
        default="_gridded_output",
        type=str,
    )

    parser.add_argument(
        "-of",
        "--outfile",
        nargs="?",
        help="Output maximum value file name. Placed in -wdir",
        default="maxval.tif",
        type=str,
    )

    parser.add_argument(
        "-cd",
        "--check_done",
        action="store_true",
        default=True,
        help="Check if file processing has already occured and only reprocess new or updated files.",
    )

    parser.add_argument(
        "-nc",
        "--num_cores",
        nargs="?",
        help="Number of cores for fortconvert",
        default=8,
        type=int,
    )

    parser.add_argument(
        "-epsg",
        "--epsg",
        nargs="?",
        help="EPSG code to use for output .tif and .shp files",
        default=None,
    )
    parser.add_argument(
        "-b", "--bilinear",
        help="use bilinear interpolation (default =False) in fortconvert",
        default = False,
        action="store_false",
        )
    parser.add_argument(
        "-w",
        "--west",
        nargs="?",
        help="West extent of bounding box.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "-e",
        "--east",
        nargs="?",
        help="East extent of bounding box.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "-n",
        "--north",
        nargs="?",
        help="North extent of bounding box.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "-s",
        "--south",
        nargs="?",
        help="South extent of bounding box.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "-r",
        "--region",
        nargs="?",
        help="Amrclaw region number to use for extent. Will owr east/south/north/west, if provided.",
        default=None,
        type=int,
    )

    parser.add_argument(
        "-ov",
        "--overwrite_level",
        help="Maxval can be overwritten so long as the level is greater than this specified level.",
        default=None,
        type=int,
    )

    parser.add_argument(
        "-wf",
        "--write_froude",
        help="Write a froude number maximum in the tif.",
        default=False,
        action="store_false",
    )

    parser.add_argument(
        "-shp",
        "--extent_shp",
        help="Write a shapefile with the extent of values greater than a specified threashold.",
        action="store_false",
    )
    parser.add_argument(
        "-sval",
        "--extent_shp_val",
        nargs="?",
        help="Quantity to use for defining the extent of values greater.",
        default="height",
        choices=["height", "momentum", "velocity"],
        type=str,
    )
    parser.add_argument(
        "-sth",
        "--extent_shp_val_thresh",
        nargs="?",
        default=0.0,
        type=float,
        help="Threshold. Default value of 0, along with 'height' yeilds inundated area",
    )
    parser.add_argument(
        "-sof",
        "--extent_shp_val_out_file",
        nargs="?",
        default="extent.shp",
        type=str,
        help="Name of output shapefile.",
    )

    args = parser.parse_args()

    # do some checking with the region.
    amrdata = get_amr2ez_data(args.wdir, args.odir)

    # full extent.
    xhi = amrdata["xupper"]
    yhi = amrdata["yupper"]
    xlow = amrdata["xlower"]
    ylow = amrdata["ylower"]

    if args.region is not None:
        region_data = get_region_data(args.wdir, args.odir)
        west = region_data[args.region]["x1"]
        east = region_data[args.region]["x2"]
        south = region_data[args.region]["y1"]
        north = region_data[args.region]["y2"]
    else:
        west = args.west
        east = args.east
        south = args.south
        north = args.north
        if west is None:
            west = xlow
        if east is None:
            east = xhi
        if south is None:
            south = ylow
        if north is None:
            north = yhi

    # get smallest cell size
    coarse_dx = (xhi - xlow) / amrdata["mx"]
    coarse_dy = (yhi - ylow) / amrdata["my"]
    fine_dx = coarse_dx
    fine_dy = coarse_dy
    for factor in amrdata["inratx"]:
        fine_dx /= factor
    for factor in amrdata["inraty"]:
        fine_dy /= factor

    mx = int((xhi - xlow) / fine_dx)
    my = int((yhi - ylow) / fine_dy)
#    print(my, mx)
#    print(xlow, xhi)
#    print(ylow, yhi)
#    print(fine_dx, fine_dy, coarse_dx, coarse_dy)

    # snap values to the grid (HERE THE GRID IS THE CELL BOUNDING BOXES, NOT
    # THE CELL CENTERS).
    if east is not None or west is not None:
        xs = np.linspace(xlow, xhi+fine_dx, mx+2)
        if east is not None:
            xhi = np.max(xs[xs <= east])
        if west is not None:
            xlow = np.min(xs[xs >= west])
        mx = int((xhi - xlow) / fine_dx)
    if north is not None or south is not None:
        ys = np.linspace(ylow, yhi+fine_dy, my+2)
        if north is not None:
            yhi = np.max(ys[ys <= north])
        if south is not None:
            ylow = np.min(ys[ys >= south])
        my = int((yhi - ylow) / fine_dy)

#    print(my, mx)
#    print(fine_dx, fine_dy)
#    print(xlow, xhi)
#    print(ylow, yhi)

    # make output dir
    if not os.path.exists(os.path.join(args.wdir, args.gdir)):
        os.mkdir(os.path.join(args.wdir, args.gdir))

    # check which files to convert. Could do a temporal filter here..
    # get all files.
    tfiles = np.sort(glob.glob(os.path.join(args.wdir, *[args.odir, "fort.t*"])))
    tfiles = [file for file in tfiles if "tck" not in file]  # remove checkpoint file.
    ntifs = np.sort(glob.glob(os.path.join(args.wdir, *[args.gdir, "fort_q*.tif"])))

    nfiles = []

    # print(args.check_done)
    for tfile in tfiles:
        file = tfile.replace("fort.t", "fort.q")
        numstr = os.path.basename(tfile)[6:]
        tifname = os.path.join(
            ".", *[os.path.join(args.wdir, args.gdir), "fort_q{}.tif".format(numstr)]
        )
        if os.path.exists(tifname) and args.check_done:
            mtime_fort = os.path.getmtime(file)
            mtime_tif = os.path.getmtime(tifname)
            if mtime_tif > mtime_fort:
                process = False
                try:
                    with rasterio.open(tifname, "r"):
                        pass  # test this more?
                except:
                    process = True
                    # if can't be opened, re-write.
            else:
                process = True
        else:
            process = True
        if process:
            nfiles.append(int(numstr))
    nfiles = np.sort(nfiles)

    if len(nfiles) == 0:
        print("check_done = True and no new files to process. No files reprocessed.")
    else:
        print("Processing {} files: {}".format(len(nfiles), nfiles))

        convertfortdir(
            "fortuniform",
            nplots=nfiles,
            outputname="fort_q",
            components="all",
            outdir=os.path.join(args.wdir, args.gdir),
            fortdir=os.path.join(args.wdir, args.odir),
            parallel=True,
            num_cores=args.num_cores,
            topotype="gtif",
            write_level=True,
            epsg=args.epsg,
            bilinear=args.bilinear,
            xlower=xlow,
            xupper=xhi,
            ylower=ylow,
            yupper=yhi,
            mx=mx,
            my=my,
        )

    dig_data = get_dig_data(args.wdir, args.odir)
    rho_f = dig_data["rho_f"]
    rho_s = dig_data["rho_s"]

    dclaw2maxval_withlev(
        wdir=args.wdir,
        odir=args.odir,
        gdir=args.gdir,
        out_file="maxval.tif",
        overwrite_level=args.overwrite_level,
        write_froude=args.write_froude,
        epsg=args.epsg,
        rho_f=rho_f,
        rho_s=rho_s,
        extent_shp=args.extent_shp,
        extent_shp_val=args.extent_shp_val,
        extent_shp_val_thresh=args.extent_shp_val_thresh,
        extent_shp_val_out_file=args.extent_shp_val_out_file,
    )


def dclaw2maxval_withlev(
    wdir=".",
    odir="_output",
    gdir="_gridded_output",
    nplots=None,
    out_file=None,
    rho_f=1000,
    rho_s=2700,
    epsg=None,
    overwrite_level=None,
    write_froude=False,
    extent_shp=True,
    extent_shp_val="height",
    extent_shp_val_thresh=0.0,
    extent_shp_val_out_file=None,
):
    owr_level = overwrite_level
    # get drytol:
    tsudata = get_tsunami_data(wdir, odir)
    drytolerance = tsudata["drytolerance"]
    wavetolerance=tsudata["wavetolerance"]
    sealevel=tsudata["sealevel"]
    # do some checking with the region.
    amrdata = get_amr2ez_data(wdir, odir)
    mxnest = amrdata["mxnest"]  # max number of levels.
    maxlevel = np.abs(mxnest)
    print("****** maxlevel is ", maxlevel)

    owr_level = owr_level or maxlevel
    print("****** owr_level is ", owr_level)

    out_file = out_file or os.path.join(wdir, "maxval.tif")
    extent_shp_val_out_file = extent_shp_val_out_file or os.path.join(
        wdir, "extent.shp"
    )

    # loop through fortq files and add:

    if epsg is not None:
        crs = rasterio.crs.CRS.from_epsg(epsg)
    else:
        crs = None

    files = np.sort(glob.glob(os.path.join(wdir, gdir, "*.tif")))
    if len(files) == 0:
        raise ValueError("no files")

    # Read them all in.
    with rasterio.open(files[0], "r") as src:
        dims = (src.meta["height"], src.meta["width"])

    # constants.
    hmin_fill = 99999.0
    nodata = -9999
    time_fill = -1.

    # initialize arrays for output values.
    eta_max = nodata * np.ones(dims, dtype="float32")
    h_max = np.zeros(dims, dtype="float32")
    h_min = np.ones(dims, dtype="float32")
    h_min[:] = hmin_fill
    m_max = np.zeros(dims, dtype="float32")
    vel_max = np.zeros(dims, dtype="float32")
    mom_max = np.zeros(dims, dtype="float32")
    lev_max = np.zeros(dims, dtype="float32")
    fr_max = np.zeros(dims, dtype="float32")

    # initialize arrays for owr_levels (e.g., what level needs to be seen to owr.)
    first_geq_owr = np.zeros(
        dims, dtype=bool
    )  # has overwrite previously been exceeded.

    wave_all = np.zeros(dims, dtype=bool)
    eta_owr_lev = np.zeros(dims, dtype=int)
    h_owr_lev = np.zeros(dims, dtype=int)
    h_min_owr_lev = np.zeros(dims, dtype=int)
    m_owr_lev = np.zeros(dims, dtype=int)
    vel_owr_lev = np.zeros(dims, dtype=int)
    mom_owr_lev = np.zeros(dims, dtype=int)
    fr_owr_lev = np.zeros(dims, dtype=int)
    arrival_lev = np.zeros(dims, dtype=int)

    arrival_time = time_fill * np.ones(dims, dtype="float32")
    eta_max_time = time_fill * np.ones(dims, dtype="float32")
    vel_max_time = time_fill * np.ones(dims, dtype="float32")

    # is thickness present?
    # if never present, its equal to 0
    # if present, its equal to the largest level seen.
    h_level_masked = np.zeros(dims, dtype=int)

    for file in files:
        fortt = file[:-4].replace("fort_q", "fort.t")
        # print(file)
        frameno = int(file[-8:-4])

        calc = True
        if nplots is not None:
            if frameno not in nplots:
                calc = False
        if calc:
            # print(file[-8:-4], frameno)
            with open(fortt, "r") as f:
                lines = f.readlines()
            time = float(lines[0].split()[0])
            with rasterio.open(file, "r") as src:
                profile = src.profile
                transform = src.transform

                dx = transform[0]

                h = src.read(1)
                hu = src.read(2)
                hv = src.read(3)
                hm = src.read(4)
                eta = src.read(8)
                level = src.read(9)

                with np.errstate(divide="ignore", invalid="ignore"):
                    m = hm / h
                    m[np.isnan(m)] = 0
                    vel = ((hu / h) ** 2 + (hv / h) ** 2) ** 0.5
                    vel[np.isnan(vel)] = 0

                    fr = vel / np.sqrt(9.81 * h)

                density = (1.0 - m) * rho_f + (m * rho_s)
                mom = (h * dx * dx) * density * vel

                # May 12, 2022, maxval algorithm updated in light of issues
                # identified using tsunami-style refinement. That is, late in
                # a simulation, after the wave has passed, a cell might get
                # refined to the highest (shoreline) level, but only because it
                # was needed to make an efficient bounding box around the
                # flagged cells. This meant that there were cells in just-deeper
                # that shoreline water which would get set to maximum eta values
                # of zero only because the refinement level went from second
                # highest to highest.
                # this issue was addressed in the following way.
                #   creation of an overwrite_level input, to permit a user to
                #   indicate that any level greater than or equal to this level
                #   should be considered equivalent in the eyes of the maxval
                #   algorithm. That is, it is only when the refinement level
                #   trips over the boundary between less than overwrite_level
                #   and geq overwrite level THE FIRST TIME that increasing the
                #   triggers resetting the maximum value to the the current
                #   value.

                # keep track of where level increased and max level.
                level_increased = level > lev_max
                lev_max[level_increased] = level[level_increased]

                # determine where h is located at this timestep.
                h_present = h > drytolerance

                # determine where h is present and the level was refined to
                # a higher level.
                h_present_and_level_higher = h_present & (level > h_level_masked)

                # determine where refinement resulted in a cell becoming dry
                # because we use level>h_level_masked this should only occur
                # the first time the cell refines to the next highest level.
                refined_to_dry = (h_present == False) & (level > h_level_masked)

                # set values of h, hmin, m, eta, vel, froude to nodata or
                # zero where refined to dry occured.
                eta_max[refined_to_dry] = nodata
                h_max[refined_to_dry] = 0
                h_min[refined_to_dry] = hmin_fill
                m_max[refined_to_dry] = 0
                vel_max[refined_to_dry] = 0
                mom_max[refined_to_dry] = 0
                fr_max[refined_to_dry] = 0
                arrival_time[refined_to_dry] = time_fill
                eta_max_time[refined_to_dry] = time_fill
                vel_max_time[refined_to_dry] = time_fill

                # if refinement is to sea-level, eta, hmax, and hmin should be
                # reset. this implies incoming refinement.
                # set values of h, hmin to sea level (or current value)
                refined_to_sea_level = (eta == sealevel) & (level > h_level_masked)

                eta_max[refined_to_sea_level] = sealevel
                h_max[refined_to_sea_level] = h[refined_to_sea_level]
                h_min[refined_to_sea_level] = h[refined_to_sea_level]
                arrival_time[refined_to_sea_level] = time_fill
                eta_max_time[refined_to_sea_level] = time_fill
                vel_max_time[refined_to_sea_level] = time_fill

                # use definition of a wave defined in tsunami refinement.
                # reset wave based on
                wave_now = (np.abs(eta - sealevel) > wavetolerance) & h_present
                wave_all[wave_now] = True

                # reset where refined to dry or sea level.
                wave_all[refined_to_sea_level] = False
                wave_all[refined_to_dry] = False

                # determine if it is the first time greater than the overwrite
                # level (and h_present)
                geq_ovr = lev_max >= owr_level
                first_time_owr = geq_ovr & (first_geq_owr == False) & h_present & wave_now
                first_geq_owr[first_time_owr] = True

                # update h_level_masked
                h_level_masked[h_present_and_level_higher] = level[
                    h_present_and_level_higher
                ]

                # set values of h, hmin, m, eta, vel, froude to value the first time
                # overwrite.
                eta_max[first_time_owr] = eta[first_time_owr]
                h_max[first_time_owr] = h[first_time_owr]
                h_min[first_time_owr] = h[first_time_owr]
                m_max[first_time_owr] = m[first_time_owr]
                vel_max[first_time_owr] = vel[first_time_owr]
                mom_max[first_time_owr] = mom[first_time_owr]
                fr_max[first_time_owr] = fr[first_time_owr]

                # determine whether to update. condition is:
                #    level greater than or equal to the overwrite level and
                #    value greater than existing value
                # OR
                #    first time greater than the overwrite level.
                #
                update_eta_max = (level >= eta_owr_lev) & (eta > eta_max) & wave_now
                update_h_max = (level >= h_owr_lev) & (h > h_max) & wave_now
                update_h_min = (level >= h_min_owr_lev) & (h < h_min) & wave_now
                update_m = (level >= m_owr_lev) & (m > m_max)& wave_now
                update_vel = (level >= vel_owr_lev) & (vel > vel_max) & wave_now
                update_mom = (level >= mom_owr_lev) & (mom > mom_max) & wave_now
                update_fr = (level >= fr_owr_lev) & (fr > fr_max) & wave_now

                # ensure owr_level arrays do not exceed owr_level
                # first update to level seen,
                eta_owr_lev[update_eta_max] = level[update_eta_max]
                h_owr_lev[update_h_max] = level[update_h_max]
                h_min_owr_lev[update_h_min] = level[update_h_min]
                m_owr_lev[update_m] = level[update_m]
                vel_owr_lev[update_vel] = level[update_vel]
                mom_owr_lev[update_mom] = level[update_mom]
                fr_owr_lev[update_fr] = level[update_fr]

                # second, ensure it doesn't exceed the owr level.
                eta_owr_lev[eta_owr_lev > owr_level] = owr_level
                h_owr_lev[h_owr_lev > owr_level] = owr_level
                h_min_owr_lev[h_min_owr_lev > owr_level] = owr_level
                m_owr_lev[m_owr_lev > owr_level] = owr_level
                vel_owr_lev[vel_owr_lev > owr_level] = owr_level
                mom_owr_lev[mom_owr_lev > owr_level] = owr_level
                fr_owr_lev[fr_owr_lev > owr_level] = owr_level

                # update max values.
                eta_max[update_eta_max] = eta[update_eta_max]
                h_max[update_h_max] = h[update_h_max]
                h_min[update_h_min] = h[update_h_min]
                m_max[update_m] = m[update_m]
                vel_max[update_vel] = vel[update_vel]
                mom_max[update_mom] = mom[update_mom]
                fr_max[update_fr] = fr[update_fr]

                # update arrival time,
                # set arrival time to the first timestep that has eta>0.01 and highest level seen.
                # here
                owr_arrival = (np.abs(eta - sealevel) > wavetolerance) & (arrival_time < 0) & (level > arrival_lev) & h_present
                arrival_lev[owr_arrival] = level[owr_arrival]
                arrival_time[owr_arrival] = time

                # set other times to arrival time, to indicate the wave has
                # arrived there and thus its valid to set a max.
                eta_max_time[owr_arrival] = time
                vel_max_time[owr_arrival] = time

                # we want the first peak
                not_super_late = ((time - eta_max_time) < (1 * 60)) & (
                    arrival_time >= 0
                )
                # use 10 minutes
                # presuming time has been set.
                # presuming the arrival time has passed.
                # and presuming that this timestep eta gets updated.
                update_eta_time = update_eta_max & not_super_late
                update_vel_time = update_vel & not_super_late

                eta_max_time[update_eta_time] = time
                vel_max_time[update_vel_time] = time

    never_inundated = h_max < drytolerance
    never_wave = wave_all == False

    zero_out = never_inundated | never_wave

    eta_max[zero_out] = nodata
    h_max[zero_out] = nodata
    h_min[h_min == 0] = nodata
    m_max[zero_out] = nodata
    eta_max[zero_out] = nodata
    vel_max[zero_out] = nodata
    mom_max[zero_out] = nodata
    eta_max_time[zero_out] = nodata
    vel_max_time[zero_out] = nodata
    arrival_time[zero_out] = nodata
    fr_max[zero_out] = nodata

    # where less than zero, wave never reached.
    eta_max_time[eta_max_time < 0] = nodata
    vel_max_time[vel_max_time < 0] = nodata
    arrival_time[arrival_time < 0] = nodata

    # then read in grided output into t, h, m, vel, mom, etc.
    # clip to required extent.

    # calculate output bands and write out as .tif
    # this is also where we figure out time of h max and time of vel max.

    out_profile = profile
    out_profile["height"], out_profile["width"] = m_max.shape
    out_profile["dtype"] = "float32"
    out_profile["count"] = 10
    out_profile["transform"] = transform
    out_profile["nodata"] = nodata

    if write_froude:
        out_profile["count"] = 11
    with rasterio.open(out_file, "w", **out_profile) as dst:
        dst.write(h_max, 1)
        dst.write(vel_max, 2)
        dst.write(mom_max, 3)
        dst.write(m_max, 4)
        dst.write(eta_max_time, 5)
        dst.write(vel_max_time, 6)
        dst.write(eta_max, 7)
        dst.write(lev_max, 8)
        dst.write(arrival_time, 9)
        dst.write(h_min, 10)
        if write_froude:
            dst.write(fr, 11)

    if extent_shp:
        if extent_shp_val == "height":
            extent = h_max > extent_shp_val_thresh
        if extent_shp_val == "momentum":
            extent = mom_max > extent_shp_val_thresh
        if extent_shp_val == "velocity":
            extent = vel_max > extent_shp_val_thresh
        transform = out_profile["transform"]
        geoms = []
        for i, (s, v) in enumerate(
            features.shapes(extent.astype(np.int16), mask=extent, transform=transform)
        ):
            shp = shape(s).buffer(0)
            if shp.is_valid == True:
                geoms.append(shp)
        out_shp = unary_union(geoms).buffer(0)

        assert out_shp.is_valid == True
        with fiona.open(
            extent_shp_val_out_file,
            "w",
            driver="Shapefile",
            crs=crs,
            schema={"properties": [], "geometry": "Polygon"},
        ) as dst:
            dst.writerecords([{"properties": {}, "geometry": mapping(out_shp)}])
