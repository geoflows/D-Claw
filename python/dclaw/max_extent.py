import glob
import os

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union

from dclaw.fortconvert import fort2list


def calc_max_extent(
    working_directory=".", output_directory="_output", minimum_level=None
):
    """
    Calculate the maximum extent seen at or above a minimum level
    """
    wd = working_directory
    od = output_directory

    # determine the number of files.
    files = glob.glob(os.path.join(wd, *[od, "fort.t[0-9][0-9][0-9][0-9]"]))
    nfiles = []
    for file in files:
        numstr = os.path.basename(file)[6:]
        nfiles.append(int(numstr))
    nfiles = np.sort(nfiles)

    # create a dictionary with all the grids across all time at
    grid_dict = {}

    for file in nfiles:
        numstring = str(10000 + file)
        framenostr = numstring[1:]
        forttname = os.path.join(wd, od, "fort.t") + framenostr
        fortqname = os.path.join(wd, od, "fort.q") + framenostr
        solutionlist = fort2list(fortqname, forttname)

        for grid in solutionlist:
            level = grid["AMR_level"]
            extent = Polygon.from_bounds(
                grid["xlow"], grid["ylow"], grid["xupper"], grid["yupper"]
            )
            if level not in grid_dict.keys():
                grid_dict[level] = []
            grid_dict[level].append(extent)

    # unary union each grid level.
    max_level_used = max(list(grid_dict.keys()))
    merge_grids = []
    minimum_level = minimum_level or max_level_used
    for level, grid_list in grid_dict.items():
        print(level, len(grid_list))
        if level >= minimum_level:
            print("using level", level)
            merge_grids.append(unary_union(grid_list))
    full_extent = unary_union(merge_grids)

    return full_extent.bounds
