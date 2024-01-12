import os

from pyclaw.data import Data

"""

Lightweight functions to get dictionaries of attributes from .data files.
KRB April 2022

"""


def get_tsunami_data(project_path, output="_output", file="settsunami.data"):
    data = Data(os.path.join(project_path, output, file))
    return {key: data.__dict__[key] for key in data.attributes}


def get_dig_data(project_path, output="_output", file="setdig.data"):
    data = Data(os.path.join(project_path, output, file))
    return {key: data.__dict__[key] for key in data.attributes}


def get_amr2ez_data(project_path, output="_output", file="amr2ez.data"):
    try:
        data = Data(os.path.join(project_path, output, file))
    except:
        data = Data(os.path.join(project_path, output, "amr.data"))
    return {key: data.__dict__[key] for key in data.attributes}


def get_gauge_data(project_path, output="_output", file="setgauges.data"):
    setgaugefile = os.path.join(project_path, output, file)

    gauge_dict = {}

    with open(setgaugefile, "r") as fid:
        inp = "#"
        while inp == "#":
            inpl = fid.readline()
            inp = inpl[0]

        inp = fid.readline()
        mgauges = int(inp.split()[0])
        linesread = 0

        while linesread < mgauges:
            row = fid.readline().split()
            if row != []:
                linesread = linesread + 1
                gaugeno = int(row[0])
                gauge_dict[gaugeno] = {}

                gauge_dict[gaugeno]["x"] = float(row[1])
                gauge_dict[gaugeno]["y"] = float(row[2])
                gauge_dict[gaugeno]["tmin"] = float(row[3])
                gauge_dict[gaugeno]["tmax"] = float(row[4])

    return gauge_dict


def get_region_data(project_path, output="_output", file="setregions.data"):
    setregionfile = os.path.join(project_path, output, file)

    region_dict = {}

    with open(setregionfile, "r") as fid:
        inp = "#"
        while inp == "#":
            inpl = fid.readline()
            inp = inpl[0]

        inp = fid.readline()
        mregions = int(inp.split()[0])
        linesread = 0

        while linesread < mregions:
            row = fid.readline().split()
            if row != []:
                linesread = linesread + 1
                regionno = (
                    len(region_dict) + 1
                )  # not officially numbered. # go in order, from 1 onwards.
                # order of .data file.
                region_dict[regionno] = {}
                region_dict[regionno]["minlevel"] = float(row[0])
                region_dict[regionno]["maxlevel"] = float(row[1])
                region_dict[regionno]["t1"] = float(row[2])
                region_dict[regionno]["t2"] = float(row[3])
                region_dict[regionno]["x1"] = float(row[4])
                region_dict[regionno]["x2"] = float(row[5])
                region_dict[regionno]["y1"] = float(row[6])
                region_dict[regionno]["y2"] = float(row[7])

        return region_dict
