#!/usr/bin/env python
# encoding: utf-8
r"""
Data Module

Contains the general class definition and the subclasses of the Clawpack data
objects.

:Authors:
    Kyle T. Mandli and Randall J. LeVeque (2008-08-07) Initial version

    Randall J. LeVeque (2008-08-07) Plotting data objects

    Alan McIntyre (2009-01-01) Speed ups and rebuilding of Data

    Kyle T. Mandli (2009-04-01) Stripped down and improved version
"""
# ============================================================================
#      Copyright (C) 2008 Kyle T. Mandli <mandli@amath.washington.edu>
#      Copyright (C) 2008 Randall J. LeVeque <rjl@amath.washington.edu>
#      Copyright (C) 2009 Alan McIntyre <amcin001@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD)
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import copy
import logging
import os
import re
import shutil

# place DIG attributes and description in one place (since this is documents )
_DIG_ATTRS = {
    "rho_s": "solid grain density (kg/m^3)",
    "rho_s": "solid grain density (kg/m^3)",
    "rho_f": "pore-fluid density  (kg/m^3)",
    "phi_bed": "basal friction angle (degrees)",
    "theta_input": "slope angle (degrees)",
    "delta": "characteristic grain diameter (m)",
    "kappita": "permeability at m=setdig.m0 (m^2), k0 in G&I eq 2.7",
    "mu": "viscosity of pore-fluid (Pa-s)",
    "alpha_c": "debris compressibility constant (#)",
    "m_crit": "critical state value of m (#)",
    "c1": "dilation coefficient 1 (#)",
    "m0": "initial solid volume fraction (#)",
    "sigma_0": "baseline stress for definition of compressibility",
    "alpha_seg": "coefficient of segregation velocity profile",
    "bed_normal": "bed_normal = 1 requires theta in aux for slope in one direction",
    "phi_seg_coeff": "adjustment to friction coefficient based on segregation",
    "entrainment": "flag for entrainment, 0 = no entrainment",
    "entrainment_rate": "rate of entrainment parameter 0-1",
    "mom_autostop": "flag for momentum autostop F = no autostop, T = autostop",
    "mom_perc": "percentage of max momentum for autostop, default is 0.05 (5%)",
    "src_ftn_num_sr": "number of in-domain sources",
    "fric_offset_val": "start/stop friction offset in degrees. if this value is >0, then hysteretic friction is used.",
    "fric_star_val": "deep friction offset in degrees",
    "chi_init_val": "initial mixture of species size 0-1",
    "kappita_diff": "permeability multiplier for different size species, kappita is used for species1, kappita*kappita_diff used for species2",
    "outaux": "flag for writing aux to output F = not written, T = written",
    "init_ptype": "0 = hydro, 1 = failure, 2= p(t)",
    "init_pmax_ratio": "p(init_ptf2)= hydro*init_pmax_ratio",
    "init_ptf": " p(init_ptf) = failure",
    "init_ptf2": "p(init_ptf2)= hydro*init_pmax_ratio",
}


# ========== Parse Value Utility Function ====================================
def _parse_value(value):
    r"""
    Attempt to make sense of a value string from a config file.  If the
    value is not obviously an integer, float, or boolean, it is returned as
    a string stripped of leading and trailing whitespace.

    :Input:
        - *value* - (string) Value string to be parsed

    :Output:
        - (id) - Appropriate object based on *value*
    """
    value = value.strip()
    if not value:
        return None

    # assume that values containing spaces are lists of values
    if len(value.split()) > 1:
        return [_parse_value(vv) for vv in value.split()]

    try:
        # see if it's an integer
        value = int(value)
    except ValueError:
        try:
            # see if it's a float
            value = float(value)
        except ValueError:
            # see if it's a bool
            if value[0] == "T":
                value = True
            elif value[0] == "F":
                value = False

    return value


# ============================================================================
#  General Data Class
# ============================================================================
class Data(object):
    r"""
    Generalized clawpack data object

    Generalized class for Clawpack data.  Contains generic methods for reading
    and writing data to and from a data file.

    :Initialization:
        Input:
         - *data_files* - (List of strings) Paths to data files to be read in,
           an empty data object can be created by providing no data files.
         - *attributes* - (List of strings) List of required attribute names
           which will be initialized to None.

    :Version: 1.2 (2009-04-01)
    """

    __name__ = "Data"

    def __init__(self, data_files=[], attributes=None):
        """
        Initialize a Data object

        See :class:`Data` for more info.
        """

        # Internal bookkeeping variables
        self.__attributes = []
        self.__owners = {}

        # Setup data logger
        self.logger = logging.getLogger("data")

        # Initialize from attribute list provided
        if attributes:
            for attr in attributes:
                self.add_attribute(attr, None, None)

        # Read data files from data_files list
        if isinstance(data_files, str):
            data_files = [data_files]
        elif not isinstance(data_files, list):
            raise Exception("data_files must be a list of strings")

        if len(data_files) > 0:
            self.read(data_files)

    # ========== Return string representation of this Data Object ========
    def __str__(self):
        output = "%s%s%s\n" % ("Name".ljust(25), "Value".ljust(12), "Owner".ljust(12))
        for (k, v) in list(self.items()):
            output += "%s%s%s\n" % (
                str(k).ljust(25),
                str(v).ljust(12),
                str(self.__owners[k]).ljust(12),
            )
        return output

    # ========== Access Methods ==============================================
    def add_attribute(self, name, value=None, owner=None):
        r"""
        Adds an attribute called name to the data object

        If an attribute needs to be added to the object, this routine must be
        called or the attribute will not be written out.

        :Input:
         - *name* - (string) Name of the data attribute
         - *value* - (id) Value to set *name* to, defaults to None
         - *owner* - (id) Owner of this particular attribute
        """
        setattr(self, name, value)
        self.__owners[name] = owner
        if name not in self.__attributes:
            self.__attributes.append(name)

    def remove_attributes(self, arg_list):
        r"""
        Remove the listed attributes.
        """

        # Convert to list if args is not already a list
        if not isinstance(arg_list, list):
            arg_list = [arg_list]

        for arg in arg_list:
            self.__owners.pop(arg)
            self.__attributes.remove(arg)
            delattr(self, arg)

    def attributes():
        def fget(self):
            return self.__attributes

        return locals()

    attributes = property(**attributes())

    def has_attribute(self, name):
        r"""
        Check if this data object has the given attributes

        :Input:
         - *name* - (string) Name of attribute

        :Output:
         - (bool) - True if data object contains a data attribute name
        """
        return name in self.__attributes

    def set_owner(self, name, owner):
        r"""
        Sets the owner of the given data

        :Input:
         - *name* - (string) Name of attribute
         - *owner* - (id) Owner of the attribute
        """
        if name not in self.__attributes:
            raise KeyError("No attribute named %s" % name)
        self.__owners[name] = owner

    def get_owner(self, name):
        r"""
        Returns the owner of the data attribute name

        :Input:
         - *name* - (string) Name of attribute

        :Output:
         - (id) - Owner of attribute
        """
        return self.__owners[name]

    def get_owners(self, supplementary_file=None):
        r"""
        Returns a list of owners excluding the owner None

        If supplementary_file is provided, None is replace by that owner.

        :Input:
         - *supplementary_file* - (string) Supplementary file, defaults to
           None.

        :Output:
         - (list) - Returns a list of owners
        """
        owners = []
        for (key, owner) in list(self.__owners.items()):
            if owner is None:
                self.__owners[key] = supplementary_file
            # This simultaneously finds one instance of an owner and tests
            # to see if the supplementary_file is not None
            owner = self.__owners[key]
            if owner not in owners and owner is not None:
                owners.append(owner)
        return owners

    def iteritems(self):
        r"""
        Returns an iterator of keys and values from this object

        :Output:
         - (Iterator) Iterator over keys and values
        """
        return [(k, getattr(self, k)) for k in self.__attributes]

    # ========== Read in a collection of data files ==========================
    def read(self, data_paths):
        r"""
        Read data in from a clawpack style data file.

        Any lines of the form::

           values  =:  name

        is used to set an attribute of self.

        INPUT

        data_paths : Path to a data file to be read in, can also be a list
                    of files to be read in.
        """
        if isinstance(data_paths, str):
            data_paths = [data_paths]

        for filename in data_paths:
            filename = os.path.abspath(filename)
            if not os.path.exists(filename):
                raise Exception("No such data file: %s" % filename)

            print("Reading from %s" % filename)
            with open(filename, "r") as f:
                lines = f.readlines()

                for lineno, line in enumerate(lines):
                    if "=:" not in line:
                        continue

                    value, tail = line.split("=:")
                    varname = tail.split()[0]

                    oldval = getattr(self, varname, None)
                    newval = _parse_value(value)
                    if oldval is not None:
                        vals = "(old=%r,new=%r)" % (oldval, newval)
                        # self.logger.debug("Overwriting %s %s" % (varname, vals))

                    # if newval is None:
                    # self.logger.warning("Empty value for %s" % varname)

                    self.add_attribute(varname, newval, filename)

    # ========== Write out the data from this object =========================
    def write(self, data_files=None, supplementary_file=None):
        r"""
        Write out the contents of the data object

        This method writes out the current required attributes of the data
        object to a file or list of files.  The format for the output will be
        in the form::

            values =: name

        The order is either retained from the files in the data list or
        written in the order they are in the list of attributes

        The behavior is determined by the arguments passed into the routine:

        No arguments:

            The contents of the object will be written out to the files listed
            in data_files only if each attribute is contained in the file.
            This implies that if an attribute is not located in any of the
            files in data_files then it will not be written out to file.

        If data_files is provided and data_files is a valid owner:

            Write out the attributes appropriate to that file

        If data_files is provided and not(data_files in owner):

            Write all attributes to the data_file given

        If supplementary_file is provided:

            Write out any attributes without an owner to this file, all owned
            attributes will be written out to the approriate files
        """

        # Expand supplementary_file if it is not None
        if supplementary_file is not None:
            supplementary_file = os.path.abspath(supplementary_file)

        # Create list of owners
        owners = self.get_owners(supplementary_file=supplementary_file)

        # print 'in write: data_files = ',data_files
        # Write to the entire owner list
        if data_files is None:
            file_list = owners
        elif isinstance(data_files, str):
            # print '<p>in write'; import sys; sys.exit(0)
            path = os.path.abspath(data_files)
            if path not in owners:
                # Create temporary data file to store all data in the object
                try:
                    temp_file = open(path, "w")
                    for key in self.__attributes:
                        temp_file.write("1 =: %s\n" % key)
                    temp_file.close()
                except IOError as xxx_todo_changeme:
                    (errno, strerror) = xxx_todo_changeme.args
                    print(("I/O error(%s): %s" % (errno, strerror)))
                    raise
                except:
                    raise
            # Add the path to the file_list
            file_list = [path]
        # Check to make sure all the paths are in the owner list
        elif isinstance(data_files, list):
            for path in data_files:
                path = os.path.abspath(path)
                if not (path in owners):
                    print("%s is not a registered owner!")
                    return
            file_list = data_files
        else:
            raise Exception("Invalid argument list given to write().")

        # Create temporary supplementary file if requested
        if supplementary_file is not None:
            try:
                sup_file = open(supplementary_file, "w")
                for attr in self.__attributes:
                    if self.__owners[attr] is supplementary_file:
                        sup_file.write("-1 =: %s\n" % attr)
                        self.__owners[attr] = supplementary_file
                sup_file.close()
            except:
                raise

        # Regular expression for searching each file
        regexp = re.compile(r"(?P<values>.*)=:(?P<name>.*)")
        # Loop over each file
        for data_path in file_list:
            # Open the data file and temporary file
            try:
                data_file = open(data_path, "r")
            except (IOError):
                raise
            try:
                temp_path = os.path.join(
                    os.path.dirname(data_path), "temp." + os.path.basename(data_path)
                )
                temp_file = open(temp_path, "w")
            except (IOError):
                print("IOERROR")
                raise

            try:
                for line in data_file:
                    result = regexp.search(line)
                    if result:
                        name = re.split(r"\s+", result.group("name").strip())[0]
                        values = re.split(r"\s+", result.group("values").strip())

                        if len(values) == 0:
                            line = ""
                        elif (
                            self.__owners[name] == data_path
                            or data_path not in self.__owners
                        ):
                            newvalues = getattr(self, name)

                            # Convert newvalues to an appropriate string repr
                            if isinstance(newvalues, tuple) | isinstance(
                                newvalues, list
                            ):
                                # Remove [], (), and ','
                                newstring = repr(newvalues)[1:-1]
                                newstring = newstring.replace(",", "")
                            elif isinstance(newvalues, bool):
                                if newvalues:
                                    newstring = "T"
                                else:
                                    newstring = "F"
                            else:
                                newstring = repr(newvalues)

                            newstart = newstring.ljust(25)
                            line = line.replace(
                                result.group("values") + "=:", newstart + " =:"
                            )
                        else:
                            print(("Error writing out %s" % name))
                            raise AttributeError(name)

                    # Write the new line
                    temp_file.write(line)
            except:
                raise

            # Close files
            data_file.close()
            temp_file.close()

            # Rename the temporary file to the data_path name
            try:
                shutil.move(temp_path, data_path)
            except:
                raise


# -----------------------------------------------------------------------

# New classes and functions for dealing with data in setrun function.


class ClawInputData(Data):
    r"""
    Object that will be written out to claw.data.
    """

    def __init__(self, ndim):
        super(ClawInputData, self).__init__()
        self.add_attribute("ndim", ndim)

        # Set default values:
        if ndim == 1:
            self.add_attribute("mx", 100)
            self.add_attribute("nout", 5)
            self.add_attribute("outstyle", 1)
            self.add_attribute("tfinal", 1.0)
            self.add_attribute("dt_initial", 1.0e-5)
            self.add_attribute("dt_max", 1.0e99)
            self.add_attribute("cfl_desired", 0.9)
            self.add_attribute("cfl_max", 1.0)
            self.add_attribute("max_steps", 5000)
            self.add_attribute("dt_variable", 1)
            self.add_attribute("order", 2)
            self.add_attribute("order_trans", 0)
            self.add_attribute("verbosity", 0)
            self.add_attribute("src_split", 0)
            self.add_attribute("mcapa", 0)
            self.add_attribute("maux", 0)
            self.add_attribute("meqn", 1)
            self.add_attribute("mwaves", 1)
            self.add_attribute("mthlim", [4])
            self.add_attribute("t0", 0.0)
            self.add_attribute("xlower", 0.0)
            self.add_attribute("xupper", 1.0)
            self.add_attribute("mbc", 2)
            self.add_attribute("mthbc_xlower", 1)
            self.add_attribute("mthbc_xupper", 1)
            self.add_attribute("restart", 0)
            self.add_attribute("N_restart", 0)

        elif ndim == 2:
            self.add_attribute("mx", 100)
            self.add_attribute("my", 100)
            self.add_attribute("nout", 5)
            self.add_attribute("outstyle", 1)
            self.add_attribute("tfinal", 1.0)
            self.add_attribute("dt_initial", 1.0e-5)
            self.add_attribute("dt_max", 1.0e99)
            self.add_attribute("cfl_desired", 0.9)
            self.add_attribute("cfl_max", 1.0)
            self.add_attribute("max_steps", 5000)
            self.add_attribute("dt_variable", 1)
            self.add_attribute("order", 2)
            self.add_attribute("order_trans", 2)
            self.add_attribute("verbosity", 0)
            self.add_attribute("src_split", 0)
            self.add_attribute("mcapa", 0)
            self.add_attribute("maux", 0)
            self.add_attribute("meqn", 1)
            self.add_attribute("mwaves", 1)
            self.add_attribute("mthlim", [4])
            self.add_attribute("t0", 0.0)
            self.add_attribute("xlower", 0.0)
            self.add_attribute("xupper", 1.0)
            self.add_attribute("ylower", 0.0)
            self.add_attribute("yupper", 1.0)
            self.add_attribute("mbc", 2)
            self.add_attribute("mthbc_xlower", 1)
            self.add_attribute("mthbc_xupper", 1)
            self.add_attribute("mthbc_ylower", 1)
            self.add_attribute("mthbc_yupper", 1)
            self.add_attribute("restart", 0)
            self.add_attribute("N_restart", 0)

        else:
            raise AttributeError("Only ndim=1 or 2 supported so far")

    def write(self):
        print("Creating data file claw.data for use with xclaw")
        make_clawdatafile(self)


class AmrclawInputData(Data):
    r"""
    Object that will be written out to amr2ez.data.
    """

    def __init__(self, ndim):
        super(AmrclawInputData, self).__init__()
        self.add_attribute("ndim", ndim)

        # Set default values:
        if ndim == 1:
            self.add_attribute("mx", 100)
            self.add_attribute("nout", 5)
            self.add_attribute("outstyle", 1)
            self.add_attribute("tfinal", 1.0)
            self.add_attribute("dt_initial", 1.0e-5)
            self.add_attribute("dt_max", 1.0e99)
            self.add_attribute("cfl_desired", 0.9)
            self.add_attribute("cfl_max", 1.0)
            self.add_attribute("max_steps", 5000)
            self.add_attribute("dt_variable", 1)
            self.add_attribute("order", 2)
            self.add_attribute("order_trans", 0)
            self.add_attribute("verbosity", 0)
            self.add_attribute("src_split", 0)
            self.add_attribute("mcapa", 0)
            self.add_attribute("maux", 0)
            self.add_attribute("meqn", 1)
            self.add_attribute("mwaves", 1)
            self.add_attribute("mthlim", [4])
            self.add_attribute("t0", 0.0)
            self.add_attribute("xlower", 0.0)
            self.add_attribute("xupper", 1.0)
            self.add_attribute("mbc", 2)
            self.add_attribute("mthbc_xlower", 1)
            self.add_attribute("mthbc_xupper", 1)
            self.add_attribute("restart", 0)
            self.add_attribute("N_restart", 0)

            # attributes need only since AMR is done using 2d amrclaw:
            self.add_attribute("my", 1)
            self.add_attribute("ylower", 0.0)
            self.add_attribute("yupper", 1.0)
            self.add_attribute("mthbc_ylower", 1)
            self.add_attribute("mthbc_yupper", 1)
            self.add_attribute("inraty", [1, 1, 1, 1, 1, 1])

        elif ndim == 2:
            self.add_attribute("mx", 100)
            self.add_attribute("my", 100)
            self.add_attribute("nout", 5)
            self.add_attribute("outstyle", 1)
            self.add_attribute("tfinal", 1.0)
            self.add_attribute("dt_initial", 1.0e-5)
            self.add_attribute("dt_max", 1.0e99)
            self.add_attribute("cfl_desired", 0.9)
            self.add_attribute("cfl_max", 1.0)
            self.add_attribute("max_steps", 5000)
            self.add_attribute("dt_variable", 1)
            self.add_attribute("order", 2)
            self.add_attribute("order_trans", 2)
            self.add_attribute("verbosity", 0)
            self.add_attribute("src_split", 0)
            self.add_attribute("mcapa", 0)
            self.add_attribute("maux", 0)
            self.add_attribute("meqn", 1)
            self.add_attribute("mwaves", 1)
            self.add_attribute("mthlim", [4])
            self.add_attribute("t0", 0.0)
            self.add_attribute("xlower", 0.0)
            self.add_attribute("xupper", 1.0)
            self.add_attribute("ylower", 0.0)
            self.add_attribute("yupper", 1.0)
            self.add_attribute("mbc", 2)
            self.add_attribute("mthbc_xlower", 1)
            self.add_attribute("mthbc_xupper", 1)
            self.add_attribute("mthbc_ylower", 1)
            self.add_attribute("mthbc_yupper", 1)
            self.add_attribute("restart", 0)
            self.add_attribute("N_restart", 0)
            self.add_attribute("inraty", [1])

        if ndim <= 2:
            # AMR parameters:
            self.add_attribute("mxnest", -1)
            self.add_attribute("inratx", [1])
            self.add_attribute("inratt", [1])
            self.add_attribute("auxtype", [])
            self.add_attribute("restart", False)
            self.add_attribute("checkpt_iousr", 1000)
            self.add_attribute("tchk", [])
            self.add_attribute("tol", -1.0)
            self.add_attribute("tolsp", 0.05)
            self.add_attribute("kcheck", 2)
            self.add_attribute("ibuff", 3)
            self.add_attribute("cutoff", 0.7)
            self.add_attribute("PRINT", False)
            self.add_attribute("NCAR", False)
            self.add_attribute("fortq", True)
            self.add_attribute("dprint", False)
            self.add_attribute("eprint", False)
            self.add_attribute("edebug", False)
            self.add_attribute("gprint", False)
            self.add_attribute("nprint", False)
            self.add_attribute("pprint", False)
            self.add_attribute("rprint", False)
            self.add_attribute("sprint", False)
            self.add_attribute("tprint", False)
            self.add_attribute("uprint", False)
        else:
            print("*** Error: only ndim=1 or 2 supported so far ***")
            raise AttributeError("Only ndim=1 or 2 supported so far")

    def write(self):
        print("Creating data file amr2ez.data for use with xamr")
        make_amrclawdatafile(self)
        make_setgauges_datafile(self)


class SharpclawInputData(Data):
    r"""
    Object that will be written out to claw.data.
    """

    def __init__(self, ndim):
        super(SharpclawInputData, self).__init__()
        self.add_attribute("ndim", ndim)

        # Set default values:
        if ndim == 1:
            self.add_attribute("mx", 100)
            self.add_attribute("nout", 5)
            self.add_attribute("outstyle", 1)
            self.add_attribute("tfinal", 1.0)
            self.add_attribute("dt_initial", 1.0e-5)
            self.add_attribute("dt_max", 1.0e99)
            self.add_attribute("cfl_desired", 0.9)
            self.add_attribute("cfl_max", 1.0)
            self.add_attribute("max_steps", 5000)
            self.add_attribute("dt_variable", 1)
            self.add_attribute("verbosity", 0)
            self.add_attribute("mcapa", 0)
            self.add_attribute("maux", 0)
            self.add_attribute("meqn", 1)
            self.add_attribute("mwaves", 1)
            self.add_attribute("mthlim", [5])
            self.add_attribute("t0", 0.0)
            self.add_attribute("xlower", 0.0)
            self.add_attribute("xupper", 1.0)
            self.add_attribute("mbc", 3)
            self.add_attribute("mthbc_xlower", 1)
            self.add_attribute("mthbc_xupper", 1)
            self.add_attribute("restart", 0)
            self.add_attribute("N_restart", 0)
            self.add_attribute("time_integrator", 2)
            self.add_attribute("tfluct_solver", 0)
            self.add_attribute("char_decomp", 0)
            self.add_attribute("lim_type", 2)
            self.add_attribute("src_term", 0)

        elif ndim == 2:
            self.add_attribute("mx", 100)
            self.add_attribute("my", 100)
            self.add_attribute("nout", 5)
            self.add_attribute("outstyle", 1)
            self.add_attribute("tfinal", 1.0)
            self.add_attribute("dt_initial", 1.0e-5)
            self.add_attribute("dt_max", 1.0e99)
            self.add_attribute("cfl_desired", 0.9)
            self.add_attribute("cfl_max", 1.0)
            self.add_attribute("max_steps", 5000)
            self.add_attribute("dt_variable", 1)
            self.add_attribute("verbosity", 0)
            self.add_attribute("mcapa", 0)
            self.add_attribute("maux", 0)
            self.add_attribute("meqn", 1)
            self.add_attribute("mwaves", 1)
            self.add_attribute("mthlim", [5])
            self.add_attribute("t0", 0.0)
            self.add_attribute("xlower", 0.0)
            self.add_attribute("xupper", 1.0)
            self.add_attribute("ylower", 0.0)
            self.add_attribute("yupper", 1.0)
            self.add_attribute("mbc", 3)
            self.add_attribute("mthbc_xlower", 1)
            self.add_attribute("mthbc_xupper", 1)
            self.add_attribute("mthbc_ylower", 1)
            self.add_attribute("mthbc_yupper", 1)
            self.add_attribute("restart", 0)
            self.add_attribute("N_restart", 0)
            self.add_attribute("time_integrator", 2)
            self.add_attribute("tfluct_solver", 0)
            self.add_attribute("char_decomp", 0)
            self.add_attribute("lim_type", 2)
            self.add_attribute("src_term", 0)

        else:
            raise AttributeError("Only ndim=1 or 2 supported so far")

    def write(self):
        print("Creating data file sharpclaw.data for use with xsclaw")
        make_sharpclawdatafile(self)


def open_datafile(name, datasource="setrun.py"):
    """
    Open a data file and write a warning header.
    Warning header starts with '#' character.  These lines are skipped if
    data file is opened using the library routine opendatafile.

    :Input:
     - *name* - (string) Name of data file
     - *datasource* - (string) Source for the data

    :Output:
     - (file) - file object
    """

    import string

    source = datasource.ljust(25)
    file = open(name, "w")
    file.write("########################################################\n")
    file.write("### DO NOT EDIT THIS FILE:  GENERATED AUTOMATICALLY ####\n")
    file.write("### To modify data, edit  %s ####\n" % source)
    file.write('###    and then "make .data"                        ####\n')
    file.write("########################################################\n\n")

    return file


def data_write(file, dataobj, name=None, descr=""):
    r"""
    Write out value to data file, in the form ::

       value =: name  descr

    Remove brackets and commas from lists, and replace booleans by T/F.
    Also convert numpy array to a list first.

    :Input:
     - *name* - (string) normally a string defining the variable,
       ``if name==None``, write a blank line.
     - *descr* - (string) A short description to appear on the line
    """

    import string

    if name is None:
        file.write("\n")
    else:
        try:
            value = getattr(dataobj, name)
        except:
            print(("Variable missing: ", name))
            print(("  from dataobj = ", dataobj))
            raise
        # Convert value to an appropriate string repr
        import numpy

        if isinstance(value, numpy.ndarray):
            value = list(value)
        if isinstance(value, tuple) | isinstance(value, list):
            # Remove [], (), and ','
            string_value = repr(value)[1:-1]
            string_value = string_value.replace(",", "")
        elif isinstance(value, bool):
            if value:
                string_value = "T"
            else:
                string_value = "F"
        else:
            string_value = repr(value)
        padded_value = string_value.ljust(25)
        padded_name = name.ljust(12)
        file.write("%s =: %s %s\n" % (padded_value, padded_name, descr))


def make_clawdatafile(clawdata):
    r"""
    Take the data specified in clawdata and write it to claw.data in the
    form required by the Fortran code lib/main.f95.
    """

    # open file and write a warning header:
    file = open_datafile("claw.data")

    ndim = clawdata.ndim
    data_write(file, clawdata, "ndim", "(number of dimensions)")
    data_write(file, clawdata, "mx", "(cells in x direction)")
    if ndim > 1:
        data_write(file, clawdata, "my", "(cells in y direction)")
    if ndim == 3:
        data_write(file, clawdata, "mz", "(cells in z direction)")
    data_write(file, clawdata, None)  # writes blank line

    data_write(file, clawdata, "nout", "(number of output times)")
    data_write(file, clawdata, "outstyle", "(style of specifying output times)")
    if clawdata.outstyle == 1:
        data_write(file, clawdata, "tfinal", "(final time)")
    elif clawdata.outstyle == 2:
        data_write(file, clawdata, "tout", "(output times)")
    elif clawdata.outstyle == 3:
        data_write(file, clawdata, "iout", "(output every iout steps)")
    elif clawdata.outstyle == 4:
        data_write(file, clawdata, "output_time_interval", "(between outputs)")
        data_write(file, clawdata, "tfinal", "(final time)")

    else:
        print("*** Error: unrecognized outstyle")
        raise
        return

    data_write(file, clawdata, None)
    data_write(file, clawdata, "dt_initial", "(initial time step dt)")
    data_write(file, clawdata, "dt_max", "(max allowable dt)")
    data_write(file, clawdata, "cfl_max", "(max allowable Courant number)")
    data_write(file, clawdata, "cfl_desired", "(desired Courant number)")
    data_write(file, clawdata, "max_steps", "(max time steps per call to claw)")
    data_write(file, clawdata, None)
    data_write(file, clawdata, "dt_variable", "(1 for variable dt, 0 for fixed)")
    data_write(file, clawdata, "order", "(1 or 2)")
    if ndim == 1:
        data_write(file, clawdata, "order_trans", "(not used in 1d)")
    else:
        data_write(file, clawdata, "order_trans", "(transverse order)")
    data_write(file, clawdata, "verbosity", "(verbosity of output)")
    data_write(file, clawdata, "src_split", "(source term splitting)")
    data_write(file, clawdata, "mcapa", "(aux index for capacity fcn)")
    data_write(file, clawdata, "maux", "(number of aux variables)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "meqn", "(number of equations)")
    data_write(file, clawdata, "mwaves", "(number of waves)")
    data_write(file, clawdata, "mthlim", "(limiter choice for each wave)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "t0", "(initial time)")
    data_write(file, clawdata, "xlower", "(xlower)")
    data_write(file, clawdata, "xupper", "(xupper)")
    if ndim > 1:
        data_write(file, clawdata, "ylower", "(ylower)")
        data_write(file, clawdata, "yupper", "(yupper)")
    if ndim == 3:
        data_write(file, clawdata, "zlower", "(zlower)")
        data_write(file, clawdata, "zupper", "(zupper)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "mbc", "(number of ghost cells)")
    data_write(file, clawdata, "mthbc_xlower", "(type of BC at xlower)")
    data_write(file, clawdata, "mthbc_xupper", "(type of BC at xupper)")
    if ndim > 1:
        data_write(file, clawdata, "mthbc_ylower", "(type of BC at ylower)")
        data_write(file, clawdata, "mthbc_yupper", "(type of BC at yupper)")
    if ndim == 3:
        data_write(file, clawdata, "mthbc_zlower", "(type of BC at zlower)")
        data_write(file, clawdata, "mthbc_zupper", "(type of BC at zupper)")

    data_write(file, clawdata, "restart", "(1 to restart from a past run)")
    data_write(file, clawdata, "N_restart", "(which frame to restart from)")
    data_write(file, clawdata, None)

    file.close()


def make_amrclawdatafile(clawdata):
    r"""
    Take the data specified in clawdata and write it to claw.data in the
    form required by the Fortran code lib/main.f95.
    """

    # open file and write a warning header:
    file = open_datafile("amr2ez.data")

    ndim = clawdata.ndim
    # data_write(file, clawdata, 'ndim', '(number of dimensions)')
    data_write(file, clawdata, "mx", "(cells in x direction)")
    data_write(file, clawdata, "my", "(cells in y direction)")
    if ndim == 3:
        data_write(file, clawdata, "mz", "(cells in z direction)")

    data_write(file, clawdata, "mxnest", "(max number of grid levels)")
    if len(clawdata.inratx) < max(abs(clawdata.mxnest) - 1, 1):
        raise ValueError(
            "*** Error in data parameter: require len(inratx) >= %s "
            % max(abs(clawdata.mxnest) - 1, 1)
        )
    data_write(file, clawdata, "inratx", "(refinement ratios)")
    if clawdata.mxnest < 0:
        # negative mxnest indicates anisotropic refinement
        if len(clawdata.inraty) < max(abs(clawdata.mxnest) - 1, 1):
            raise ValueError(
                "*** Error in data parameter: require len(inraty) >= %s "
                % max(abs(clawdata.mxnest) - 1, 1)
            )
        data_write(file, clawdata, "inraty", "(refinement ratios)")
        if ndim == 3:
            if len(clawdata.inratz) < max(abs(clawdata.mxnest) - 1, 1):
                raise ValueError(
                    "*** Error in data parameter: require len(inratz) >= %s "
                    % max(abs(clawdata.mxnest) - 1, 1)
                )
            data_write(file, clawdata, "inratz", "(refinement ratios)")
        if len(clawdata.inratt) < max(abs(clawdata.mxnest) - 1, 1):
            raise ValueError(
                "*** Error in data parameter: require len(inratt) >= %s "
                % max(abs(clawdata.mxnest) - 1, 1)
            )
        data_write(file, clawdata, "inratt", "(refinement ratios)")

    data_write(file, clawdata, None)  # writes blank line

    data_write(file, clawdata, "nout", "(number of output times)")
    data_write(file, clawdata, "outstyle", "(style of specifying output times)")
    if clawdata.outstyle == 1:
        data_write(file, clawdata, "tfinal", "(final time)")
    elif clawdata.outstyle == 2:
        data_write(file, clawdata, "tout", "(output times)")
    elif clawdata.outstyle == 3:
        data_write(file, clawdata, "iout", "(output every iout steps)")
    elif clawdata.outstyle == 4:
        data_write(file, clawdata, "output_time_interval", "(between outputs)")
        data_write(file, clawdata, "tfinal", "(final time)")
    else:
        print("*** Error: unrecognized outstyle")
        raise
        return

    data_write(file, clawdata, None)
    data_write(file, clawdata, "dt_initial", "(initial time step dt)")
    data_write(file, clawdata, "dt_max", "(max allowable dt)")
    data_write(file, clawdata, "cfl_max", "(max allowable Courant number)")
    data_write(file, clawdata, "cfl_desired", "(desired Courant number)")
    data_write(file, clawdata, "max_steps", "(max time steps per call to claw)")
    data_write(file, clawdata, None)
    data_write(file, clawdata, "dt_variable", "(1 for variable dt, 0 for fixed)")
    data_write(file, clawdata, "order", "(1 or 2)")
    if ndim == 1:
        data_write(file, clawdata, "order_trans", "(not used in 1d)")
    else:
        data_write(file, clawdata, "order_trans", "(transverse order)")
    data_write(file, clawdata, "verbosity", "(verbosity of output)")
    data_write(file, clawdata, "src_split", "(source term splitting)")
    data_write(file, clawdata, "mcapa", "(aux index for capacity fcn)")
    data_write(file, clawdata, "maux", "(number of aux variables)")
    if len(clawdata.auxtype) != clawdata.maux:
        file.close()
        print("*** Error: An auxtype array must be specified of length maux")
        raise AttributeError("require len(clawdata.auxtype) == clawdata.maux")
    for i in range(clawdata.maux):
        file.write("'%s'\n" % clawdata.auxtype[i])
    data_write(file, clawdata, None)

    data_write(file, clawdata, "meqn", "(number of equations)")
    data_write(file, clawdata, "mwaves", "(number of waves)")
    data_write(file, clawdata, "mthlim", "(limiter choice for each wave)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "t0", "(initial time)")
    data_write(file, clawdata, "xlower", "(xlower)")
    data_write(file, clawdata, "xupper", "(xupper)")
    data_write(file, clawdata, "ylower", "(ylower)")
    data_write(file, clawdata, "yupper", "(yupper)")
    if ndim == 3:
        data_write(file, clawdata, "zlower", "(zlower)")
        data_write(file, clawdata, "zupper", "(zupper)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "mbc", "(number of ghost cells)")
    data_write(file, clawdata, "mthbc_xlower", "(type of BC at xlower)")
    data_write(file, clawdata, "mthbc_xupper", "(type of BC at xupper)")
    data_write(file, clawdata, "mthbc_ylower", "(type of BC at ylower)")
    data_write(file, clawdata, "mthbc_yupper", "(type of BC at yupper)")
    if ndim == 3:
        data_write(file, clawdata, "mthbc_zlower", "(type of BC at zlower)")
        data_write(file, clawdata, "mthbc_zupper", "(type of BC at zupper)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "restart", "(1 to restart from a past run)")
    data_write(file, clawdata, "checkpt_iousr", "(how often to checkpoint)")
    if clawdata.checkpt_iousr < 0:
        data_write(file, clawdata, "tchk", "(checkpoint times)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "tol", "(tolerance for Richardson extrap)")
    data_write(file, clawdata, "tolsp", "(tolerance used in flag2refine)")
    data_write(file, clawdata, "kcheck", "(how often to regrid)")
    data_write(file, clawdata, "ibuff", "(buffer zone around flagged pts)")
    data_write(file, clawdata, "cutoff", "(efficiency cutoff for grid gen.)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "PRINT", "(print to fort.amr)")
    data_write(file, clawdata, "NCAR", "(obsolete!)")
    data_write(file, clawdata, "fortq", "(Output to fort.q* files)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "dprint", "(print domain flags)")
    data_write(file, clawdata, "eprint", "(print err est flags)")
    data_write(file, clawdata, "edebug", "(even more err est flags)")
    data_write(file, clawdata, "gprint", "(grid bisection/clustering)")
    data_write(file, clawdata, "nprint", "(proper nesting output)")
    data_write(file, clawdata, "pprint", "(proj. of tagged points)")
    data_write(file, clawdata, "rprint", "(print regridding summary)")
    data_write(file, clawdata, "sprint", "(space/memory output)")
    data_write(file, clawdata, "tprint", "(time step reporting each level)")
    data_write(file, clawdata, "uprint", "(update/upbnd reporting)")

    file.close()


def make_sharpclawdatafile(clawdata):
    r"""
    Take the data specified in clawdata and write it to sharpclaw.data in the
    form required by the Fortran code lib/main.f95.
    """

    # open file and write a warning header:
    file = open_datafile("sharpclaw.data")

    ndim = clawdata.ndim
    data_write(file, clawdata, "ndim", "(number of dimensions)")
    data_write(file, clawdata, "mx", "(cells in x direction)")
    if ndim > 1:
        data_write(file, clawdata, "my", "(cells in y direction)")
    if ndim == 3:
        data_write(file, clawdata, "mz", "(cells in z direction)")
    data_write(file, clawdata, None)  # writes blank line

    data_write(file, clawdata, "nout", "(number of output times)")
    data_write(file, clawdata, "outstyle", "(style of specifying output times)")
    if clawdata.outstyle == 1:
        data_write(file, clawdata, "tfinal", "(final time)")
    elif clawdata.outstyle == 2:
        data_write(file, clawdata, "tout", "(output times)")
    elif clawdata.outstyle == 3:
        data_write(file, clawdata, "iout", "(output every iout steps)")
    else:
        print("*** Error: unrecognized outstyle")
        raise
        return

    data_write(file, clawdata, None)
    data_write(file, clawdata, "dt_initial", "(initial time step dt)")
    data_write(file, clawdata, "dt_max", "(max allowable dt)")
    data_write(file, clawdata, "cfl_max", "(max allowable Courant number)")
    data_write(file, clawdata, "cfl_desired", "(desired Courant number)")
    data_write(file, clawdata, "max_steps", "(max time steps per call to claw)")
    data_write(file, clawdata, None)
    data_write(file, clawdata, "dt_variable", "(1 for variable dt, 0 for fixed)")
    data_write(file, clawdata, "time_integrator", "(time stepping scheme)")
    data_write(file, clawdata, "verbosity", "(verbosity of output)")
    data_write(file, clawdata, "src_term", "(source term present)")
    data_write(file, clawdata, "mcapa", "(aux index for capacity fcn)")
    data_write(file, clawdata, "maux", "(number of aux variables)")
    data_write(file, clawdata, "tfluct_solver", "(total fluctuation solver)")
    data_write(file, clawdata, "char_decomp", "(characteristic decomposition)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "meqn", "(number of equations)")
    data_write(file, clawdata, "mwaves", "(number of waves)")
    data_write(file, clawdata, "lim_type", "(0=None, 1=TVD, 2=WENO)")
    data_write(file, clawdata, "mthlim", "(limiter choice for each wave)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "t0", "(initial time)")
    data_write(file, clawdata, "xlower", "(xlower)")
    data_write(file, clawdata, "xupper", "(xupper)")
    if ndim > 1:
        data_write(file, clawdata, "ylower", "(ylower)")
        data_write(file, clawdata, "yupper", "(yupper)")
    if ndim == 3:
        data_write(file, clawdata, "zlower", "(zlower)")
        data_write(file, clawdata, "zupper", "(zupper)")
    data_write(file, clawdata, None)

    data_write(file, clawdata, "mbc", "(number of ghost cells)")
    data_write(file, clawdata, "mthbc_xlower", "(type of BC at xlower)")
    data_write(file, clawdata, "mthbc_xupper", "(type of BC at xupper)")
    if ndim > 1:
        data_write(file, clawdata, "mthbc_ylower", "(type of BC at ylower)")
        data_write(file, clawdata, "mthbc_yupper", "(type of BC at yupper)")
    if ndim == 3:
        data_write(file, clawdata, "mthbc_zlower", "(type of BC at zlower)")
        data_write(file, clawdata, "mthbc_zupper", "(type of BC at zupper)")

    data_write(file, clawdata, "restart", "(1 to restart from a past run)")
    data_write(file, clawdata, "N_restart", "(which frame to restart from)")
    data_write(file, clawdata, None)

    file.close()


def make_userdatafile(userdata):
    r"""
    Create the data file using the parameters in userdata.
    The parameters will be written to this file in the same order they were
    specified using userdata.add_attribute.
    Presumably the user will read these in using a Fortran routine, such as
    setprob.f95, and the order is important.
    """

    # open file and write a warning header:
    file = open_datafile(userdata._UserData__fname)

    # write all the parameters:
    for param in userdata.attributes:
        data_write(file, userdata, param, userdata._UserData__descr[param])

    file.close()


def make_setgauges_datafile(clawdata):
    """
    Create setgauges.data using gauges attribute of clawdata.
    """
    gauges = getattr(clawdata, "gauges", [])
    ngauges = len(gauges)

    print("Creating data file setgauges.data")
    # open file and write a warning header:
    file = open_datafile("setgauges.data")
    file.write("%4i   =: ngauges\n" % ngauges)
    gaugeno_used = []
    for gauge in gauges:
        gaugeno = gauge[0]
        if gaugeno in gaugeno_used:
            print(("*** Gauge number %s used more than once! " % gaugeno))
            raise Exception("Repeated gauge number")
        else:
            gaugeno_used.append(gauge[0])
        file.write("%4i %19.10e  %17.10e  %13.6e  %13.6e\n" % tuple(gauge))
        # or use this variant with =:
        # gauge.append(gaugeno)
        # file.write("%4i %19.10e  %17.10e  %13.6e  %13.6e  =: gauge%s\n" % tuple(gauge))
    file.close()


# -----------------------------------------------------
# New version 6/30/09


class ClawRunData(Data):
    r"""
    Object that will be written out to claw.data.
    """

    def __init__(self, pkg, ndim):
        super(ClawRunData, self).__init__()
        self.add_attribute("pkg", pkg)
        self.add_attribute("ndim", ndim)
        self.add_attribute("datalist", [])

        if pkg.lower() in ["classic", "classicclaw"]:
            self.add_attribute("xclawcmd", "xclaw")

            # Required data set for basic run parameters:
            clawdata = ClawInputData(ndim)
            self.add_attribute("clawdata", clawdata)
            self.datalist.append(clawdata)

        elif pkg.lower() in ["amrclaw", "amr"]:
            self.add_attribute("xclawcmd", "xamr")

            # Required data set for basic run parameters:
            clawdata = AmrclawInputData(ndim)
            self.add_attribute("clawdata", clawdata)
            self.datalist.append(clawdata)

        elif pkg.lower() in ["geoclaw"]:
            self.add_attribute("xclawcmd", "xgeoclaw")

            # Required data set for basic run parameters:
            clawdata = AmrclawInputData(ndim)
            self.add_attribute("clawdata", clawdata)
            self.datalist.append(clawdata)
            geodata = GeoclawInputData(ndim)
            self.add_attribute("geodata", geodata)
            self.datalist.append(geodata)

        elif pkg.lower() in ["digclaw"]:
            self.add_attribute("xclawcmd", "digclaw")

            # Required data set for basic run parameters:
            clawdata = AmrclawInputData(ndim)
            self.add_attribute("clawdata", clawdata)
            self.datalist.append(clawdata)
            geodata = GeoclawInputData(ndim)
            self.add_attribute("geodata", geodata)
            self.datalist.append(geodata)
            digdata = DigclawInputData(ndim)
            self.add_attribute("digdata", digdata)
            self.datalist.append(digdata)

        elif pkg.lower() in ["sharpclaw"]:
            self.add_attribute("xclawcmd", "xsclaw")

            # Required data set for basic run parameters:
            clawdata = SharpclawInputData(ndim)
            self.add_attribute("clawdata", clawdata)
            self.datalist.append(clawdata)

        else:
            raise AttributeError("Unrecognized Clawpack pkg = %s" % pkg)

    def new_UserData(self, name, fname):
        r"""
        Create a new attribute called name
        for application specific data to be written
        to the data file fname.
        """
        userdata = UserData(fname)
        self.datalist.append(userdata)
        exec("self.%s = userdata" % name)
        return userdata

    def add_GaugeData(self):
        r"""
        Create a gaugedata attribute for writing to gauges.data.
        """
        gaugedata = GaugeData(self.ndim)
        self.datalist.append(gaugedata)
        self.gaugedata = gaugedata
        return gaugedata

    def write(self):
        for d in self.datalist:
            d.write()


class UserData(Data):
    r"""
    Object that will be written out to user file such as setprob.data, as
    determined by the fname attribute.
    """

    def __init__(self, fname):
        super(UserData, self).__init__()
        self.__fname = fname  # file to be read by Fortran for this data
        self.__descr = {}  # dictionary to hold descriptions

    def add_param(self, name, value, descr=""):
        self.add_attribute(name, value)
        self.__descr[name] = descr

    def write(self):
        print(("Creating data file %s" % self.__fname))
        make_userdatafile(self)


class GaugeData(Data):
    r"""
    Data to be written out to gauge.data specifying gauges.
    DEPRECATED:  Use GeoclawInputData  instead.
    """

    def __init__(self, ndim):
        super(GaugeData, self).__init__()
        self.add_attribute("ndim", ndim)
        self.add_attribute("ngauges", 0)
        self.__gauge_dict = {}

    def add_gauge(self, gaugeno, location, time_interval):
        self.__gauge_dict[gaugeno] = (gaugeno, location, time_interval)
        self.ngauges = len(self.__gauge_dict)

    def write(self):
        print("Creating data file gauges.data")

        # open file and write a warning header:
        file = open_datafile("gauges.data")

        data_write(file, self, "ngauges", "Number of gauges")
        data_write(file, self, None)

        ndim = self.ndim

        # write a line for each gauge:
        for (gaugeno, gdata) in list(self.__gauge_dict.items()):
            tmin = gdata[2][0]
            tmax = gdata[2][1]
            if isinstance(gdata[1], (list, tuple)):
                xyz = gdata[1]
                x = xyz[0]
                if ndim > 1:
                    y = xyz[1]
                if ndim > 2:
                    z = xyz[2]
            else:
                x = gdata[1]

            if ndim == 1:
                file.write("%i   %e   %e   %e" % (gdata[0], x, tmin, tmax))
            elif ndim == 2:
                file.write("%i   %e   %e   %e   %e" % (gdata[0], x, y, tmin, tmax))
            elif ndim == 3:
                file.write(
                    "%i   %e   %e   %e   %e   %e" % (gdata[0], x, y, z, tmin, tmax)
                )

        printxyz = {1: "x  ", 2: "x  y  ", 3: "x  y  z"}
        file.write(
            "\n\n# Format of each line: \n#   gaugeno  %s tmin tmax" % printxyz[ndim]
        )
        file.close()


class GeoclawInputData(Data):
    r"""
    Object that will be written out to the various GeoClaw data files.
    """

    def __init__(self, ndim):
        super(GeoclawInputData, self).__init__()

        # Set default values:
        self.add_attribute("igravity", 1)
        self.add_attribute("iqinit", 0)
        self.add_attribute("icoriolis", 1)
        self.add_attribute("Rearth", 6367500.0)
        self.add_attribute("variable_dt_refinement_ratios", False)
        # NEED TO CONTINUE!

    def write(self):

        print("Creating data file setgeo.data")
        # open file and write a warning header:
        file = open_datafile("setgeo.data")
        data_write(file, self, "igravity")
        data_write(file, self, "gravity")
        data_write(file, self, "icoordsys")
        data_write(file, self, "icoriolis")
        data_write(file, self, "Rearth")
        data_write(file, self, "variable_dt_refinement_ratios")
        file.close()

        print("Creating data file settsunami.data")
        # open file and write a warning header:
        file = open_datafile("settsunami.data")
        data_write(file, self, "sealevel")
        data_write(file, self, "drytolerance")
        data_write(file, self, "wavetolerance")
        data_write(file, self, "depthdeep")
        data_write(file, self, "maxleveldeep")
        data_write(file, self, "ifriction")
        data_write(file, self, "coeffmanning")
        data_write(file, self, "frictiondepth")
        file.close()

        print("Creating data file settopo.data")
        # open file and write a warning header:
        file = open_datafile("settopo.data")
        self.ntopofiles = len(self.topofiles)
        data_write(file, self, "ntopofiles")
        for tfile in self.topofiles:
            try:
                fname = os.path.abspath(tfile[-1])
            except:
                print(("*** Error: file not valid string ", tfile[-1]))
                raise ValueError("file not valid string")
            # todo could also test for length of maximum string (as defined)
            if not os.path.exists(fname):
                print(("*** Error: file not found: ", tfile[-1]))
                raise ValueError("file not found")

            file.write("\n'%s' \n " % fname)
            file.write("%3i %3i %3i %20.10e %20.10e \n" % tuple(tfile[:-1]))
        file.close()

        print("Creating data file setdtopo.data")
        # open file and write a warning header:
        file = open_datafile("setdtopo.data")
        self.mdtopofiles = len(self.dtopofiles)
        data_write(file, self, "mdtopofiles")
        data_write(file, self, None)
        for tfile in self.dtopofiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                print(("*** Error: file not valid string ", tfile[-1]))
                raise ValueError("file not valid string")
            # todo could also test for length of maximum string (as defined)
            if not os.path.exists(tfile[-1]):
                print(("*** Error: file not found: ", tfile[-1]))
                raise ValueError("file not found")
            file.write("\n%s \n" % fname)
            file.write("%3i %3i %3i\n" % tuple(tfile[:-1]))
        file.close()

        print("Creating data file setqinit.data")
        # open file and write a warning header:
        file = open_datafile("setqinit.data")
        # self.iqinit tells which component of q is perturbed!
        self.nqinits = len(self.qinitfiles)
        data_write(file, self, "nqinits")
        data_write(file, self, None)
        for tfile in self.qinitfiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                print(("*** Error: file not valid string ", tfile[-1]))
                raise ValueError("file not valid string")
            # todo could also test for length of maximum string (as defined)
            if not os.path.exists(tfile[-1]):
                print(("*** Error: file not found: ", tfile[-1]))
                raise ValueError("file not found")
            file.write("\n%s  \n" % fname)
            file.write("%3i %3i %3i %3i \n" % tuple(tfile[:-1]))
        file.close()

        print("Creating data file setauxinit.data")
        # open file and write a warning header:
        file = open_datafile("setauxinit.data")
        # self.iauxinit tells which component of q is perturbed!
        self.nauxinits = len(self.auxinitfiles)
        data_write(file, self, "nauxinits")
        data_write(file, self, None)
        for tfile in self.auxinitfiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                print(("*** Error: file not valid string ", tfile[-1]))
                raise ValueError("file not valid string")
            # todo could also test for length of maximum string (as defined)
            if not os.path.exists(tfile[-1]):
                print(("*** Error: file not found: ", tfile[-1]))
                raise ValueError("file not found")
            file.write("\n%s  \n" % fname)
            file.write("%3i %3i %3i %3i \n" % tuple(tfile[:-1]))
        file.close()

        make_setgauges_datafile(self)

        #        print 'Creating data file setgauges.data'
        #        # open file and write a warning header:
        #        file = open_datafile('setgauges.data')
        #        self.ngauges = len(self.gauges)
        #        data_write(file, self, 'ngauges')
        #        data_write(file, self, None)
        #        gaugeno_used = []
        #        for gauge in self.gauges:
        #            gaugeno = gauge[0]
        #            if gaugeno in gaugeno_used:
        #                print "*** Gauge number %s used more than once! " % gaugeno
        #                raise Exception("Repeated gauge number")
        #            else:
        #                gaugeno_used.append(gauge[0])
        #            gauge.append(gaugeno)
        #            file.write("%3i %19.10e  %19.10e  %15.6e  %15.6e  =: gauge%s\n" % tuple(gauge))
        #        file.close()

        print("Creating data file setfixedgrids.data")
        # open file and write a warning header:
        file = open_datafile("setfixedgrids.data")
        self.nfixedgrids = len(self.fixedgrids)
        data_write(file, self, "nfixedgrids")
        data_write(file, self, None)
        for fixedgrid in self.fixedgrids:
            file.write(
                "%18.8e %18.8e %i %18.8e %18.8e %18.8e %18.8e %i %i %i %i"
                % tuple(fixedgrid)
                + "\n"
            )
        file.close()

        print("Creating data file setregions.data")
        # open file and write a warning header:
        file = open_datafile("setregions.data")
        self.nregions = len(self.regions)
        data_write(file, self, "nregions")
        data_write(file, self, None)
        for regions in self.regions:
            file.write(8 * "%g  " % tuple(regions) + "\n")
        file.close()

        print("Creating data file setflowgrades.data")
        # open file and write a warning header:
        file = open_datafile("setflowgrades.data")
        self.nflowgrades = len(self.flowgrades)
        data_write(file, self, "nflowgrades")
        data_write(file, self, None)
        for flowgrade in self.flowgrades:
            file.write(4 * "%g  " % tuple(flowgrade) + "\n")
        file.close()


class DigclawInputData(Data):
    r"""
    Object that will be written out to the various GeoClaw & data files.
    """

    def __init__(self, ndim):
        super(DigclawInputData, self).__init__()

        # Set default values:
        self.add_attribute("rho_s", 2700.0, _DIG_ATTRS["rho_s"])
        self.add_attribute("rho_f", 1000.0, _DIG_ATTRS["rho_f"])
        self.add_attribute("phi_bed", 40.0, _DIG_ATTRS["phi_bed"])
        self.add_attribute("theta_input", 0.0, _DIG_ATTRS["theta_input"])
        self.add_attribute("delta", 0.01, _DIG_ATTRS["delta"])
        self.add_attribute("kappita", 0.0001, _DIG_ATTRS["kappita"])
        self.add_attribute("mu", 0.001, _DIG_ATTRS["mu"])
        self.add_attribute("alpha_c", 1.0, _DIG_ATTRS["alpha_c"])
        self.add_attribute("m_crit", 0.62, _DIG_ATTRS["m_crit"])
        self.add_attribute("c1", 1.0, _DIG_ATTRS["c1"])
        self.add_attribute("m0", 0.52, _DIG_ATTRS["m0"])
        self.add_attribute("sigma_0", 1.0e3, _DIG_ATTRS["sigma_0"])
        self.add_attribute("alpha_seg", 0.0, _DIG_ATTRS["alpha_seg"])
        self.add_attribute("init_ptype", 0, _DIG_ATTRS["init_pmax_ratio"])
        self.add_attribute("init_pmax_ratio", 1.0, _DIG_ATTRS[""])
        self.add_attribute("init_ptf", 1.0, _DIG_ATTRS["init_ptf"])
        self.add_attribute("init_ptf2", 0.0, _DIG_ATTRS["init_ptf2"])
        self.add_attribute("bed_normal", 0, _DIG_ATTRS["bed_normal"])
        self.add_attribute("phi_seg_coeff", 0.0, _DIG_ATTRS["phi_seg_coeff"])
        self.add_attribute("entrainment", 0, _DIG_ATTRS["entrainment"])
        self.add_attribute("entrainment_rate", 0.2, _DIG_ATTRS["entrainment_rate"])
        self.add_attribute("mom_autostop", False, _DIG_ATTRS["mom_autostop"])
        self.add_attribute("mom_perc", 0.05, _DIG_ATTRS["mom_perc"])
        self.add_attribute("src_ftn_num_sr", 0, _DIG_ATTRS["src_ftn_num_sr"])
        self.add_attribute("fric_offset_val", 0.0, _DIG_ATTRS["fric_offset_val"])
        self.add_attribute("fric_star_val", 0.0, _DIG_ATTRS["fric_star_val"])
        self.add_attribute("chi_init_val", 0.0, _DIG_ATTRS["chi_init_val"])
        self.add_attribute("kappita_diff", 1.0, _DIG_ATTRS["kappita_diff"])
        self.add_attribute("outaux", False, _DIG_ATTRS["outaux"])
        # self.add_attribute('m_crit2', 0.62, 'critical state value of m (#) for different size species')
        # self.add_attribute('rho_s2', 2700.0, 'solid grain density (kg/m^3) different size species')
        # self.add_attribute('fric_offset_val2', 0.0, 'start/stop friction offset in degrees for different size species')
        # self.add_attribute('fric_star_val2', 0.0, 'deep friction offset in degrees for different size species')

    def write(self):

        print("Creating data file setdig.data")
        # open file and write a warning header:
        file = open_datafile("setdig.data")
        data_write(file, self, "rho_s", _DIG_ATTRS["rho_s"])
        data_write(file, self, "rho_f", _DIG_ATTRS["rho_f"])
        data_write(file, self, "phi_bed", _DIG_ATTRS["phi_bed"])
        data_write(file, self, "theta_input", _DIG_ATTRS["theta_input"])
        data_write(file, self, "delta", _DIG_ATTRS["delta"])
        data_write(file, self, "kappita", _DIG_ATTRS["kappita"])
        data_write(file, self, "mu", _DIG_ATTRS["mu"])
        data_write(file, self, "alpha_c", _DIG_ATTRS["alpha_c"])
        data_write(file, self, "m_crit", _DIG_ATTRS["m_crit"])
        data_write(file, self, "c1", _DIG_ATTRS["c1"])
        data_write(file, self, "m0", _DIG_ATTRS["m0"])
        data_write(file, self, "sigma_0", _DIG_ATTRS["sigma_0"])
        data_write(file, self, "alpha_seg", _DIG_ATTRS["alpha_seg"])
        data_write(file, self, "bed_normal", _DIG_ATTRS["bed_normal"])
        data_write(file, self, "phi_seg_coeff", _DIG_ATTRS["phi_seg_coeff"])
        data_write(file, self, "entrainment", _DIG_ATTRS["entrainment"])
        data_write(file, self, "entrainment_rate", _DIG_ATTRS["entrainment_rate"])
        data_write(file, self, "mom_autostop", _DIG_ATTRS["mom_autostop"])
        data_write(file, self, "mom_perc", _DIG_ATTRS["mom_perc"])
        data_write(file, self, "src_ftn_num_sr", _DIG_ATTRS["src_ftn_num_sr"])
        data_write(file, self, "fric_offset_val", _DIG_ATTRS["fric_offset_val"])
        data_write(file, self, "fric_star_val", _DIG_ATTRS["fric_star_val"])
        data_write(file, self, "chi_init_val", _DIG_ATTRS["chi_init_val"])
        data_write(file, self, "kappita_diff", _DIG_ATTRS["kappita_diff"])
        data_write(file, self, "outaux", _DIG_ATTRS["outaux"])

        # data_write(file, self, 'm_crit2', 'critical state value of m (#) for different size species')
        # data_write(file, self, 'rho_s2', 'solid grain density (kg/m^3) different size species')
        # data_write(file, self, 'fric_offset_val2', 'start/stop friction offset in degrees for different size species')
        # data_write(file, self, 'fric_star_val2', 'deep friction offset in degrees for different size species')

        file.close()

        print("Creating data file setpinit.data")
        # open file and write a warning header:
        file = open_datafile("setpinit.data")
        data_write(file, self, "init_ptype", _DIG_ATTRS["init_ptype"])
        data_write(file, self, "init_pmax_ratio", _DIG_ATTRS["init_pmax_ratio"])
        data_write(file, self, "init_ptf", _DIG_ATTRS["init_ptf"])
        data_write(file, self, "init_ptf2", _DIG_ATTRS["init_ptf2"])
        file.close()
