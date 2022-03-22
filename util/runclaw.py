"""
Generic code for running the fortran version of Clawpack and sending the
results to subdirectory output of the directory from which this is executed.
Execute via
    $ python $CLAW/python/pyclaw/runclaw.py
from a directory that contains a claw.data file and a Clawpack executable.
"""


def runclaw(xclawcmd=None, outdir=None, overwrite=True, restart=False, rundir=None):
    """
    Run the Fortran version of Clawpack using executable xclawcmd, which is
    typically set to 'xclaw', 'xamr', etc.

    If it is not set by the call, get it from the environment variable
    CLAW_EXE.  Default to 'xclaw' if that's not set.
    
    If rundir is None, all *.data is copied from current directory, if a path 
    is given, data files are copied from there instead.
    """

    import os

    if type(overwrite) is str:
        # convert to boolean
        overwrite = overwrite.lower() in ["true", "t"]

    if type(restart) is str:
        # convert to boolean
        restart = restart.lower() in ["true", "t"]

    # importing these modules requires $CLAW/python in PYTHONPATH:
    from pyclaw.controller import Controller

    if xclawcmd is None:
        # Determine what executable to use from environment variable CLAW_EXE
        # Default to 'xclaw' if it's not set:
        xclawcmd = os.environ.get("CLAW_EXE", "xclaw")

    if outdir is None:
        outdir = "."

    if rundir is None:
        rundir = os.getcwd()
    rundir = os.path.abspath(rundir)
    print("Will take data from ", rundir)

    # directory for fort.* files:
    outdir = os.path.abspath(outdir)
    print("== runclaw: Will write output to ", outdir)

    clawjob = Controller()
    clawjob.xdir = os.getcwd()
    clawjob.rundir = rundir  # use data files from current directory
    clawjob.outdir = outdir  # write fort files to outdir
    clawjob.xclawcmd = xclawcmd  # Clawpack executable
    clawjob.overwrite = overwrite  # Ok to overwrite outdir and plotdir?
    clawjob.restart = restart  # Restarting a previous run?

    returncode = clawjob.runxclaw()
    if returncode != 0:
        print("== runclaw: *** fortran returncode = ", returncode, "   aborting")
    print("== runclaw: Done executing %s via pyclaw.runclaw.py" % xclawcmd)
    print("== runclaw: Output is in ", outdir)


# ----------------------------------------------------------

if __name__ == "__main__":
    """
    If executed at command line prompt, simply call the function, with
    any argument used as setplot:
    """
    import sys

    args = sys.argv[1:]  # any command line arguments
    runclaw(*args)
