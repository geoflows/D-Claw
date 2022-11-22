"""
Generic code for plotting Clawpack results.

Execute from unix by
    $ python $CLAW/python/pyclaw/plotclaw.py
from a directory that contains a file setplot.py that sets other
parameters to specify the desired plots, location of data, etc.

If a different file is to be used to define the function setplot, this can
be given as an argument, e.g.
    $ python $CLAW/python/pyclaw/plotclaw.py setplot_alternative.py

From most Clawpack applications directories the command
    $ make .plots
will call the plotclaw function from this module.

"""

import os
import glob
import sys

import matplotlib

# Use the Agg backend for plotting -- this doesn't require opening windows
# when plotting remotely.


matplotlib.use("Agg")
import matplotlib.pyplot as plt

if sys.platform in ["win32", "cygwin"]:
    pypath = "C:/cygwin" + os.environ["CLAW"] + "/python"
    sys.path.append(pypath)


def plotclaw(outdir=".", plotdir="_plots", setplot="setplot.py"):
    """
    Create html and/or latex versions of plots.

    INPUT:
        setplot is a module containing a function setplot that will be called
                to set various plotting parameters.
    """

    from pyclaw.plotters import plotpages
    from pyclaw.plotters.data import ClawPlotData

    data_file_search = os.path.abspath(os.path.join(outdir, "*.data"))
    data_files = glob.glob(data_file_search)
    plotdata = ClawPlotData(data_files=data_files)
    plotdata.outdir = outdir
    plotdata.plotdir = plotdir
    plotdata.setplot = setplot

    plotpages.plotclaw_driver(plotdata, verbose=False)


# ----------------------------------------------------------

if __name__ == "__main__":
    """
    If executed at command line prompt, simply call the function, with
    any arguments passed in.
    """
    import sys

    if len(sys.argv) == 4:
        plotclaw(sys.argv[1], sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 3:
        plotclaw(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2:
        plotclaw(sys.argv[1])
    else:
        plotclaw()
