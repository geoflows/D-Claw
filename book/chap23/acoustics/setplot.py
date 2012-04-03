
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from mapc2p import mapc2p

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from pyclaw.plotters import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotdata.MappedGrid = True

    plotdata.mapc2p = mapc2p
    

    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='pressure', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'pressure'
    plotaxes.afteraxes = "pylab.axis('scaled')" 

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = -2
    plotitem.pcolor_cmax = 2
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    

    # Figure for u
    plotfigure = plotdata.new_plotfigure(name='u', figno=1)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'u'
    plotaxes.afteraxes = "pylab.axis('scaled')" 

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 1
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    

    # Figure for v
    plotfigure = plotdata.new_plotfigure(name='v', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'v'
    plotaxes.afteraxes = "pylab.axis('scaled')" 

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    
    # Figure for grid
    plotfigure = plotdata.new_plotfigure(name='grid', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'grid'
    plotaxes.afteraxes = "pylab.axis('scaled')" 

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_grid')


    # Scatter plot of q[0]
    plotfigure = plotdata.new_plotfigure(name='Scatter', figno=10)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0., 2.0]
    plotaxes.ylimits = [-4.0, 5.0]
    plotaxes.title = 'Scatter plot of pressure'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = 0
    def q_vs_radius(current_data):
        from numpy import sqrt
        xc = current_data.x
        yc = current_data.y
        x,y = mapc2p(xc,yc)
        r = sqrt(x**2 + y**2)
        q = current_data.q[:,:,0]
        return r,q
    plotitem.map_2d_to_1d = q_vs_radius
    plotitem.MappedGrid = False
    plotitem.plotstyle = 'o'

    # Plot the 1drad solution on scatter plot:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.MappedGrid = False
    plotitem.plot_var = 0
    import os
    plotitem.outdir = os.path.abspath('1drad/_output')
    plotitem.plotstyle = 'r-'


    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
