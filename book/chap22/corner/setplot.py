
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from pyclaw.plotters import colormaps
    from numpy import linspace

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    def plot_corner(current_data):
        from pylab import plot
        plot([.0,.0],[-1,0],'r',linewidth=2)
        plot([.0,1],[0,.55],'r',linewidth=2)

    def sigmatr(current_data):
        # trace of sigma
        q = current_data.q
        return q[:,:,0] + q[:,:,1]



    # Figure for trace of sigma
    plotfigure = plotdata.new_plotfigure(name='trace(sigma)', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'trace(sigma)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_corner

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = -1.
    plotitem.pcolor_cmax = 1.
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    


    # Figure for shear stress
    plotfigure = plotdata.new_plotfigure(name='shear', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'shear stress'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_corner

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 2  # sigma_12
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = -0.2
    plotitem.pcolor_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    

    # Figure for contours 
    plotfigure = plotdata.new_plotfigure(name='contours', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'pressure(black) and shear(green)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_corner

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = sigmatr
    plotitem.contour_levels = linspace(-2,8,50)
    plotitem.contour_colors = 'k'
    plotitem.show = True       # show on plot?

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 2  # sigma_12
    plotitem.contour_levels = linspace(-0.4,0.4,30)
    plotitem.contour_colors = 'g'
    plotitem.show = True       # show on plot?
    
    

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

    
