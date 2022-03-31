def fix_long_tick_labels(xlimits="auto", ylimits="auto", kwargs={}):
    """
    Utility to stop matplotlib from only displaying part of long
    tickmark labels.

    xlimits, ylimits can be set to lists of 2 elements each since
    limits have to be set before grabbing the ticks and reformatting.

    kwargs can be a dictionary of keyword arguments to be passed
    to xticks and yticks.  For example: kwargs={'fontsize' : 15}
    for larger labels.

    Problem with this way of doing it:  Ticks are fixed if you zoom in on plot.
    """

    from pylab import figure, xlim, xticks, ylim, yticks

    if xlimits != "auto":
        xlim(xlimits)
    if ylimits != "auto":
        ylim(ylimits)
    locs, labels = xticks()
    xticks(locs, [str(loc) for loc in locs], rotation=20, **kwargs)
    locs, labels = yticks()
    yticks(locs, [str(loc) for loc in locs], **kwargs)
