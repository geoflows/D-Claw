def mapc2p(xc,yc):
    """
    specifies the mapping to curvilinear coordinates -- should be consistent
    with mapc2p.f
    """
    from numpy import abs
    xp = xc + (abs(yc+.2)+ .8)/2
    yp = yc
    return xp,yp
