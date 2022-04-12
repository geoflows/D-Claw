"""
find log spiral given two points and two angles
call with xhi, zhi, xlo, zlo, alpha1, alpha2
"""

import numpy as np
import scipy.optimize as scio

beta1 = 1.0e0
beta2 = 1.0e0
deg2rad = np.pi / 180.0


def fitnorm(x, x1, y1, x2, y2, alpha1, alpha2):
    """
    least squares fit log spiral intersecting x1,y1,x2,y2
    with tangents alpha1 and alpha2 (angles with horizontal)
    ... -90.0 < alpha1 < 0.0 (scarp angle)
    ... -90.0 < alpha2 < 90.0 (toe angle)
    ...points will be on the same arm of the logspiral
    ...above contraint means tangents not guaranteed to be satisfied
    """
    import numpy as np

    pi = np.pi
    deg2rad = pi / 180.0
    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]
    theta1 = x[4]
    theta2 = x[5]

    f1x = (x1 - xc) - a * np.exp(b * theta1) * np.cos(theta1)
    f1y = (y1 - yc) - a * np.exp(b * theta1) * np.sin(theta1)

    f2x = (x2 - xc) - a * np.exp(b * theta2) * np.cos(theta2)
    f2y = (y2 - yc) - a * np.exp(b * theta2) * np.sin(theta2)

    psi1 = (alpha1) * deg2rad
    psi2 = (alpha2) * deg2rad

    falpha1 = beta1 * (
        np.cos(psi1) * (b * np.cos(theta1) - np.sin(theta1))
        + np.sin(psi1) * (b * np.sin(theta1) + np.cos(theta1))
        - np.sqrt(b ** 2 + 1.0)
    )

    falpha2 = beta2 * (
        np.cos(psi2) * (b * np.cos(theta2) - np.sin(theta2))
        + np.sin(psi2) * (b * np.sin(theta2) + np.cos(theta2))
        - np.sqrt(b ** 2 + 1.0)
    )

    f = f1x ** 2 + f1y ** 2 + f2x ** 2 + f2y ** 2 + falpha1 ** 2 + falpha2 ** 2

    g = np.ones(6)
    g[0] = -2.0 * f1x - 2.0 * f2x
    g[1] = -2.0 * f1y - 2.0 * f2y

    g[2] = (
        -2.0 * f1x * np.exp(b * theta1) * np.cos(theta1)
        - 2.0 * f1y * np.exp(b * theta1) * np.sin(theta1)
        - 2.0 * f2x * np.exp(b * theta2) * np.cos(theta2)
        - 2.0 * f2y * np.exp(b * theta2) * np.sin(theta2)
    )
    g[3] = (
        -2.0 * f1x * a * theta1 * np.exp(b * theta1) * np.cos(theta1)
        - 2.0 * f1y * a * theta1 * np.exp(b * theta1) * np.sin(theta1)
        - 2.0 * f2x * a * theta2 * np.exp(b * theta2) * np.cos(theta2)
        - 2.0 * f2y * a * theta2 * np.exp(b * theta2) * np.sin(theta2)
        + 2.0
        * falpha1
        * (np.cos(psi1) * np.cos(theta1) + np.sin(psi1) * np.sin(theta1))
        - (b / np.sqrt(b ** 2 + 1.0))
        + 2.0
        * falpha2
        * (np.cos(psi2) * np.cos(theta2) + np.sin(psi2) * np.sin(theta2))
        - (b / np.sqrt(b ** 2 + 1.0))
    )
    g[4] = (
        2.0
        * f1x
        * (
            -a * b * np.exp(b * theta1) * np.cos(theta1)
            + a * np.exp(b * theta1) * np.sin(theta1)
        )
        + 2.0
        * f1y
        * (
            -a * b * np.exp(b * theta1) * np.sin(theta1)
            - a * np.exp(b * theta1) * np.cos(theta1)
        )
        + 2.0
        * falpha1
        * (
            np.cos(psi1) * (-b * np.sin(theta1) - np.cos(theta1))
            + np.sin(psi1) * (b * np.cos(theta1) - np.sin(theta1))
        )
    )

    g[5] = (
        2.0
        * f2x
        * (
            -a * b * np.exp(b * theta2) * np.cos(theta2)
            + a * np.exp(b * theta2) * np.sin(theta2)
        )
        + 2.0
        * f2y
        * (
            -a * b * np.exp(b * theta2) * np.sin(theta2)
            - a * np.exp(b * theta2) * np.cos(theta2)
        )
        + 2.0
        * falpha2
        * (
            np.cos(psi2) * (-b * np.sin(theta2) - np.cos(theta2))
            + np.sin(psi2) * (b * np.cos(theta2) - np.sin(theta2))
        )
    )

    # return (f,g)

    return f


def fitprime(x, x1, y1, x2, y2, alpha1, alpha2):
    """
    least squares fit log spiral intersecting x1,y1,x2,y2
    with tangents alpha1 and alpha2 (angles with horizontal)
    ... -90.0 < alpha1 < 0.0 (scarp angle)
    ... -90.0 < alpha2 < 90.0 (toe angle)
    ...points will be on the same arm of the logspiral
    ...above contraint means tangents not guaranteed to be satisfied
    """
    import numpy as np

    pi = np.pi
    deg2rad = pi / 180.0
    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]
    theta1 = x[4]
    theta2 = x[5]

    f1x = (x1 - xc) - a * np.exp(b * theta1) * np.cos(theta1)
    f1y = (y1 - yc) - a * np.exp(b * theta1) * np.sin(theta1)

    f2x = (x2 - xc) - a * np.exp(b * theta2) * np.cos(theta2)
    f2y = (y2 - yc) - a * np.exp(b * theta2) * np.sin(theta2)

    psi1 = (alpha1) * deg2rad
    psi2 = (alpha2) * deg2rad

    falpha1 = beta1 * (
        np.cos(psi1) * (b * np.cos(theta1) - np.sin(theta1))
        + np.sin(psi1) * (b * np.sin(theta1) + np.cos(theta1))
        - np.sqrt(b ** 2 + 1.0)
    )

    falpha2 = beta2 * (
        np.cos(psi2) * (b * np.cos(theta2) - np.sin(theta2))
        + np.sin(psi2) * (b * np.sin(theta2) + np.cos(theta2))
        - np.sqrt(b ** 2 + 1.0)
    )

    f = f1x ** 2 + f1y ** 2 + f2x ** 2 + f2y ** 2 + falpha1 ** 2 + falpha2 ** 2

    g = np.ones(6)
    g[0] = -2.0 * f1x - 2.0 * f2x
    g[1] = -2.0 * f1y - 2.0 * f2y

    g[2] = (
        -2.0 * f1x * np.exp(b * theta1) * np.cos(theta1)
        - 2.0 * f1y * np.exp(b * theta1) * np.sin(theta1)
        - 2.0 * f2x * np.exp(b * theta2) * np.cos(theta2)
        - 2.0 * f2y * np.exp(b * theta2) * np.sin(theta2)
    )
    g[3] = (
        -2.0 * f1x * a * theta1 * np.exp(b * theta1) * np.cos(theta1)
        - 2.0 * f1y * a * theta1 * np.exp(b * theta1) * np.sin(theta1)
        - 2.0 * f2x * a * theta2 * np.exp(b * theta2) * np.cos(theta2)
        - 2.0 * f2y * a * theta2 * np.exp(b * theta2) * np.sin(theta2)
        + 2.0
        * falpha1
        * (np.cos(psi1) * np.cos(theta1) + np.sin(psi1) * np.sin(theta1))
        - (b / np.sqrt(b ** 2 + 1.0))
        + 2.0
        * falpha2
        * (np.cos(psi2) * np.cos(theta2) + np.sin(psi2) * np.sin(theta2))
        - (b / np.sqrt(b ** 2 + 1.0))
    )
    g[4] = (
        2.0
        * f1x
        * (
            -a * b * np.exp(b * theta1) * np.cos(theta1)
            + a * np.exp(b * theta1) * np.sin(theta1)
        )
        + 2.0
        * f1y
        * (
            -a * b * np.exp(b * theta1) * np.sin(theta1)
            - a * np.exp(b * theta1) * np.cos(theta1)
        )
        + 2.0
        * falpha1
        * (
            np.cos(psi1) * (-b * np.sin(theta1) - np.cos(theta1))
            + np.sin(psi1) * (b * np.cos(theta1) - np.sin(theta1))
        )
    )

    g[5] = (
        2.0
        * f2x
        * (
            -a * b * np.exp(b * theta2) * np.cos(theta2)
            + a * np.exp(b * theta2) * np.sin(theta2)
        )
        + 2.0
        * f2y
        * (
            -a * b * np.exp(b * theta2) * np.sin(theta2)
            - a * np.exp(b * theta2) * np.cos(theta2)
        )
        + 2.0
        * falpha2
        * (
            np.cos(psi2) * (-b * np.sin(theta2) - np.cos(theta2))
            + np.sin(psi2) * (b * np.cos(theta2) - np.sin(theta2))
        )
    )

    # return (f,g)

    return g


def spiral_def(x1, y1, x2, y2, alpha1, alpha2, npoints=100):
    """
    return parameters to uniquely define the logspiral (2d ...x along a transect)
    """
    import numpy as np

    deg2rad = np.pi / 180.0
    psi1 = (alpha1) * deg2rad
    psi2 = (alpha2) * deg2rad

    # initial guess
    xc0 = 0.5 * (x1 + x2)
    yc0 = y1 + 0.5 * (y1 - y2)
    r1 = np.sqrt((x1 - xc0) ** 2 + (y1 - yc0) ** 2)
    r2 = np.sqrt((x2 - xc0) ** 2 + (y2 - yc0) ** 2)
    a0 = 0.5 * (r1 + r2)
    b0 = 0.1
    theta1 = 1.5 * np.pi + psi1 + 2.0 * np.pi  # np.pi
    theta2 = 1.5 * np.pi + psi2 + 2.0 * np.pi

    x0 = np.ones(6)
    x0[0] = xc0
    x0[1] = yc0
    x0[2] = a0
    x0[3] = b0
    x0[4] = theta1
    x0[5] = theta2

    bounds = [(x1, None)]
    bounds.append((y2, None))
    bounds.append((0.0, None))
    bounds.append((0.0, None))
    bounds.append((0.0, None))
    bounds.append((0.0, None))
    # bounds.append((0.8*np.pi ,1.5*np.pi))
    # bounds.append((1.2*np.pi ,1.6*np.pi))

    (x, j1, j2) = scio.fmin_tnc(
        fitnorm,
        x0,
        fprime=fitprime,
        args=(x1, y1, x2, y2, alpha1, alpha2),
        approx_grad=False,
        bounds=bounds,
        maxfun=10000000,
    )
    # (x,j1,j2) = scio.fmin_l_bfgs_b(fitnorm,x0,fprime=fitprime,args=(x1,y1,x2,y2,alpha1,alpha2),approx_grad=False,bounds=bounds)
    # (x,j1,j2) = scio.fmin_tnc(fitnorm,x0,fprime=None,args=(x1,y1,x2,y2,alpha1,alpha2),approx_grad=True,bounds=bounds)
    # (x,j1,j2) = scio.differential_evolution(fitnorm,bounds=bounds,args=(x1,y1,x2,y2,alpha1,alpha2))

    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]
    theta1 = x[4]
    theta2 = x[5]
    # import pdb;pdb.set_trace()

    theta = np.linspace(theta1, theta2, npoints)
    theta = np.linspace(0.0, 4.0 * np.pi, npoints)
    xv = xc + a * np.exp(b * theta) * np.cos(theta)
    yv = yc + a * np.exp(b * theta) * np.sin(theta)

    import matplotlib.pyplot as plt

    plt.plot(
        xv,
        yv,
        [x1, x2],
        [y1, y2],
        "ro",
        [x1, x1 + 20 * np.cos(psi1)],
        [y1, y1 + 20 * np.sin(psi1)],
        "k-",
        [x2, x2 + 20.0 * np.cos(psi2)],
        [y2, y2 + 20.0 * np.sin(psi2)],
        "k-",
    )
    plt.plot(xc, yc, "ko")
    plt.axis("equal")
    plt.show()

    return x


def spiral(x1, y1, x2, y2, alpha1, alpha2, npoints=100):
    """
    return parameters to uniquely define the logspiral (2d ...x along a transect)
    """
    import numpy as np

    # initial guess
    xc0 = 0.5 * (x1 + x2)
    yc0 = y1 + 0.5 * (y1 - y2)
    r1 = np.sqrt((x1 - xc0) ** 2 + (y1 - yc0) ** 2)
    r2 = np.sqrt((x2 - xc0) ** 2 + (y2 - yc0) ** 2)
    a0 = 0.5 * (r1 + r2)
    b0 = 0.1
    theta1 = np.pi
    theta2 = 1.5 * np.pi

    x0 = np.ones(6)
    x0[0] = xc0
    x0[1] = yc0
    x0[2] = a0
    x0[3] = b0
    x0[4] = theta1
    x0[5] = theta2

    bounds = [(x1, None)]
    bounds.append((y2, None))
    bounds.append((0.0, None))
    bounds.append((0.0, None))
    bounds.append((0.8 * np.pi, 1.5 * np.pi))
    bounds.append((1.2 * np.pi, 1.5 * np.pi))

    deg2rad = np.pi / 180.0
    psi1 = (alpha1) * deg2rad
    psi2 = (alpha2) * deg2rad

    (x, j1, j2) = scio.fmin_tnc(
        fitnorm,
        x0,
        fprime=fitprime,
        args=(x1, y1, x2, y2, alpha1, alpha2),
        approx_grad=False,
        bounds=bounds,
    )
    (x, j1, j2) = scio.fmin_l_bfgs_b(
        fitnorm,
        x0,
        fprime=fitprime,
        args=(x1, y1, x2, y2, alpha1, alpha2),
        approx_grad=False,
        bounds=bounds,
    )

    # (x,j1,j2) = scio.differential_evolution(fitnorm,bounds=bounds,args=(x1,y1,x2,y2,alpha1,alpha2))

    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]
    theta1 = x[4]
    theta2 = x[5]

    theta = np.linspace(theta1, theta2, npoints)
    # theta = np.linspace(0.0,2.*np.pi)
    xv = xc + a * np.exp(b * theta) * np.cos(theta)
    yv = yc + a * np.exp(b * theta) * np.sin(theta)

    # import matplotlib.pyplot as plt
    # plt.plot(xv,yv,[x1,x2],[y1,y2],'ro',[x1,x1+20*np.cos(psi1)],[y1,y1+20*np.sin(psi1)],'k-',[x2,x2+20.*np.cos(psi2)],[y2,y2+20.*np.sin(psi2)],'k-')

    # plt.show()

    return (xv, yv)


def test_grad(x1, y1, x2, y2, alpha1, alpha2):

    import numpy as np

    # initial guess
    xc0 = 0.5 * (x1 + x2)
    yc0 = y1 + 0.5 * (y1 - y2)
    r1 = np.sqrt((x1 - xc0) ** 2 + (y1 - yc0) ** 2)
    r2 = np.sqrt((x2 - xc0) ** 2 + (y2 - yc0) ** 2)
    a0 = 0.5 * (r1 + r2)
    b0 = 0.1
    theta1 = np.pi
    theta2 = 1.5 * np.pi

    x0 = np.ones(4)
    x0[0] = xc0
    x0[1] = yc0
    x0[2] = a0
    x0[3] = b0
    # x0[4] = theta1
    # x0[5] = theta2

    # args={'x1':x1,'y1':y1,'x2':y2,'alpha1':alpha1,'alpha2':alpha2}
    xargs = (x1, y1, x2, y2, alpha1, alpha2)
    err = scio.check_grad(
        objectf_4x4, objectfprime_4x4, x0, x1, y1, x2, y2, alpha1, alpha2
    )

    print(err)


def objectf_4x4_scalar(x, x1, y1, x2, y2, alpha1, alpha2):

    """
    least squares fit log spiral intersecting x1,y1,x2,y2
    with tangents alpha1 and alpha2 (angles with horizontal)
    ... -90.0 < alpha1 < 0.0 (scarp angle)
    ... -90.0 < alpha2 < 90.0 (toe angle)
    ...points will be on the same arm of the logspiral
    ...
    """

    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]

    psi1 = (alpha1) * deg2rad + 2.0 * np.pi
    psi2 = (alpha2) * deg2rad + 2.0 * np.pi

    theta1 = np.arctan2(y1 - yc, x1 - xc) + 2.0 * np.pi
    theta2 = np.arctan2(y2 - yc, x2 - xc) + 2.0 * np.pi

    f1 = (
        (x1 - xc) ** 2.0
        + (y1 - yc) ** 2.0
        - a * np.exp(b * theta1) * a * np.exp(b * theta1)
    )
    f2 = (
        (x2 - xc) ** 2.0
        + (y2 - yc) ** 2.0
        - a * np.exp(b * theta2) * a * np.exp(b * theta2)
    )

    falpha1 = -np.sin(psi1) * (b * np.cos(theta1) - np.sin(theta1)) + np.cos(psi1) * (
        b * np.sin(theta1) + np.cos(theta1)
    )

    falpha2 = -np.sin(psi2) * (b * np.cos(theta2) - np.sin(theta2)) + np.cos(psi2) * (
        b * np.sin(theta2) + np.cos(theta2)
    )

    linelength = np.sqrt((y2 - y1) ** 2.0 + (x2 - x1) ** 2.0)
    alength = (
        a * np.sqrt(1.0 + b ** 2.0) * (np.exp(b * theta2) - np.exp(b * theta1)) / b
        - 1.0 * linelength
    )

    f = f1 ** 2.0 + f2 ** 2.0 + falpha1 ** 2.0 + falpha2 ** 2.0

    return f


def objectfprime_4x4_scalar(x, x1, y1, x2, y2, alpha1, alpha2, npoints=100):

    """
    gradient
    """

    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]

    psi1 = (alpha1) * deg2rad + 2.0 * np.pi
    psi2 = (alpha2) * deg2rad + 2.0 * np.pi

    theta1 = np.arctan2(y1 - yc, x1 - xc) + 2.0 * np.pi
    theta2 = np.arctan2(y2 - yc, x2 - xc) + 2.0 * np.pi

    f1 = (
        (x1 - xc) ** 2.0
        + (y1 - yc) ** 2.0
        - a * np.exp(b * theta1) * a * np.exp(b * theta1)
    )
    f2 = (
        (x2 - xc) ** 2.0
        + (y2 - yc) ** 2.0
        - a * np.exp(b * theta2) * a * np.exp(b * theta2)
    )

    falpha1 = np.sin(psi1) * (b * np.cos(theta1) - np.sin(theta1)) - np.cos(psi1) * (
        b * np.sin(theta1) + np.cos(theta1)
    )

    falpha2 = np.sin(psi2) * (b * np.cos(theta2) - np.sin(theta2)) - np.cos(psi2) * (
        b * np.sin(theta2) + np.cos(theta2)
    )

    dtanrat1_dxc = (y1 - yc) / (x1 - xc) ** 2
    dtanrat2_dxc = (y2 - yc) / (x2 - xc) ** 2

    dtanrat1_dyc = -1.0 / (x1 - xc)
    dtanrat2_dyc = -1.0 / (x2 - xc)

    dtheta1_dxc = dtanrat1_dxc / (1.0 + ((y1 - yc) / (x1 - xc)) ** 2)
    dtheta2_dxc = dtanrat2_dxc / (1.0 + ((y2 - yc) / (x2 - xc)) ** 2)

    dtheta1_dyc = dtanrat1_dyc / (1.0 + ((y1 - yc) / (x1 - xc)) ** 2)
    dtheta2_dyc = dtanrat2_dyc / (1.0 + ((y2 - yc) / (x2 - xc)) ** 2)

    df1_dxc = (
        -2.0 * (x1 - xc) - a ** 2.0 * 2.0 * b * np.exp(2.0 * b * theta1) * dtheta1_dxc
    )
    df1_dyc = (
        -2.0 * (y1 - yc) - a ** 2.0 * 2.0 * b * np.exp(2.0 * b * theta1) * dtheta1_dyc
    )
    df1_da = -2.0 * a * np.exp(2.0 * b * theta1)
    df1_db = -2.0 * a ** 2 * theta1 * np.exp(2.0 * b * theta1)

    df2_dxc = (
        -2.0 * (x2 - xc) - a ** 2.0 * 2.0 * b * np.exp(2.0 * b * theta2) * dtheta2_dxc
    )
    df2_dyc = (
        -2.0 * (y2 - yc) - a ** 2.0 * 2.0 * b * np.exp(2.0 * b * theta2) * dtheta2_dyc
    )
    df2_da = -2.0 * a * np.exp(2.0 * b * theta2)
    df2_db = -2.0 * a ** 2 * theta2 * np.exp(2.0 * b * theta2)

    dalpha1_dxc = dtheta1_dxc * (
        -np.sin(psi1) * (b * np.sin(theta1) + np.cos(theta1))
        - np.cos(psi1) * (b * np.cos(theta1) - np.sin(theta1))
    )
    dalpha1_dyc = dtheta1_dyc * (
        -np.sin(psi1) * (b * np.sin(theta1) + np.cos(theta1))
        - np.cos(psi1) * (b * np.cos(theta1) - np.sin(theta1))
    )

    dalpha1_db = np.cos(psi1) * np.sin(theta1) - np.sin(psi1) * np.cos(theta1)
    dalpha1_da = 0.0

    dalpha2_dxc = dtheta1_dxc * (
        -np.sin(psi1) * (b * np.sin(theta1) + np.cos(theta1))
        - np.cos(psi1) * (b * np.cos(theta1) - np.sin(theta1))
    )
    dalpha2_dyc = dtheta1_dyc * (
        -np.sin(psi1) * (b * np.sin(theta1) + np.cos(theta1))
        - np.cos(psi1) * (b * np.cos(theta1) - np.sin(theta1))
    )
    dalpha2_db = np.cos(psi2) * np.sin(theta2) - np.sin(psi2) * np.cos(theta2)
    dalpha2_da = 0.0

    g = np.ones(4)

    g[0] = (
        2.0 * f1 * df1_dxc
        + 2.0 * f2 * df2_dxc
        + 2.0 * falpha1 * dalpha1_dxc
        + 2.0 * falpha2 * dalpha2_dxc
    )
    g[1] = (
        2.0 * f1 * df1_dyc
        + 2.0 * f2 * df2_dyc
        + 2.0 * falpha1 * dalpha1_dyc
        + 2.0 * falpha2 * dalpha2_dyc
    )
    g[2] = (
        2.0 * f1 * df1_da
        + 2.0 * f2 * df2_da
        + 2.0 * falpha1 * dalpha1_da
        + 2.0 * falpha2 * dalpha2_da
    )
    g[3] = (
        2.0 * f1 * df1_db
        + 2.0 * f2 * df2_db
        + 2.0 * falpha1 * dalpha1_db
        + 2.0 * falpha2 * dalpha2_db
    )

    return g


def objectf_4x4(x, x1, y1, x2, y2, alpha1, alpha2):

    """
    least squares fit log spiral intersecting x1,y1,x2,y2
    with tangents alpha1 and alpha2 (angles with horizontal)
    ... -90.0 < alpha1 < 0.0 (scarp angle)
    ... -90.0 < alpha2 < 90.0 (toe angle)
    ...points will be on the same arm of the logspiral
    ...
    """

    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]

    psi1 = (alpha1) * deg2rad + 2.0 * np.pi
    psi2 = (alpha2) * deg2rad + 2.0 * np.pi

    theta1 = np.arctan2(y1 - yc, x1 - xc) + 2.0 * np.pi
    theta2 = np.arctan2(y2 - yc, x2 - xc) + 2.0 * np.pi

    f1 = np.sqrt((x1 - xc) ** 2.0 + (y1 - yc) ** 2.0) - a * np.exp(
        b * theta1
    )  # *a*np.exp(b*theta1)
    f2 = np.sqrt((x2 - xc) ** 2.0 + (y2 - yc) ** 2.0) - a * np.exp(
        b * theta2
    )  # *a*np.exp(b*theta2)

    falpha1 = -np.sin(psi1) * (b * np.cos(theta1) - np.sin(theta1)) + np.cos(psi1) * (
        b * np.sin(theta1) + np.cos(theta1)
    )

    falpha2 = -np.sin(psi2) * (b * np.cos(theta2) - np.sin(theta2)) + np.cos(psi2) * (
        b * np.sin(theta2) + np.cos(theta2)
    )

    f = np.ones(4)
    f[0] = f1
    f[1] = f2
    f[2] = falpha1
    f[3] = falpha2

    return f


def objectD_4x4(x, x1, y1, x2, y2, alpha1, alpha2):

    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]

    psi1 = (alpha1) * deg2rad + 2.0 * np.pi
    psi2 = (alpha2) * deg2rad + 2.0 * np.pi

    theta1 = np.arctan2(y1 - yc, x1 - xc) + 2.0 * np.pi
    theta2 = np.arctan2(y2 - yc, x2 - xc) + 2.0 * np.pi

    rs1 = np.sqrt((x1 - xc) ** 2.0 + (y1 - yc) ** 2.0)
    rs2 = np.sqrt((x2 - xc) ** 2.0 + (y2 - yc) ** 2.0)
    r1 = a * np.exp(b * theta1)
    r2 = a * np.exp(b * theta2)

    dtanrat1_dxc = (y1 - yc) / (x1 - xc) ** 2
    dtanrat2_dxc = (y2 - yc) / (x2 - xc) ** 2

    dtanrat1_dyc = -1.0 / (x1 - xc)
    dtanrat2_dyc = -1.0 / (x2 - xc)

    dtheta1_dxc = dtanrat1_dxc / (1.0 + ((y1 - yc) / (x1 - xc)) ** 2)
    dtheta2_dxc = dtanrat2_dxc / (1.0 + ((y2 - yc) / (x2 - xc)) ** 2)

    dtheta1_dyc = dtanrat1_dyc / (1.0 + ((y1 - yc) / (x1 - xc)) ** 2)
    dtheta2_dyc = dtanrat2_dyc / (1.0 + ((y2 - yc) / (x2 - xc)) ** 2)

    dalpha_dtheta1 = np.sin(psi1) * (b * np.sin(theta1) + np.cos(theta1)) + np.cos(
        psi1
    ) * (b * np.cos(theta1) - np.sin(theta1))
    dalpha_dtheta2 = np.sin(psi2) * (b * np.sin(theta2) + np.cos(theta2)) + np.cos(
        psi2
    ) * (b * np.cos(theta2) - np.sin(theta2))

    dalpha_dtheta2

    df1_dxc = -(x1 - xc) / rs1 - b * r1 * dtheta1_dxc
    df1_dyc = -(y1 - yc) / rs1 - b * r1 * dtheta1_dyc
    df1_da = -np.exp(b * theta1)
    df1_db = -theta1 * r1

    df2_dxc = -(x2 - xc) / rs2 - b * r2 * dtheta2_dxc
    df2_dyc = -(y2 - yc) / rs2 - b * r2 * dtheta2_dyc
    df2_da = -np.exp(b * theta2)
    df2_db = -theta2 * r2

    dalpha1_dxc = dalpha_dtheta1 * dtheta1_dxc
    dalpha1_dyc = dalpha_dtheta1 * dtheta1_dyc
    dalpha1_db = np.cos(psi1) * np.sin(theta1) - np.sin(psi1) * np.cos(theta1)
    dalpha1_da = 0.0

    dalpha2_dxc = dalpha_dtheta2 * dtheta2_dxc
    dalpha2_dyc = dalpha_dtheta2 * dtheta2_dyc
    dalpha2_db = np.cos(psi2) * np.sin(theta2) - np.sin(psi2) * np.cos(theta2)
    dalpha2_da = 0.0

    D = np.ones((4, 4))

    D[0, 0] = df1_dxc
    D[0, 1] = df1_dyc
    D[0, 2] = df1_da
    D[0, 3] = df1_db

    D[1, 0] = df2_dxc
    D[1, 1] = df2_dyc
    D[1, 2] = df2_da
    D[1, 3] = df2_db

    D[2, 0] = dalpha1_dxc
    D[2, 1] = dalpha1_dyc
    D[2, 2] = dalpha1_da
    D[2, 3] = dalpha1_db

    D[3, 0] = dalpha2_dxc
    D[3, 1] = dalpha2_dyc
    D[3, 2] = dalpha2_da
    D[3, 3] = dalpha2_db

    return D


def spiral2(
    x1, y1, x2, y2, alpha1, alpha2, npoints=100, plotspiral=False, verbose=True
):
    """
    return parameters to uniquely define the logspiral (2d ...x along a transect)
    """
    import numpy as np

    slopeangle = np.arctan2(y2 - y1, x2 - x1)

    deg2rad = np.pi / 180.0
    psi1 = (alpha1) * deg2rad
    psi2 = (alpha2) * deg2rad

    # import pdb;pdb.set_trace()
    if slopeangle < psi1:
        # print 'INPUTERROR: scarp angle (alpha1) must be less than average slope angle'
        raise Exception(
            "INPUT ERROR: scarp angle (alpha1) must be less than average slope angle: %s"
            % np.rad2deg(slopeangle)
        )

    if slopeangle > psi2:
        # print 'INPUTERROR: scarp angle (alpha1) must be less than average slope angle'
        raise Exception(
            "INPUT ERROR: toe angle (alpha2) must be greater than average slope angle: %s"
            % np.rad2deg(slopeangle)
        )

    slopeangledeg = slopeangle / deg2rad
    scarpdiff = (slopeangle - psi1) / deg2rad
    toediff = -(slopeangle - psi2) / deg2rad
    anglerat1 = scarpdiff / toediff
    anglerat2 = toediff / scarpdiff

    if scarpdiff / toediff < 0.3:
        errstr = (
            "****WARNING: scarp and toe angles differ significantly from a log-spiral \n"
            + "*************for a better fit decrease alpha1 or decrease alpha2 \n"
            + (
                "           0.3>(slopeangle-alpha1)/(alpha2-slopeangle)=%s \n"
                % anglerat1
            )
            + ("           slopeangle= %s \n" % slopeangledeg)
        )
        print(errstr)
        # raise Exception(errstr)

    if toediff / scarpdiff < 0.3:
        errstr = (
            "****WARNING: scarp and toe angles differ significantly from a log-spiral \n"
            + "*************for a better fit increase alpha1 or increase alpha2 \n"
            + (
                "           0.3>(alpha2-slopeangle)/(slopeangle-alpha1)=%s \n"
                % anglerat2
            )
            + ("           slopeangle= %s \n" % slopeangledeg)
        )
        print(errstr)
        # raise Exception(errstr)

    # initial guess
    xc0 = 0.5 * (x1 + x2)
    yc0 = y1 + 0.5 * (y1 - y2)
    r1 = np.sqrt((x1 - xc0) ** 2 + (y1 - yc0) ** 2)
    r2 = np.sqrt((x2 - xc0) ** 2 + (y2 - yc0) ** 2)
    a0 = 0.5 * (r1 + r2)
    b0 = 0.1
    # theta1 = 1.5*np.pi + psi1 + 2.*np.pi#np.pi
    # theta2 = 1.5*np.pi + psi2 + 2.*np.pi

    x0 = np.ones(4)
    x0[0] = xc0
    x0[1] = yc0
    x0[2] = a0
    x0[3] = b0

    (x, cov_x, infodict, mesg, ier) = scio.leastsq(
        objectf_4x4,
        x0,
        args=(x1, y1, x2, y2, alpha1, alpha2),
        maxfev=10000,
        col_deriv=False,
        Dfun=objectD_4x4,
        full_output=True,
    )

    xc = x[0]
    yc = x[1]
    a = x[2]
    b = x[3]

    theta1 = np.arctan2(y1 - yc, x1 - xc) + 2.0 * np.pi
    theta2 = np.arctan2(y2 - yc, x2 - xc) + 2.0 * np.pi

    if theta1 > theta2:
        errstr = "ERROR: fit failed*** try decreasing alpha2-alpha1"
        raise Exception(errstr)

    t1 = np.array(
        [b * np.cos(theta1) - np.sin(theta1), b * np.sin(theta1) + np.cos(theta1)]
    )
    t2 = np.array(
        [b * np.cos(theta2) - np.sin(theta2), b * np.sin(theta2) + np.cos(theta2)]
    )
    a1 = np.array([np.cos(psi1), np.sin(psi1)])
    a2 = np.array([np.cos(psi2), np.sin(psi2)])
    norm = np.linalg.norm

    # import pdb;pdb.set_trace()

    erroralpha1 = abs(
        np.rad2deg(np.arccos(min(1.0, np.dot(t1, a1) / (norm(a1) * norm(t1)))))
    )
    erroralpha2 = abs(
        np.rad2deg(np.arccos(min(1.0, np.dot(t2, a2) / (norm(a2) * norm(t2)))))
    )
    errorx1 = np.sqrt(
        (x1 - xc - a * np.exp(b * theta1) * np.cos(theta1)) ** 2
        + (y1 - yc - a * np.exp(b * theta1) * np.sin(theta1)) ** 2
    )
    errorx2 = np.sqrt(
        (x2 - xc - a * np.exp(b * theta2) * np.cos(theta2)) ** 2
        + (y2 - yc - a * np.exp(b * theta2) * np.sin(theta2)) ** 2
    )

    if verbose:
        print(mesg)
        print(
            ("logspiral parameters:\n xc=%s \n yc=%s \n a=%s \n b=%s " % (xc, yc, a, b))
        )

    if verbose:
        print("\n Goodness of fit:\n")
        print(("error in scarp angle: %s (degrees)\n" % erroralpha1))
        print(("error in toe angle: %s (degrees)\n" % erroralpha2))
        print(("error in (x1,y1): %s (meters)\n" % errorx1))
        print(("error in (x2,y2): %s (meters)\n" % errorx2))

    if max((erroralpha1, erroralpha2)) > 5.0:
        errstr = (
            "***WARNING: poor fit: Angle error> 5 degrees\n"
            + "***try adjusting angles so that following relative differences are closer in magnitude\n"
            + ("***headscarp: slopeangle-alpha1 = %s degrees\n" % scarpdiff)
            + ("***toe: alpha2-slopeangle = %s degrees\n" % toediff)
        )
        if max((erroralpha1, erroralpha2)) > 15.0:
            raise Exception(errstr)
        else:
            print(errstr)

    if max((errorx1, errorx2)) > 1.0:
        errstr = (
            "***WARNING: poor fit: error in enpoints > 1.0 meter\n"
            + "***try adjusting angles so that following relative differences are closer in magnitude\n"
            + ("***headscarp: slopeangle-alpha1 = %s degrees\n" % scarpdiff)
            + ("***toe:       alpha2-slopeangle = %s degrees\n" % toediff)
        )
        if max((errorx1, errorx2)) > 2.0:
            raise Exception(errstr)
        else:
            print(errstr)

    # import pdb;pdb.set_trace()

    if plotspiral == True:
        thetap = np.linspace(theta1 - 0.1 * np.pi, theta2 + 0.1 * np.pi, npoints)
        xvp = xc + a * np.exp(b * thetap) * np.cos(thetap)
        yvp = yc + a * np.exp(b * thetap) * np.sin(thetap)
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        plt.plot(
            xvp,
            yvp,
            [x1, x2],
            [y1, y2],
            "ro",
            [x1, x1 + 20 * np.cos(psi1)],
            [y1, y1 + 20 * np.sin(psi1)],
            "k-",
            [x2, x2 + 20.0 * np.cos(psi2)],
            [y2, y2 + 20.0 * np.sin(psi2)],
            "k-",
        )
        plt.plot(xc, yc, "ko")
        plt.axis("equal")
        plt.show()

    # import matplotlib.pyplot as plt
    # plt.plot(xv,yv,[x1,x2],[y1,y2],'ro',[x1,x1+20*np.cos(psi1)],[y1,y1+20*np.sin(psi1)],'k-',[x2,x2+20.*np.cos(psi2)],[y2,y2+20.*np.sin(psi2)],'k-')

    # plt.show()

    theta = np.linspace(theta1, theta2, npoints)
    # theta = np.linspace(0.0,4.0*np.pi,npoints)
    xv = xc + a * np.exp(b * theta) * np.cos(theta)
    yv = yc + a * np.exp(b * theta) * np.sin(theta)

    return (xv, yv)
