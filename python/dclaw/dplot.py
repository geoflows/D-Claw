"""
Useful things for plotting GeoClaw results.
"""

import numpy as np
from matplotlib.colors import Normalize
from numpy import ma as ma
from pyclaw.geotools import topotools
from pyclaw.plotters import colormaps

i_eta = 7

i_dig = 3
i_phi = i_dig
i_theta = i_dig + 1
i_fs = i_dig + 2
i_cohesion = i_dig + 3
i_taudir_x = i_dig + 4
i_taudir_y = i_dig + 5


# Colormaps from geoclaw
# Color attributes, single instance per run
# Colors
black = [0.0, 0.0, 0.0]
white = [1.0, 1.0, 1.0]
red = [1.0, 0.0, 0.0]
green = [0.0, 1.0, 0.0]
dark_green = [0.1, 0.4, 0.0]
light_green = [0.8, 1.0, 0.5]
blue = [0.0, 0.0, 1.0]
dark_blue = [0.2, 0.2, 0.7]
light_blue = [0.5, 0.5, 1.0]
blue_green = [0.0, 1.0, 1.0]
tan = [0.9, 0.8, 0.2]
tan = [0.8, 0.5, 0.2]
tan0 = [0.95, 0.95, 0.2]
tan1 = [0.9, 0.9, 0.2]
tan2 = [0.85, 0.7, 0.2]
tan3 = [0.8, 0.5, 0.2]
brown = [0.9, 0.8, 0.2]
gray5 = [0.5, 0.5, 0.5]
gray6 = [0.6, 0.6, 0.6]
gray7 = [0.7, 0.7, 0.7]
gray8 = [0.8, 0.8, 0.8]
gray9 = [0.9, 0.9, 0.9]
purple = [0.8, 0.3, 0.8]

green0 = [0.05, 0.2, 0.013]
green1 = [0.056, 0.24, 0.01]
green2 = [0.068, 0.3, 0.015]
green3 = [0.075, 0.35, 0.014]
green4 = [0.083, 0.4, 0.012]
green5 = [0.09, 0.45, 0.01]
green6 = [0.1, 0.5, 0.01]
green7 = [0.11, 0.56, 0.017]
green8 = [0.12, 0.63, 0.01]
green9 = [0.136, 0.7, 0.007]
green10 = [0.16, 0.78, 0.025]
green11 = [0.17, 0.86, 0.017]
green12 = [0.19, 0.93, 0.02]
green13 = [0.188, 0.98, 0.01]
green14 = [0.69, 0.98, 0.01]

yellow1 = [0.98, 0.96, 0.01]
tan1 = [0.98, 0.83, 0.01]
tan2 = [0.71, 0.62, 0.085]
tan3 = [0.65, 0.56, 0.07]
tan4 = [0.58, 0.50, 0.064]
tan5 = [0.48, 0.415, 0.043]
tan6 = [0.42, 0.36, 0.038]
tan7 = [0.36, 0.233, 0.029]

brown1 = [0.42, 0.167, 0.03]


# Colormaps
TSUNAMI_MAX_AMPLITUDE = 0.6
tsunami_colormap = colormaps.make_colormap(
    {-TSUNAMI_MAX_AMPLITUDE: blue, 0.0: blue_green, TSUNAMI_MAX_AMPLITUDE: red}
)


flume_colormap = colormaps.make_colormap(
    {0.0: tan0, 0.06: tan1, 0.12: tan2, 0.18: tan3}
)

oso_debris_colormap = colormaps.make_colormap(
    {0.0: white, 10.0: light_blue, 20.0: blue, 30.0: dark_blue}
)

oso_debris_colormap_invert = colormaps.make_colormap(
    {
        0.0: dark_blue,
        # 10.0: blue,
        # 20.0: light_blue,
        30.0: white,
    }
)

oso_land_colormap_gray = colormaps.make_colormap({0.0: gray8, 1.0: gray5})

oso_fs_colormap = colormaps.make_colormap({0.0: blue, 30.0: red})

oso_phi_colormap = colormaps.make_colormap({0.0: blue, 0.75: red})

oso_taudir_colormap = colormaps.make_colormap({0.0: blue, 1.0: red})

oso_cohesion_colormap = colormaps.make_colormap({0.0: blue, 1.0e6: red})


oso_land_colormap = colormaps.make_colormap(
    {
        0: green0,
        0.01: green1,
        0.02: green2,
        0.5 * 0.025: green3,
        0.5 * 0.05: green4,
        0.5 * 0.075: green5,
        0.5 * 0.100: green6,
        0.5 * 0.125: green7,
        0.5 * 0.150: green8,
        0.5 * 0.175: green9,
        0.5 * 0.200: green10,
        0.5 * 0.225: green11,
        0.5 * 0.250: green12,
        0.5 * 0.275: green13,
        0.5 * 0.300: green14,
        0.325: yellow1,
        0.70: tan1,
        0.75: tan2,
        0.80: tan3,
        0.85: tan4,
        0.90: tan5,
        0.95: tan6,
        1.0: tan7,
    }
)

oso_land_colormap2 = colormaps.make_colormap(
    {
        0: green0,
        0.01: green1,
        0.1: green14,
        0.2: yellow1,
        0.5: tan5,
        0.6: tan7,
        1.0: brown1,
    }
)

land1_colormap = colormaps.make_colormap(
    {0.0: dark_green, 1000.0: green, 2000.0: light_green, 4000.0: tan}
)

land2_colormap = colormaps.make_colormap(
    {0: dark_green, 50: green, 100: light_green, 200: tan}
)

runoutpad_colormap = colormaps.make_colormap({0.0: gray8, 1.0: gray5})

water_land_colormap = colormaps.make_colormap(
    {
        -1000: dark_blue,
        -500: blue,
        0: light_blue,
        0.1: tan,
        5: tan,
        6: dark_green,
        1000: green,
        2000: light_green,
        4000: tan,
    }
)

bathy1_colormap = colormaps.make_colormap(
    {-1000: brown, 0: tan, 0.1: dark_green, 1000: green, 2000: light_green}
)

bathy2_colormap = colormaps.make_colormap(
    {
        -1000: brown,
        -100: tan,
        0: dark_green,
        0.1: dark_green,
        1000: green,
        2000: light_green,
    }
)

bathy3_colormap = colormaps.make_colormap(
    {
        -1: [0.3, 0.2, 0.1],
        -0.01: [0.95, 0.9, 0.7],
        0.01: [0.5, 0.7, 0],
        1: [0.2, 0.5, 0.2],
    }
)

seafloor_colormap = colormaps.make_colormap({-1: [0.3, 0.2, 0.1], 0: [0.95, 0.9, 0.7]})

land_colormap = colormaps.make_colormap({0: [0.95, 0.9, 0.7], 1: [0.2, 0.5, 0.2]})


colormaps_list = {
    "tsunami": tsunami_colormap,
    "land1": land1_colormap,
    "land2": land2_colormap,
    "water_land": water_land_colormap,
    "bathy1": bathy1_colormap,
    "bathy2": bathy2_colormap,
    "runoutpad": runoutpad_colormap,
    "oso_land": oso_land_colormap,
}


def plot_colormaps():
    r"""Plots all colormaps avaiable or the ones specified"""

    import matplotlib.pyplot as plt
    import numpy as np

    a = np.linspace(0, 1, 256).reshape(1, -1)
    a = np.vstack((a, a))

    nmaps = len(colormaps_list) + 1

    fig = plt.figure(figsize=(5, 10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)

    for i, name in enumerate(colormaps_list):
        ax = plt.subplot(nmaps, 1, i + 1)
        plt.axis("off")
        plt.imshow(a, aspect="auto", cmap=colormaps_list[name], origin="lower")
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], name, fontsize=10, horizontalalignment="right")

    # plt.show()


land_colors = colormaps.make_colormap({0: [0.5, 0.7, 0], 1: [0.2, 0.5, 0.2]})
# water_colors = colormaps.make_colormap({-1.:'r', 0.:[0, .8, .8], 1.: 'b'})
# land_colors = land2_colormap
water_colors = tsunami_colormap

# Plotting functions

# The drytol parameter is used in masking land and water and
# affects what color map is used for cells with small water depth h.
# The best value to use often depends on the application and can
# be set for an application by setting current_data.user.drytol in
# a beforeframe function, for example.  If it's not set by the user,
# the following default value is used (in meters):

drytol_default = 1.0e-3
rho_f_default = 1000.0
rho_s_default = 2700.0
grav_default = 9.81
mu_default = 0.001
kappita_default = 0.0001
kappita_diff_default = 1.0
m0_default = 0.52
delta_default = 0.01
m_crit_default = 0.62
sigma_0_default = 1.0e3
alpha_c_default = 1.0
alpha_seg_default = 0.0
phi_bed_default = 40
bed_normal_default = 1
c1_default = 1.0


def gmod(current_data):
    # gravity
    grav = getattr(current_data.user, "gravity", grav_default)
    bed_normal = getattr(current_data.user, "bed_normal", bed_normal_default)

    gmod = grav

    if bed_normal==1:
        aux = current_data.aux
        theta = aux[:, :, i_theta]
        gmod = grav*np.cos(theta)

    return gmod


def solid_frac(current_data):
    # solid volume fraction
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    hm = q[:, :, 3]
    with np.errstate(divide="ignore", invalid="ignore"):
        m = ma.masked_where(h < drytol, hm / h)
    return m


def m_minus_mcrit(current_data):
    drytol = getattr(current_data.user, "drytol", drytol_default)
    m_crit = getattr(current_data.user, "m_crit", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    hm = q[:, :, 3]
    with np.errstate(divide="ignore", invalid="ignore"):
        m = ma.masked_where(h < drytol, (hm / h) - m_crit)
    return m


def solid_frac_gt03(current_data):
    # solid fraction where it is greater than 0.3
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    hm = q[:, :, 3]
    with np.errstate(divide="ignore", invalid="ignore"):
        m = ma.masked_where(h < drytol, hm / h)
        m = ma.masked_where(m < 0.3, m)
    return m


def density(current_data):
    # new segregation might modify.
    m = solid_frac(current_data)
    rho_f = getattr(current_data.user, "rho_f", rho_f_default)
    rho_s = getattr(current_data.user, "rho_s", rho_s_default)
    rho = (1.0 - m) * rho_f + m * rho_s
    return rho


def basalP(current_data):
    # basal pressure.
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    basalP = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 4])
    return basalP


def lithostaticP(current_data):
    # lithostatic pressure
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    rho = density(current_data)
    var = ma.masked_where(h < drytol, gmod(current_data) * rho * h)
    return var


def hydrostaticP(current_data):
    # hydrostatic pressure
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    rho_f = getattr(current_data.user, "rho_f", rho_f_default)
    var = ma.masked_where(h < drytol, gmod(current_data) * rho_f * h)
    return var


def sigma_e(current_data):
    # effective basal pressure.
    drytol = getattr(current_data.user, "drytol", drytol_default)
    rho = density(current_data)
    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    p = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 4])
    sigma_e = rho * gmod(current_data) * h - p
    sigma_e[sigma_e < 0.0] = 0.0
    return sigma_e


def sigma_e_over_hydrostatic(current_data):
    # effective basal pressure.
    return sigma_e(current_data)/hydrostaticP(current_data)


def N(current_data):  # dimensionless state parameter
    mu = getattr(current_data.user, "mu", mu_default)
    rho_s = getattr(current_data.user, "rho_s", rho_s_default)
    delta = getattr(current_data.user, "delta", delta_default)
    gamma = shear(current_data)
    sigbedc = rho_s * (gamma * delta) ** 2.0 + sigma_e(current_data)
    N = mu * gamma / (sigbedc)
    N[sigbedc < 0.0] = 0.0
    return N


def m_eqn(current_data):
    # equilibrium value for m.
    drytol = getattr(current_data.user, "drytol", drytol_default)
    rho_f = getattr(current_data.user, "rho_f", rho_f_default)
    sigma_0 = getattr(current_data.user, "sigma_0", sigma_0_default)
    alpha_c = getattr(
        current_data.user, "alpha_c", alpha_c_default
    )  # this is a in part 2, equation 2.8
    m_crit = getattr(current_data.user, "m_crit", m_crit_default)
    alpha_seg = getattr(current_data.user, "alpha_seg", alpha_seg_default)

    alpha_seg = 1.0 - alpha_seg  # digclaw.mod.f90 line 121

    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    hm = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 3])
    with np.errstate(divide="ignore", invalid="ignore"):
        m = hm / h

    pm = species1_fraction(current_data)

    # if segregation occurs, then need to reduce
    # mcrit by m_crit_pm
    # TODO. this part of code may change as segregation use is changed.

    if alpha_seg - 1.0 < 1.0e-6:
        seg = 0.0
        rho_fp = rho_f
        pmtanh01 = 0.0
    else:
        seg = 1.0
        pmtanh01 = seg * (0.5 * (np.tanh(40.0 * (pm - 0.90)) + 1.0))
        rho_fp = (1.0 - pmtanh01) * rho_f

    m_crit_pm = pmtanh01 * 0.09
    m_crit_m = m_crit - m_crit_pm

    m_eqn = m_crit_m / (1.0 + np.sqrt(N(current_data)))

    return m_eqn


def meqn_over_mcrit(current_data):
    mcrit = getattr(current_data.user, "mcrit", drytol_default)
    return m_eqn(current_data) / mcrit


def m_minus_meqn(current_data):
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    hm = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 3])
    with np.errstate(divide="ignore", invalid="ignore"):
        m = hm / h
    return m - m_eqn(current_data)


def shear(current_data):
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    hu = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 1])
    hv = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 2])
    with np.errstate(divide="ignore", invalid="ignore"):
        u = hu / h
        v = hv / h
    vnorm = np.sqrt(u ** 2.0 + v ** 2.0)
    return 2.0 * vnorm / h  # in code refered to as hbounded (defined as h)


def kperm(current_data):
    # permeability
    drytol = getattr(current_data.user, "drytol", drytol_default)
    kappita = getattr(current_data.user, "kappita", kappita_default)
    m0 = getattr(current_data.user, "m0", m0_default)
    kappita_diff = getattr(current_data.user, "kappita_diff", kappita_diff_default)
    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    hm = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 3])
    with np.errstate(divide="ignore", invalid="ignore"):
        m = hm / h
    pm = species1_fraction(current_data)
    kappita2 = kappita * kappita_diff
    kequiv = kappita2 * pm + kappita * (1 - pm)
    return kequiv * np.exp(-(m - m0) / (0.04))


def dilatency(current_data):
    # depth averaged dilatency rate
    mu = getattr(current_data.user, "mu", mu_default)
    rho_f = getattr(current_data.user, "rho_f", rho_f_default)
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    p = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 4])
    # Royal Society, Part 2, Eq 2.6
    D = 2.0 * (kperm(current_data) / (mu * h)) * (rho_f * gmod(current_data) * h - p)
    return D


def tanpsi(current_data):
    # dilation angle.
    c1 = getattr(current_data.user, "c1", c1_default)
    gamma = shear(current_data)

    # in code, m-meqn is regularized based on shear.
    # c1*(m-m_eqn)*tanh(shear/0.1)
    return c1 * m_minus_meqn(current_data) * np.tanh(shear/0.1)


def psi(current_data):
    return np.arctan(m_minus_meqn(current_data)) # maybe this should be arctan of tanpsi (with the regularization)


def Fgravitational(current_data):
    # gravitational driving force.

    drytol = getattr(current_data.user, "drytol", drytol_default)
    bed_normal = getattr(current_data.user, "bed_normal", bed_normal)
    gmod = gmod(current_data)

    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    eta = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 7])

    if bed_normal == 1:
        theta = q[:, :, i_theta]
        sintheta = np.sin(theta)
    else:
        sintheta = 0.

    dx = current_data.dx
    dy = current_data.dy

    hL = np.roll(h.copy(), 1, axis = 1) # roll right on columns so that value at (i, j-1) is at (i,j)
    hL[:, 0] = np.nan # first column has undefined values
    hR = np.roll(h.copy(), -1, axis = 1)
    hR[:, -1] = np.nan
    hT = np.roll(h.copy, -1, axis = 0)
    hT[-1, :] = np.nan
    hB = np.roll(h.copy, 1, axis = 0)
    hB[0, :] = np.nan

    etaL = np.roll(eta.copy(), 1, axis = 1)
    etaL[:, 0] = np.nan
    etaR = np.roll(eta.copy(), -1, axis = 1)
    etaR[:, -1] = np.nan
    etaT = np.roll(eta.copy, -1, axis = 0)
    etaT[-1, :] = np.nan
    etaB = np.roll(eta.copy, 1, axis = 0)
    etaB[0, :] = np.nan

    FxL = -gmod*0.5*(h+hL)*(eta-etaL)/(dx) + gmod*0.5*(h+hL)*np.sin(theta)
    FyL = -gmod*0.5*(h+hB)*(eta-etaB)/(dy)

    FxR = -gmod*0.5*(h+hR)*(etaR-eta)/(dx) + gmod*0.5*(h+hR)*np.sin(theta)
    FyR = -gmod*0.5*(h+hT)*(etaT-eta)/(dy)

    Fx = np.min(np.abs(FxR), np.abs(FxL))
    Fy = np.min(np.abs(FyR), np.abs(FyL))

    F = np.sqrt(Fx**2, Fy**2)
    # Royal Proceedings, Part 2, equation 2.4b,c (momentum source terms) first term on RHS
    # deta/dx portion taken from calc_taudir
    return np.sqrt(Fx**2, Fy**2)


def Fdrag(current_data):  # units of force per unit area
# Royal Proceedings, Part 2, equation 2.4b,c (momentum source terms) second term on RHS
# what is this term?
# katy asks: driving force due to longitudinal stress gradients? (based on text right before part 2 eq 2.15
# dave says: I don't know if there's a simple interpretation...it sort of drops out from the derivation and then rearrangement of the equations into conservative form (ie. derivatives on hu not u). I think it might be similar to a drag term that appears on fully two-phase equations for solid and fluid velocity fields. ...I'll revisit the derivation and see if I can shed more light on it.
    drytol = getattr(current_data.user, "drytol", drytol_default)
    rho_f = getattr(current_data.user, "rho_f", rho_f_default)

    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    hu = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 1])
    hv = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 2])
    with np.errstate(divide="ignore", invalid="ignore"):
        u = hu / h
        v = hv / h
    vnorm = np.sqrt(u ** 2.0 + v ** 2.0)
    D = dilatency(current_data)
    rho = density(current_data)

    return vnorm * D * (rho - rho_f) / rho


def Ffluid(current_data):  # units of force per unit area
    # Resisting force due to fluid.
    drytol = getattr(current_data.user, "drytol", drytol_default)
    mu = getattr(current_data.user, "mu", drytol_default)
    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    hu = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 1])
    hv = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 2])
    hm = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 3])
    with np.errstate(divide="ignore", invalid="ignore"):
        u = hu / h
        v = hv / h
        m = hm / h
    vnorm = np.sqrt(u ** 2.0 + v ** 2.0)
    tauf = 2.0 * mu * (1.0 - m) * vnorm / h
    return tauf


def Fsolid(current_data):  # units of force per unit area
    # resisting force due to solid.
    drytol = getattr(current_data.user, "drytol", drytol_default)
    phi_bed = getattr(current_data.user, "phi_bed", phi_bed_default)

    q = current_data.q
    h = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 0])
    p = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 4])

    rho_f = getattr(current_data.user, "rho_f", rho_f_default)
    hydrostatic = ma.masked_where(h < drytol, gmod(current_data) * rho_f * h)
    sigma_e = hydrostatic - p
    tau_s = sigma_e * np.tan(phi_bed + psi(current_data))

    # where tau_s>F, then material is static. Make taus = F driving.
    F = Fdriving(current_data)
    tau_s[tau_s > F] = F[tau_s > F]
    return tau_s


def Fdriving(current_data):  # units of force per unit area
    return Fgravitational(current_data) + Fdrag(current_data)


def Fresisting(current_data):  # units of force per unit area
    return Ffluid(current_data) + Fsolid(current_data)


def Fnet(current_data):  # units of force per unit area
    return Fresisting(current_data) - Fdriving(current_data)


def liquefaction_ratio(current_data):
    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    p = ma.masked_where(q[:, :, 0] < drytol, q[:, :, 4])
    litho = lithostaticP(current_data)
    ratio = ma.masked_where(q[:, :, 0] < drytol, p / litho)
    return ratio


def eta(current_data):
    """
    Return eta
    """
    q = current_data.q
    eta = q[:, :, i_eta]
    return eta


def topo(current_data):
    """
    Return topography = eta - h.
    """
    q = current_data.q
    h = q[:, :, 0]
    eta = q[:, :, i_eta]
    topo = eta - h
    return topo


def land(current_data):
    """
    Return a masked array containing the surface elevation only in dry cells.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    # drytol = 5.e-2
    q = current_data.q
    h = q[:, :, 0]
    eta = q[:, :, i_eta]
    land = ma.masked_where(h > drytol, eta)
    return land


def water(current_data):
    """Deprecated: use surface instead."""
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    eta = q[:, :, i_eta]
    water = ma.masked_where(h <= drytol, eta)
    return water


def depth(current_data):
    """
    Return a masked array containing the depth of fluid only in wet cells.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    # drytol = 5.e-2
    q = current_data.q
    h = q[:, :, 0]
    depth = ma.masked_where(h <= drytol, h)
    return depth


def surface(current_data):
    """
    Return a masked array containing the surface elevation only in wet cells.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    eta = q[:, :, i_eta]
    water = ma.masked_where(h <= drytol, eta)
    return water


def surface_solid_frac_lt03(current_data):
    """
    Return a masked array containing the surface elevation only in wet cells.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    eta = q[:, :, i_eta]
    hm = q[:, :, 3]

    with np.errstate(divide="ignore", invalid="ignore"):
        m = hm / h

    water = ma.masked_where(h <= drytol, eta)
    water = ma.masked_where(m > 0.3, water)

    return water


def surface_or_depth(current_data):
    """
    Return a masked array containing the surface elevation where the topo is
    below sea level or the water depth where the topo is above sea level.
    Mask out dry cells.  Assumes sea level is at topo=0.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    from numpy import ma, where

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    eta = q[:, :, i_eta]
    topo = eta - h
    surface = ma.masked_where(h <= drytol, eta)
    depth = ma.masked_where(h <= drytol, h)
    surface_or_depth = where(topo < 0, surface, depth)
    return surface_or_depth


def velocity_u(current_data):
    """
    Return a masked array containing velocity u in wet cells.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    drytol = 1.0
    q = current_data.q
    h = q[:, :, 0]
    hu = q[:, :, 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = ma.masked_where(h <= drytol, hu / h)
    return u


def velocity_v(current_data):
    """
    Return a masked array containing velocity v in wet cells.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    drytol = 1.0
    q = current_data.q
    h = q[:, :, 0]
    hv = q[:, :, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        v = ma.masked_where(h <= drytol, hv / h)
    return v


def species1_fraction(current_data):
    """
    Return a masked array containing the fraction of species 1 in wet cells.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    # drytol = 1.0
    q = current_data.q
    h = q[:, :, 0]
    hchi = q[:, :, 5]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi1 = ma.masked_where(h <= drytol, hchi / h)
    return chi1


def species2_fraction(current_data):
    """
    Return a masked array containing the fraction of species 2 in wet cells.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    # drytol = 1.0
    q = current_data.q
    h = q[:, :, 0]
    hchi = q[:, :, 5]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = ma.masked_where(h <= drytol, 1.0 - (hchi / h))
    return chi2


def velocity(current_data):
    """
    Return a masked array containing a tuple of x and y directed velocity (u,v)
    in wet cells.

    velocity defined as sqrt(u**2 + v**2)
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    hu = q[:, :, 1]
    hv = q[:, :, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = ma.masked_where(h <= drytol, hu / h)
        v = ma.masked_where(h <= drytol, hv / h)
    return (u, v)


def velocity_magnitude(current_data):
    """
    Return a masked array of the magnitude of velocity at wet cells.

    velocity defined as sqrt(u**2 + v**2)
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    hu = q[:, :, 1]
    hv = q[:, :, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = hu / h
        v = hv / h
    vel = np.sqrt(u ** 2 + v ** 2)

    vel = ma.masked_where(h <= drytol, vel)

    return vel


def fs(current_data):
    """
    Return a masked array containing factor of safety.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    aux = current_data.aux
    fs = aux[:, :, i_fs]
    fs = ma.masked_where(h <= drytol, fs)
    return fs


def phi(current_data):
    """
    Return a masked array containing factor of safety.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    aux = current_data.aux
    fs = aux[:, :, i_phi]
    fs = ma.masked_where(h <= drytol, fs)

    return fs


def cohesion(current_data):
    """
    Return a masked array containing factor of safety.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    aux = current_data.aux
    c = aux[:, :, i_cohesion]
    c = ma.masked_where(h <= drytol, c)

    return c


def taudir(current_data):
    """
    Return a masked array containing factor of safety.
    """
    from numpy import ma

    drytol = getattr(current_data.user, "drytol", drytol_default)
    q = current_data.q
    h = q[:, :, 0]
    aux = current_data.aux
    c = aux[:, :, i_taudir_x]
    c = ma.masked_where(h <= drytol, c)

    return c


class TopoPlotData(object):
    def __init__(self, fname):
        self.fname = fname
        self.topotype = 3
        self.neg_cmap = None
        self.pos_cmap = None
        self.cmax = 100.0
        self.cmin = -4000.0
        self.climits = None
        self.figno = 200
        self.addcolorbar = False
        self.addcontour = False
        self.contour_levels = [0, 0]
        self.xlimits = None
        self.ylimits = None
        self.coarsen = 1
        self.imshow = True
        self.gridedges_show = True
        self.print_fname = True

    def plot(self):
        plot_topo_file(self)


def plot_topo_file(topoplotdata):
    """
    Read in a topo or bathy file and produce a pcolor map.
    """

    import os

    import pylab
    from pyclaw.data import Data

    fname = topoplotdata.fname
    topotype = topoplotdata.topotype
    if topoplotdata.climits:
        # deprecated option
        cmin = topoplotdata.climits[0]
        cmax = topoplotdata.climits[1]
    else:
        cmin = topoplotdata.cmin
        cmax = topoplotdata.cmax
    figno = topoplotdata.figno
    addcolorbar = topoplotdata.addcolorbar
    addcontour = topoplotdata.addcontour
    contour_levels = topoplotdata.contour_levels
    xlimits = topoplotdata.xlimits
    ylimits = topoplotdata.ylimits
    coarsen = topoplotdata.coarsen
    imshow = topoplotdata.imshow
    gridedges_show = topoplotdata.gridedges_show
    neg_cmap = topoplotdata.neg_cmap
    pos_cmap = topoplotdata.pos_cmap
    print_fname = topoplotdata.print_fname

    if neg_cmap is None:
        neg_cmap = colormaps.make_colormap({cmin: [0.3, 0.2, 0.1], 0: [0.95, 0.9, 0.7]})
    if pos_cmap is None:
        pos_cmap = colormaps.make_colormap({0: [0.5, 0.7, 0], cmax: [0.2, 0.5, 0.2]})

    if abs(topotype) == 1:

        X, Y, topo = topotools.topofile2griddata(fname, topotype)
        topo = pylab.flipud(topo)
        Y = pylab.flipud(Y)
        x = X[0, :]
        y = Y[:, 0]
        xllcorner = x[0]
        yllcorner = y[0]
        cellsize = x[1] - x[0]

    elif abs(topotype) == 3:

        file = open(fname, "r")
        lines = file.readlines()
        ncols = int(lines[0].split()[0])
        nrows = int(lines[1].split()[0])
        xllcorner = float(lines[2].split()[0])
        yllcorner = float(lines[3].split()[0])
        cellsize = float(lines[4].split()[0])
        NODATA_value = int(lines[5].split()[0])

        print(("Loading file ", fname))
        print(("   nrows = %i, ncols = %i" % (nrows, ncols)))
        topo = pylab.loadtxt(fname, skiprows=6, dtype=float)
        print("   Done loading")

        if 0:
            topo = []
            for i in range(nrows):
                topo.append(
                    pylab.array(
                        lines[6 + i],
                    )
                )
            print(("+++ topo = ", topo))
            topo = pylab.array(topo)

        topo = pylab.flipud(topo)

        x = pylab.linspace(xllcorner, xllcorner + ncols * cellsize, ncols)
        y = pylab.linspace(yllcorner, yllcorner + nrows * cellsize, nrows)
        print(("Shape of x, y, topo: ", x.shape, y.shape, topo.shape))

    else:
        raise Exception("*** Only topotypes 1 and 3 supported so far")

    if coarsen > 1:
        topo = topo[slice(0, nrows, coarsen), slice(0, ncols, coarsen)]
        x = x[slice(0, ncols, coarsen)]
        y = y[slice(0, nrows, coarsen)]
        print(("Shapes after coarsening: ", x.shape, y.shape, topo.shape))

    if topotype < 0:
        topo = -topo

    if figno:
        pylab.figure(figno)

    if topoplotdata.imshow:
        color_norm = Normalize(cmin, cmax, clip=True)
        xylimits = (x[0], x[-1], y[0], y[-1])
        # pylab.imshow(pylab.flipud(topo.T), extent=xylimits, \
        pylab.imshow(
            pylab.flipud(topo),
            extent=xylimits,
            cmap=cmap,
            interpolation="nearest",
            norm=color_norm,
        )
    else:
        neg_topo = ma.masked_where(topo > 0.0, topo)
        all_masked = ma.count(neg_topo) == 0
        if not all_masked:
            pylab.pcolormesh(x, y, neg_topo, cmap=neg_cmap)
            pylab.clim([cmin, 0])
            if addcolorbar:
                pylab.colorbar()

        pos_topo = ma.masked_where(topo < 0.0, topo)
        all_masked = ma.count(pos_topo) == 0
        if not all_masked:
            pylab.pcolormesh(x, y, pos_topo, cmap=pos_cmap)
            pylab.clim([0, cmax])
            if addcolorbar:
                pylab.colorbar()

    pylab.axis("scaled")

    if addcontour:
        pylab.contour(x, y, topo, levels=contour_levels, colors="k")

    if gridedges_show:
        pylab.plot([x[0], x[-1]], [y[0], y[0]], "k")
        pylab.plot([x[0], x[-1]], [y[-1], y[-1]], "k")
        pylab.plot([x[0], x[0]], [y[0], y[-1]], "k")
        pylab.plot([x[-1], x[-1]], [y[0], y[-1]], "k")

    if print_fname:
        fname2 = os.path.splitext(fname)[0]
        pylab.text(xllcorner + cellsize, yllcorner + cellsize, fname2, color="m")

    topodata = Data()
    topodata.x = x
    topodata.y = y
    topodata.topo = topo

    return topodata
