"""
explore steady-state dynamics
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scio

from mpl_toolkits.mplot3d import Axes3D

g=9.81
theta=np.deg2rad(30.0)
gx = g*np.sin(theta)
gz = g*np.cos(theta)
rho_s = 2700.0
rho_f = 1000.0
rho_a = 2000.0
#rho = 2000.0
delta = 0.05
mu = 5.0*1.e-3
alpha = rho_a*delta*np.sqrt(gx)
k0 = 80.0*1.e-9
a0 = 0.005
sig0 = 1000.0
compress = .01/1000.
m_crit = 0.64
phi = np.deg2rad(43.0)

def plot_meqn():

    m_critv = np.linspace(0.0,1.0)

    a = 1.0 +0.0J
    b = alpha/mu + 0.0J
    c = -1.0 + 0.0J
    d = (m_crit-1.0)*alpha/mu + 0.0J

    q = (3.*c - b**2)/9.0
    r = (c*b - 3.0*d)/6.0 - b**3/27.0
    s1 = (r + np.sqrt(q**3 + r**2))**(1./3.)
    s2 = (r - np.sqrt(q**3 + r**2))**(1./3.)

    #q = (3.*c - b**2)/9.0
    #r = -27*d + b*(9.0*c-2.0*b**2)
    #discriminant = q**3 + r**2
    #s = r + np.sqrt(discriminant) 
    #t = r - np.sqrt(discriminant)
    #term1 = np.sqrt(3.0)*((-t + s)/2.0)
    #r13 = 2.0*np.sqrt(q)
    #x1 = (-term1 + r13*np.cos(q**3/3.0))
    x1 = (s1+s2) - b/3.0
    meqn = (1.0 - x1**2)

    #import pdb;pdb.set_trace()
    plt.plot(m_crit,m_crit,'r')
    plt.plot(m_crit,meqn,'b')
    plt.show()

def quiver_plot():

    def meqn04mcrit(m_crit):
        a = 1.0 +0.0J
        b = alpha/mu + 0.0J
        c = -1.0 + 0.0J
        d = (m_crit-1.0)*alpha/mu + 0.0J

        q = (3.*c - b**2)/9.0
        r = (c*b - 3.0*d)/6.0 - b**3/27.0
        s1 = (r + np.sqrt(q**3 + r**2))**(1./3.)
        s2 = (r - np.sqrt(q**3 + r**2))**(1./3.)

        x1 = (s1+s2) - b/3.0
        meqn = np.real(1.0 - x1**2)

        return meqn

    mcrit = 0.64
    meqn0 = meqn04mcrit(mcrit)

    A = np.sqrt(g)*mu/(rho*g*k)
    C = k/(compress*mu*np.sqrt(g))
    B = rho*np.sqrt(g)/mu

    m = np.linspace(0,1)
    p = np.linspace(-1.,1.)

    M,P = np.meshgrid(m,p)

    u0 = B*gx/(2.0*g*(1-meqn0))
    m0 = meqn0
    p0 = 0.0
    import pdb; pdb.set_trace()

    fM = (2./A)*M*P
    fP = -3.*C*P - 3.*A*C*(M-meqn)*u

    plt.quiver(M,P,fM,fP)
    plt.show()

def integrate(m0,h0,u0,pe0,npoints=100,tend=10.0):


    def rhs(mp,hp,up,pep):

        pp = pep + rho_f*gz*hp
        rhop = mp*rho_s + (1.0-mp)*rho_f
        shear = 2.*up/hp
        sigbed = max(0.0,rhop*gz*hp - pp)
        sigbedc = rho_s*(shear*delta)**2.0 + sigbed
        N = shear*mu/(sigbedc)
        meqn = m_crit/(1.0+np.sqrt(N))
        compress = a0/(mp*(sigbed+sig0))
        k = k0*np.exp((0.6-mp)/0.04)

        f = np.ones(4)

        rhorat = (rhop-rho_f)/rhop
        
        f[0] = 2.*k*mp*pep/(mu*hp**2)
        f[1] = rhorat*2.0*k*pep/(mu*hp)
        f[2] =  gx - max(0.,gz*rhorat - pep/(rhop*hp))*max(0.0,np.tan(phi+np.arctan(mp-meqn))) - (1.-mp)*2.0*mu*up/(rhop*hp**2)
        f[3] = -(3.*k/(compress*mu*hp**2))*(1.0 - 0.5*compress*rho_f*gz*hp*rhorat)*pep - 3*up*np.tan(mp-meqn)/(compress*hp)
        #import pdb;pdb.set_trace()
        tt=5.
        return (f,meqn)


    t = np.linspace(0.0,tend,npoints)
    dt = t[1]-t[0]

    m = np.zeros(npoints)
    h = np.zeros(npoints)
    u = np.zeros(npoints)
    pe = np.zeros(npoints)
    m_eqn = np.zeros(npoints)

    m[0] = m0
    h[0] = h0
    u[0] = u0
    pe[0]= pe0
    m_eqn[0]=m_crit

    for n in range(1,npoints):
        (f,meqn) = rhs(m[n-1],h[n-1],u[n-1],pe[n-1])
        m[n] = max(0.0,m[n-1] + dt*f[0])
        m[n] = min(1.0,m[n])
        h[n] = h[n-1] + dt*f[1]
        u[n] = max(0.0,u[n-1] + dt*f[2])
        pe[n] = pe[n-1] + dt*f[3]
        m_eqn[n] = meqn

    #import pdb;pdb.set_trace()
    p = pe + rho_f*gz*h
    rho = m*rho_s + (1.0-m)*rho_f
    pa = pe/(gz*h*rho-rho_f*(gz*h))

    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(m[0],u[0],pa[0],c='r',marker='o')

    ax.scatter(m,u,pe,c='b',marker='o')
    ax.set_xlabel('m')
    ax.set_ylabel('u')
    ax.set_zlabel('pa')

    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(u,pa,'bo-')
    ax.set_xlabel('u')
    ax.set_ylabel('pa')

    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(m,pa,'bo-')
    ax.set_xlabel('m')
    ax.set_ylabel('pa')

    fig = plt.figure(4)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(t,pa,'b-')
    ax.set_xlabel('t')
    ax.set_ylabel('pa')

    fig = plt.figure(5)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(t,u,'b-')
    ax.set_xlabel('t')
    ax.set_ylabel('u')

    fig = plt.figure(6)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(t,m,'b-',t,m_eqn,'r-')
    ax.set_xlabel('t')
    ax.set_ylabel('m, m_eqn')

    fig = plt.figure(7)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(t,h,'b-')
    ax.set_xlabel('t')
    ax.set_ylabel('h')

    plt.show()


def integrate2(m0,h0,u0,pe0,npoints=100,tend=10.0):


    def fobject(q,mn,hn,un,pen,dt,mu,compress,k0,a0,gz,rho_f,rho_s,delta,m_crit):

        mp=q[0] 
        hp=q[1] 
        up=q[2] 
        pep=q[3] 
        pp = pep + rho_f*gz*hp
        rhop = mp*rho_s + (1.0-mp)*rho_f
        shear = 2.*up/hp
        sigbed = max(0.0,rhop*gz*hp - pp)
        sigbedc = rho_s*(shear*delta)**2.0 + sigbed
        N = shear*mu/(sigbedc)
        meqn = m_crit/(1.0+np.sqrt(N))
        compress = a0/(mp*(sigbed+sig0))
        k = k0*np.exp((0.6-mp)/0.04)
        rhorat = (rhop-rho_f)/rhop

        psi0 = 2.*k*mp*pep/(mu*hp**2)
        psi1 = rhorat*2.0*k*pep/(mu*hp)
        psi2 = gx - max(0.,gz*rhorat - pep/(rhop*hp))*max(0.0,np.tan(phi+np.arctan(mp-meqn))) - (1.-mp)*2.0*mu*up/(rhop*hp**2)
        psi3 = -(3.*k/(compress*mu*hp**2))*(1.0 - 0.5*compress*rho_f*gz*hp*rhorat)*pep - 3.*up*np.tan(mp-meqn)/(compress*hp)
    
        f = np.ones(4)
        f[0] = mn + dt*psi0 - mp
        f[1] = hn + dt*psi1 - hp
        f[2] = un + dt*psi2 - up
        f[3] = pen + dt*psi3 - pep

        return f


    t = np.linspace(0.0,tend,npoints)
    dt = t[1]-t[0]

    m = np.zeros(npoints)
    h = np.zeros(npoints)
    u = np.zeros(npoints)
    pe = np.zeros(npoints)
    m_eqn = np.zeros(npoints)

    m[0] = m0
    h[0] = h0
    u[0] = u0
    pe[0]= pe0
    m_eqn[0]=m_crit

    for n in range(0,npoints-1):
        q0 = np.array([m[n],h[n],u[n],pe[n]])
        argstup = (m[n],h[n],u[n],pe[n],dt,mu,compress,k0,a0,gz,rho_f,rho_s,delta,m_crit)
        (q,cov_x,infodict,mesg,ier) = scio.leastsq(fobject,q0,args=argstup,maxfev=10000,col_deriv=False,Dfun=None,full_output=True)
        #import pdb;pdb.set_trace()
        m[n+1] = q[0]
        h[n+1] = q[1]
        u[n+1] = max(0.0,q[2])
        pe[n+1] = q[3]

    p = pe + rho_f*gz*h
    rho = m*rho_s + (1.0-m)*rho_f
    pa = pe/(gz*h*rho-rho_f*(gz*h))

    #import pdb;pdb.set_trace()
    p = pe + rho_f*gz*h
    rho = m*rho_s + (1.0-m)*rho_f
    pa = pe/(gz*h*rho-rho_f*(gz*h))

    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(m[0],u[0],pa[0],c='r',marker='o')

    ax.scatter(m,u,pe,c='r',marker='o')
    ax.set_xlabel('m')
    ax.set_ylabel('u')
    ax.set_zlabel('pa')

    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(u,pa,'r-')
    ax.set_xlabel('u')
    ax.set_ylabel('pa')

    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(m,pa,'r-')
    ax.set_xlabel('m')
    ax.set_ylabel('pa')

    fig = plt.figure(4)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(t,pa,'r-')
    ax.set_xlabel('t')
    ax.set_ylabel('pa')

    fig = plt.figure(5)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(t,u,'r-')
    ax.set_xlabel('t')
    ax.set_ylabel('u')

    fig = plt.figure(6)
    ax = fig.add_subplot(111)
    #plt.plot(u[0],pe[0],'ro')
    plt.plot(t,m,'k-',t,m_eqn,'r-')
    ax.set_xlabel('t')
    ax.set_ylabel('m, m_eqn')

    plt.show()