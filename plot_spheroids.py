"""GDit: code for plotting a grid of spheroids"""

__version__ = '0'
__author__ = 'Aaron Dotter'

from pylab import *
from scipy.optimize import root
from scipy.integrate import simps, trapz
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from GDit import *

#plot one spheroid with Teff variation across the surface
def plot_one_spheroid(omega):
    n=100
    
    def get_one_T(omega,theta):
        r,g,T,F=solve_ELR(omega,theta)
        return T
    
    Re=1.0
    Rp=Rp_div_Re(omega)*Re
   
    # Set of all spherical angles:
    phi = linspace(-pi, pi, n)
    nu = linspace(-pi/2, pi/2, n)

    angles=outer(ones_like(phi),nu)
    T=outer(ones_like(phi),ones_like(nu))
    for j in range(n):
        theta=angles[0,j]
        T[:,j]=get_one_T(omega,theta)
        
    if(Rp==Re): #spherical
        a=Re
        x=a*outer(cos(phi),cos(nu))
        y=a*outer(sin(phi),cos(nu))
        z=a*outer(ones_like(phi), sin(nu))
        
        T=outer(ones_like(phi), ones_like(nu))
        
    else: #spheroidal
        mu=arctanh(Rp/Re)
        a=sqrt(Re*Re - Rp*Rp)

        # Cartesian coordinates that correspond to the spherical angles:
        # (this is the equation of an ellipsoid):
        x = a*cosh(mu) * outer(cos(phi), cos(nu))
        y = a*cosh(mu) * outer(sin(phi), cos(nu))
        z = a*sinh(mu) * outer(ones_like(phi), sin(nu))

        Tmin=T.min()
        Tmax=T.max()
        T=(T-Tmin)/(Tmax - Tmin)

        
    fig = figure(figsize=figaspect(1))  # Square figure
    ax = fig.add_subplot(111, projection='3d', aspect='equal')
    ax.plot_surface(x,y,z,rstride=4,cstride=4,facecolors=cm.plasma_r(T) )
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(-1,1)
    title(r'$\omega$='+str(omega))
    show()




#plot several spheroids with Teff variation across the surface
def plot_more_spheroids():
    n=100
    m=7
    
    def get_one_T(omega,theta):
        r,T,F=solve_ELR(omega,theta)
        return T

    fig=figure(figsize=(15,9))
    
    w=linspace(0,0.9,m)

    for i,omega in enumerate(w):
        Re=1.0
        Rp=Rp_div_Re(omega)*Re

        # Set of all spherical angles:
        phi = linspace(-pi, pi, n)
        nu = linspace(-pi/2, pi/2, n)

        angles=outer(ones_like(phi),nu)
        T=outer(ones_like(phi),ones_like(nu))
        for j in range(n):
            theta=angles[0,j]
            T[:,j]=get_one_T(omega,theta)

        if(Rp==Re): #spherical
            a=Re
            x=a*outer(cos(phi),cos(nu))
            y=a*outer(sin(phi),cos(nu))
            z=a*outer(ones_like(phi), sin(nu))

            T=outer(ones_like(phi), ones_like(nu))

        else: #spheroidal
            mu=arctanh(Rp/Re)
            a=sqrt(Re*Re - Rp*Rp)

            # Cartesian coordinates that correspond to the spherical angles:
            # (this is the equation of an ellipsoid):
            x = a*cosh(mu) * outer(cos(phi), cos(nu))
            y = a*cosh(mu) * outer(sin(phi), cos(nu))
            z = a*sinh(mu) * outer(ones_like(phi), sin(nu))

            #print omega, T.min(), T.max()

        Tmax = 1.25
        Tmin = 0.87
        T = (Tmax-T)/(Tmax-Tmin)

        for q in range(3):
            ax = fig.add_subplot(3,m,q*m+i+1, projection='3d', aspect='equal')
            ax.view_init(elev=q*45.,azim=0)
            s=ax.plot_surface(x,y,z,rstride=3,cstride=3,facecolors=cm.magma(T) )
            ax.set_xlim(-1,1)
            ax.set_ylim(-1,1)
            ax.set_zlim(-1,1)
            ax.set_axis_off()
            title(r'$\omega$='+str(omega)+' $, i=$'+str(q*45))
            
    tight_layout()
    show()
    return s
