from pylab import *
from scipy.optimize import root
from scipy.integrate import simps, trapz
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

""" Implements Espinosa Lara & Rieutord 2011 equations for gravity darkening """

"""xticks([0, pi/8, pi/4, 3*pi/8, pi/2], ['0', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2'], fontsize=18"""

#constants
G=6.67428e-8
sigma=5.67040e-5
f23=2.0/3.0
Lsun=3.8418e33
Msun=1.989e33
Rsun=6.96e10

#equations
#gives the value of phi
def eq24(phi,theta,omega,rtw):
    tau = (pow(omega,2) * pow(rtw*cos(theta),3) )/3.0 + cos(theta) + log(tan(0.5*theta))
    return cos(phi) + log(tan(0.5*phi)) - tau


#solve for rtw given omega
def eq30(rtw,theta,omega):
    w2=omega*omega
    return (1./w2)*(1./rtw - 1.0) + 0.5*(pow(rtw*sin(theta),2) - 1.0)


#ratio of equatorial to polar Teff
def eq32(omega):
    w2=omega*omega
    return sqrt(2./(2.+w2))*pow(1.-w2, 1./12.)*exp(-(4./3.)*w2/pow(2+w2, 3))



def stuff(omega,theta): #eq.26, 27, 28; returns rtw, Fw, Teff_ratio
    """calculates r~, factors multiplying geff and Teff, and Fw"""

    if theta > pi/2:
        theta -= pi/2
    if theta < 0:
        theta = abs(theta)

    if omega==0.0: #common sense
        return ones(4)
    
    else: 
        q = root(fun=eq30,args=(theta, omega), x0=1.0)
        rtw = asscalar(q['x'])

        w2r3=pow(omega,2)*pow(rtw,3)

        if theta==0.0: #pole, eq. 27
            Fw = exp( f23 * w2r3 )

        elif theta==0.5*pi: #equator, eq. 28
            Fw = pow(1.0 - w2r3, -f23)

        else: #general case
            q = root(fun=eq24,args=(theta, omega, rtw), x0=theta)
            phi = asscalar(q['x'])
            
            Fw = pow(tan(phi)/tan(theta), 2)

        term1 = pow(rtw,-4)
        term2 = pow(omega,4)*pow(rtw*sin(theta),2)
        term3 = -2*pow(omega*sin(theta),2)/rtw
        gterm = term1+term2+term3
        geff_ratio = sqrt( gterm )
        Teff_ratio = pow(Fw*geff_ratio,0.25)
        
        return rtw, geff_ratio, Teff_ratio, Fw

def star_stuff(omega,theta,L,M,Re):
    rtw, geff_ratio, Teff_ratio, Fw = stuff(omega,theta)
    eta = L/(4*pi*G*M)
    R=rtw*Re
    geff = G*M*pow(Re,-2)*geff_ratio
    Teff = pow( L/(4*pi*sigma*Re*Re), 0.25) * Teff_ratio
    Flux = eta*Fw*geff
    return R, geff, Teff, Flux


def ellipse(Re,Rp,i): #from B&H 2015
    Re2=Re*Re
    Rp2=Rp*Rp
    b=sqrt((Re2*Rp2)/(Re2*pow(cos(i),2)+Rp2*pow(sin(i),2)))
    a=Re

    return a, b, pi*a*b




#convenience function for calculating the ratio of Rp to Re
def Rp_div_Re(omega,type='numerical'):
    if type=='numerical':
        y= 0.7761742 +  0.22298167*cos(  2.07345414 *omega)
    elif type=='analytical':
        y= (7. + 2.*cos(2.*pi*omega/3.))/9.0
    return y


#plot one spheroid with Teff variation across the surface
def plot_one_spheroid(omega):
    n=100
    
    def get_one_T(omega,theta):
        r,g,T,F=stuff(omega,theta)
        return T
    
    Re=1.0
    Rp=Rp_div_Re(omega,type='analytical')*Re
   
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
        r,g,T,F=stuff(omega,theta)
        return T

    fig=figure(figsize=(15,9))
    
    w=linspace(0,0.9,m)

    for i,omega in enumerate(w):
        Re=1.0
        Rp=Rp_div_Re(omega,type='analytical')*Re

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

            print omega, T.min(), T.max()

        Tmax = 1.25
        Tmin = 0.87
        T = (Tmax-T)/(Tmax-Tmin)

        for q in range(3):
            ax = fig.add_subplot(3,m,q*m+i+1, projection='3d', aspect='equal')
            ax.view_init(elev=q*45.,azim=0)
            ax.plot_surface(x,y,z,rstride=3,cstride=3,facecolors=cm.viridis_r(T))
            ax.set_xlim(-1,1)
            ax.set_ylim(-1,1)
            ax.set_zlim(-1,1)
            ax.set_axis_off()
            title(r'$\omega$='+str(omega)+', $i$='+str(q*45)+'$^{\circ}$')

    tight_layout()
    show()


    
