"""Gravity Darkening (GD) module:
   Implements the Espinosa Lara & Rieutord GD model and applies it to an
   oblate spheroidal surface. Includes a process to compute the observed
   Teff and luminosity projected along the line of sight."""

__version__ = '0'
__author__ = 'Aaron Dotter'

from pylab import *
from scipy.optimize import root
from scipy.integrate import simps, trapz
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.interpolate import RectBivariateSpline

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

def stuff(omega,theta): #eq.26, 27, 28; solve the ELR11 equations
    """calculates r~, Teff_ratio, and Flux_ratio"""
    #theta is the polar angle.
    #this routine calculates values for 0 <= theta <= pi/2
    #everything else is mapped into this interval by symmetry
    # theta = 0 at the pole(s)
    # theta = pi/2 at the equator
    # -pi/2 < theta < 0: theta -> abs(theta)
    #  pi/2 > theta > pi: theta -> pi - theta
    if pi/2 < theta <= pi:
        theta = pi - theta
    if -pi/2 <= theta < 0:
        theta = abs(theta)

    if omega==0.0: #common sense
        return ones(3)
    
    else:
        #first we solve equation 30 for rtw
        q = root(fun=eq30,args=(theta, omega), x0=1.0)
        rtw = asscalar(q['x'])

        #the following are special solutions for extreme values of theta
        w2r3=pow(omega,2)*pow(rtw,3)

        if theta==0.0: #pole, eq. 27
            Fw = exp( f23 * w2r3 )

        elif theta==0.5*pi: #equator, eq. 28
            Fw = pow(1.0 - w2r3, -f23)

        else: #general case for Fw
            q = root(fun=eq24,args=(theta, omega, rtw), x0=theta)
            phi = asscalar(q['x'])
            
            Fw = pow(tan(phi)/tan(theta), 2)

        #equation 31 and similar for Fw
        term1 = pow(rtw,-4)
        term2 = pow(omega,4)*pow(rtw*sin(theta),2)
        term3 = -2*pow(omega*sin(theta),2)/rtw
        gterm = sqrt(term1+term2+term3)
        Flux_ratio = Fw*gterm
        Teff_ratio = pow(Flux_ratio,0.25)
        
        return rtw, Teff_ratio, Flux_ratio

#convenience functions
def Rp_div_Re(omega):
    rtw, Teff_ratio, Flux_ratio = stuff(omega, theta=0)
    return rtw

def R_div_Re(omega):
    return cbrt(Rp_div_Re(omega))

def Re_div_R(omega):
    return 1.0/R_div_Re(omega)

def Rp_div_R(omega):
    return pow(Rp_div_Re(omega), f23)

def R_div_Rp(omega):
    return 1.0/Rp_div_R(omega)

#analytical formulas for the surface area and volume of oblate spheroid
def spheroid_surface_area(Rp,Re):
    #c=polar, a=equatorial; c/a=Rp/Re
    c_div_a=Rp / Re
    a=Re
    if c_div_a == 1.0 : #spherical
        extra = 1.0
    else:
        e=sqrt(1-pow(c_div_a,2))
        extra = (1-e*e)*arctanh(e)/e
    return 2*pi*a*a*(1+extra)

def spheroid_volume(Rp,Re):
    return 4*pi*Rp*Re*Re/3

#from Binggeli et al. (1980) the following calculates the ratio of semi-major and -minor axes
#of the ellipse resulting from the projection of an oblate spheroid along the line of sight
def beta(theta,q): #q is ratio Rp/Re; Binggeli et al. 1980
    j= pow(q*sin(theta),2) + pow(cos(theta),2)
    l=1.0
    top=j+l-sqrt( pow(j-l,2))
    bottom=j+l+sqrt( pow(j-l,2))
    return sqrt(top/bottom)

#the ares of the projected ellipse; inclication angle is defined differently in this case
def ellipse(Rp, Re, i): 
    b=beta(theta=pi/2-i,q=Rp/Re)
    return pi*Re*Re*b

def geometric_factors(omega, i, n_nu=50, n_phi=50, do_checks=False):
    """solves for geometric factors C_T and C_L for arbitrary omega, and inclination"""
    if omega==0: #perfect sphere
        return ones(2)
    
    #line of sight, inclination angle i
    LOS = array([cos(i), 0, sin(i)])

    #spherical angles
    phi_array = linspace(-pi, pi, n_phi)
    dphi=diff(phi_array)[0]
    nu_array = linspace(-pi/2, pi/2, n_nu)
    dnu=diff(nu_array)[0]

    #find the polar radius
    Re=1.
    Rp=Rp_div_Re(omega)

    mu=arctanh(Rp/Re)

    a=sqrt(Re*Re - Rp*Rp)

    sinh_mu = sinh(mu)
    cosh_mu = cosh(mu)

    #now do area integrals
    nu_int1=empty(n_nu)
    nu_int2=empty(n_nu)
    nu_int3=empty(n_nu)
    nu_int4=empty(n_nu)

    #loop over polar angle nu
    for j,nu in enumerate(nu_array):

        theta = pi/2 - nu
        r,T,F=stuff(omega=omega,theta=theta)
        
        cos_nu = cos(nu)
        sin_nu = sin(nu)

        #scale factors
        K = sqrt(pow(sinh_mu,2) + pow(sin_nu,2))
        h_nu = a*K
        h_phi = a*cosh_mu*cos_nu
        
        phi_int1 = empty(n_phi)
        phi_int2 = empty(n_phi)

        #loop over azimuth angle phi
        for k,phi in enumerate(phi_array):
            e_mu_x = sinh_mu * cos_nu * cos(phi)
            e_mu_y = sinh_mu * cos_nu * sin(phi)
            e_mu_z = cosh_mu * sin_nu
            e_mu = array([e_mu_x, e_mu_y, e_mu_z]) / K

            dArea = h_nu * h_phi
            proj = max(0,dot(e_mu,LOS)) #only projected toward observer
            phi_int1[k] = dArea
            phi_int2[k] = dArea*proj

        #phi integral at constant nu
        nu_int1[j] = simps(y=phi_int1,dx=dphi) #total surface area
        nu_int2[j] = simps(y=phi_int2,dx=dphi) #projected area
        nu_int3[j] = F*nu_int2[j] #integral of Flux_ratio, projected
        nu_int4[j] = F*nu_int1[j] #integral of Flux_ratio, total

    #nu integrals
    surface_area = simps(y=nu_int1, dx=dnu)
    projected_area = simps(y=nu_int2, dx=dnu)
    total_flux = simps(y=nu_int4, dx=dnu)
    
    C_L = 4 * simps(y=nu_int3, dx=dnu) / total_flux
    C_T = pow( C_L * surface_area / projected_area / 4, 0.25 )

    if do_checks:
        surface_area_check = abs( surface_area / spheroid_surface_area(Rp,Re) - 1.0)
        projected_area_check = abs( projected_area / ellipse(Rp, Re, i) - 1.0)

        if surface_area_check > 0.001 or projected_area_check > 0.001:
            print(' omega, i = {0:n}, {1:n}'.format(omega,i))
            print(' surface area ratio: {0:n}'.format(surface_area_ratio))
            print(' projected area ratio: {0:n}'.format(projected_area_ratio))
        
    return C_T, C_L

#compute the coefficients C_T and C_L on an nxn matrix
def save_coefficients(n,output='GD.npz'):
    omega=linspace(0,1,n)
    inclination=linspace(0,pi/2,n)
    C_T=empty((n,n))
    C_L=empty((n,n))
    for i,omega in enumerate(omega):
        for j,incl in enumerate(inclination):
            C_T[j,i], C_L[j,i] = geometric_factors(omega, incl)
    savez(output, C_T=C_T, C_L=C_L)

def print_coefficients(n,output='data.txt'):
    omega=linspace(0,1,n)
    inclination=linspace(0,pi/2,n)
    with open(output,'w') as f:
        for i,w in enumerate(omega):
            for j,incl in enumerate(inclination):
                C_T, C_L = geometric_factors(w, incl)
                f.write('{0:12.8f} {1:12.8f} {2:12.8f} {3:12.8f}\n'.format(w,incl,C_T,C_L))

    
#returns instance of the BivariateRect
def create_interpolants(npz):
    data=load(npz)
    C_L=data['C_L']
    C_T=data['C_T']
    omega=data['omega']
    inclination=data['inclination']
    f_T=RectBivariateSpline(y=omega, x=inclination, z=C_T)
    f_L=RectBivariateSpline(y=omega, x=inclination, z=C_L)
    return f_T, f_L

    
