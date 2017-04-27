from pylab import *
from scipy.interpolate import RectBivariateSpline
from read_mist_models import *

#this creates two 2-D interpolants for C_T and C_L
def create_interpolants(npz):
    data=load(npz)
    C_L=data['C_L']
    C_T=data['C_T']
    omega=data['omega']
    inclination=data['inclination']
    f_T=RectBivariateSpline(y=omega, x=inclination, z=C_T)
    f_L=RectBivariateSpline(y=omega, x=inclination, z=C_L)
    return f_T, f_L

#read in an isochrone set
x=ISO('MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_full.iso')

f_T, f_L = create_interpolants('GD.npz')

#choose one isochrone for demonstration
iso=x.isos[x.age_index(7.5)]
i=squeeze(where( (iso['EEP']<808) & (iso['EEP']>202) ))
iso=iso[i]
omega=iso['surf_avg_omega_div_omega_crit']
T=pow(10,iso['log_Teff'])
L=pow(10,iso['log_L'])

figure(1)
#plot the intrinsic model quantities
plot(T,log10(L), color='Lime', label='Intrinsic')

#plot the gravity-darkened models at i=0, "edge-on"
C_T=f_T.ev(yi=omega, xi=zeros(len(omega)))
C_L=f_L.ev(yi=omega, xi=zeros(len(omega)))
plot(C_T*T, log10(C_L*L), color='Blue', label=r'$i=0$')

#plot the gravity-darkened models at i=90 deg, "face-on"
C_T=f_T.ev(yi=omega, xi=(pi/2)*ones(len(omega)))
C_L=f_L.ev(yi=omega, xi=(pi/2)*ones(len(omega)))
plot(C_T*T, log10(C_L*L), color='Red', label=r'$i=\pi/2$')

xlim(22000,2000)
ylim(3,4.75)

xlabel(r'$T_{eff} (K)$', fontsize=18)
ylabel(r'$L/L_{\odot}$', fontsize=18)

legend(loc='lower right', fontsize=18, fancybox=True)

show()
