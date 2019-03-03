'''
pablo's roche model for gravity darkening

working in dimensionless variables

re=omega^(2/3)
x=r*sin(theta)
y=r*cos(theta)




'''

from pylab import *

tiny=1.0e-6
f23=2./3.

def y(x,w): #from equation 5
    a=pow(w,-2./3.)
    b=0.5*pow(w,4./3.)
    x2=x*x
    c=-0.5*x*x
    q=pow(a+b+c,-2)
    return sqrt(q-x2)

def ellipse(x,w):
    re=pow(w,f23)
    re_div_rp = 1 + 0.5*w*w
    return sqrt(re*re - x*x)/re_div_rp


close(1)
figure(1,figsize=(10,10))
for w in linspace(0,1,11):
    X=linspace(0,pow(w,f23)-tiny,1000)
    Y=y(X,w)
    Z=ellipse(X,w)
    plot(X,Y,color='Red')
    plot(X,Z,color='Blue',ls='--')
xlabel(r'${\rm x^{\prime}}$',fontsize=32)
ylabel(r'${\rm y^{\prime}}$',fontsize=32)
xticks(fontsize=32)
yticks(fontsize=32)
xlim(0,1)
ylim(0,1)
text(0.6,0.8,'Spheroid',fontsize=32,color='Blue')
text(0.6,0.9,'Roche',fontsize=32,color='Red')
show()
savefig('Roche_vs_Spheroid.pdf')
