import numpy as np

Tint = 400.0
Rd = 3568.0
p = 1e5
cp = 12039.0
g = 10.0**(4.5)/100.0

print(g)

gam = 0.1
sb = 5.670374419e-8
F = sb * Tint**4
nf = 40.0

tau_rad = p/g * (cp/(4.0*sb*Tint**3)) 

dzdH = np.cbrt(F/(Rd*p)) * (cp/Tint)**(1.0/6.0) * (1.0/np.sqrt(gam))

famp = dzdH * Tint/(tau_rad*np.sqrt(nf))

print(tau_rad,dzdH,famp)
