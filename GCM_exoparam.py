import numpy as np
import configparser
from scipy import interpolate

def read_parameters(file_path):
    parameters = {}
    current_section = None

    with open(file_path, 'r') as file:
        for line in file:
            stripped_line = line.strip()
            
            # Check if the line is a section header
            if stripped_line.startswith("#"):
                current_section = stripped_line[1:].strip()
                parameters[current_section] = {}
            else:
                # Split the line into key and value
                if '=' in stripped_line:
                    key, value = stripped_line.split('=')
                    key = key.strip()
                    value = float(value.strip())
                    # Add the key-value pair to the current section
                    if current_section is not None:
                        parameters[current_section][key] = value
                    else:
                        parameters[key] = value

    return parameters

# Physical constants
GM = 1.2668653e17
G =  6.67430e-11
au = 1.495978707e11
Rj = 7.1492e7
Re = 6.3781e6
Me = 5.9722e24
Mj = GM/G
Rsun = 6.95700e8
Msun = 1.98847e30
sb = 5.670374419e-8
daysec = 86400.0


# Read input file
file_path = 'input.txt'
parameters = read_parameters(file_path)
print('Input dictionary: ')
print(parameters)
print('----')

## Calculate bas

##System parameters
# Stellar effective temperature
Teff = parameters['Stellar parameters']['Teff']
# Radius of star
Rs = parameters['Stellar parameters']['Rs'] * Rsun
# Semi-Major axis
a = parameters['Planetary parameters']['a'] * au
# Mass of planet
Mp =  parameters['Planetary parameters']['Mp'] * Mj
# Radius of planet
Rp = parameters['Planetary parameters']['Rp'] * Rj
# Orbital period (days)
P = parameters['Planetary parameters']['P']  * daysec 


print('Basic values for exoplanet GCMs: ')

print('Radius of planet [m], [R_jup]: ', Rp, Rp/Rj)

omega = (2.0 * np.pi)/P
print('rotation rate (omega) [rad s-1]: ', omega)

Tirr = Teff * np.sqrt(Rs/a)
print('irradition temperature (T_irr) [K]: ', Tirr)

Teq = Teff * np.sqrt(Rs/(2.0*a))
print('equilibrium temperature (T_eq) [K]: ', Teq)

gsurf = (G * Mp) / Rp**2
print('Surface gravity from Mass, Radius relation (g_surf) [m s-2]: ', gsurf)

sb_G = 5.670374419e-14
F_Thorn = sb_G * Tirr**4
Teq_Thorn = (F_Thorn/4.0/sb_G)**0.25
Tint_Thorn = 0.39 * Teq_Thorn * np.exp(-((np.log10(F_Thorn) - 0.14)**2 / 1.0952))
Fint_Thorn = sb * Tint_Thorn**4
print('Thorngren expression internal temperature (T_int) [K]: ', Tint_Thorn)

print('----')

# Calculate specific gas constants Rd, cp and kappa

print('Finding Rd, cp and kappa for input values ')

T_in = parameters['Gas']['T']
p_in = parameters['Gas']['p']
met = parameters['Gas']['met']

# Read the CE table
ifile = 'CE_tables/mini_chem_IC_FastChem_'+str(int(met))+'x.txt'

data = np.loadtxt(ifile,skiprows=2,max_rows=1)

Tf = data[:]

data = np.loadtxt(ifile,skiprows=3,max_rows=1)

P = data[:]

lP = np.log10(P)
lT = np.log10(Tf)

nP = len(lP)
nT = len(lT)

data =  np.loadtxt(ifile,skiprows=4)

# Mean molecular weight interpolator
# Convert 1D to 3D arrays
mul = data[:,0]
mu = np.reshape(mul, (nT, nP))
# Create a scipy interpolation function for the mean molecular weight
f_mu = []
f_mu.append(interpolate.interp2d(lP[:], lT[:], mu[:,:], kind='linear'))

# VMR interpolator 
sp = ['OH', 'H2', 'H2O' , 'H',  'CO',  'CO2', 'O', 'CH4', 'C2H2', 'NH3', 'N2', 'HCN', 'He']
nsp = len(sp)
f_VMR = []
for i in range(nsp):
  # Convert 1D to 3D arrays
  VMRl = data[:,1+i]
  VMR = np.log10(np.reshape(VMRl, (nT, nP)))
  # Create a scipy interpolation function for the VMR for each species
  f_VMR.append(interpolate.interp2d(lP[:], lT[:], VMR[:,:], kind='linear'))


VMR = np.zeros(nsp)
for i in range(nsp):
  VMR[i] = 10.0**f_VMR[i](np.log10(p_in),np.log10(T_in))

mu = f_mu[0](np.log10(p_in),np.log10(T_in))

print('----')
print('T_in [K]: ', T_in)
print('p_in [bar]: ', p_in )
print('species and VMR: ')
for i in range(nsp):
  print(i, sp[i], VMR[i])


print('----')
print('mean molecular weight: ')
print(mu[0])
print('----')

# Read in NASA polynomial

fname = 'NASA9/NASA_9_poly.txt'

# Arrays for data
mw = np.zeros(nsp)
a_l = np.zeros((nsp,9))
a_h = np.zeros((nsp,9))

# Read the file line by line
f = open(fname,'r')

# Skip first few rows
for i in range(8):
  line = f.readline()

for i in range(nsp):
  line1 = f.readline().split()
  mw[i] = float(line1[1])
  line2 = f.readline().split()
  a_l[i,:] = line2[:]
  line3 = f.readline().split()
  a_h[i,:] = line3[:]
  #print(mw[i])
  #print(a_l[i,:])
  #print(a_h[i,:])


# # Calculate Rd, cp and kappa
T = T_in
T2 = T**2
T3 = T**3
T4 = T**4

R = 8.31446261815324

#Main loop calculations
R_bar = 0.0
cp_bar = 0.0
for i in range(nsp):

  # Mass mixing ratio = VMR * molecular weight / mean molecular weight
  mmr = VMR[i] * mw[i]/mu[0]

  # Contribution to R_bar = mmr * specific gas constant
  R_bar = R_bar + mmr * R/mw[i]

  if T_in <= 1000.0:
    cp_val = a_l[i,0]/T2 + a_l[i,1]/T + a_l[i,2] + a_l[i,3]*T + a_l[i,4]*T2 + a_l[i,5]*T3 + a_l[i,6]*T4
  else:
    cp_val = a_h[i,0]/T2 + a_h[i,1]/T + a_h[i,2] + a_h[i,3]*T + a_h[i,4]*T2 + a_h[i,5]*T3 + a_h[i,6]*T4

  cp_val = cp_val*R

  # Contribution to cp_bar = mmr * cp_val / molecular weight
  cp_bar = cp_bar + mmr * cp_val/mw[i]


# Convert R_bar [J g-1 K-1] to SI units [J kg-1 K-1]
R_bar = R_bar * 1000.0

# Convert cp_bar [J g-1 K-1] to SI units [J kg-1 K-1]
cp_bar = cp_bar * 1000.0

# kappa_prime evaluation
k_prime = R_bar / cp_bar

print('specific gas constant Rd [J kg-1 K-1]: ', R_bar)
print('specific gas heat capacity cp [J kg-1 K-1]: ', cp_bar)
print('adibatic coefficent kappa [-]: ', k_prime)

print('----')
