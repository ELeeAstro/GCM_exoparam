# GCM_exoparam
Python script for finding useful parameters for exoplanet GCMs

To run, enter 

python GCM_exoparam.py

This generates useful values for GCMs and other exoplanet models. And then interpolates chemical equilibrium tables and uses NASA polynomials to calculate the specific gas constant, specific heat capacity and adibiatic constant at the gas temperature, pressure and metallicity. The mean molecular weight is also calculated.

# Inputs

Edit the input.txt file to change the input.
To get the Rd, cp and kappa, edit the T, p and met values in the 'Gas' subheading

### Stellar parameters
Teff = Stellar effective temperature [K] \
Rs = Stellar radius [Rsun] 

### Planetary parameters
Rp = Radius of planet [Rj] \
Mp = Mass of planet [Mj] \
P = Orbital period [days] \
a = Semi-major axis [au]

### Gas
T = Temperature of gas [K] \
p = Pressure of gas [bar] \
met = Metallicity [1,10,100,500]

Can only do 1x, 10x, 100x, 500x solar metallicity currently

TO DO: Change depreciating scipy interpolation functions. Add more useful auxiluary values. Allow interpolation to arbitary metalliacity.
