#
# water_states.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2015-02-12 -- created
# 2015-11-13 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script provides various functions of state equations for calculating
# the density and viscosity of water.
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import numpy
import matplotlib.pyplot as plt


###############################################################################
## FUNCTIONS:
###############################################################################
def tumlirz_density(tc, p):
    """
    Name:     tumlirz_density
    Input:    - float, ambient temperature (tc), degrees C
              - float, ambient pressure (p), Pa
    Output:   float, density of water (rho), kg m^3
    Features: Calculates density of water at a given temperature and pressure.
    Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of pure
              water and sea water, Tech. Rept., Marine Physical Laboratory,
              San Diego, CA.
    """
    # Calculate lambda, (bar cm^3)/g:
    my_lambda = 1788.316
    my_lambda += 21.55053*tc
    my_lambda += -0.4695911*tc*tc
    my_lambda += (3.096363e-3)*tc*tc*tc
    my_lambda += -(7.341182e-6)*tc*tc*tc*tc
    #
    # Calculate pure water pressure (po), bar
    po = 5918.499
    po += 58.05267*tc
    po += -1.1253317*tc*tc
    po += (6.6123869e-3)*tc*tc*tc
    po += -(1.4661625e-5)*tc*tc*tc*tc
    #
    # Calculate vinf, cm^3/g
    vinf = 0.6980547
    vinf += -(7.435626e-4)*tc
    vinf += (3.704258e-5)*tc*tc
    vinf += -(6.315724e-7)*tc*tc*tc
    vinf += (9.829576e-9)*tc*tc*tc*tc
    vinf += -(1.197269e-10)*tc*tc*tc*tc*tc
    vinf += (1.005461e-12)*tc*tc*tc*tc*tc*tc
    vinf += -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc
    vinf += (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc
    vinf += -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc
    #
    # Convert pressure to bars (1 bar = 100000 Pa)
    pbar = (1e-5)*p
    #
    # Calculate the specific volume (cm^3 g^-1):
    v = vinf + my_lambda/(po + pbar)
    #
    # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
    rho = (1e3/v)
    #
    return rho


def chen_density(tc, p):
    """
    Name:     chen_density
    Input:    - float, air temperature (tc), degrees C
              - float, atmospheric pressure (p), Pa
    Output:   float, density of water, kg/m^3
    Features: Calculates density of water at a given temperature and pressure.
    Ref:      Chen, C.-T., R. A. Fine, and F. J. Millero (1977) The equation of
              state of pure water determined from sound speeds, J. Chem. Phys.,
              66(5), 2142--2144.
    """
    # Calculate density at 1 atm:
    po = 0.99983952
    po += (6.788260e-5)*tc
    po += -(9.08659e-6)*tc*tc
    po += (1.022130e-7)*tc*tc*tc
    po += -(1.35439e-9)*tc*tc*tc*tc
    po += (1.471150e-11)*tc*tc*tc*tc*tc
    po += -(1.11663e-13)*tc*tc*tc*tc*tc*tc
    po += (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc
    po += -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
    #
    # Calculate bulk modulus at 1 atm:
    ko = 19652.17
    ko += 148.1830*tc
    ko += -2.29995*tc*tc
    ko += 0.01281*tc*tc*tc
    ko += -(4.91564e-5)*tc*tc*tc*tc
    ko += (1.035530e-7)*tc*tc*tc*tc*tc
    #
    # Calculate temperature dependent coefficients:
    ca = 3.26138
    ca += (5.223e-4)*tc
    ca += (1.324e-4)*tc*tc
    ca += -(7.655e-7)*tc*tc*tc
    ca += (8.584e-10)*tc*tc*tc*tc
    #
    cb = (7.2061e-5)
    cb += -(5.8948e-6)*tc
    cb += (8.69900e-8)*tc*tc
    cb += -(1.0100e-9)*tc*tc*tc
    cb += (4.3220e-12)*tc*tc*tc*tc
    #
    # Convert atmospheric pressure to bar (1 bar = 100000 Pa)
    pbar = (1e-5)*p
    #
    pw = (ko + ca*pbar + cb*pbar**2.0)
    pw /= (ko + ca*pbar + cb*pbar**2.0 - pbar)
    pw *= (1e3)*po
    #
    return pw


def huber_viscosity(tc, p):
    """
    Name:     huber_viscosity
    Input:    - float, ambient temperature (tc), degrees C
              - float, ambient pressure (p), Pa
    Return:   float, viscosity of water (mu), Pa s
    Features: Calculates water viscosity at a given temperature and pressure.
    Depends:  tumlirz_density
    Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V.
              Sengers, M. J. Assael, ..., K. Miyagawa (2009) New
              international formulation for the viscosity of H2O, J. Phys.
              Chem. Ref. Data, Vol. 38(2), pp. 101-125.
    """
    # Define reference temperature, density, and pressure values:
    tk_ast = 647.096      # Kelvin
    rho_ast = 322.0       # kg/m^3
    mu_ast = (1e-6)       # Pa s
    #
    # Get the density of water, kg/m^3
    rho = tumlirz_density(tc, p)
    #
    # Calculate dimensionless parameters:
    tbar = (tc + 273.15)/tk_ast
    tbarx = tbar**(0.5)
    tbar2 = tbar**2
    tbar3 = tbar**3
    rbar = rho/rho_ast
    #
    # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
    mu0 = 1.67752
    mu0 += 2.20462/tbar
    mu0 += 0.6366564/tbar2
    mu0 += -0.241605/tbar3
    mu0 = 1e2*tbarx/mu0
    #
    # Create Table 3, Huber et al. (2009):
    hj0 = (0.520094, 0.0850895, -1.08374, -0.289555, 0., 0.)
    hj1 = (0.222531, 0.999115, 1.88797, 1.26613, 0., 0.120573)
    hj2 = (-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.)
    hj3 = (0.161913,  0.257399, 0., 0., 0., 0.)
    hj4 = (-0.0325372, 0., 0., 0.0698452, 0., 0.)
    hj5 = (0., 0., 0., 0., 0.00872102, 0.)
    hj6 = (0., 0., 0., -0.00435673, 0., -0.000593264)
    h = hj0 + hj1 + hj2 + hj3 + hj4 + hj5 + hj6
    h_array = numpy.reshape(numpy.array(h), (7, 6))
    #
    # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
    mu1 = 0.
    ctbar = (1./tbar) - 1.
    for i in xrange(6):
        coef1 = numpy.power(ctbar, i)
        coef2 = 0.
        for j in xrange(7):
            coef2 += h_array[j][i]*numpy.power((rbar - 1.), j)
        mu1 += coef1*coef2
    mu1 = numpy.exp(rbar*mu1)
    #
    # Calculate mu_bar (Eq. 2, Huber et al., 2009)
    #   assumes mu2 = 1
    mu_bar = mu0*mu1
    #
    # Calculate mu (Eq. 1, Huber et al., 2009)
    mu = mu_bar*mu_ast    # Pa s
    #
    return mu


def vogel_viscosity(tc):
    """
    Name:     vogel_viscosity
    Input:    float, ambient temperature (tc), degrees C
    Return:   float, viscosity of water, Pa s
    Features: Calculates the viscosity of water at a given temperature.
    Ref:      Vogel Coefficients for Water (mPa s): n = exp[A + B/(C + T)]
                 where: A = -3.7188;   B = 578.919;   C = -137.546
    """
    mu = 578.919
    mu /= ((tc + 273.15) - 137.546)
    mu += -3.7188
    mu = numpy.exp(mu)
    mu *= 1e-3
    return mu

###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':
    n = 61
    my_temps = [i - 20. for i in xrange(n)]
    my_chen = numpy.zeros((n,))
    my_tumlirz = numpy.zeros((n,))
    my_huber = numpy.zeros((n,))
    my_vogel = numpy.zeros((n,))
    my_press = 101325.
    for i in xrange(n):
        my_chen[i] = chen_density(my_temps[i], my_press)
        my_tumlirz[i] = tumlirz_density(my_temps[i], my_press)
        my_huber[i] = huber_viscosity(my_temps[i], my_press)
        my_vogel[i] = vogel_viscosity(my_temps[i])

    # Density:
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.plot(my_temps, my_chen, 'r-', label='Chen et al. (1977)')
    ax1.plot(my_temps, my_tumlirz, 'k--', label='Tumlirz Equation')
    ax1.set_ylabel('Density of water, kg m$^{-3}$', fontsize=16)
    ax1.set_xlabel('Temperature, $^{\circ}$C', fontsize=16)
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)
    plt.show()

    # Viscosity
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.plot(my_temps, 1e3*my_huber, 'r-', label='Huber et al. (2009)')
    ax1.plot(my_temps, 1e3*my_vogel, 'k--', label='Vogel Equation')
    ax1.set_ylabel('Viscosity of water, mPa s', fontsize=16)
    ax1.set_xlabel('Temperature, $^{\circ}$C', fontsize=16)
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)
    plt.show()
