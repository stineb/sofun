# RStudio v 0.98
#
# water_states.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-01-30 -- created
# 2014-10-30 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script provides various functions of state equations for calculating 
# the density and viscosity of water
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### FUNCTIONS ################################################################
# /////////////////////////////////////////////////////////////////////////////

# ************************************************************************
# * Name: tumlirz_density
# *
# * Input: ambient temperature, degrees C (tc)
# *        ambient pressure, Pa (P)
# *
# * Return: density of water, rho, in kg m^3
# *
# * Features: This function calculates the density of water based on 
# *           temperature and pressure relationships defined in the 
# *           Tumlirz Equation:
# *
# *           V = Vinf - K1 S + lamda (Po + K2 S + P)^-1
# *           where:
# *             V :: specific volume (1/density), cm^3 g^-1
# *             S :: salinity (assumed zero for pure water), ppm
# *             P :: pressure, bar
# *
# * Ref:  F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of pure
# *         water and sea water, Tech. Rept., Marine Physical Laboratory,
# *         San Diego, CA.
# ************************************************************************
tumlirz_density <- function(tc, P){
  # Calculate lamda, (bar c^3)/g:
  L1 <- 1788.316
  L2 <- 21.55053
  L3 <- -0.4695911
  L4 <- (3.096363e-3)
  L5 <- -(7.341182e-6)
  #
  lamda <- L1 +
    L2*tc +
    L3*tc*tc +
    L4*tc*tc*tc +
    L5*tc*tc*tc*tc
  #
  # Calculate Po, bar
  P1 <- 5918.499
  P2 <- 58.05267
  P3 <- -1.1253317
  P4 <- (6.6123869e-3)
  P5 <- -(1.4661625e-5)
  #
  Po <- P1 +
    P2*tc +
    P3*tc*tc +
    P4*tc*tc*tc +
    P5*tc*tc*tc*tc
  #
  # Calculate Vinf, c^3/g
  V1 <- 0.6980547
  V2 <- -(7.435626e-4)
  V3 <- (3.704258e-5)
  V4 <- -(6.315724e-7)
  V5 <- (9.829576e-9)
  V6 <- -(1.197269e-10)
  V7 <- (1.005461e-12)
  V8 <- -(5.437898e-15)
  V9 <- (1.69946e-17)
  V10 <- -(2.295063e-20)
  # 
  Vinf <- V1 +
    V2*tc +
    V3*tc*tc +
    V4*tc*tc*tc +
    V5*tc*tc*tc*tc +
    V6*tc*tc*tc*tc*tc +
    V7*tc*tc*tc*tc*tc*tc +
    V8*tc*tc*tc*tc*tc*tc*tc +
    V9*tc*tc*tc*tc*tc*tc*tc*tc +
    V10*tc*tc*tc*tc*tc*tc*tc*tc*tc
  #
  # Convert pressure to bars (1 bar = 100000 Pa)
  Pbar <- (1e-5)*P
  #
  # Calculate the specific volume (cm^3 g^-1):
  V <- Vinf + lamda/(Po + Pbar)
  #
  # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  rho <- (1000/V)
  #
  return(rho)
}

# ************************************************************************
# * Name: chen_density
# *
# * Input: ambient temperature, degrees C (tc)
# *        ambient pressure, Pa (P)
# *
# * Return: density of water, rho, in kg m^3
# *
# * Features: This function calculates the density of water based on 
# *           temperature and pressure relationships defined in the 
# *           Equation of State:
# *
# *           V = Vo - Vo P / (Ko + A P + B P^2)
# *           where:
# *             V ..... specific volume, cm^3/g
# *             P ..... pressure, bar
# *             Vo .... specific volume at 1 atm
# *             Ko .... secant bulk moduli at 1 atm
# *             A, B .. temperature-dependent parameters
# *
# * Ref: Chen, C.-T., R. A. Fine, and F. J. Millero (1977) The equation 
# *        of state of pure water determined from sound speeds, Journal 
# *        of Chemical Physics, vol. 66, pp. 2142-2144.
# *      Kell, G. S. (1975) J. Chem. Eng. Data, vol. 20, p. 97.
# *
# ************************************************************************
chen_density <- function(tc, P){
  # Vo can be found by using do (1/Vo), g/cm^3, (Kell, 1975):
  D1 <- 0.99983952
  D2 <- (6.78826e-5)
  D3 <- -(9.08659e-6)
  D4 <- (1.02213e-7)
  D5 <- -(1.35439e-9)
  D6 <- (1.47115e-11)
  D7 <- -(1.11663e-13)
  D8 <- (5.04407e-16)
  D9 <- -(1.00659e-18)
  #
  do <- D1 +
    D2*tc +
    D3*tc*tc +
    D4*tc*tc*tc +
    D5*tc*tc*tc*tc +
    D6*tc*tc*tc*tc*tc +
    D7*tc*tc*tc*tc*tc*tc +
    D8*tc*tc*tc*tc*tc*tc*tc +
    D9*tc*tc*tc*tc*tc*tc*tc*tc
  #
  # Calculate Ko, atm
  K1 <- 19652.17
  K2 <- 148.183
  K3 <- -2.29995
  K4 <- 0.01281
  K5 <- -(4.91564e-5)
  K6 <- (1.03553e-7)
  # 
  Ko <- K1 +
    K2*tc +
    K3*tc*tc +
    K4*tc*tc*tc +
    K5*tc*tc*tc*tc +
    K6*tc*tc*tc*tc*tc
  #
  # Calculate A
  A1 <- 3.26138
  A2 <- (5.223e-4)
  A3 <- (1.324e-4)
  A4 <- -(7.655e-7)
  A5 <- (8.584e-10)
  #
  A <- A1 +
    A2*tc +
    A3*tc*tc +
    A4*tc*tc*tc +
    A5*tc*tc*tc*tc 
  #
  # Calculat B
  B1 <- (7.2061e-5)
  B2 <- -(5.8948e-6)
  B3 <- (8.699e-8)
  B4 <- -(1.01e-9)
  B5 <- (4.322e-12)
  #
  B <- B1 +
    B2*tc +
    B3*tc*tc +
    B4*tc*tc*tc +
    B5*tc*tc*tc*tc 
  #
  # Convert pressure to bar (1 bar = 100000 Pa)
  Pbar <- (1e-5)*P
  #
  # Calculate the density (i.e., 1/V), kg/m^3:
  rho <- 1000*do*(Ko + A*Pbar + B*Pbar^2)/(Ko + A*Pbar + B*Pbar^2 - Pbar)
  return(rho)
}

# ************************************************************************
# * Name: huber_visc
# *
# * Input: numeric (vector), ambient temperature, degrees C (tc)
# *        numeric, ambient pressure, Pa (P)
# *
# * Return: numeric (vector), viscosity of water in mPa s (mu)
# *
# * Features: This function calculates the viscosity of water based on 
# *           temperature and pressure relationships defined in the 
# *           new international equation for viscosity:
# *
# *           mu = mu_0(Tbar) * mu_1(rbar,Tbar) * mu_2(rbar,Tbar)
# *           where:
# *             mu_0 :: coefficient function
# *             mu_1 :: ceofficient function
# *             rbar :: rho / rho_ast
# *             Tbar :: TK / TK_ast
# *
# * Depends: chen_density (may be changed with other density function)
# *
# * Ref: Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. 
# *        Sengers, M. J. Assael, ..., K. Miyagawa (2009) New 
# *        international formulation for the viscosity of H2O, J. Phys. 
# *        Chem. Ref. Data, Vol. 38(2), pp. 101-125.
# *
# ************************************************************************
huber_visc <- function(tc, P){
  # Define reference temperature, density, and pressure values:
  TK_ast <- 647.096      # Kelvin
  rho_ast <- 322.0       # kg/m^3
  P_ast <- (2.2064e7)    # Pa
  mu_ast <- (1e-6)       # Pa s
  #
  # Get the density of water:
  rho <- chen_density(tc, P) # kg/m^3
  #
  # Calculate dimensionless parameters:
  Tbar <- (tc + 273.15)/TK_ast
  rbar <- rho/rho_ast
  #
  # Calculate mu_0:
  mu_0 <- 100.0*sqrt(Tbar)/(
    1.67752 + 
      2.20462*Tbar^{-1} +
      0.6366564*Tbar^{-2} +
      -0.241605*Tbar^{-3})
  #
  # Define function for mu_1:
  mu_1 <- function(Tbar, rbar){
    # Table 3 (Huber et al., 2009)
    H <- matrix(
      data=c(
        # i=0      i=1         i=2        i=3         i=4         i=5
         0.520094,  0.0850895, -1.08374,  -0.289555,   0,          0,          # j=0
         0.222531,  0.999115,   1.88797,   1.26613,    0,          0.120573,   # j=1
        -0.281378, -0.906851,  -0.772479, -0.489837,  -0.257040,   0,          # j=2
         0.161913,  0.257399,   0,         0,          0,          0,          # j=3
        -0.0325372, 0,          0,         0.0698452,  0,          0,          # j=4
         0,         0,          0,         0,          0.00872102, 0,          # j=5
         0,         0,          0,        -0.00435673, 0,         -0.000593264 # j=6
      ),
      nrow = 7,
      ncol = 6,
      byrow = TRUE
    )
    #
    my_mu <- 0
    for(i in seq(6)){
      coef1 <- (Tbar^{-1} - 1)^{i-1}
      coef2 <- 0
      for (j in seq(7)){
        coef2 <- coef2 + H[j,i]*(rbar-1)^{j-1}
      }
      my_mu <- my_mu + coef1*coef2
    }
    my_mu <- exp(rbar * my_mu)
    my_mu
  }
  #
  mu_bar <- mu_0*mu_1(Tbar, rbar)
  mu <- mu_bar*mu_ast    # Pa s
  mu <- (1e3)*mu         # mPa s
  return(mu)
}

# ************************************************************************
# * Name: vogel_visc
# *
# * Input: numeric (vector), ambient temperature, degrees C (tc)
# *
# * Return: numeric (vector), viscosity of water in mPa s
# *
# * Features: This function calculates the viscosity of water based on 
# *           the Vogel Equation.
# ************************************************************************
vogel_visc <- function(tc){
  0.024263*exp(578.919/((tc + 273.15) - 137.546))
}
