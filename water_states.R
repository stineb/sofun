# RStudio v 0.98
#
# water_states.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-01-30 -- created
# 2015-03-10 -- last updated
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
# *        ambient pressure, Pa (p)
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
tumlirz_density <- function(tc, p){
  # Calculate lamda, (bar c^3)/g:
  c1 <- 1788.316
  c2 <- 21.55053
  c3 <- -0.4695911
  c4 <- (3.096363e-3)
  c5 <- -(7.341182e-6)
  #
  my_lambda <- c1 +
    c2*tc +
    c3*tc*tc +
    c4*tc*tc*tc +
    c5*tc*tc*tc*tc
  #
  # Calculate Po, bar
  p1 <- 5918.499
  p2 <- 58.05267
  p3 <- -1.1253317
  p4 <- (6.6123869e-3)
  p5 <- -(1.4661625e-5)
  #
  po <- p1 +
    p2*tc +
    p3*tc*tc +
    p4*tc*tc*tc +
    p5*tc*tc*tc*tc
  #
  # Calculate Vinf, c^3/g
  v1 <- 0.6980547
  v2 <- -(7.435626e-4)
  v3 <- (3.704258e-5)
  v4 <- -(6.315724e-7)
  v5 <- (9.829576e-9)
  v6 <- -(1.197269e-10)
  v7 <- (1.005461e-12)
  v8 <- -(5.437898e-15)
  v9 <- (1.69946e-17)
  v10 <- -(2.295063e-20)
  # 
  vinf <- v1 +
    v2*tc +
    v3*tc*tc +
    v4*tc*tc*tc +
    v5*tc*tc*tc*tc +
    v6*tc*tc*tc*tc*tc +
    v7*tc*tc*tc*tc*tc*tc +
    v8*tc*tc*tc*tc*tc*tc*tc +
    v9*tc*tc*tc*tc*tc*tc*tc*tc +
    v10*tc*tc*tc*tc*tc*tc*tc*tc*tc
  #
  # Convert pressure to bars (1 bar = 100000 Pa)
  pbar <- (1e-5)*p
  #
  # Calculate the specific volume (cm^3 g^-1):
  v <- vinf + my_lambda/(po + pbar)
  #
  # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  rho <- (1e3/v)
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
chen_density <- function(tc, p){
  # Vo can be found by using do (1/Vo), g/cm^3, (Kell, 1975):
  d1 <- 0.99983952
  d2 <- (6.78826e-5)
  d3 <- -(9.08659e-6)
  d4 <- (1.02213e-7)
  d5 <- -(1.35439e-9)
  d6 <- (1.47115e-11)
  d7 <- -(1.11663e-13)
  d8 <- (5.04407e-16)
  d9 <- -(1.00659e-18)
  #
  do <- d1 +
    d2*tc +
    d3*tc*tc +
    d4*tc*tc*tc +
    d5*tc*tc*tc*tc +
    d6*tc*tc*tc*tc*tc +
    d7*tc*tc*tc*tc*tc*tc +
    d8*tc*tc*tc*tc*tc*tc*tc +
    d9*tc*tc*tc*tc*tc*tc*tc*tc
  #
  # Calculate Ko, atm
  k1 <- 19652.17
  k2 <- 148.183
  k3 <- -2.29995
  k4 <- 0.01281
  k5 <- -(4.91564e-5)
  k6 <- (1.03553e-7)
  # 
  ko <- k1 +
    k2*tc +
    k3*tc*tc +
    k4*tc*tc*tc +
    k5*tc*tc*tc*tc +
    k6*tc*tc*tc*tc*tc
  #
  # Calculate A
  a1 <- 3.26138
  a2 <- (5.223e-4)
  a3 <- (1.324e-4)
  a4 <- -(7.655e-7)
  a5 <- (8.584e-10)
  #
  a <- a1 +
    a2*tc +
    a3*tc*tc +
    a4*tc*tc*tc +
    a5*tc*tc*tc*tc 
  #
  # Calculat B
  b1 <- (7.2061e-5)
  b2 <- -(5.8948e-6)
  b3 <- (8.699e-8)
  b4 <- -(1.01e-9)
  b5 <- (4.322e-12)
  #
  b <- b1 +
    b2*tc +
    b3*tc*tc +
    b4*tc*tc*tc +
    b5*tc*tc*tc*tc 
  #
  # Convert pressure to bar (1 bar = 100000 Pa)
  pbar <- (1e-5)*p
  #
  # Calculate the density (i.e., 1/V), kg/m^3:
  rho <- (1e3)*do*(ko + a*pbar + b*pbar^2)/(ko + a*pbar + b*pbar^2 - pbar)
  #
  return(rho)
}

# ************************************************************************
# * Name: huber_visc
# *
# * Input: numeric (vector), ambient temperature, degrees C (tc)
# *        numeric, ambient pressure, Pa (P)
# *        character, the water density method (e.g., 'tumlirz' or 'chen')
# *
# * Return: numeric (vector), viscosity of water in Pa s (mu)
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
# * Depends: density function (either chen or tumlirz)
# *
# * Ref: Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. 
# *        Sengers, M. J. Assael, ..., K. Miyagawa (2009) New 
# *        international formulation for the viscosity of H2O, J. Phys. 
# *        Chem. Ref. Data, Vol. 38(2), pp. 101-125.
# *
# ************************************************************************
huber_visc <- function(tc, p, method){
  # Define reference temperature, density, and pressure values:
  tk_ast <- 647.096      # Kelvin
  rho_ast <- 322.0       # kg/m^3
  mu_ast <- (1e-6)       # Pa s
  #
  # Get the density of water:
  if (method == 'chen'){
    rho <- chen_density(tc, p)         # kg/m^3
  } else if (method == 'tumlirz') {
    rho <- tumlirz_density(tc, p)     # kg/m^3
  } else {
    rho <- 1e3                         # kg/m^3
  }
  #
  # Calculate dimensionless parameters:
  tbar <- (tc + 273.15)/tk_ast
  tbarx <- tbar^(0.5)
  tbar2 <- tbar^2
  tbar3 <- tbar^3
  rbar <- rho/rho_ast
  #
  # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
  mu0 <- 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
  mu0 <- 1e2*tbarx/mu0
  #
  # Create Table 3, Huber et al. (2009)
  h_array <- matrix(
    data=c(
      # i=0       i=1         i=2        i=3         i=4         i=5
      0.520094,    0.0850895, -1.08374,  -0.289555,   0,          0,          # j=0
      0.222531,    0.999115,   1.88797,   1.26613,    0,          0.120573,   # j=1
     -0.281378,   -0.906851,  -0.772479, -0.489837,  -0.257040,   0,          # j=2
      0.161913,    0.257399,   0,         0,          0,          0,          # j=3
     -0.0325372,   0,          0,         0.0698452,  0,          0,          # j=4
      0,           0,          0,         0,          0.00872102, 0,          # j=5
      0,           0,          0,        -0.00435673, 0,         -0.000593264 # j=6
    ),
    nrow = 7,
    ncol = 6,
    byrow = TRUE
  )
  #
  # Calculate mu1:
  mu1 <- 0
  ctbar <- (1.0/tbar) - 1.0
  for(i in seq(6)){
    coef1 <- ctbar^(i - 1)
    coef2 <- 0
    for (j in seq(7)){
      coef2 <- coef2 + h_array[j,i]*(rbar - 1)^(j - 1)
    }
    mu1 <- mu1 + coef1*coef2
  }
  mu1 <- exp(rbar*mu1)
  #
  # Calculate mu_bar (Eq. 2, Huber et al., 2009)
  # * assumes mu2 = 1.0
  mu_bar <- mu0*mu1
  #
  # Calculate mu (Eq. 1, Huber et al., 2009)
  mu <- mu_bar*mu_ast    # Pa s
  #
  return(mu)
}

# ************************************************************************
# * Name: vogel_visc
# *
# * Input: numeric (vector), ambient temperature, degrees C (tc)
# *
# * Return: numeric (vector), viscosity of water in Pa s
# *
# * Features: This function calculates the viscosity of water based on 
# *           the Vogel Equation.
# ************************************************************************
vogel_visc <- function(tc){
  (2.4263e-5)*exp(578.919/((tc + 273.15) - 137.546))
}
