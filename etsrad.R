# R-Studio
#
# etsrad.R
#
# written by Tyler W. Davis
# Imperial College London
#
# created: 2014-03-30
# updated: 2014-11-27
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script calculates half-hourly extraterrestrial solar radiation based on
# STASH 2.0 methodology
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 01. updated function documentation [14.11.26]
# 02. updated solar function and added dependencies [14.11.26]
# 03. separated global constants [14.11.26]
# 04. updated variable names [14.11.27]
# 05. added figure plot [14.11.27]
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////

# ************************************************************************
# * Name: etsrad
# *
# * Input: - double, longitude, degrees (user.lon)
# *        - double, latitude, degrees (user.lat)
# *        - double, day of the year (user.day)
# *        - double, year (user.year)
# *        - double, daylight savings factor (user.ds)
# *
# * Return: list
# *
# * Features: Returns a list object with extraterrstrial solar radiation
# *           variables.
# *
# * Depends: - berger_tls() 
# *          - berger_dr()
# *          - dcos()
# *          - dsin()
# *          - julian_day()
# *          - woolf_delta()
# *          - spencer_eot()
# *
# ************************************************************************
solar <- function(lon, lat, n, y, dsf){
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  et.srad <- list()
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # Half-hourly time series:
  local_hh <- (seq(48) - 1.0)*0.5
  et.srad$local_time <- local_hh
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (lon > 180 || lon < -180){
    stop("Warning: Longitude outside range of validity (-180 to 180)!")
  } else {
    et.srad$user_lon <- lon
  }
  if (lat > 90 || lat < -90){
    stop("Warning: Latitude outside range of validity (-90 to 90)!")
  } else {
    et.srad$user_lat <- lat
  }
  if (n < 1 || n > 366){
    stop("Warning: Day outside range of validity (1 to 366)!")
  } else {
    et.srad$user_day <- n
  }
  if (dsf < 0 || dsf > 1){
    stop("Warning: Daylight savings factor outside range of validity (0 or 1)")
  } else {
    et.srad$dsf <- dsf
  }
  et.srad$user_year <- y
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate the number of days in the year
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (y == 0){
    kN = 365
  } else {
    kN = (julian_day(y+1, 1, 1) - julian_day(y, 1, 1))
  }
  et.srad$kN <- kN
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Heliocentric longitudes: nu and lambda (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  my_lons = berger_tls(n, kN)
  my_nu = my_lons[1]
  my_lambda = my_lons[2]
  #
  et.srad$lambda_deg <- my_lambda
  et.srad$nu_deg <- my_nu
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate distance factor (unitless)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dr <- berger_dr(my_nu)
  et.srad$dr <- dr
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Calculate the declination angle (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  delta <- woolf_delta(my_lambda)
  et.srad$delta_deg <- delta
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Calculate time zone hour (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Correction for local time zone; based on one hour per 15 degrees longitude
  # * positive or negative deviation from UTC/GMT
  if (lon < 0){
    temp_lon <- -1.0*lon
    temp_tzh <- floor(temp_lon/15)
    tz_hour <- -1.0*temp_tzh
  } else {
    tz_hour <- floor(lon/15)
  }
  et.srad$tz_hour <- tz_hour
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 06. Calculate the equation of time (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eot <- spencer_eot(n, kN)
  et.srad$eot_hour <- eot
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 07. Longitudinal correction factor (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Correction factor for longitude
  # ref: Eq. 3.6, Stine and Geyer (2001)
  lc <- (15*tz_hour - lon)/15
  et.srad$lc_hour <- lc
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 08. Calculate solar time (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate the solar time based on the local time
  # ref: Eq. 3.5, Stine and Geyer (2001)
  ts_hh <- local_hh + eot - lc - dsf
  et.srad$ts_hh <- ts_hh
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 09. Calculate the hour angle (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Solar zenith angle
  # ref: Eq. 3.1, Stine and Geyer (2001)
  # * solar noon is at 0 degrees (i.e., vertical line)
  # * based once again on 15 deg movement per hour
  w_hh <- (360/24)*(ts_hh - 12)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 10. Calculate the variable substitutes
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ru <- dsin(delta)*dsin(lat)
  rv <- dcos(delta)*dcos(lat)
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 11. Calculate the sunset angle (degrees)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 3.22, Stein & Geyer (2001)
  # - represents the number of hours from solar noon when the sun sets (i.e., 
  #   the Ra curve dips below zero)
  # - when using hs to integrate over daylight hours, multiply by 2
  if (ru/rv >= 1.0){
    hs <- 180  # Polar day (no sunset)
  } else if (ru/rv <= -1.0){ 
    hs <- 0 # Polar night (no sunrise)
  } else {
    hs <- (180/pi)*acos(-1.0*ru/rv)
  }
  et.srad$hs_deg <- hs
  et.srad$hs_lct <- (hs/15.0) - eot + lc + dsf + 12
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 12. Calculate daylight hours (hours)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 34, Allen et al. (1998)
  ds <- (24.0/180)*hs
  et.srad$ds_hour <- ds
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 13. Calculate solar radiation flux (W/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extraterrestrial solar radiation w.r.t. a flat horizon
  # ref: Eq. 1.10.2 in Duffie and Beckman (1991)
  io_hh <- kGsc*dr*(ru + rv*dcos(w_hh))
  io_hh[io_hh < 0 ] <- 0
  et.srad$io_wm2 <- io_hh
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 14. Calculate half-hourly PPFD (umol/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ppfd_hh = kfFEC*io_hh
  et.srad$par_umolm2 <- ppfd_hh
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 15. Calculate daily solar irradiation (J/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 1.10.3, Duffy & Beckman (1991)
  ho <- (86400.0/pi)*kGsc*dr*(ru*pi*hs/180 + rv*dsin(hs))
  et.srad$ho_jm2 <- ho
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 16. Calculate daily top of the atmosphere PPFD (mol/m^2)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  qo <- (1e-6)*kfFEC*ho
  et.srad$qo_molm2 <- qo
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  et.srad
}

# ************************************************************************
# * Name: berger_dr
# *
# * Input: double, true anomaly, degrees (lon)
# *
# * Return: double, distance factor, unitless
# *
# * Features: Returns the distance factor based on heliocentric 
# *           longitude and eccentricity:
# *           dr = (1+e cos(nu))^2/(1-e^2)^2
# *
# * Depends: - dcos()
# *          - ke
# *
# * Ref: Berger, A. L., M. F. Loutre, and C. Tricot (1993), Insolation
# *        and earth's orbital periods, J. Geophys. Res., 98, 10341-
# *        10362.
# *
# ************************************************************************
berger_dr <- function(lon){
  rho <- (1.0 - ke^2)/(1.0 + ke*dcos(lon))
  return (1/rho^2)
}

# ************************************************************************
# * Name: berger_tls
# *
# * Input: - double, day of the year (n)
# *        - double, days in the year (P)
# *
# * Return: - double, true anomaly, degrees (nu)
# *         - double, true longitude, degrees (tls)
# *
# * Features: Returns the true anomaly (nu) and true longitude (lambda) 
# *           based on the methodology in BERGER78 fortran code 
# *           (version 2014)
# *
# * Depends: - ke ...... eccentricity, unitless
# *          - komega .. longitude of perihelion, degrees
# *
# * Ref: Berger (1978) Long-term variations of daily insolation and 
# *      quarternary climatic changes, Journal of Atmospheric Sciences, 
# *      vol. 35, pp. 2362--2367.
# *
# ************************************************************************
berger_tls <- function(n, P=365){
  # Variable substitutes:
  pir <- pi/180.0
  xse <- sqrt(1 - ke^2)
  xee <- ke^2
  xec <- ke^3
  #
  # Mean longitude for vernal equinox:
  xlam <- (ke/2.0 + xec/8.0)*(1 + xse)*sin(komega*pir) - 
    xee/4.0*(0.5 + xse)*sin(2.0*komega*pir) + 
    xec/8.0*(1.0/3.0 + xse)*sin(3.0*komega*pir)
  xlam <- 2.0*xlam/pir
  #
  # Mean longitude for day of year:
  dlamm <- xlam + (n - 80.0)*(360.0/P)
  #
  # Mean anomaly for day of year:
  anm <- dlamm - komega
  ranm <- anm*pir
  #
  # True anomaly for day of year:
  ranv <- ranm + (2.0*ke - xec/4.0)*sin(ranm) +
    5.0/4.0*xee*sin(2.0*ranm) + 
    13.0/12.0*xec*sin(3.0*ranm)
  anv <- ranv/pir
  #
  # True longitude for day of year:
  tls <- anv + komega
  if (tls < 0){
    tls <- tls + 360
  } else if (tls > 360) {
    tls <- tls - 360
  }
  #
  # True anomaly:
  nu <- tls - komega
  if (nu < 0){
    nu <- nu + 360.
  } 
  #
  return (c(nu, tls))
}

# ************************************************************************
# * Name: dcos
# *
# * Input: double, angle, degrees (d)
# *
# * Return: double, cosine of angle, unitless
# *
# * Features: This function calculates the cosine of an angle (d) given
# *           in degrees.
# *
# ************************************************************************
dcos <- function(d) {
  cos(d*pi/180)
}

# ************************************************************************
# * Name: dsin
# *
# * Input: double, angle, degrees (d)
# *
# * Return: double, sine of angle, unitless
# *
# * Features: This function calculates the sine of an angle (d) given
# *           in degrees.
# *
# ************************************************************************
dsin <- function(d) {
  sin(d*pi/180)
}

# ************************************************************************
# * Name: julian_day
# *
# * Input: int (year) 
# *        int (month) 
# *        int (day)
# *
# * Return: int, Julian day
# *
# * Features: This function converts a date in the Gregorian calendar
# *           to a Julian day number (i.e., a method of consecutative 
# *           numbering of days---does not have anything to do with 
# *           the Julian calendar!)
# *
# * Ref: Eq. 7.1 J. Meeus (1991), Chapter 7 "Julian Day", 
# *      Astronomical Algorithms
# *
# ************************************************************************
julian_day <- function(Y, M, D){
  if(M<=2){
    Y <- Y-1
    M <- M+12
  }
  A <- floor(Y/100)
  B <- 2 - A + floor(A/4)  # modified leap year def. for Gegorian calendar
  # for Julian calendar, B = 0
  # Eq. 7.1, J. Meeus (1991), Astronomical Algorithms
  JDE <- 1.0*floor(365.25*(Y+4716)) + 
    1.0*floor(30.6001*(M+1)) + 1.0*D + 1.0*B - 1524.5
  JDE
}

# ************************************************************************
# * Name: spencer_eot
# *
# * Input: double, day of the year
# *
# * Return: double, equation of time, hours
# *
# * Features: This function calculates the equation of time.
# *
# * Ref: Spencer, J.W. (1971), Fourier series representation of the 
# *      position of the sun, Search, 2 (5), p. 172.
# *
# * Note: Equation represented is the corrected version by Oglesby, M. 
# *       (1998), Fourier Paper. Available: 
# *       https://www.mail-archive.com/sundial@uni-koeln.de/msg01050.html
# *
# ************************************************************************
spencer_eot <- function(n, P=365){
  B <- 2*pi*(n - 1)/P
  eot <- 12/pi*(
    (7.5e-6) + (1.868e-3)*cos(B) - (3.2077e-2)*sin(B) - 
      (1.4615e-2)*cos(2*B) - (4.0849e-2)*sin(2*B)
  )
  return(eot)
}

# ************************************************************************
# * Name: woolf_delta
# *
# * Input: double, true longitude, degrees (lambda)
# *
# * Return: double, declination angle, degrees
# *
# * Features: Returns the declination angle based on heliocentric 
# *           longitude and obliquity:
# *           delta = asin(sin(lambda)*sin(epsilon))
# *
# * Depends: - dsin()
# *          - keps
# *
# * Ref: Woolf, H. M. (1968), On the computation of solar evaluation 
# *        angles and the determination of sunrise and sunset times, 
# *        Tech. rep. NASA-TM-X-164, National Aeronautics and Space 
# *        Administration, Washington, DC.
# *
# ************************************************************************
woolf_delta <- function(lon){
  return (asin(dsin(lon)*dsin(keps))*180.0/pi)
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define constants #########################################################
# /////////////////////////////////////////////////////////////////////////////
ke <- 0.0167   # eccentricity for 2000 CE (Berger, 1978)
keps <- 23.44  # obliquity for 2000 CE, degrees (Berger, 1978)
kfFEC <- 2.04  # from flux to energy conversion, umol/J (Meek et al., 1984)
kGsc <- 1360.8 # solar constant, W/m^2 (Kopp & Lean, 2011)
komega <- 283  # longitude of perihelion for 2000 CE, degrees (Berger, 1978)

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Main program #############################################################
# /////////////////////////////////////////////////////////////////////////////
my_lat <- 51.408
my_lon <- -0.64
my_day <- 83
my_year <- 2001

my_solar <- solar(my_lon, my_lat, my_day, my_year, 0)

par(mar=c(4.5, 4.5, 1, 1))
plot(my_solar$local_time, my_solar$io_wm2, type='l', col='red', lwd=2,
     xlab='Local Time (hours)', ylab=expression(italic(I[o])~(W%.%m^{-2})))

