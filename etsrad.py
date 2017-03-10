#!/usr/bin/python
#
# etsrad.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-09-10 -- created
# 2015-11-13 -- last updated
#
# ------------
# description:
# ------------
# This script calculates extraterrestrial solar radiation.
# * half-hour time series and daily totals
# * total daylight hours
#
# ----------
# changelog:
# ----------
# 01. change lat/lon order in __init__() to lon, lat, day [13.09.13]
# 02. added srad_to_ppfd conversion factor [13.09.13]
# 03. added local_sec time series [13.09.16]
# 04. calculates daylight hours based on 1-sec time series [13.09.16]
# 05. added daylight radiation integrals & their calculation [13.09.16]
# 06. updated calc_hours() [13.09.16]
# --> daily integral of PPFD (units of umol m-2); instead of daylight average
# --> implemented Simpson's rule instead of Trapezoidal rule
# 07. added functions [14.07.16]
# --> earth_period()
# --> earth_velocity()
# --> correct_kepler()
# --> simplified_kepler()
# --> julian_day()
# --> equinox()
# --> get_lambda()
# 08. moved functions into Solar class [14.10.15]
# 09. added daylight savings factor to class input [14.10.17]
# 10. added datetime module [14.10.18]
# 11. added local_time & solar_time as SOLAR class variables [14.10.18]
# 12. added check for ds input [14.10.18]
# 13. updated variable names [14.11.27]
# 14. added figure plot [14.11.27]
# 15. PEP8 style fixes [15.11.13]
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import datetime
import matplotlib.pyplot as plt
import numpy
from sys import exit


###############################################################################
## CLASSES:
###############################################################################
class SOLAR:
    """
    Name:     SOLAR
    Features: This class calculates the half-hourly extraterrestrial PPFD
              [umol m-2 s-1], and the daily extraterrestrial PPFD [umol m-2]
              based on the STASH 2.0 methodology
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    ke = 0.0167    # eccentricity for 2000 CE (Berger, 1978)
    keps = 23.44   # obliquity for 2000 CE, degrees (Berger, 1978)
    kfFEC = 2.04   # from flux to energy conver., umol/J (Meek et al., 1984)
    kGsc = 1360.8  # Solar constant, W/m^2 (Kopp & Lean, 2011)
    komega = 283.  # longitude of perihelion for 2000 CE, deg (Berger, 1978)
    #
    # List of local time at half-hourly time step:
    local_hh = numpy.array([0.5*i for i in xrange(48)])

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, lon, lat, n, y, dsf=0):
        """
        Name:     SOLAR.__init__
        Input:    - float, longitude, degrees (lon)
                  - float, latitude, degrees (lat)
                  - int, day of year (n)
                  - int, year
                  - int, daylight savings (ds)
        """
        # Error handle and assign required public variables:
        if y == 0:
            self.user_year = 1950
        else:
            self.user_year = y
        if lat > 90.0 or lat < -90.0:
            print "Latitude outside range of validity (-90 to 90)!"
            exit(1)
        else:
            self.user_lat = lat
        if lon > 180.0 or lon < -180.0:
            print "Longitude outside range of validity (-180 to 180)!"
            exit(1)
        else:
            self.user_lon = lon
        if n < 1 or n > 366:
            print "Day outside range of validity (1 to 366)!"
            exit(1)
        else:
            self.user_day = n
        if dsf < 0 or dsf > 1:
            print "Set daylight savings time to 0 or 1"
            exit(1)
        else:
            self.dsf = dsf
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Create datetime series
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.local_time = numpy.array([
            datetime.datetime(self.user_year, 1, 1, 0, 0, 0) +
            datetime.timedelta(days=(n-1)) +
            datetime.timedelta(hours=i) for i in self.local_hh
        ])
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year (kN), days
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if y == 0:
            self.kN = 365
        else:
            self.kN = (
                self.julian_day((y + 1), 1, 1) -
                self.julian_day(y, 1, 1)
            )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lamb), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        my_nu, my_lambda = self.berger_tls(n)
        self.nu_deg = my_nu
        self.lambda_deg = my_lambda
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        dr = self.berger_dr(my_nu)
        self.dr = dr
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        delta = self.woolf_delta(my_lambda)
        self.delta_deg = delta
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate time zone hour, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if lon < 0:
            # Swap to positive to "round down" negative numbers:
            temp_lon = -1.0*lon
            temp_tzh = int(temp_lon/15)
            tz_hour = -1.0*temp_tzh
        else:
            tz_hour = int(lon/15)
        self.tz_hour = tz_hour
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the equation of time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eot = self.spencer_eot(n, self.kN)
        self.eot_hour = eot
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate the longitude correction, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        lc = (15.0*tz_hour - lon)/15.0
        self.lc_hour = lc
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8a. Calculate the solar time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ts_hh = self.local_hh + eot - lc - dsf
        self.ts_hh = ts_hh
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8b. Create solar datetime series
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.solar_time = numpy.array([
            datetime.datetime(self.user_year, 1, 1, 0, 0, 0) +
            datetime.timedelta(days=(n-1)) +
            datetime.timedelta(hours=i) for i in self.ts_hh
        ])
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate the hour angle, degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        w_hh = (360./24.)*(ts_hh - 12.0)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(delta)*self.dsin(lat)
        rv = self.dcos(delta)*self.dcos(lat)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate the sunset hour angle (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Note: ru/rv == tan(delta)*tan(lat)
        # Eq. 3.22, Stine & Geyer (2001)
        if (ru/rv) >= 1.0:
            hs = 180.0   # Polar day (no sunset)
        elif (ru/rv) <= -1.0:
            hs = 0.0     # Polar night (no sunrise)
        else:
            hs = (180.0/numpy.pi)*numpy.arccos(-1.0*ru/rv)
        self.hs_deg = hs
        self.hs_lct = (hs/15.0) - eot + lc + dsf + 12
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate daylight hours, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ds = (24.0/180.0)*hs
        self.ds_hour = ds
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate the solar radiation flux, W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        io_hh = self.kGsc*dr*(ru + rv*self.dcos(w_hh))
        io_hh = io_hh.clip(min=0)
        self.io_wm2 = io_hh
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate the half-hourly PPFD, umol/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ppfd_hh = self.kfFEC*io_hh
        self.par_umolm2 = ppfd_hh
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 15. Calculate the daily solar irradiation, J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ho = (86400.0/numpy.pi)*self.kGsc*dr*(
            ru*(numpy.pi/180.0)*hs + rv*self.dsin(hs)
        )
        self.ho_jm2 = ho
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 16. Calculate the daily top of atm. PPFD, mol/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        qo = (1.0e-6)*self.kfFEC*ho
        self.qo_molm2 = qo

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def berger_dr(self, lon):
        """
        Name:     SOLAR.berger_dr
        Input:    float, true anomaly, degrees (lon)
        Output:   float, distance factor, unitless
        Features: Returns distance factor for a given true anomaly
        Ref:      Berger, A. L., M. F. Loutre, and C. Tricot (1993), Insolation
                  and earth's orbital periods, J. Geophys. Res., 98, 10341-
                  10362.
        """
        my_rho = (1.0 - self.ke**2)/(1.0 + self.ke*self.dcos(lon))
        my_dr = (1.0/my_rho)**2
        return(my_dr)

    def berger_tls(self, n):
        """
        Name:     SOLAR.berger_tls
        Input:    int, day of year
        Output:   tuple,
                  - true anomaly, degrees
                  - true longitude, degrees
        Features: Returns true anomaly and true longitude for a given day
        Ref:      Berger, A. L. (1978), Long term variations of daily
                  insolation and quaternary climatic changes, J. Atmos. Sci.,
                  35, 2362-2367.
        """
        # Variable substitutes:
        xee = self.ke**2
        xec = self.ke**3
        xse = numpy.sqrt(1.0 - xee)
        #
        # Mean longitude for vernal equinox:
        xlam = (
            (self.ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(self.komega) -
            xee/4.0*(0.5 + xse)*self.dsin(2.0*self.komega) +
            xec/8.0*(1.0/3.0 + xse)*self.dsin(3.0*self.komega)
            )
        xlam = numpy.degrees(2.0*xlam)
        #
        # Mean longitude for day of year:
        dlamm = xlam + (n - 80.0)*(360.0/self.kN)
        #
        # Mean anomaly:
        anm = dlamm - self.komega
        ranm = numpy.radians(anm)
        #
        #  True anomaly:
        ranv = (ranm + (2.0*self.ke - xec/4.0)*numpy.sin(ranm) +
                5.0/4.0*xee*numpy.sin(2.0*ranm) +
                13.0/12.0*xec*numpy.sin(3.0*ranm))
        anv = numpy.degrees(ranv)
        #
        # True longitude:
        my_tls = anv + self.komega
        if my_tls < 0:
            my_tls += 360.0
        elif my_tls > 360:
            my_tls -= 360.0
        #
        # True anomaly:
        my_nu = (my_tls - self.komega)
        if my_nu < 0:
            my_nu += 360.0
        #
        return(my_nu, my_tls)

    def dcos(self, x):
        """
        Name:     SOLAR.dcos
        Input:    float, angle, degrees (x)
        Output:   float, cos(x*pi/180)
        Features: Calculates the cosine of an angle given in degrees
        """
        return numpy.cos(x*numpy.pi/180.0)

    def dsin(self, x):
        """
        Name:     SOLAR.dsin
        Input:    float, angle, degrees (x)
        Output:   float, sin(x*pi/180)
        Features: Calculates the sine of an angle given in degrees
        """
        return numpy.sin(x*numpy.pi/180.0)

    def julian_day(self, y, m, i):
        """
        Name:     SOLAR.julian_day
        Input:    - int, year (y)
                  - int, month (m)
                  - int, day of month (i)
        Output:   float, Julian Ephemeris Day
        Features: Converts Gregorian date (year, month, day) to Julian
                  Ephemeris Day
        Ref:      Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day," Astronomical
                  Algorithms
        """
        if m <= 2.0:
            y -= 1.0
            m += 12.0
        #
        a = int(y/100)
        b = 2 - a + int(a/4)
        #
        jde = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5
        return jde

    def spencer_eot(self, n, P):
        """
        Name:     SOLAR.spencer_eot
        Input:    - int, day of the year (n)
                  - int, number of days per year (P)
        Output:   float, equation of time, hours
        Features: Returns the equation of time
        Ref:      Spencer, J.W. (1971), Fourier series representation of the
                  position of the sun, Search, 2 (5), p. 172.
        """
        B = 2.0*numpy.pi*(n - 1.0)/P
        my_eot = 12.0/(numpy.pi)*(
            (7.5e-6) + (1.868e-3)*self.dcos(B) - (3.2077e-2)*self.dsin(B) -
            (1.4615e-2)*self.dcos(2.0*B) - (4.0849e-2)*self.dsin(2.0*B)
        )
        return(my_eot)

    def woolf_delta(self, lon):
        """
        Name:     SOLAR.woolf_delta
        Input:    float, true longitude, degrees (lon)
        Output:   float, declination angle, degrees
        Features: Returns the declination angle for a given true longitude
        Ref:      Woolf, H. M. (1968), On the computation of solar evaluation
                  angles and the determination of sunrise and sunset times,
                  Tech. rep. NASA-TM-X-164, National Aeronautics and Space
                  Administration, Washington, DC.
        """
        my_delta = numpy.arcsin(self.dsin(lon)*self.dsin(self.keps))
        my_delta *= (180.0/numpy.pi)
        return(my_delta)

###############################################################################
## MAIN:
###############################################################################
# User options
my_lat = 51.408
my_lon = -0.64
my_day = 83
my_year = 2001

my_solar = SOLAR(my_lon, my_lat, my_day, my_year, 0)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(my_solar.local_time, my_solar.io_wm2, 'r-', label="Io")
ax1.set_ylabel('$I_o$ (W m$^{-2}$)')
ax1.set_xlabel('Local Time (hours)')
plt.show()
