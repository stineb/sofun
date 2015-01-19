  module waterbal_stash
  !////////////////////////////////////////////////////////////////
  ! STASH WATERBALANCE MODULE
  ! Contains the "main" subroutine 'waterbal' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'waterbal' must contain this list 
  ! of subroutines (names that way).
  !   - getpar_modl_waterbal
  !   - initio_waterbal
  !   - initoutput_waterbal
  !   - getout_daily_waterbal
  !   - getout_monthly_waterbal
  !   - writeout_ascii_waterbal
  !   - waterbal
  ! Required module-independent model state variables (necessarily 
  ! updated by 'waterbal') are:
  !   - daytime net radiation ('rn_d')
  !   - soil moisture ('dwn')
  !   - runoff ('dro')
  !----------------------------------------------------------------
  use params_core, only: ndayyear, nmonth, nlu, maxgrid

  implicit none

  ! STASH PARAMETERS
  ! read from file, therefore not declared as parameter
  real :: kA                ! constant for Rnl (Monteith & Unsworth, 1990)
  real :: kalb_sw           ! shortwave albedo (Federer, 1968)
  real :: kalb_vis          ! visible light albedo (Sellers, 1985)
  real :: kb                ! constant for Rnl (Linacre, 1968)
  real :: kc                ! cloudy transmittivity (Linacre, 1968)
  real :: kCw               ! supply constant, mm/hr (Federer, 1982)
  real :: kd                ! angular coefficient of transmittivity (Linacre, 1968)
  real :: ke                ! eccentricity for 2000 CE (Berger, 1978)
  real :: keps              ! obliquity for 2000 CE, degrees (Berger, 1978)
  real :: kfFEC             ! from flux to energy conversion, umol/J (Meek et al., 1984)
  real :: kG                ! gravitational acceleration, m/s^2 (Allen, 1973)
  real :: kGsc              ! solar constant, W/m^2 (Kopp & Lean, 2011)
  real :: kL                ! temperature lapse rate, K/m (Cavcar, 2000)
  real :: kMa               ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  real :: kMv               ! molecular weight of water vapor, kg/mol (Tsilingiris, 2008)
  real :: kPo               ! standard atmosphere, Pa (Allen, 1973)
  real :: kR                ! universal gas constant, J/mol/K (Allen, 1973)
  real :: kTo               ! base temperature, K (Prentice, unpublished)
  real :: kWm               ! soil moisture capacity, mm (Cramer & Prentice, 1988)
  real :: kw                ! entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
  real :: komega            ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)

  ! MODULE-SPECIFIC VARIABLES
  ! These are not required outside module STASH
  ! adjust these and make global if required by GePiSaT
  real :: ra_d                   ! daily extraterrestrial solar radiation (J/m^2)
  real :: ppfd_d                 ! daily PPFD (mol/m^2)
  real :: eet_d                  ! daily EET (mm)
  real :: pet_d                  ! daily PET (mm)
  real, dimension(nlu) :: aet_d  ! daily AET (mm)
  real, dimension(nlu) :: wc     ! daily condensation (mm)
  real, dimension(nlu) :: sw     ! evaporative supply rate (mm/h)

  ! DAILY OUTPUT VARIABLES
  ! xxx add jpngr dim
  real, dimension(nlu,ndayyear,maxgrid) :: outdwn             ! daily soil moisture, mm
  real, dimension(ndayyear,maxgrid)     :: outdho             ! daily solar irradiation, J/m2
  real, dimension(ndayyear,maxgrid)     :: outdhn             ! daily net radiation, J/m2
  real, dimension(ndayyear,maxgrid)     :: outdqn             ! daily PPFD, mol/m2
  real, dimension(nlu,ndayyear,maxgrid) :: outdcn             ! daily condensation water, mm
  real, dimension(nlu,ndayyear,maxgrid) :: outdro             ! daily runoff, mm
  real, dimension(ndayyear,maxgrid)     :: outdeq_n           ! daily equilibrium ET, mm
  real, dimension(ndayyear,maxgrid)     :: outdep_n           ! daily potential ET, mm
  real, dimension(nlu,ndayyear,maxgrid) :: outdea_n           ! daily actual ET, mm



  ! MONTHLY OUTPUT VARIABLES
  ! xxx add jpngr dim
  real, dimension(nmonth,maxgrid) :: outmeq_m
  real, dimension(nmonth,maxgrid) :: outmep_m
  real, dimension(nlu,nmonth,maxgrid) :: outmea_m
  real, dimension(nlu,nmonth,maxgrid) :: outmcpa
  real, dimension(nlu,nmonth,maxgrid) :: outmcwd
  real, dimension(nmonth,maxgrid) :: outmqm


contains

  subroutine waterbal( doy, lon, lat, dprec, dtemp, dsun, elv )
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates daily and monthly quantities for one year
    !-------------------------------------------------------------------------
    use params_core, only: ndayyear, ndaymonth, nlu
    use fluxes, only: dwn, dro

    ! ARGUMENTS
    integer, intent(in) :: doy                       ! day of year
    real, intent(in) :: lon                          ! longitude (degrees)
    real, intent(in) :: lat                          ! latitude (degrees)
    real, intent(in) :: dprec       ! daily precip (mm) 
    real, intent(in) :: dtemp        ! mean monthly temperature (deg C)
    real, intent(in) :: dsun        ! mean monthly sunshine fraction (unitless)
    !real, dimension(nmonth), intent(in) :: ppt       ! monthly precip (mm) 
    !real, dimension(nmonth), intent(in) :: tc        ! mean monthly temperature (deg C)
    !real, dimension(nmonth), intent(in) :: sf        ! mean monthly sunshine fraction (unitless)
    real, intent(in) :: elv                          ! altitude (m)

    ! LOCAL VARIABLES
    integer :: lu                        ! land unit (gridcell tile)
    integer :: moy                       ! month of year
    integer :: idx                       ! day of year corresponding to yesterday
    integer :: dm                        ! day of month

    ! Reset monthly totals
    !call initmonthly

    do lu=1,nlu

      ! Calculate evaporative supply rate, mm/h
      sw(lu) = kCw * dwn(lu) / kWm

      ! Calculate radiation and evaporation quantities
      !print*,'calling evap for doy', doy, '...'
      call evap( lu, lon, lat, doy, elv=elv, sf=dsun, tc=dtemp )
      !print*,'... done'

      ! Update soil moisture
      dwn(lu) = dwn(lu) + dprec + wc(lu) - aet_d(lu)

      ! Bucket model for runoff generation
      if (dwn(lu)>kWm) then
        ! Bucket is full 
        ! * set soil moisture to capacity
        ! * add remaining water to monthly runoff total
        dro(lu) = dwn(lu)
        dro(lu) = dro(lu) - kWm
        dwn(lu) = kWm

      elseif (dwn(lu)<0.0) then
        ! Bucket is empty
        ! * set soil moisture to zero
        dwn(lu) = 0.0
        dro(lu) = 0.0

      else
        ! No runoff occurrs
        dro(lu) = 0.0

      end if

    end do

    return

  end subroutine waterbal


  subroutine evap( lu, lon, lat, doy, elv, sf, tc )
    !/////////////////////////////////////////////////////////////////////////
    ! This subroutine calculates daily radiation and evapotranspiration
    ! quantities
    ! - daily PPFD (ppfd_d), mol/m^2
    ! - daily EET (eet_d), mm
    ! - daily PET (pet_d), mm
    ! - daily AET (aet_d), mm
    ! - daily condensation (wc), mm
    !-------------------------------------------------------------------------  
    use params_core, only: ndayyear, pi
    use fluxes, only: rn_d

    implicit none

    ! ARGUMENTS
    integer, intent(in) :: lu            ! land unit
    real, intent(in) :: lon              ! longitude, degrees
    real, intent(in) :: lat              ! latitude, degrees
    integer, intent(in) :: doy           ! day of the year (formerly 'n')
    real,    intent(in) :: elv ! elevation, metres
    real,    intent(in) :: sf  ! fraction of sunshine hours 
    real,    intent(in) :: tc  ! mean daily air temperature, C


    ! LOCAL VARIABLES
    real, dimension(2) :: out_berger     ! temporary var to store output of berger_tls
    real :: my_nu                        ! heliocentric longitudes
    real :: my_lambda                    ! heliocentric longitudes
    real :: my_rho
    real :: dr                           ! distance factor
    real :: delta                        ! declination angle 
    real :: ru                           ! variable substitute for u
    real :: rv                           ! variable substitute for v
    real :: hs                           ! sunset hour angle
    real :: tau_o, tau                   ! transmittivity (unitless)
    real :: rnl                          ! net longwave radiation (W/m^2)
    real :: rnn_d                        ! nighttime net radiation (J/m^2)
    real :: rw                           ! variable substitute (W/m^2)
    real :: hn                           ! net radiation cross-over hour angle
    real :: s                            ! slope of saturation vap press temp curve, Pa/K
    real :: pw                           ! density of water, kg/m^3
    real :: lv                           ! enthalpy of vaporization, J/kg
    real :: g                            ! psychrometric constant, Pa/K
    real :: econ                         ! Eq. 58, STASH 2.0 Documentation
    real :: rx                           ! variable substitute (mm/hr)/(W/m^2)
    real :: hi, cos_hi                   ! intersection hour angle, degrees


    ! Set default values for optional arguments
    ! xxx really keep them as optional inputs?

    ! Error handle and assign required public variables:
    if (lon>180.0 .or. lon<(-180.0)) then
      print*,"Longitude outside range of validity (-180 to 180)!"
      stop
    end if

    if (lat>90.0 .or. lat<(-90.0)) then
      print*,"Latitude outside range of validity (-90 to 90)!"
      stop
    end if

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 1. Calculate number of days in year (kN), days
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! (is provided within SOFUN)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 2. Calculate heliocentric longitudes (nu and lambda), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Berger (1978)
    out_berger = berger_tls( doy )
    
    ! xxx debug
    !print*,'doy, out_berger', doy, out_berger

    ! xxx my_nu is not perfectly identical to Python version
    my_nu      = out_berger(1)
    my_lambda  = out_berger(2)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3. Calculate distance factor (dr), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Berger et al. (1993)
    my_rho = (1.0 - ke**2)/(1.0 + ke * dgcos( my_nu ))
    dr = (1.0/my_rho)**2

    ! xxx debug
    !print*,'dr', dr

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 4. Calculate declination angle (delta), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Woolf (1968)
    ! delta = dasin( dgsin( my_lambda ) * dgsin( keps ) )   ! xxx arcsin is asin in Fortran?
    delta = asin( dgsin( my_lambda ) * dgsin( keps ) )   ! xxx use asin with single-precision compilation
    delta = degrees( delta )

    ! xxx debug
    !print*,'delta', delta

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 5. Calculate variable substitutes (u and v), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ru = dgsin(delta) * dgsin(lat)
    rv = dgcos(delta) * dgcos(lat)

    ! xxx debug
    !print*,'dgsin(lat)',dgsin(lat)
    !stop
    !print*,'ru,rv',ru,rv

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 6. Calculate the sunset hour angle (hs), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Note: u/v == tan(delta)*tan(lat)
    ! Eq. 3.22, Stine & Geyer (2001)
    if ((ru/rv) >= 1.0) then
      ! Polar day (no sunset)
      hs = 180.0 
    elseif ((ru/rv) <= -1.0) then
      ! Polar night (no sunrise)
      hs = 0.0
    else
      !hs = dacos(-1.0*ru/rv)
      hs = acos(-1.0*ru/rv)  ! xxx use acos with single-precision compilation
      hs = degrees(hs)
    end if

    ! xxx debug
    !print*,'hs',hs

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 7. Calculate daily extraterrestrial solar radiation (ra_d), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 1.10.3, Duffy & Beckman (1993)
    ra_d = (86400.0/pi)*kGsc*dr*(radians(ru)*hs + rv * dgsin(hs))

    ! xxx debug
    !print*,'pi, kGsc, dr, radians(ru), hs, rv dgsin(hs)', pi, kGsc, dr, radians(ru), hs, rv, dgsin(hs)
    !print*,'doy, ra_d', doy, ra_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 8. Calculate transmittivity (tau), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 11, Linacre (1968)
    tau_o = (kc + kd*sf)

    ! Eq. 2, Allen (1996)
    tau = tau_o*(1.0 + (2.67e-5)*elv)

    ! xxx debug
    !print*,'tau',tau

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 9. Calculate daily PPFD (ppfd_d), mol/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 57, STASH 2.0 Documentation
    ppfd_d = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*ra_d

    ! xxx debug
    !print*,'ppfd_d',ppfd_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 10. Estimate net longwave radiation (rnl), W/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
    rnl = (kb + (1.0 - kb)*sf)*(kA - tc)

    ! xxx debug
    !print*,'rnl',rnl

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 11. Calculate variable substitute (rw), W/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rw = (1.0-kalb_sw)*tau*kGsc*dr

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 12. Calculate net radiation cross-over hour angle (hn), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ((rnl - rw*ru)/(rw*rv) >= 1.0) then
      ! Net radiation negative all day
      hn = 0.0
    else if ((rnl - rw*ru)/(rw*rv) <= -1.0) then
      ! Net radiation positive all day
      hn = 180.0
    else
      !hn = degrees( dacos((rnl - rw*ru)/(rw*rv)) )
      hn = degrees( acos((rnl - rw*ru)/(rw*rv)) )   ! use acos with single precision compilation
    end if

    ! xxx debug
    !print*,'hn',hn

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 13. Calculate daytime net radiation (rn_d), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 53, STASH 2.0 Documentation
    rn_d = (86400.0/pi) * (hn*(pi/180.0)*(rw*ru - rnl) + rw*rv*dgsin(hn))

    ! xxx debug: slight deviation from python code (-6.156683e-05 %)
    !print*,'rn_d',rn_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 14. Calculate nighttime net radiation (rnn_d), J/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 56, STASH 2.0 Documentation
    rnn_d = (86400.0/pi)*(radians(rw*ru*(hs-hn)) + rw*rv*(dgsin(hs)-dgsin(hn)) + rnl*(pi - 2.0*radians(hs) + radians(hn)))

    ! xxx debug
    !print*,'rnn_d',rnn_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 15. Calculate water-to-energy conversion (econ), m^3/J
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Slope of saturation vap press temp curve, Pa/K
    s = sat_slope(tc)
    ! Enthalpy of vaporization, J/kg
    lv = enthalpy_vap(tc)
    ! Density of water, kg/m^3
    pw = density_h2o(tc, elv2pres(elv))
    ! Psychrometric constant, Pa/K
    g = psychro(tc, elv2pres(elv))

    ! Eq. 58, STASH 2.0 Documentation
    econ = s/(lv*pw*(s + g))

    ! xxx debug
    !print*,'econ',econ

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 16. Calculate daily condensation (wc), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 68, STASH 2.0 Documentation
    wc(lu) = 1000.0 * econ * abs(rnn_d)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 17. Estimate daily EET (eet_d), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 70, STASH 2.0 Documentation
    eet_d = 1000.0*econ*(rn_d)

    ! xxx debug
    !print*,'eet_d',eet_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 18. Estimate daily PET (pet_d), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 72, STASH 2.0 Documentation
    pet_d = (1.0+kw)*eet_d

    ! xxx debug
    !print*,'pet_d',pet_d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rx = 1000.0*3600.0*(1.0+kw)*econ

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 20. Calculate the intersection hour angle (hi), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cos_hi = sw(lu)/(rw*rv*rx) + rnl/(rw*rv) - ru/rv
    if (cos_hi >= 1.0) then
      ! Supply exceeds demand:
      hi = 0.0
    elseif (cos_hi <= -1.0) then
      ! Supply limits demand everywhere:
      hi = 180.0
    else
      hi = degrees(acos(cos_hi))
    end if

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 21. Estimate daily AET (aet_d), mm
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 81, STASH 2.0 Documentation
    aet_d(lu) = (24.0/pi)*(radians(sw(lu)*hi) + rx*rw*rv*(dgsin(hn) - dgsin(hi)) + radians((rx*rw*ru - rx*rnl)*(hn - hi)))

    ! xxx debug
    !print*,'aet_d',aet_d

    return

    !-------------------------------------------------------------   
    ! Refs: Allen, R.G. (1996), Assessing integrity of weather data for 
    !         reference evapotranspiration estimation, Journal of Irrigation
    !         and Drainage Engineering, vol. 122, pp. 97--106.
    !       Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998), 
    !         'Meteorological data,' Crop evapotranspiration - Guidelines for 
    !         computing crop water requirements - FAO Irrigation and drainage 
    !         paper 56, Food and Agriculture Organization of the United 
    !         Nations, online: http://www.fao.org/docrep/x0490e/x0490e07.htm
    !       Berger, A.L. (1978), Long-term variations of daily insolation and 
    !         quarternary climatic changes, Journal of Atmospheric Sciences, 
    !         vol. 35, pp. 2362--2367.
    !       Berger, A.L., M.F. Loutre, and C. Tricot (1993), Insolation and 
    !         Earth's orbital periods, J. Geophys. Res., 98, 10341--10362.
    !       Duffie, J. A. and W. A. Beckman (1991). Solar engineering of 
    !         thermal processes. 4th ed. New Jersey: John Wiley and Sons
    !       Federer (1982), Transpirational supply and demand: plant, soil, 
    !         and atmospheric effects evaluated by simulation, Water 
    !         Resources Research, vol. 18, no. 2, pp. 355--362.
    !       Ge, S., R.G. Smith, C.P. Jacovides, M.G. Kramer, R.I. Carruthers 
    !         (2011), Dynamics of photosynthetic photon flux density (PPFD) 
    !         and estimates in coastal northern California, Theoretical and 
    !         Applied Climatology, vol. 105, pp. 107--118.
    !       Henderson-Sellers, B. (1984), A new formula for latent heat of 
    !         vaporization of water as a function of temperature, Quarterly 
    !         Journal of the Royal Meteorological Society 110, pp. 1186–1190
    !       Linacre (1968), Estimating the net-radiation flux, Agricultural 
    !         Meteorology, vol. 5, pp. 49--63.
    !       Prentice, I.C., M.T. Sykes, W. Cramer (1993), A simulation model 
    !         for the transient effects of climate change on forest 
    !         landscapes, Ecological Modelling, vol. 65, pp. 51--70.
    !       Priestley, C.H.B. and R.J. Taylor (1972), On the assessment of 
    !         surface heat flux and evaporation using large-scale parameters, 
    !         Monthly Weather Review, vol. 100 (2), pp. 81--92.
    !       Spencer, J. W. (1971), Fourier series representation of the 
    !         position of the sun, Search, vol. 2, p. 172.
    !       Stine, W. B. and M. Geyer (2001). “Power from the Sun”. 
    !         online: http://www.powerfromthesun.net/Book/chapter03/chapter03
    !       Wetherald, R.T., S. Manabe (1972), Response to joint ocean-
    !         atmosphere model to the seasonal variation of the solar 
    !         radiation, Monthly Weather Review, vol. 100 (1), pp. 42--59.
    !       Woolf, H. M. (1968). On the computation of solar evaluation 
    !         angles and the determination of sunrise and sunset times. 
    !         Tech. rep. NASA-TM-X-164. National Aeronautics and Space 
    !         Administration (NASA).
    !-------------------------------------------------------------   

  end subroutine evap


  subroutine getpar_modl_waterbal( paramfilnam )
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads waterbalance module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    implicit none

    ! ARGUMENTS
    character(len=*) paramfilnam

    ! LOCAL VARIABLES
    integer, parameter :: npar = 21
    real, dimension(npar) :: params_array


    open(unit=01,file=trim(paramfilnam),status='OLD')      
    read (01,*) params_array
    close (01)

    kA = params_array(1)         ! constant for Rnl (Monteith & Unsworth, 1990)
    kalb_sw = params_array(2)    ! shortwave albedo (Federer, 1968)
    kalb_vis = params_array(3)   ! visible light albedo (Sellers, 1985)
    kb = params_array(4)         ! constant for Rnl (Linacre, 1968)
    kc = params_array(5)         ! cloudy transmittivity (Linacre, 1968)
    kCw = params_array(6)        ! supply constant, mm/hr (Federer, 1982)
    kd = params_array(7)         ! angular coefficient of transmittivity (Linacre, 1968)
    ke = params_array(8)         ! eccentricity for 2000 CE (Berger, 1978)
    keps = params_array(9)       ! obliquity for 2000 CE, degrees (Berger, 1978)
    kfFEC = params_array(10)     ! from flux to energy conversion, umol/J (Meek et al., 1984)
    kG = params_array(11)        ! gravitational acceleration, m/s^2 (Allen, 1973)
    kGsc = params_array(12)      ! solar constant, W/m^2 (Kopp & Lean, 2011)
    kL = params_array(13)        ! temperature lapse rate, K/m (Cavcar, 2000)
    kMa = params_array(14)       ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)
    kMv = params_array(15)       ! molecular weight of water vapor, kg/mol (Tsilingiris, 2008)
    kPo = params_array(16)       ! standard atmosphere, Pa (Allen, 1973)
    kR = params_array(17)        ! universal gas constant, J/mol/K (Allen, 1973)
    kTo = params_array(18)       ! base temperature, K (Prentice, unpublished)
    kWm = params_array(19)       ! soil moisture capacity, mm (Cramer & Prentice, 1988)
    kw = params_array(20)        ! entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
    komega = params_array(21)    ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)

    return

  end subroutine getpar_modl_waterbal


  subroutine initio_waterbal( prefix )
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    implicit none

    ! ARGUMENTS
    character(len=*) :: prefix

    ! LOCAL VARIABLES
    character(len=256) :: filnam


    !----------------------------------------------------------------
    ! DAILY OUTPUT
    !----------------------------------------------------------------

    ! HO: daily solar irradiation, J/m2
    filnam=trim(prefix)//'.d.ho.out'
    open(151,file=filnam,err=888,status='unknown')

    ! HN: daily net radiation, J/m2
    filnam=trim(prefix)//'.d.hn.out'
    open(152,file=filnam,err=888,status='unknown')

    ! QN: daily PPFD, mol/m2
    filnam=trim(prefix)//'.d.qn.out'
    open(153,file=filnam,err=888,status='unknown')

    ! CN: daily condensation water, mm
    filnam=trim(prefix)//'.d.cn.out'
    open(154,file=filnam,err=888,status='unknown')

    ! WN: daily soil moisture, mm
    filnam=trim(prefix)//'.d.wn.out'
    open(155,file=filnam,err=888,status='unknown')

    ! PN: daily precipitation, mm
    filnam=trim(prefix)//'.d.pn.out'
    open(156,file=filnam,err=888,status='unknown')

    ! RO: daily runoff, mm
    filnam=trim(prefix)//'.d.ro.out'
    open(157,file=filnam,err=888,status='unknown')

    ! EQ_N: daily equilibrium ET, mm
    filnam=trim(prefix)//'.d.eq_n.out'
    open(158,file=filnam,err=888,status='unknown')

    ! EP_N: daily potential ET, mm
    filnam=trim(prefix)//'.d.ep_n.out'
    open(159,file=filnam,err=888,status='unknown')

    ! EA_N: daily actual ET, mm
    filnam=trim(prefix)//'.d.ea_n.out'
    open(160,file=filnam,err=888,status='unknown')


    !----------------------------------------------------------------
    ! MONTHLY OUTPUT
    !----------------------------------------------------------------

    ! eq_m
    filnam=trim(prefix)//'.m.eq_m.out'
    open(251,file=filnam,err=888,status='unknown')

    ! ep_m
    filnam=trim(prefix)//'.m.ep_m.out'
    open(252,file=filnam,err=888,status='unknown')

    ! ea_m
    filnam=trim(prefix)//'.m.ea_m.out'
    open(253,file=filnam,err=888,status='unknown')

    ! cpa
    filnam=trim(prefix)//'.m.cpa.out'
    open(254,file=filnam,err=888,status='unknown')

    ! cwd
    filnam=trim(prefix)//'.m.cwd.out'
    open(255,file=filnam,err=888,status='unknown')

    ! qm
    filnam=trim(prefix)//'.m.qm.out'
    open(256,file=filnam,err=888,status='unknown')

    return

  888  stop 'INITIO_WATERBAL: error opening output files'

  end subroutine initio_waterbal


  subroutine initoutput_waterbal
    !////////////////////////////////////////////////////////////////
    !  Initialises waterbalance-specific output variables
    !----------------------------------------------------------------
    implicit none

    ! xxx remove their day-dimension
    outdwn(:,:,:) = 0.0
    outdho(:,:) = 0.0
    outdhn(:,:) = 0.0
    outdqn(:,:) = 0.0
    outdcn(:,:,:) = 0.0
    outdro(:,:,:) = 0.0
    outdeq_n(:,:) = 0.0
    outdep_n(:,:) = 0.0
    outdea_n(:,:,:) = 0.0

    ! xxx remove their month-dimension
    outmeq_m(:,:) = 0.0
    outmep_m(:,:) = 0.0
    outmea_m(:,:,:) = 0.0
    outmcpa(:,:,:) = 0.0
    outmcwd(:,:,:) = 0.0
    outmqm(:,:) = 0.0

    return

  end subroutine initoutput_waterbal


  subroutine getout_daily_waterbal( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    use fluxes, only: rn_d, dwn

    implicit none

    ! ARGUMENTS
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy    
    integer, intent(in) :: doy    

    ! Save the daily totals:
    ! xxx add lu-dimension and jpngr-dimension
    outdwn(:,doy,jpngr)   = dwn(:)
    outdho(doy,jpngr)     = ra_d  
    outdhn(doy,jpngr)     = rn_d
    outdqn(doy,jpngr)     = ppfd_d
    outdcn(:,doy,jpngr)   = wc(:)
    outdeq_n(doy,jpngr)   = eet_d
    outdep_n(doy,jpngr)   = pet_d
    outdea_n(:,doy,jpngr) = aet_d(:)

    ! Update monthly totals:
    outmeq_m(moy,jpngr)   = outmeq_m(moy,jpngr) + eet_d
    outmep_m(moy,jpngr)   = outmep_m(moy,jpngr) + pet_d
    outmea_m(:,moy,jpngr) = outmea_m(:,moy,jpngr) + aet_d(:)
    outmqm(moy,jpngr)     = outmqm(moy,jpngr) + ppfd_d

    return  

  end subroutine getout_daily_waterbal


  subroutine getout_monthly_waterbal( jpngr, moy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    implicit none

    ! ARGUMENTS
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy

    ! Calculate other monthly totals:
    outmcpa(:,moy,jpngr) = outmea_m(:,moy,jpngr)
    outmcpa(:,moy,jpngr) = outmcpa(:,moy,jpngr) / outmep_m(moy,jpngr)
    outmcwd(:,moy,jpngr) = outmep_m(moy,jpngr)
    outmcwd(:,moy,jpngr) = outmcwd(:,moy,jpngr) - outmea_m(:,moy,jpngr)

    return  

  end subroutine getout_monthly_waterbal


  subroutine writeout_ascii_waterbal( year, spinup )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use params_core, only: ndayyear, nmonth
    use params_siml, only: firstyeartrend, spinupyears

    implicit none

    ! Arguments
    integer, intent(in) :: year       ! simulation year
    logical, intent(in) :: spinup     ! true during spinup years

    ! Local variables
    real :: itime
    integer :: day, moy, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1


    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    !-------------------------------------------------------------------------
    if (.not.spinup) then
      ! Write daily output only during transient simulation
      do day=1,ndayyear

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

        if (nlu>1) stop 'writeout_ascii_waterbal: write out lu-area weighted sum'

        ! xxx lu-area weighted sum if nlu>0
        write(151,999) itime, outdho(day,jpngr)
        write(152,999) itime, outdhn(day,jpngr)
        write(153,999) itime, outdqn(day,jpngr)
        write(154,999) itime, outdcn(1,day,jpngr)
        write(155,999) itime, outdwn(1,day,jpngr)
        write(157,999) itime, outdro(1,day,jpngr)
        write(158,999) itime, outdeq_n(day,jpngr)
        write(159,999) itime, outdep_n(day,jpngr)
        write(160,999) itime, outdea_n(1,day,jpngr)

      end do
    end if


    !-------------------------------------------------------------------------
    ! MONTHLY OUTPUT
    ! Write monthly value
    !-------------------------------------------------------------------------
    do moy=1,nmonth

      itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(moy-1)/real(nmonth)

      ! STASH variables
      write(111,999) itime, outmeq_m(moy,jpngr)
      write(112,999) itime, outmep_m(moy,jpngr)
      write(113,999) itime, outmea_m(1,moy,jpngr)
      write(114,999) itime, outmcpa(1,moy,jpngr)
      write(115,999) itime, outmcwd(1,moy,jpngr)
      write(116,999) itime, outmqm(moy,jpngr)

    end do

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_waterbal


  function dgcos( x )
    !----------------------------------------------------------------   
    ! Calculates the cosine of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    use params_core, only: pi
    implicit none

    ! arguments
    real, intent(in) :: x  ! angle, degrees (0-360)

    ! function return value
    real, intent(out) :: dgcos ! cosine value of x when x is in degrees

    !dgcos = dcos(x*pi/180.0)
    dgcos = cos(x*pi/180.0)  ! xxx use cos with single-precision compilation

  end function dgcos


  function dgsin( x )
    !----------------------------------------------------------------   
    ! Calculates the sinus of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    use params_core, only: pi
    implicit none

    ! arguments
    real, intent(in) :: x  ! angle, degrees (0-360)

    ! function return value
    real, intent(out) :: dgsin ! sinus value of x when x is in degrees

    !dgsin = dsin(x*pi/180.0)
    dgsin = sin(x*pi/180.0)   ! xxx use cos with single-precision compilation

  end function dgsin


  function degrees( x )
    !----------------------------------------------------------------   
    ! Returns corresponding degrees if x is given in radians
    !----------------------------------------------------------------   
    use params_core, only: pi
    implicit none

    ! arguments
    real, intent(in) :: x  ! angle, radians

    ! function return value
    real, intent(out) :: degrees

    degrees = x*180.0/pi

  end function degrees


  function radians( x )
    !----------------------------------------------------------------   
    ! Returns corresponding radians if x is given in degrees
    !----------------------------------------------------------------   
    use params_core, only: pi
    implicit none

    ! arguments
    real, intent(in) :: x  ! angle, radians

    ! function return value
    real, intent(out) :: radians

    radians = x*pi/180.0

  end function radians


  function berger_tls( day )
    !----------------------------------------------------------------   
    ! Returns true anomaly and true longitude for a given day
    ! Reference: Berger, A. L. (1978), Long term variations of daily 
    ! insolation and quaternary climatic changes, J. Atmos. Sci., 35, 
    ! 2362-2367.
    !----------------------------------------------------------------   
    implicit none

    ! arguments
    integer, intent(in) :: day   ! day of the year

    ! local variables
    real :: anm, ranm, anv, ranv
    real :: dlamm                ! Mean longitude for day of year
    real :: my_nu
    real :: my_tls
    real :: xee, xec, xse        ! variable substitutes
    real :: xlam                 ! Mean longitude for vernal equinox
    real :: tmp1, tmp2, tmp3     ! variable substitutes

    ! function return value
    real, dimension(2) :: berger_tls ! return "tuple": (my_nu = true anomaly, degrees; my_tls = true longitude, degrees)

    ! Variable substitutes:
    xee = ke**2 
    xec = ke**3
    !xse = dsqrt(1.0 - xee)
    xse = sqrt(1.0 - xee)  ! xxx use sqrt with single-precision compilation

    ! Mean longitude for vernal equinox:
    tmp1 = (ke/2.0 + xec/8.0)*(1.0 + xse)*dgsin(komega)
    tmp2 = xee/4.0*(0.5 + xse)*dgsin(2.0*komega)
    tmp3 = xec/8.0*(1.0/3.0 + xse)*dgsin(3.0*komega)
    xlam = tmp1 - tmp2 + tmp3
    xlam = degrees(2.0*xlam)

    ! xxx debug
    !print*,'xlam', xlam

    ! Mean longitude for day of year:
    dlamm = xlam + (day - 80.0)*(360.0/ndayyear)

    ! Mean anomaly:
    anm = dlamm - komega
    ranm = radians(anm)

    ! True anomaly:
    !ranv = (ranm + (2.0*ke - xec/4.0)*dsin(ranm) + 5.0/4.0*xee*dsin(2.0*ranm) + 13.0/12.0*xec*dsin(3.0*ranm))
    ranv = (ranm + (2.0*ke - xec/4.0)*sin(ranm) + 5.0/4.0*xee*sin(2.0*ranm) + 13.0/12.0*xec*sin(3.0*ranm))  ! xxx use sin with single-precision compilation
    anv = degrees(ranv)

    ! True longitude:
    my_tls = anv + komega
    if (my_tls < 0.0) then
        my_tls = my_tls + 360.0
    else if (my_tls > 360.0) then
        my_tls = my_tls - 360.0
    end if

    ! True anomaly:
    my_nu = (my_tls - komega)
    if (my_nu < 0.0) then
        my_nu = my_nu + 360.0
    end if

    berger_tls = (/my_nu, my_tls/)

  end function berger_tls


  function sat_slope( tc )
    !----------------------------------------------------------------   
    ! Calculates the slope of the sat pressure temp curve, Pa/K
    ! Ref:      Eq. 13, Allen et al. (1998)
    !----------------------------------------------------------------   
    implicit none

    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C

    ! function return value
    real :: sat_slope  ! slope of the sat pressure temp curve, Pa/K

    sat_slope = (17.269)*(237.3)*(610.78)*(exp(tc*17.269/(tc + 237.3))/((tc + 237.3)**2))

  end function sat_slope


  function enthalpy_vap( tc )
    !----------------------------------------------------------------   
    ! Calculates the enthalpy of vaporization, J/kg
    ! Ref:      Eq. 8, Henderson-Sellers (1984)
    !----------------------------------------------------------------   
    implicit none

    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C

    ! function return value
    real ::  enthalpy_vap ! enthalpy of vaporization, J/kg

    enthalpy_vap = 1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2

  end function enthalpy_vap


  function elv2pres( alt )
    !----------------------------------------------------------------   
    ! Calculates atm. pressure for a given elevation
    ! Ref:      Allen et al. (1998)
    !----------------------------------------------------------------   
    implicit none

    ! arguments
    real, intent(in) :: alt ! elevation above sea level, m

    ! function return value
    real ::  elv2pres ! atm. pressure for a given elevation

    elv2pres = kPo*(1.0 - kL*alt/kTo)**(kG*kMa/(kR*kL))

  end function elv2pres


  function density_h2o( tc, press )
    !----------------------------------------------------------------   
    ! Calculates density of water at a given temperature and pressure
    ! Ref: Chen et al. (1977)
    !----------------------------------------------------------------   
    implicit none

    ! arguments
    real, intent(in) :: tc     ! air temperature (degrees C)
    real, intent(in) :: press  ! atmospheric pressure (Pa)

    ! local variables
    real :: po, ko, ca, cb
    real :: pbar               ! atmospheric pressure (bar)

    ! function return value
    real :: density_h2o  ! density of water, kg/m^3

    ! Calculate density at 1 atm:
    po = (&
             0.99983952&
             + 6.788260e-5  *tc&
             - 9.08659e-6   *tc*tc&
             + 1.022130e-7  *tc*tc*tc  &
             - 1.35439e-9   *tc*tc*tc*tc&
             + 1.471150e-11 *tc*tc*tc*tc*tc&
             - 1.11663e-13  *tc*tc*tc*tc*tc*tc&
             + 5.044070e-16 *tc*tc*tc*tc*tc*tc*tc&
             - 1.00659e-18  *tc*tc*tc*tc*tc*tc*tc*tc&
         )

    ! Calculate bulk modulus at 1 atm:
    ko = (&
             19652.17&
             + 148.1830   *tc&
             - 2.29995    *tc*tc&
             + 0.01281    *tc*tc*tc&
             - 4.91564e-5 *tc*tc*tc*tc&
             + 1.035530e-7*tc*tc*tc*tc*tc&
         )

    ! Calculate temperature dependent coefficients:
    ca = (&
             3.26138&
             + 5.223e-4  *tc&
             + 1.324e-4  *tc*tc&
             - 7.655e-7  *tc*tc*tc&
             + 8.584e-10 *tc*tc*tc*tc&
         )
    cb = (&
             7.2061e-5&
             - 5.8948e-6  *tc&
             + 8.69900e-8 *tc*tc&
             - 1.0100e-9  *tc*tc*tc&
             + 4.3220e-12 *tc*tc*tc*tc&
         )

    ! Convert atmospheric pressure to bar (1 bar = 100000 Pa)
    pbar = (1.0e-5)*press

    density_h2o = 1000.0*po*(ko + ca*pbar + cb*pbar**2.0) &
      /(ko + ca*pbar + cb*pbar**2.0 - pbar)

  end function density_h2o


  function psychro( tc, press )
    !----------------------------------------------------------------   
    ! Calculates the psychrometric constant for a given temperature and pressure
    ! Ref: Allen et al. (1998); Tsilingiris (2008) 
    !----------------------------------------------------------------   
    implicit none

    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C
    real, intent(in) :: press  ! atmospheric pressure, Pa

    ! local variables
    real :: lv  ! latent heat of vaporization (J/kg)
    real :: cp

    ! function return value
    real :: psychro  ! psychrometric constant, Pa/K

    ! Calculate the specific heat capacity of water, J/kg/K
    ! Eq. 47, Tsilingiris (2008)
    cp = 1.0e3*(&
               1.0045714270&
             + 2.050632750e-3  *tc&
             - 1.631537093e-4  *tc*tc&
             + 6.212300300e-6  *tc*tc*tc&
             - 8.830478888e-8  *tc*tc*tc*tc&
             + 5.071307038e-10 *tc*tc*tc*tc*tc&
            )

    ! Calculate latent heat of vaporization, J/kg
    lv = enthalpy_vap(tc)

    ! Calculate psychrometric constant, Pa/K
    ! Eq. 8, Allen et al. (1998)
    psychro = cp*kMa*press/(kMv*lv)

  end function psychro


end module waterbal_stash