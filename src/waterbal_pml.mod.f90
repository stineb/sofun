module md_waterbal
  !////////////////////////////////////////////////////////////////
  ! SPLASH WATERBALANCE MODULE
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
  !   - daytime net radiation ('rn')
  !   - soil water conent ('psoilphys%wcont')
  !   - runoff ('soilphys%dro')
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  ! ...
  !----------------------------------------------------------------
  use md_params_core, only: ndayyear, nmonth, nlu, maxgrid, kTo, kR, &
    kMv, kMa, kfFEC, secs_per_day, pi, dummy
  use md_params_core, only: ndayyear, ndaymonth, npft, kTkelvin
  use md_tile, only: tile_type, tile_fluxes_type
  use md_forcing, only: climate_type
  use md_grid, only: gridtype

  implicit none

  private
  public waterbal, solar,                           &
    getout_daily_waterbal, initoutput_waterbal,     &
    getpar_modl_waterbal, initio_nc_waterbal,       &
    writeout_nc_waterbal, get_rlm_waterbal,         & !, init_rlm_waterbal
    getrlm_daily_waterbal

  !-----------------------------------------------------------------------
  ! Uncertain (unknown) parameters. Runtime read-in
  !-----------------------------------------------------------------------
  real :: maxmeltrate       ! maximum snow melting rate (mm d-1) (Orth et al., 2013) 
  real :: kA                ! constant for dRnl (Monteith & Unsworth, 1990)
  real :: kalb_sw           ! shortwave albedo (Federer, 1968)
  real :: kalb_vis          ! visible light albedo (Sellers, 1985)
  real :: kb                ! constant for dRnl (Linacre, 1968)
  real :: kc                ! cloudy transmittivity (Linacre, 1968)
  real :: kCw               ! supply constant, mm/hr (Federer, 1982)
  real :: kd                ! angular coefficient of transmittivity (Linacre, 1968)
  real :: ke                ! eccentricity for 2000 CE (Berger, 1978)
  real :: keps              ! obliquity for 2000 CE, degrees (Berger, 1978)
  real :: kw                ! entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
  real :: komega            ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)
  real :: vpdstress_par_a   ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())
  real :: vpdstress_par_b   ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())
  real :: vpdstress_par_m   ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC, PRIVATE VARIABLES
  !----------------------------------------------------------------
  ! Orbit parameters 
  type outtype_berger
    real :: nu
    real :: lambda
  end type outtype_berger

  ! ! Radiation variables. aet, sw, and cpa are affected by soil moisture.
  ! type evaptype
  !   real :: rho_air      ! density of air (g m-3)
  !   real :: sat_slope    ! slope of saturation vap press temp curve, Pa/K 
  !   real :: lv           ! enthalpy of vaporization, J/kg
  !   real :: rho_water    ! density of water (g m-3)
  !   real :: econ         ! water-to-energy conversion factor (econ), m^3/J
  !   real :: rn           ! daily net radiation (J/m2/d)
  !   real :: rnn          ! nighttime net radiation (J/m^2/d)
  !   real :: rnl          ! net longwave radiation (W/m^2)
  !   real :: cn           
  !   real :: aet          
  !   real :: aet_e        
  ! end type evaptype

  ! holds return variables of function get_snow_rain()
  type outtype_snow_rain
    real :: snow_updated     ! snow depth in water equivalents (mm)
    real :: liquid_to_soil   ! water 
  end type outtype_snow_rain

  ! Holds output of function calc_et 
  type outtype_et
    real :: cn             ! daily condensation (mm d-1)
    real :: pet            ! daily potential evapotranspiration (mm d-1)
    real :: pet_e          ! daily potential evapotranspiration (J m-2 d-1)
    real :: aet            ! daily actual evapotranspiration (mm d-1)
    real :: aet_e          ! daily actual evapotranspiration (J m-2 d-1)
    real :: aet_soil       ! daily soil evaporation (mm d-1)
    real :: aet_e_soil     ! daily soil evaporation (J m-2 d-1)
    real :: aet_canop      ! daily canopy evaporation (mm d-1)
    real :: aet_e_canop    ! daily canopy evaporation (J m-2 d-1)
  end type outtype_et

  type outtype_netrad
    real :: rnl
    real :: rn
    real :: rnn
  end type outtype_netrad

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC, KNOWN PARAMETERS
  !----------------------------------------------------------------
  logical :: outenergy = .true.

  !----------------------------------------------------------------
  ! Module-specific rolling mean variables
  !----------------------------------------------------------------
  real, dimension(nlu), save :: aalpha                      ! annual mean alpha (AET/PET)
  integer, parameter :: nyrs_rlmalpha = 5                   ! number of years for rolling mean (=width of sliding window)

  !----------------------------------------------------------------
  ! Daily module-specific output variables
  !----------------------------------------------------------------
  real, allocatable, dimension(:,:)   :: outdpet            ! daily potential ET, mm r J/m2/d depending on 'outenergy'
  real, allocatable, dimension(:,:,:) :: outdaet            ! daily actual ET, mm or J/m2/d depending on 'outenergy'
  real, allocatable, dimension(:,:,:) :: outdwbal           ! daily water balance, mm
  real, allocatable, dimension(:,:,:) :: outdwcont          ! daily water content = soil moisture, mm
  real, allocatable, dimension(:,:,:) :: outdrn             ! daily net radiation, J/m2
  ! real, allocatable, dimension(:,:)   :: outdra             ! daily solar irradiation, J/m2
  ! real, allocatable, dimension(:,:)   :: outdayl            ! daily day length, h
  ! real, allocatable, dimension(:,:)   :: outdcn             ! daily condensation water, mm
  ! real, allocatable, dimension(:,:,:) :: outdro             ! daily runoff, mm
  ! real, allocatable, dimension(:,:,:) :: outdfleach         ! daily NO3 leaching fraction, (unitless)
  ! real, allocatable, dimension(:,:)   :: outdeet            ! daily equilibrium ET, mm
  ! real, allocatable, dimension(:,:)   :: outdecon           ! daily water-to-energy conversion factor m TJ-1 = mm GJ-1

  !----------------------------------------------------------------
  ! Annual module-specific output variables
  !----------------------------------------------------------------
  real, allocatable, dimension(:)     :: outapet            ! annual total potential ET, mm r J/m2/yr depending on 'outenergy'
  real, allocatable, dimension(:,:)   :: outaaet            ! annual total actual ET, mm or J/m2/yr depending on 'outenergy'
  real, allocatable, dimension(:,:)   :: outaalpha          ! annual mean AET/PET (of daily values!), unitless

  !----------------------------------------------------------------
  ! Module-specific NetCDF output file and variable names
  !----------------------------------------------------------------
  ! Annual output files
  character(len=256) :: ncoutfilnam_apet
  character(len=256) :: ncoutfilnam_aaet
  character(len=256) :: ncoutfilnam_aalpha

  ! Daily output files
  character(len=256) :: ncoutfilnam_dwcont
  character(len=256) :: ncoutfilnam_dwbal
  character(len=256) :: ncoutfilnam_dpet
  character(len=256) :: ncoutfilnam_daet
  character(len=256) :: ncoutfilnam_drn

  character(len =*), parameter :: WCONT_NAME = "wcont"
  character(len =*), parameter :: PET_NAME = "pet"
  character(len =*), parameter :: AET_NAME = "aet"
  character(len =*), parameter :: RN_NAME = "netrad"
  character(len =*), parameter :: ALPHA_NAME = "alpha"
  character(len =*), parameter :: WBAL_NAME = "wbal"

  character(len=7) :: in_ppfd       ! information whether PPFD is prescribed from meteo file for global attribute in NetCDF file


  logical, parameter :: splashtest = .false.

contains

  subroutine waterbal( tile, tile_fluxes, grid, climate, doy ) !, lai, fapar, h_canopy, g_stomata )
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates soil water balance
    !-------------------------------------------------------------------------
    ! arguments
    type(tile_type), dimension(nlu), intent(inout)        :: tile
    type(tile_fluxes_type), dimension(nlu), intent(inout) :: tile_fluxes
    type(gridtype), intent(in)                            :: grid
    type(climate_type), intent(in)                        :: climate
    integer, intent(in) :: doy          ! day of year

    ! real, intent(in)    :: pr           ! daily precip as rain (liquid) (mm) 
    ! real, intent(in)    :: sn           ! daily precip as snow (mm water equivalent) 
    ! real, intent(in)    :: tc           ! mean monthly temperature (deg C)
    ! real, intent(in)    :: sf           ! mean monthly sunshine fraction (unitless)
    ! real, intent(in)    :: netrad       ! net radiation (J m-2 d-1), may be dummy (in which case this is not used)
    ! real, intent(in)    :: vpd          ! vapour pressure deficit (Pa)
    ! real, intent(in)    :: patm         ! atmospheric pressure (Pa)
    ! real, intent(in)    :: vwind        ! wind speed (m s-1)

    ! local variables
    ! type(canopy_type), pointer :: canopy
    type(outtype_netrad)    :: out_netrad
    type(outtype_et)        :: out_et
    type(outtype_snow_rain) :: out_snow_rain
    real                    :: g_aero         ! aerodynamic conductance (m s-1)
    real                    :: g_canopy       ! canopy conductance (m s-1)
    integer                 :: lu             ! land unit (gridcell tile)
    real                    :: tk


    tk = climate%dtemp + kTkelvin

    ! Loop over gricell tiles
    do lu=1,nlu

      ! canopy => tile(lu)%canopy

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculate canopy conductance (m s-1)
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      g_aero = calc_g_aero( climate%dvwind, tile(lu)%canopy%height )

      ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ! Net radiation (J m-2 d-1)
      ! ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! call calc_netrad( tile_fluxes(lu), grid, climate, doy )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Update soil moisture and snow pack
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      out_snow_rain = get_snow_rain( climate%dprec + out_et%cn, climate%dsnow, climate%dtemp, tile(lu)%soil%phy%snow )
      tile(lu)%soil%phy%snow = out_snow_rain%snow_updated 

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Canopy transpiration and soil evaporation
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call calc_et( tile(lu), tile_fluxes(lu), climate, out_snow_rain%liquid_to_soil, g_aero )

      ! Update soil moisture
      tile(lu)%soil%phy%wcont = tile(lu)%soil%phy%wcont + out_snow_rain%liquid_to_soil - out_et%aet

      ! Bucket model for runoff generation
      if (tile(lu)%soil%phy%wcont > tile(lu)%soil%params%whc) then
        ! -----------------------------------
        ! Bucket is full 
        ! -----------------------------------
        ! determine NO3 leaching fraction 
        tile_fluxes(lu)%canopy%dfleach = 1.0 - tile(lu)%soil%params%whc / tile(lu)%soil%phy%wcont

        ! add remaining water to monthly runoff total
        tile_fluxes(lu)%canopy%dro = tile(lu)%soil%phy%wcont - tile(lu)%soil%params%whc

        ! set soil moisture to capacity
        tile(lu)%soil%phy%wcont = tile(lu)%soil%params%whc

      else if (tile(lu)%soil%phy%wcont < 0.0) then
        ! -----------------------------------
        ! Bucket is empty
        ! -----------------------------------
        ! set soil moisture to zero
        out_et%aet = out_et%aet + tile(lu)%soil%phy%wcont
        tile(lu)%soil%phy%wcont        = 0.0
        tile_fluxes(lu)%canopy%dro     = 0.0
        tile_fluxes(lu)%canopy%dfleach = 0.0

      else
        ! No runoff occurrs
        tile_fluxes(lu)%canopy%dro     = 0.0
        tile_fluxes(lu)%canopy%dfleach = 0.0

      end if

      ! water scalar (fraction of plant-available water holding capacity; water storage at wilting point is already accounted for in tile(lu)%soil%params%whc)
      tile(lu)%soil%phy%wscal = tile(lu)%soil%phy%wcont / tile(lu)%soil%params%whc

      ! save daily water fluxes
      tile_fluxes(lu)%canopy%dwbal        = out_snow_rain%liquid_to_soil
      tile_fluxes(lu)%canopy%dpet         = out_et%pet
      tile_fluxes(lu)%canopy%dpet         = out_et%pet
      tile_fluxes(lu)%canopy%daet         = out_et%aet
      tile_fluxes(lu)%canopy%daet_e       = out_et%aet_e
      tile_fluxes(lu)%canopy%daet_soil    = out_et%aet_soil
      tile_fluxes(lu)%canopy%daet_e_soil  = out_et%aet_e_soil
      tile_fluxes(lu)%canopy%daet_canop   = out_et%aet_canop
      tile_fluxes(lu)%canopy%daet_e_canop = out_et%aet_e_canop

    end do

  end subroutine waterbal


  subroutine solar( tile_fluxes, grid, sf, doy)
    !/////////////////////////////////////////////////////////////////////////
    ! This subroutine calculates daily PPFD. Code is an extract of the subroutine
    ! 'evap', adopted from the evap() function in GePiSaT (Python version). 
    ! This subroutine ('get_solar') is called before the daily loop.
    ! Output:
    ! - daily extraterrestrial solar radiation (dra), J/m^2
    ! - daily PPFD (dppfd), mol/m^2
    !-------------------------------------------------------------------------  
    use md_params_core, only: ndayyear, pi, dummy, kGsc
    use md_sofunutils, only: daily2monthly
    use md_tile, only: tile_fluxes_type
    use md_grid, only: gridtype

    ! arguments
    type(tile_fluxes_type), dimension(nlu), intent(inout) :: tile_fluxes
    type(gridtype), intent(inout) :: grid
    real, intent(in) :: sf            ! fraction of sunshine hours 
    integer, intent(in) :: doy        ! day of year

    ! local variables
    real    :: dr                    ! distance factor
    real    :: delta                 ! declination angle 
    real    :: ru                    ! variable substitute for u
    real    :: rv                    ! variable substitute for v
    real    :: hs                    ! sunset hour angle
    real    :: tau                   ! transmittivity (unitless)
    real    :: rw                    ! variable substitute (W/m^2)
    real, dimension(2) :: out_ru_rv  ! function return variable containing 'ru' and 'rv'.

    ! real :: daysecs ! daylight seconds for each DOY
    ! real, dimension(nmonth)   :: monsecs ! daylight seconds for each MOY

    type(outtype_berger) :: out_berger    ! stores output of function berger_tls

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 2. Calculate heliocentric longitudes (nu and lambda), degrees
    ! Store daily return values for later use in subroutine 'evap'.
    ! Other variables defined and over-written below may be stored
    ! for later use in 'evap' the same way. However function 
    ! 'out_berger' is by far the most expensive one. This is there-
    ! fore a compromise.
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Berger (1978)
    out_berger  = get_berger_tls( doy )
    grid%nu     = out_berger%nu
    grid%lambda = out_berger%lambda

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3. Calculate distance factor (dr), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dr = calc_dr( grid%nu )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 4. Calculate declination angle (delta), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    delta = calc_delta( grid%lambda )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 5. Calculate variable substitutes (u and v), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    out_ru_rv = calc_ru_rv( delta, lat )
    ru = out_ru_rv(1)
    rv = out_ru_rv(2)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 6. Calculate the sunset hour angle (hs), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    hs = calc_hs( ru, rv )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 6.a Calculate day length from sunset hour angle, h
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tile_fluxes(:)%canopy%dayl = 24.0 * hs / 180.0  ! hs is in degrees (pi = 180 deg)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 7. Calculate daily extraterrestrial solar radiation (dra), J/m^2/d
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 1.10.3, Duffy & Beckman (1993)
    tile_fluxes(:)%canopy%dra = ( secs_per_day / pi ) * kGsc * dr * ( radians(ru) * hs + rv * dgsin(hs) )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 8. Calculate transmittivity (tau), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tau = calc_tau( sf, elv )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 9. Calculate daily PPFD (dppfd), mol/m^2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 57, SPLASH 2.0 Documentation
    tile_fluxes(:)%canopy%ppfd_splash = (1.0e-6) * kfFEC * ( 1.0 - kalb_vis ) * tau * tile_fluxes(:)%canopy%dra

    if (splashtest) then
      print*,'transmittivity, tau: ', tau
      print*,'daily TOA radiation: ', (1.0e-6)*tile_fluxes(:)%canopy%dra
      print*,'sunset angle, hs: ', hs
      print*,'true anomaly, nu: ', tile_fluxes(:)%canopy%nu
      print*,'true longitude, lambda: ', tile_fluxes(:)%canopy%lambda
      print*,'distance factor, dr: ', dr
      print*,'declination, delta: ', delta
      print*,'variable substitute, ru: ', ru
      print*,'variable substitute, rv: ', rv
      print*,'daily PPFD: ', tile_fluxes(:)%canopy%ppfd_splash
    end if

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
  end subroutine solar


  subroutine calc_et( tile, tile_fluxes, climate, liquid_to_soil, g_aero )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns ecosystem-level evapotranspiration as a function of LAI following
    ! Eq. 1 in Zhang et al., 2017 JGR (doi:10.1002/2017JD027025)
    !-------------------------------------------------------------------------  
    use md_params_core, only: cp, kTkelvin
    use md_sofunutils, only: dampen_variability

    ! arguments
    type(tile_type), intent(inout) :: tile
    type(tile_fluxes_type), intent(inout) :: tile_fluxes
    type(climate_type), intent(in) :: climate
    real, intent(in) :: liquid_to_soil  ! liquid water infiltrating the soil (mm d-1)
    real, intent(in) :: g_aero          ! aerodynamic conductance (m s-1)

    ! real, intent(in) :: tc              ! air temperature (deg C)
    ! real, intent(in) :: pr              ! precipitation (mm d-1)
    ! real, intent(in) :: patm            ! atmospheric pressure (Pa)
    ! real, intent(in) :: lai             ! leaf area index (m2-leaf m-2-ground)
    ! real, intent(in) :: fapar           ! fraction of absorbed photosynthetically active radiation (unitless)
    ! real, intent(in) :: netrad          ! net radiation (J m-2 d-1)   xxx in Zhang described as (MJ m-2 d-1)
    ! real, intent(in) :: vpd             ! vapour pressure deficit of air (Pa)
    ! real, intent(in) :: g_canopy        ! canopy conductance (m s-1)
    ! real, intent(in) :: g_aero          ! aerodynamic conductance (m s-1)

    ! function return variable
    type(outtype_et) :: out_et

    ! local variables
    real :: e_avl_soil                      ! available energy of soil surface (corresponds to net radiation minus ground heat flux) (MJ m-2 d-1)
    real :: e_avl_canopy                    ! available energy of canopy (corresponds to net radiation minus ground heat flux (MJ m-2 d-1)
    real :: epsilon                         ! = s/gamma; s: slope of the  
    real :: gamma                           ! psychrometric constant (Pa K-1) ! xxx Zhang et al. use it in units of (kPa K-1), probably they use sat_slope in kPa/K, too.
    real :: sat_slope                       ! slope of saturation vapour pressure vs. temperature curve, Pa K-1
    real :: rho_air                         ! density of air (g m-3)
    real :: lv                              ! enthalpy of vaporization, J/kg
    real :: rho_water                       ! density of water (g m-3)
    real :: econ                            ! water-to-energy conversion factor (econ), m^3/J
    real :: p_over_pet                      ! precipitation divided by potential evapotranspiration
    real, save :: p_over_pet_memory = 1.0   ! precipitation over equilibrium evapotranspiration, damped variability

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Calculate temperature and pressure-dependent quantities
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Density of (dry?) air
    rho_air = calc_density_air(climate%dtemp, climate%dpatm)

    ! Psychrometric constant, Pa/K
    gamma = psychro(climate%dtemp, climate%dpatm)

    ! Slope of saturation vap press temp curve, Pa/K
    sat_slope = calc_sat_slope(climate%dtemp)

    ! Enthalpy of vaporization, J/kg
    lv = calc_enthalpy_vap(climate%dtemp)
    
    ! Density of water, kg/m^3
    rho_water = density_h2o(climate%dtemp, climate%dpatm)

    ! Eq. 51, SPLASH 2.0 Documentation
    ! out_evap%econ = 1.0 / ( lv * rho_water ) ! this is to convert energy into mass (water)
    econ = sat_slope / (lv * rho_water * (sat_slope + gamma)) ! MORE PRECISELY - this is to convert energy into mass (water)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Daily condensation, mm d-1
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (netrad < 0.0) then
      tile_fluxes%canopy%dcn = 1000.0 * econ * abs(netrad)
    else
      tile_fluxes%canopy%dcn = 0.0
    end if

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Implementation of PML model (Zhang et al., 2017)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! variable substitute as used in Zhang et al. 2017 JGR    
    epsilon = sat_slope / gamma

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Energy partitioning to soil and canopy
    ! assuming zero ground heat flux 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    e_avl_soil = netrad * (1.0 - fapar)
    e_avl_canopy = netrad - e_avl_soil

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Calculate soil evaporation
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tile_fluxes%canopy%deet_e = epsilon * e_avl_soil / (epsilon + 1.0) 
    tile_fluxes%canopy%deet = tile_fluxes%canopy%deet_e * econ * 1000.0
    p_over_pet = liquid_to_soil / tile_fluxes%canopy%deet

    p_over_pet_memory = dampen_variability( p_over_pet, 30.0, p_over_pet_memory )   ! corresponds to f in Zhang et al., 2017 Eq. 9
    
    tile_fluxes%canopy%daet_soil = min(p_over_pet_memory, 1.0) * tile_fluxes%canopy%deet
    tile_fluxes%canopy%daet_e_soil = tile_fluxes%canopy%daet_soil / (econ * 1000.0)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Penman-Monteith Equation for canopy transpiration; Eq. 1 in Zhang et al., 2017
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tile_fluxes%canopy%aet_e_canop = (epsilon * e_avl_canopy + (rho_air * cp / gamma) * vpd * g_aero) / (1.0 + epsilon + g_aero / g_canopy)  ! g_canopy from P-model!

    ! Total evapotranspiration 
    tile_fluxes%canopy%aet_e     = tile_fluxes%canopy%aet_e_canop + tile_fluxes%canopy%daet_e_soil

    ! quantities in mm water (equivalent to kg m-2)
    tile_fluxes%canopy%aet       = tile_fluxes%canopy%aet_e       * econ * 1000.0
    tile_fluxes%canopy%daet_soil  = tile_fluxes%canopy%daet_e_soil  * econ * 1000.0
    tile_fluxes%canopy%aet_canop = tile_fluxes%canopy%aet_e_canop * econ * 1000.0

  end subroutine calc_et


  subroutine calc_netrad( tile_fluxes, grid, climate, doy ) result( out_netrad )
    !/////////////////////////////////////////////////////////////////////////
    ! This subroutine calculates daily quantities required for evapotranspiration
    ! calculations. 
    ! Implements the SPLASH model (Davis et al., 2017 GMD; doi:10.5194/gmd-10-689-2017)
    !-------------------------------------------------------------------------  
    use md_params_core, only: ndayyear, pi, dummy, kGsc
    use md_sofunutils, only: calc_patm

    ! arguments
    type(tile_fluxes_type), intent(inout) :: tile_fluxes
    type(gridtype), intent(in)            :: grid
    type(climate_type), intent(in)        :: climate
    integer, intent(in)                   :: doy          ! day of year

    ! function return variable
    type(outtype_netrad)  :: out_netrad

    ! local variables
    real :: dr                           ! distance factor
    real :: delta                        ! declination angle 
    real :: ru                           ! variable substitute for u
    real :: rv                           ! variable substitute for v
    real :: hs                           ! sunset hour angle
    real :: tau                          ! transmittivity (unitless)
    real :: rw                           ! variable substitute (W/m^2)
    real :: hn                           ! net radiation cross-over hour angle
    real, dimension(2) :: out_ru_rv      ! function return variable containing 'ru' and 'rv'.

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3. Calculate distance factor (dr), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dr = calc_dr( nu )
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 4. Calculate declination angle (delta), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    delta = calc_delta( lambda )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 5. Calculate variable substitutes (u and v), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    out_ru_rv = calc_ru_rv( delta, lat )
    ru = out_ru_rv(1)
    rv = out_ru_rv(2)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 6. Calculate the sunset hour angle (hs), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    hs = calc_hs( ru, rv )
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 8. Calculate transmittivity (tau), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tau = calc_tau( sf, elv )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 10. Estimate net longwave radiation (out_evap%rnl), W m-2
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
    out_netrad%rnl = ( kb + (1.0 - kb ) * sf ) * ( kA - tc )
  
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 11. Calculate variable substitute (rw), W m-2 -- shortwave radiation?
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rw = ( 1.0 - kalb_sw ) * tau * kGsc * dr

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 12. Calculate net radiation cross-over hour angle (hn), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ((out_netrad%rnl - rw*ru)/(rw*rv) >= 1.0) then
      ! Net radiation negative all day
      hn = 0.0
    else if ((out_netrad%rnl - rw*ru)/(rw*rv) <= -1.0) then
      ! Net radiation positive all day
      hn = 180.0
    else
      !hn = degrees( dacos((out_netrad%rnl - rw*ru)/(rw*rv)) )
      hn = degrees( acos((out_netrad%rnl - rw*ru)/(rw*rv)) )   ! use acos with single precision compilation
    end if

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 13. Calculate daytime total net radiation (out_netrad%rn), J m-2 d-1
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 53, SPLASH 2.0 Documentation
    out_netrad%rn = (secs_per_day/pi) * (hn*(pi/180.0)*(rw*ru - out_netrad%rnl) + rw*rv*dgsin(hn))

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 14. Calculate nighttime total net radiation (out_netrad%rnn), J m-2 d-1
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Eq. 56, SPLASH 2.0 Documentation
    ! adopted bugfix from Python version (iss#13)
    out_netrad%rnn = (86400.0/pi)*(radians(rw*ru*(hs-hn)) + rw*rv*(dgsin(hs)-dgsin(hn)) - out_netrad%rnl * (pi - radians(hn)))

  end function calc_netrad


  function calc_g_aero(v_wind, h_canopy) result( g_aero )
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates aerodynamic conductance as a function of wind speed and
    ! vegetation height. 
    ! Following Eq. 7 in Zhang et al., 2017 JGR (doi:10.1002/2017JD027025)
    !-------------------------------------------------------------------------  
    use md_params_core, only: k_karman

    ! arguments
    real, intent(in) :: v_wind    ! wind speed (m s-1)
    real, intent(in) :: h_canopy  ! canopy height (m) 

    ! function return variable
    real :: g_aero                ! aerodynamic conductance (m s-1)

    ! local variables
    real :: zom, zov, zm, d_canopy

    zom = 0.123 * h_canopy
    zov = 0.1 * zom
    zm = 10.0                     ! measurement height for wind speed (m), must be consistent with input data information
    d_canopy = (2.0/3.0) * h_canopy
    g_aero = (k_karman**2 * v_wind) / (log((zm - d_canopy)/zom) * log((zm - d_canopy)/zov))

  end function calc_g_aero


  function calc_density_air(tc_air, patm) result( rho )
    !/////////////////////////////////////////////////////////////////////////
    ! Air density of moist air from air temperature and pressure.
    !-------------------------------------------------------------------------
    use md_params_core, only: kTkelvin, Rd

    ! arguments
    real, intent(in) :: tc_air  ! air temperature (degrees Celsius)
    real, intent(in) :: patm    ! atmospheric pressure (Pa)

    ! local 
    real :: tk

    !function return variable
    real :: rho

    tk = tc_air + kTkelvin
    rho = patm / (Rd * tk)

  end function calc_density_air

  
  function get_snow_rain( pr, sn, tc, snow ) result( out_snow_rain )
    !/////////////////////////////////////////////////////////////////////////
    ! Translates precipitation into change in snow depth and liquid water
    ! input to soil.
    !-------------------------------------------------------------------------
    ! arguments
    real, intent(in) :: pr     ! daily precip (mm), includes condensation
    real, intent(in) :: sn     ! daily precip as snow (mm water equivalent) 
    real, intent(in) :: tc     ! mean monthly temperature (deg C)
    real, intent(in) :: snow   ! snow depth, water equivalents (mm)

    ! function return variable
    type(outtype_snow_rain) :: out_snow_rain

    ! local variables
    real :: fsnow                             ! fraction of precipitation as snow (temperature dependent)
    real :: melt                              ! snow melting rate (mm d-1)
    real, parameter :: temp_threshold = 1.0   ! deg C

    if ( snow > 0.0 .and. tc > temp_threshold ) then
      melt  = min( snow, maxmeltrate * ( tc - temp_threshold ) )
    else
      melt = 0.0
    end if 

    if (sn==dummy) then
      fsnow = max( min( 1.0,  1.0 - ( 1.0 / 2.0 ) * tc ), 0.0 )
      out_snow_rain%snow_updated   = snow + fsnow * pr - melt
      out_snow_rain%liquid_to_soil = pr * ( 1.0 - fsnow ) + melt
    else
      out_snow_rain%snow_updated   = snow + sn - melt
      out_snow_rain%liquid_to_soil = pr + melt
    end if

  end function get_snow_rain


  function calc_vdpstress( vpd ) result( out_vpdstress )
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates a VPD stress function based on Oren et al. 2001 Eq. 4
    ! Reference: 
    ! Oren et al.: Sensitivity of mean canopy stomatal conductance
    ! to vapor pressure deficit in a flooded Taxodium distichum L. forest:
    ! hydraulic and non-hydraulic effectsOecologia (2001) 126:21–29, 
    ! DOI 10.1007/s004420000497
    !-------------------------------------------------------------------------
    ! arguments
    real, intent(in) :: vpd    ! Vapour pressure deficit (Pa)

    ! function return variable
    real :: out_vpdstress

    if (vpd<1) then
      out_vpdstress = 1.0
    else
      out_vpdstress = vpdstress_par_a * (vpdstress_par_b - vpdstress_par_m * (log(0.001) + log(vpd)))
      if (out_vpdstress > 1.0) out_vpdstress = 1.0
      if (out_vpdstress < 0.0) out_vpdstress = 0.0
    end if

  end function calc_vdpstress


  function calc_dr( nu ) result( dr )
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Calculates distance factor (dr), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! arguments
    real, intent(in) :: nu

    ! local variables
    real :: rho

    ! function return variable
    real :: dr

    ! Berger et al. (1993)
    rho = (1.0 - ke**2)/(1.0 + ke * dgcos( nu ))        
    dr = (1.0/rho)**2

  end function calc_dr


  function calc_delta( lambda ) result( delta )
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Calculates declination angle (delta), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! arguments
    real, intent(in) :: lambda

    ! function return variable
    real :: delta

    ! Woolf (1968)
    delta = asin( dgsin( lambda ) * dgsin( keps ) )   ! xxx use asin with single-precision compilation
    delta = degrees( delta )

  end function calc_delta


  function calc_ru_rv( delta, lat ) result( out_ru_rv )
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Calculates variable substitutes (u and v), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! arguments
    real, intent(in) :: delta
    real, intent(in) :: lat

    ! local variables
    real :: ru, rv

    ! function return variable
    real, dimension(2) :: out_ru_rv

    ru = dgsin(delta) * dgsin(lat)
    rv = dgcos(delta) * dgcos(lat)

    out_ru_rv(1) = ru
    out_ru_rv(2) = rv

  end function calc_ru_rv


  function calc_hs( ru, rv ) result( hs )
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Calculates the sunset hour angle (hs), degrees
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! arguments
    real, intent(in) :: ru, rv

    ! function return variable
    real :: hs

    ! Note: u/v == tan(delta)*tan(lat)
    ! Eq. 3.22, Stine & Geyer (2001)
    if ((ru/rv) >= 1.0) then
      ! Polar day (no sunset)
      hs = 180.0 
    elseif ((ru/rv) <= -1.0) then
      ! Polar night (no sunrise)
      hs = 0.0
    else
      hs = degrees(acos(-1.0*ru/rv))
    end if

  end function calc_hs


  function calc_tau( sf, elv ) result( tau )
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Calculates transmittivity (tau), unitless
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! arguments
    real, intent(in) :: sf     ! sunshine fraction
    real, intent(in) :: elv    ! elevation

    ! local variables
    real :: tau_o

    ! function return variable
    real :: tau

    ! Eq. 11, Linacre (1968)
    tau_o = (kc + kd*sf)

    ! Eq. 2, Allen (1996)
    tau = tau_o * (1.0 + (2.67e-5) * elv)

  end function calc_tau


  subroutine getpar_modl_waterbal()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads waterbalance module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal
    use md_interface, only: interface

    print*,'reading waterbal parameters ...'

    ! constant for dRnl (Monteith & Unsworth, 1990)
    kA       = getparreal( 'params/params_waterbal_splash.dat', 'kA' )
    
    ! shortwave albedo (Federer, 1968)
    kalb_sw  = getparreal( 'params/params_waterbal_splash.dat', 'kalb_sw' )
    
    ! visible light albedo (Sellers, 1985) xxx planetary albedo? xxx
    kalb_vis = getparreal( 'params/params_waterbal_splash.dat', 'kalb_vis' )
    
    ! constant for dRnl (Linacre, 1968)
    kb       = getparreal( 'params/params_waterbal_splash.dat', 'kb' )
    
    ! cloudy transmittivity (Linacre, 1968)
    kc       = getparreal( 'params/params_waterbal_splash.dat', 'kc' )
    
    ! supply constant, mm/hr (Federer, 1982)
    kCw      = getparreal( 'params/params_waterbal_splash.dat', 'kCw' )
    
    ! angular coefficient of transmittivity (Linacre, 1968)
    kd       = getparreal( 'params/params_waterbal_splash.dat', 'kd' )
    
    ! eccentricity for 2000 CE (Berger, 1978)
    ke       = getparreal( 'params/params_waterbal_splash.dat', 'ke' )
    
    ! obliquity for 2000 CE, degrees (Berger, 1978)
    keps     = getparreal( 'params/params_waterbal_splash.dat', 'keps' )

    ! entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
    kw       = getparreal( 'params/params_waterbal_splash.dat', 'kw' )
    
    ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)
    komega   = getparreal( 'params/params_waterbal_splash.dat', 'komega' )

    ! maximum snow melting rate (mm d-1) (Orth et al., 2013)
    maxmeltrate = getparreal( 'params/params_waterbal_splash.dat', 'maxmeltrate' )

    if (interface%params_siml%is_calib) then

      ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())
      vpdstress_par_a = interface%params_calib%vpdstress_par_a

      ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())
      vpdstress_par_b = interface%params_calib%vpdstress_par_b

      ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())
      vpdstress_par_m = interface%params_calib%vpdstress_par_m


    else

      ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())
      vpdstress_par_a = getparreal( 'params/params_waterbal_splash.dat', 'vpdstress_par_a' )

      ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())
      vpdstress_par_b = getparreal( 'params/params_waterbal_splash.dat', 'vpdstress_par_b' )

      ! Parameter for Oren et al.-VPD stress function (see calc_vpdstress())
      vpdstress_par_m = getparreal( 'params/params_waterbal_splash.dat', 'vpdstress_par_m' )

    end if

    print*,'... done'

  end subroutine getpar_modl_waterbal

  ! xxx put these functions into a 'contain' within calling SR?

  function dgcos( x ) result( dgcos_out )
    !----------------------------------------------------------------   
    ! Calculates the cosine of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    use md_params_core, only: pi

    ! arguments
    real, intent(in) :: x  ! angle, degrees (0-360)

    ! function return value
    real :: dgcos_out ! cosine value of x when x is in degrees

    !dgcos = dcos(x*pi/180.0)
    dgcos_out = cos(x*pi/180.0)  ! xxx use cos with single-precision compilation

  end function dgcos


  function dgsin( x ) result( dgsin_out )
    !----------------------------------------------------------------   
    ! Calculates the sinus of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    use md_params_core, only: pi

    ! arguments
    real, intent(in) :: x  ! angle, degrees (0-360)

    ! function return value
    real :: dgsin_out ! sinus value of x when x is in degrees

    !dgsin_out = dsin(x*pi/180.0)
    dgsin_out = sin(x*pi/180.0)   ! xxx use cos with single-precision compilation

  end function dgsin


  function degrees( x ) result( degrees_out )
    !----------------------------------------------------------------   
    ! Returns corresponding degrees if x is given in radians
    !----------------------------------------------------------------   
    use md_params_core, only: pi

    ! arguments
    real, intent(in) :: x  ! angle, radians

    ! function return value
    real :: degrees_out

    degrees_out = x*180.0/pi

  end function degrees


  function radians( x ) result( radians_out )
    !----------------------------------------------------------------   
    ! Returns corresponding radians if x is given in degrees
    !----------------------------------------------------------------   
    use md_params_core, only: pi

    ! arguments
    real, intent(in) :: x  ! angle, radians

    ! function return value
    real :: radians_out

    radians_out = x*pi/180.0

  end function radians


  function get_berger_tls( day ) result( out_berger )
    !----------------------------------------------------------------   
    ! Returns true anomaly and true longitude for a given day
    ! Reference: Berger, A. L. (1978), Long term variations of daily 
    ! insolation and quaternary climatic changes, J. Atmos. Sci., 35, 
    ! 2362-2367.
    !----------------------------------------------------------------   
    ! arguments
    integer, intent(in) :: day   ! day of the year

    ! function return value
    type(outtype_berger) :: out_berger  ! stores output of function berger_tls

    ! local variables
    real :: anm, ranm, anv, ranv
    real :: dlamm                ! Mean longitude for day of year
    real :: my_nu
    real :: my_tls
    real :: xee, xec, xse        ! variable substitutes
    real :: xlam                 ! Mean longitude for vernal equinox
    real :: tmp1, tmp2, tmp3     ! variable substitutes

    ! Variable substitutes:
    xee = ke**2 
    xec = ke**3
    xse = sqrt(1.0 - xee)

    ! Mean longitude for vernal equinox:
    tmp1 = (ke/2.0 + xec/8.0)*(1.0 + xse)*dgsin(komega)
    tmp2 = xee/4.0*(0.5 + xse)*dgsin(2.0*komega)
    tmp3 = xec/8.0*(1.0/3.0 + xse)*dgsin(3.0*komega)
    xlam = tmp1 - tmp2 + tmp3
    xlam = degrees(2.0*xlam)

    ! Mean longitude for day of year:
    dlamm = xlam + (day - 80.0)*(360.0/ndayyear)

    ! Mean anomaly:
    anm = dlamm - komega
    ranm = radians(anm)

    ! True anomaly:
    ranv = (ranm + (2.0*ke - xec/4.0)*sin(ranm) + 5.0/4.0*xee*sin(2.0*ranm) + 13.0/12.0*xec*sin(3.0*ranm))  ! xxx use sin with single-precision compilation
    anv = degrees(ranv)

    ! True longitude:
    out_berger%lambda = anv + komega
    if (out_berger%lambda < 0.0) then
      out_berger%lambda = out_berger%lambda + 360.0
    else if (out_berger%lambda > 360.0) then
      out_berger%lambda = out_berger%lambda - 360.0
    endif

    ! True anomaly:
    out_berger%nu = (out_berger%lambda - komega)
    if (out_berger%nu < 0.0) then
      out_berger%nu = out_berger%nu + 360.0
    endif

  end function get_berger_tls


  function calc_sat_slope( tc ) result( sat_slope )
    !----------------------------------------------------------------   
    ! Calculates the slope of the sat pressure temp curve, Pa/K
    ! Ref:      Eq. 13, Allen et al. (1998)
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C

    ! function return value
    real :: sat_slope  ! slope of the sat pressure temp curve, Pa/K

    sat_slope = (17.269)*(237.3)*(610.78)*(exp(tc*17.269/(tc + 237.3))/((tc + 237.3)**2))

  end function calc_sat_slope


  function calc_enthalpy_vap( tc ) result( enthalpy_vap )
    !----------------------------------------------------------------   
    ! Calculates the enthalpy of vaporization, J/kg
    ! Ref:      Eq. 8, Henderson-Sellers (1984)
    !----------------------------------------------------------------   

    ! arguments
    real, intent(in) :: tc ! air temperature, degrees C

    ! function return value
    real ::  enthalpy_vap ! enthalpy of vaporization, J/kg

    enthalpy_vap = 1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))**2

  end function calc_enthalpy_vap


  function density_h2o( tc, press )
    !----------------------------------------------------------------   
    ! Calculates density of water at a given temperature and pressure
    ! Ref: Chen et al. (1977)
    !----------------------------------------------------------------   

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
    ! arguments
    real, intent(in) :: tc     ! air temperature, degrees C
    real, intent(in) :: press  ! atmospheric pressure, Pa

    ! local variables
    real :: lv  ! latent heat of vaporization (J/kg)
    real :: cp

    ! function return value
    real :: psychro  ! psychrometric constant, Pa/K

    ! local variables
    real :: my_tc    ! adjusted temperature to avoid numerical blow-up 

    ! Adopted temperature adjustment from SPLASH, Python version
    my_tc = tc
    if (my_tc < 0) then
      my_tc = 0.0
    else if (my_tc > 100) then
      my_tc = 100.0
    end if

    ! Calculate the specific heat capacity of water, J/kg/K
    ! Eq. 47, Tsilingiris (2008)
    cp = 1.0e3*(&
               1.0045714270&
             + 2.050632750e-3  *my_tc&
             - 1.631537093e-4  *my_tc*my_tc&
             + 6.212300300e-6  *my_tc*my_tc*my_tc&
             - 8.830478888e-8  *my_tc*my_tc*my_tc*my_tc&
             + 5.071307038e-10 *my_tc*my_tc*my_tc*my_tc*my_tc&
            )

    ! Calculate latent heat of vaporization, J/kg
    lv = calc_enthalpy_vap(tc)

    ! Calculate psychrometric constant, Pa/K
    ! Eq. 8, Allen et al. (1998)
    psychro = cp * kMa * press / (kMv * lv)

  end function psychro


  subroutine initio_nc_waterbal()
    !////////////////////////////////////////////////////////////////
    ! Opens NetCDF output files.
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: init_nc_3D_time, check
    use md_interface, only: interface

    ! local variables
    character(len=256) :: prefix

    character(len=*), parameter :: TITLE = "SOFUN GP-model output, module md_waterbal (SPLASH)"
    character(len=4) :: year_char

    integer :: jpngr, doy

    write(year_char,999) interface%steering%outyear

    prefix = "./output_nc/"//trim(interface%params_siml%runname)

    if ( .not. interface%steering%spinup ) then
      !----------------------------------------------------------------
      ! Annual NetCDF output
      !----------------------------------------------------------------
      if (interface%params_siml%loutwaterbal) then
        !----------------------------------------------------------------
        ! Annual PET output file 
        !----------------------------------------------------------------
        ncoutfilnam_apet = trim(prefix)//'.'//year_char//".a.pet.nc"
        print*,'initialising ', trim(ncoutfilnam_apet), '...'
        call init_nc_3D_time( filnam  = trim(ncoutfilnam_apet), &
                        nlon     = interface%domaininfo%nlon, &
                        nlat     = interface%domaininfo%nlat, &
                        lon      = interface%domaininfo%lon, &
                        lat      = interface%domaininfo%lat, &
                        outyear  = interface%steering%outyear, &
                        outdt    = 365, &
                        outnt    = 1, &
                        varnam   = PET_NAME, &
                        varunits = "mm yr-1", &
                        longnam  = "potential evapotranspiration", &
                        title    = TITLE, &
                        globatt2_nam = "in_ppfd",   globatt2_val = trim(in_ppfd)   &
                        )

        !----------------------------------------------------------------
        ! Annual AET output file 
        !----------------------------------------------------------------
        ncoutfilnam_aaet = trim(prefix)//'.'//year_char//".a.aet.nc"
        print*,'initialising ', trim(ncoutfilnam_aaet), '...'
        call init_nc_3D_time( filnam  = trim(ncoutfilnam_aaet), &
                        nlon     = interface%domaininfo%nlon, &
                        nlat     = interface%domaininfo%nlat, &
                        lon      = interface%domaininfo%lon, &
                        lat      = interface%domaininfo%lat, &
                        outyear  = interface%steering%outyear, &
                        outdt    = 365, &
                        outnt    = 1, &
                        varnam   = AET_NAME, &
                        varunits = "mm yr-1", &
                        longnam  = "actual evapotranspiration", &
                        title    = TITLE, &
                        globatt2_nam = "in_ppfd",   globatt2_val = trim(in_ppfd)   &
                        )

        !----------------------------------------------------------------
        ! Annual ALPHA (AET/PET) output file 
        !----------------------------------------------------------------
        ncoutfilnam_aalpha = trim(prefix)//'.'//year_char//".a.alpha.nc"
        print*,'initialising ', trim(ncoutfilnam_aalpha), '...'
        call init_nc_3D_time( filnam  = trim(ncoutfilnam_aalpha), &
                        nlon     = interface%domaininfo%nlon, &
                        nlat     = interface%domaininfo%nlat, &
                        lon      = interface%domaininfo%lon, &
                        lat      = interface%domaininfo%lat, &
                        outyear  = interface%steering%outyear, &
                        outdt    = 365, &
                        outnt    = 1, &
                        varnam   = ALPHA_NAME, &
                        varunits = "", &
                        longnam  = "AET/PET, mean of daily values", &
                        title    = TITLE, &
                        globatt2_nam = "in_ppfd",   globatt2_val = trim(in_ppfd)   &
                        )
      end if


      if (       interface%steering%outyear>=interface%params_siml%daily_out_startyr &
           .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then      
        !----------------------------------------------------------------
        ! Daily NetCDF output
        !----------------------------------------------------------------

        !----------------------------------------------------------------
        ! Daily WCONT output file 
        !----------------------------------------------------------------
        if (interface%params_siml%loutdwcont) then
          ncoutfilnam_dwcont = trim(prefix)//'.'//year_char//".d.wcont.nc"
          print*,'initialising ', trim(ncoutfilnam_dwcont), '...'
          call init_nc_3D_time( filnam  = trim(ncoutfilnam_dwcont), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = interface%params_siml%outdt, &
                          outnt    = interface%params_siml%outnt, &
                          varnam   = WCONT_NAME, &
                          varunits = "mm", &
                          longnam  = "soil water content", &
                          title    = TITLE, &
                          globatt2_nam = "in_ppfd",   globatt2_val = trim(in_ppfd)   &
                          )
        end if

        !----------------------------------------------------------------
        ! Daily water balance
        !----------------------------------------------------------------
        if (interface%params_siml%loutdwbal) then
          ncoutfilnam_dwbal = trim(prefix)//'.'//year_char//".d.wbal.nc"
          print*,'initialising ', trim(ncoutfilnam_dwbal), '...'
          call init_nc_3D_time( filnam  = trim(ncoutfilnam_dwbal), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = interface%params_siml%outdt, &
                          outnt    = interface%params_siml%outnt, &
                          varnam   = WBAL_NAME, &
                          varunits = "mm d-1", &
                          longnam  = "water balance as precipitation and snow melt minus runoff and evapotranspiration", &
                          title    = TITLE, &
                          globatt2_nam = "in_ppfd",   globatt2_val = trim(in_ppfd)   &
                          )
        end if

        !----------------------------------------------------------------
        ! Daily PET output file 
        !----------------------------------------------------------------
        if (interface%params_siml%loutdpet) then
          ncoutfilnam_dpet = trim(prefix)//'.'//year_char//".d.pet.nc"
          print*,'initialising ', trim(ncoutfilnam_dpet), '...'
          call init_nc_3D_time( filnam  = trim(ncoutfilnam_dpet), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = interface%params_siml%outdt, &
                          outnt    = interface%params_siml%outnt, &
                          varnam   = PET_NAME, &
                          varunits = "mm d-1", &
                          longnam  = "potential evapotranspiration", &
                          title    = TITLE, &
                          globatt2_nam = "in_ppfd",   globatt2_val = trim(in_ppfd)   &
                          )
        end if

        !----------------------------------------------------------------
        ! Daily net radiation output file 
        !----------------------------------------------------------------
        if (interface%params_siml%loutdnetrad) then
          ncoutfilnam_drn = trim(prefix)//'.'//year_char//".d.netrad.nc"
          print*,'initialising ', trim(ncoutfilnam_drn), '...'
          call init_nc_3D_time( filnam  = trim(ncoutfilnam_drn), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = interface%params_siml%outdt, &
                          outnt    = interface%params_siml%outnt, &
                          varnam   = RN_NAME, &
                          varunits = "J m-2 d-1", &
                          longnam  = "net radiation", &
                          title    = TITLE, &
                          globatt2_nam = "in_ppfd",   globatt2_val = trim(in_ppfd)   &
                          )
        end if        

        !----------------------------------------------------------------
        ! Daily AET output file 
        !----------------------------------------------------------------
        if (interface%params_siml%loutdaet) then
          ncoutfilnam_daet = trim(prefix)//'.'//year_char//".d.aet.nc"
          print*,'initialising ', trim(ncoutfilnam_daet), '...'
          call init_nc_3D_time( filnam  = trim(ncoutfilnam_daet), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = interface%params_siml%outdt, &
                          outnt    = interface%params_siml%outnt, &
                          varnam   = AET_NAME, &
                          varunits = "mm d-1", &
                          longnam  = "actual evapotranspiration", &
                          title    = TITLE, &
                          globatt2_nam = "in_ppfd",   globatt2_val = trim(in_ppfd)   &
                          )        
        end if

      end if

    end if

    888  format (F12.6)
    999  format (I4.4)
    
  end subroutine initio_nc_waterbal


  subroutine initoutput_waterbal( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises waterbalance-specific output variables
    ! The same subroutine is used here for initialising rolling mean variables
    !----------------------------------------------------------------
    use md_interface, only: interface

    ! arguments
    integer, intent(in) :: ngridcells

    ! Annual output variables
    if (interface%params_siml%loutwaterbal) then

      if (interface%steering%init) then
        allocate( outapet(ngridcells) )
        allocate( outaaet(nlu,ngridcells) )
        allocate( outaalpha(nlu,ngridcells) )
      end if
      outapet(:)     = 0.0
      outaaet(:,:)   = 0.0
      outaalpha(:,:) = 0.0
    end if

    ! Daily output variables
    if (interface%steering%init) then
      if (interface%params_siml%loutdpet)   allocate( outdpet(   interface%params_siml%outnt,ngridcells)     )     ! daily potential ET, mm
      if (interface%params_siml%loutdaet)   allocate( outdaet(   nlu,interface%params_siml%outnt,ngridcells) )     ! daily actual ET, mm
      if (interface%params_siml%loutdwcont) allocate( outdwcont( nlu,interface%params_siml%outnt,ngridcells) )     ! daily soil moisture, mm
      if (interface%params_siml%loutdwbal)  allocate( outdwbal(  nlu,interface%params_siml%outnt,ngridcells) )     ! daily Cramer-Prentice-Alpha, (unitless)
      if (interface%params_siml%loutdnetrad)allocate( outdrn(    nlu,interface%params_siml%outnt,ngridcells) )     ! daily net radiation J m-2 d-1
    end if
    if (interface%params_siml%loutdwcont)  outdwcont(:,:,:)  = 0.0
    if (interface%params_siml%loutdpet)    outdpet(:,:)      = 0.0
    if (interface%params_siml%loutdaet)    outdaet(:,:,:)    = 0.0
    if (interface%params_siml%loutdwbal)   outdwbal(:,:,:)   = 0.0
    if (interface%params_siml%loutdnetrad) outdrn(:,:,:)     = 0.0

  end subroutine initoutput_waterbal


  ! subroutine init_rlm_waterbal( ngridcells )
  !   !////////////////////////////////////////////////////////////////
  !   ! Initialises waterbalance-specific output variables
  !   ! The same subroutine is used here for initialising rolling mean variables
  !   !----------------------------------------------------------------
  !   use md_interface, only: interface

  !   ! arguments
  !   integer, intent(in) :: ngridcells

  !   ! Rolling mean variables
  !   if (interface%steering%init) then
  !     if (.not.allocated(rlmalpha)) allocate( rlmalpha(nlu,ngridcells) )
  !   end if
  !   rlmalpha(:,:) = 0.0

  ! end subroutine init_rlm_waterbal


  subroutine getout_daily_waterbal( tile, tile_fluxes, jpngr, doy )
    !////////////////////////////////////////////////////////////////
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    use md_interface, only: interface
    use md_tile, only: tile_type, tile_fluxes_type

    ! argument
    type(tile_type), dimension(nlu), intent(in)        :: tile
    type(tile_fluxes_type), dimension(nlu), intent(in) :: tile_fluxes
    integer, intent(in)                                :: jpngr
    integer, intent(in)                                :: doy    

    ! local variables
    integer :: it

    it = floor( real( doy - 1 ) / real( interface%params_siml%outdt ) ) + 1

    ! Annual output variables
    if (interface%params_siml%loutwaterbal) then
      if (outenergy) then
        outapet(jpngr)    = outapet(jpngr)   + tile_fluxes(1)%canopy%dpet_e
        outaaet(:,jpngr)  = outaaet(:,jpngr) + tile_fluxes(:)%canopy%daet_e
      else 
        outapet(jpngr)    = outapet(jpngr)   + tile_fluxes(1)%canopy%dpet
        outaaet(:,jpngr)  = outaaet(:,jpngr) + tile_fluxes(:)%canopy%daet
      end if
      if (tile_fluxes(1)%canopy%dpet > 0.0) then
        outaalpha(:,jpngr)  = outaalpha(:,jpngr) + (tile_fluxes(:)%canopy%daet / tile_fluxes(:)%canopy%dpet) / ndayyear
      else
        outaalpha(:,jpngr)  = outaalpha(:,jpngr) + 1.0 / ndayyear
      end if
    end if

    ! Daily output variables
    if (interface%params_siml%loutdwcont)  outdwcont(:,it,jpngr)  = outdwcont(:,it,jpngr) + tile(:)%soil%phy%wcont / real( interface%params_siml%outdt )
    if (interface%params_siml%loutdwbal)   outdwbal(:,it,jpngr)   = outdwbal(:,it,jpngr)  + tile_fluxes(:)%canopy%dwbal / real( interface%params_siml%outdt )
    if (interface%params_siml%loutdnetrad) outdrn(:,it,jpngr)     = outdrn(:,it,jpngr)    + (tile_fluxes(:)%canopy%drn + tile_fluxes(:)%canopy%drnn) / real( interface%params_siml%outdt ) ! daytime plus nighttime

    if (outenergy) then
      if (interface%params_siml%loutdpet) outdpet(it,jpngr)    = outdpet(it,jpngr)   + tile_fluxes(1)%canopy%dpet_e / real( interface%params_siml%outdt )
      if (interface%params_siml%loutdaet) outdaet(:,it,jpngr)  = outdaet(:,it,jpngr) + tile_fluxes(:)%canopy%daet_e / real( interface%params_siml%outdt )
    else 
      if (interface%params_siml%loutdpet) outdpet(it,jpngr)    = outdpet(it,jpngr)   + tile_fluxes(1)%canopy%dpet / real( interface%params_siml%outdt )
      if (interface%params_siml%loutdaet) outdaet(:,it,jpngr)  = outdaet(:,it,jpngr) + tile_fluxes(:)%canopy%daet / real( interface%params_siml%outdt )
    end if
    
    ! outdecon(it,jpngr)     = outdecon(it,jpngr)    + tile_fluxes(lu)%canopy%econ * 1.0e12 / real( interface%params_siml%outdt ) ! converting from m J-1 to mm GJ-1 = m TJ-1

  end subroutine getout_daily_waterbal


  subroutine getrlm_daily_waterbal( tile_fluxes, jpngr, doy )
    !////////////////////////////////////////////////////////////////
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    use md_interface, only: interface
    use md_tile, only: tile_fluxes_type

    ! argument
    type(tile_fluxes_type), dimension(nlu), intent(in) :: tile_fluxes
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy    

    if (tile_fluxes(1)%canopy%dpet > 0.0) then
      aalpha(:) = aalpha(:) + (tile_fluxes(:)%canopy%daet / tile_fluxes(:)%canopy%dpet) / ndayyear
    else
      aalpha(:) = aalpha(:) + 1.0 / ndayyear
    end if

  end subroutine getrlm_daily_waterbal


  subroutine get_rlm_waterbal( tile, init )
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates the rolling mean of relevant variables
    ! This requires the full arrays (all gridcells) to be stored.
    !-------------------------------------------------------------------------
    use md_params_core, only: nlu
    use md_tile, only: tile_type

    ! arguments
    type(tile_type), dimension(nlu), intent(inout) :: tile
    logical :: init

    ! local variables
    integer, save :: ncalls
    integer :: nyrs_uptonow
    integer :: lu

    if (init) ncalls = 0
    ncalls = ncalls + 1
    nyrs_uptonow = min( ncalls, nyrs_rlmalpha )

    tile(:)%rlmalpha = ( tile(:)%rlmalpha * (nyrs_uptonow - 1) + aalpha(:) ) / nyrs_uptonow

  end subroutine get_rlm_waterbal


  ! subroutine writeout_ascii_waterbal()
  !   !/////////////////////////////////////////////////////////////////////////
  !   ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
  !   !-------------------------------------------------------------------------
  !   use md_params_core, only: ndayyear, nmonth
  !   use md_interface, only: interface

  !   ! Local variables
  !   real :: itime
  !   integer :: it, jpngr
    
  !   ! xxx implement this: sum over gridcells? single output per gridcell?
  !   if (maxgrid>1) stop 'writeout_ascii_waterbal: think of something ...'
  !   jpngr = 20000

  !   !-------------------------------------------------------------------------
  !   ! DAILY OUTPUT
  !   !-------------------------------------------------------------------------
  !   if (interface%params_siml%loutwaterbal) then

  !     if ( .not. interface%steering%spinup &
  !          .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
  !          .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

  !       ! Write daily output only during transient simulation
  !       do it=1,interface%params_siml%outnt

  !         ! Define 'itime' as a decimal number corresponding to day in the year + year
  !         itime = real(interface%steering%outyear) + real( it - 1 ) * interface%params_siml%outdt / real( ndayyear )

  !         if (nlu>1) stop 'writeout_ascii_waterbal: write out lu-area weighted sum'

  !         ! xxx lu-area weighted sum if nlu>0
  !         write(255,999) itime, outdwcont(1,it,jpngr)
  !         write(259,999) itime, outdpet(it,jpngr)
  !         write(260,999) itime, outdaet(1,it,jpngr)
  !         ! write(253,999) itime, outdwbal(it,jpngr)
  !         ! write(251,999) itime, outdra(it,jpngr)
  !         ! write(252,999) itime, outdrn(it,jpngr)
  !         ! write(254,999) itime, outdcn(it,jpngr)
  !         ! write(257,999) itime, outdro(1,it,jpngr)
  !         ! write(263,999) itime, outdfleach(1,it,jpngr)
  !         ! write(258,999) itime, outdeet(it,jpngr)
  !         ! write(261,999) itime, outdayl(it,jpngr)
  !         ! write(262,999) itime, outdwbal(1,it,jpngr)
  !         ! write(264,999) itime, outdecon(it,jpngr)

  !       end do
  !     end if
  !   end if

  !   return
    
  !   888 format (F20.8,E20.8)
  !   999 format (F20.8,F20.8)

  ! end subroutine writeout_ascii_waterbal


  subroutine writeout_nc_waterbal()
    !/////////////////////////////////////////////////////////////////////////
    ! Write NetCDF output
    !-------------------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: write_nc_2D, write_nc_3D_time, check
    use md_interface, only: interface

    if (nlu>1) stop 'writeout_nc_waterbal(): nlu > 1. Think of something...'

    if ( .not. interface%steering%spinup ) then
      !-------------------------------------------------------------------------
      ! Annual output
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutwaterbal) then
        !-------------------------------------------------------------------------
        ! PET
        !-------------------------------------------------------------------------
        print*,'writing ', trim(ncoutfilnam_apet), '...'
        call write_nc_2D( trim(ncoutfilnam_apet), &
                          PET_NAME, &
                          interface%domaininfo%maxgrid, &
                          interface%domaininfo%nlon, &
                          interface%domaininfo%nlat, &
                          interface%grid(:)%ilon, &
                          interface%grid(:)%ilat, &
                          interface%grid(:)%dogridcell, &
                          outapet(:) &
                          )

        !-------------------------------------------------------------------------
        ! AET
        !-------------------------------------------------------------------------
        if (nlu>1) stop 'writeout_nc_waterbal: nlu>1. Think of something clever!'
        print*,'writing ', trim(ncoutfilnam_aaet), '...'
        call write_nc_2D( trim(ncoutfilnam_aaet), &
                          AET_NAME, &
                          interface%domaininfo%maxgrid, &
                          interface%domaininfo%nlon, &
                          interface%domaininfo%nlat, &
                          interface%grid(:)%ilon, &
                          interface%grid(:)%ilat, &
                          interface%grid(:)%dogridcell, &
                          outaaet(1,:) &
                          )

        !-------------------------------------------------------------------------
        ! ALPHA (AET/PET)
        !-------------------------------------------------------------------------
        print*,'writing ', trim(ncoutfilnam_aalpha), '...'
        call write_nc_2D( trim(ncoutfilnam_aalpha), &
                          ALPHA_NAME, &
                          interface%domaininfo%maxgrid, &
                          interface%domaininfo%nlon, &
                          interface%domaininfo%nlat, &
                          interface%grid(:)%ilon, &
                          interface%grid(:)%ilat, &
                          interface%grid(:)%dogridcell, &
                          outaalpha(1,:) &
                          )        

      end if

      if (       interface%steering%outyear>=interface%params_siml%daily_out_startyr &
           .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then
        !-------------------------------------------------------------------------
        ! Daily output
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        ! soil water content
        !-------------------------------------------------------------------------
        if (interface%params_siml%loutdwcont) then
          print*,'writing ', trim(ncoutfilnam_dwcont), '...'
          call write_nc_3D_time( trim(ncoutfilnam_dwcont), &
                            WCONT_NAME, &
                            interface%domaininfo%maxgrid, &
                            interface%domaininfo%nlon, &
                            interface%domaininfo%nlat, &
                            interface%grid(:)%ilon, &
                            interface%grid(:)%ilat, &
                            interface%params_siml%outnt, &
                            interface%grid(:)%dogridcell, &
                            outdwcont(1,:,:) &
                            )
        end if

        !-------------------------------------------------------------------------
        ! Daily water balance
        !-------------------------------------------------------------------------
        if (interface%params_siml%loutdwbal) then
          print*,'writing ', trim(ncoutfilnam_dwbal), '...'
          call write_nc_3D_time( trim(ncoutfilnam_dwbal), &
                            WBAL_NAME, &
                            interface%domaininfo%maxgrid, &
                            interface%domaininfo%nlon, &
                            interface%domaininfo%nlat, &
                            interface%grid(:)%ilon, &
                            interface%grid(:)%ilat, &
                            interface%params_siml%outnt, &
                            interface%grid(:)%dogridcell, &
                            outdwbal(1,:,:) &
                            )
        end if

        !-------------------------------------------------------------------------
        ! PET
        !-------------------------------------------------------------------------
        if (interface%params_siml%loutdpet) then
          print*,'writing ', trim(ncoutfilnam_dpet), '...'
          call write_nc_3D_time( trim(ncoutfilnam_dpet), &
                            PET_NAME, &
                            interface%domaininfo%maxgrid, &
                            interface%domaininfo%nlon, &
                            interface%domaininfo%nlat, &
                            interface%grid(:)%ilon, &
                            interface%grid(:)%ilat, &
                            interface%params_siml%outnt, &
                            interface%grid(:)%dogridcell, &
                            outdpet(:,:) &
                            )
        end if

        !-------------------------------------------------------------------------
        ! AET
        !-------------------------------------------------------------------------
        if (interface%params_siml%loutdaet) then
          if (nlu>1) stop 'writeout_nc_waterbal: nlu>1. Think of something clever!'
          print*,'writing ', trim(ncoutfilnam_daet), '...'
          call write_nc_3D_time( trim(ncoutfilnam_daet), &
                            AET_NAME, &
                            interface%domaininfo%maxgrid, &
                            interface%domaininfo%nlon, &
                            interface%domaininfo%nlat, &
                            interface%grid(:)%ilon, &
                            interface%grid(:)%ilat, &
                            interface%params_siml%outnt, &
                            interface%grid(:)%dogridcell, &
                            outdaet(1,:,:) &
                            )
        end if

        !-------------------------------------------------------------------------
        ! Net radiation
        !-------------------------------------------------------------------------
        if (interface%params_siml%loutdnetrad) then
          if (nlu>1) stop 'writeout_nc_waterbal: nlu>1. Think of something clever!'
          print*,'writing ', trim(ncoutfilnam_drn), '...'
          call write_nc_3D_time( trim(ncoutfilnam_drn), &
                            RN_NAME, &
                            interface%domaininfo%maxgrid, &
                            interface%domaininfo%nlon, &
                            interface%domaininfo%nlat, &
                            interface%grid(:)%ilon, &
                            interface%grid(:)%ilat, &
                            interface%params_siml%outnt, &
                            interface%grid(:)%dogridcell, &
                            outdrn(1,:,:) &
                            )
        end if

      end if

    end if

  end subroutine writeout_nc_waterbal

end module md_waterbal
