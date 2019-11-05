module md_gpp
  !////////////////////////////////////////////////////////////////
  ! P-MODEL MODULE
  ! 
  ! Contains P-model functions for calculating ci, vcmax, and light
  ! use efficiency (LUE) as a function of ambient conditions (temperature,
  ! vapour pressure difference, CO2, and elevation), following 
  ! Prentice et al., 2014, and Wang et al., 2017. This is implemented
  ! by function pmodel().
  ! 
  ! Outputs of pmodel() can used within a LUE model to calculate gross 
  ! primary production by multiplying with aborbed light. This is 
  ! implemented by function calc_dgpp().
  !
  ! High-level wrapper functions and subroutines (getlue(), and 
  ! gpp()) are provided for use within SOFN.
  !
  ! An empirical soil moisture stress function is implemented 
  ! following Stocker et al., 2018.
  !
  ! Note that getpar_modl_gpp() must be invoked prior to first
  ! call of any other function/subroutine. An example is given 
  ! by demo_pmodel.f90.
  !
  ! Copyright (C) 2015, see LICENSE, Benjamin Stocker
  ! Written by Benjamin Stocker, partly based on Python code by
  ! Tyler Davis.
  !----------------------------------------------------------------
  ! load core parameters
  use md_params_core, only: nmonth, npft, nlu, c_molmass, h2o_molmass, maxgrid, ndayyear, dummy
  use md_sofunutils, only: calc_patm

  implicit none

  private
  public params_pft_gpp, getpar_modl_gpp, initoutput_gpp, &
    gpp, getlue, pmodel, getout_daily_gpp, getout_annual_gpp, &
    outtype_pmodel, calc_ftemp_kphio, calc_dgpp, calc_drd, &
    initio_nc_gpp, writeout_nc_gpp
    
  !----------------------------------------------------------------
  ! Module-specific state variables
  !----------------------------------------------------------------
  real, dimension(npft) :: dassim           ! daily leaf-level assimilation rate (per unit leaf area) [gC/m2/d]
  real, dimension(npft) :: dgs              ! stomatal conductance (per unit leaf area, average daily) [mol H2O m-2 s-1]  
  real, dimension(npft) :: dvcmax_canop     ! canopy-level Vcmax [gCO2/m2-ground/s]
  real, dimension(npft) :: dvcmax_leaf      ! leaf-level Vcmax [gCO2/m2-leaf/s]

  !-----------------------------------------------------------------------
  ! Known parameters, therefore hard-wired.
  !-----------------------------------------------------------------------
  real, parameter :: kPo = 101325.0        ! standard atmosphere, Pa (Allen, 1973)
  real, parameter :: kTo = 25.0            ! base temperature, deg C (Prentice, unpublished)
  real, parameter :: temp0 = 0.0           ! temperature below which all quantities are zero (deg C)

  !-----------------------------------------------------------------------
  ! Metabolic N ratio (N per unit Vcmax)
  ! Reference: Harrison et al., 2009, Plant, Cell and Environment; Eq. 3
  !-----------------------------------------------------------------------
  real, parameter :: mol_weight_rubisco    = 5.5e5    ! molecular weight of Rubisco, (g R)(mol R)-1
  real, parameter :: n_conc_rubisco        = 1.14e-2  ! N concentration in rubisco, (mol N)(g R)-1
  real, parameter :: cat_turnover_per_site = 2.33     ! catalytic turnover rate per site at 25 deg C, (mol CO2)(mol R sites)-1; use 2.33 instead of (3.5) as not all Rubisco is active (see Harrison et al., 2009)  
  real, parameter :: cat_sites_per_mol_R   = 8.0      ! number of catalytic sites per mol R, (mol R sites)(mol R)-1

  ! Metabolic N ratio (= 336.3734 mol N s (mol CO2)-1 )
  real, parameter :: n_v = mol_weight_rubisco * n_conc_rubisco / ( cat_turnover_per_site * cat_sites_per_mol_R )

  !-----------------------------------------------------------------------
  ! Uncertain (unknown) parameters. Runtime read-in
  !-----------------------------------------------------------------------
  type paramstype_gpp
    real :: beta         ! Unit cost of carboxylation (dimensionless)
    real :: soilm_par_a
    real :: soilm_par_b
    real :: rd_to_vcmax  ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  end type paramstype_gpp

  type(paramstype_gpp) :: params_gpp

  ! PFT-DEPENDENT PARAMETERS
  type pftparamstype_gpp
    real :: kphio        ! quantum efficiency (Long et al., 1993)  
  end type pftparamstype_gpp

  type(pftparamstype_gpp), dimension(npft) :: params_pft_gpp

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC, PRIVATE VARIABLES
  !----------------------------------------------------------------
  ! Function return variables as derived types
  type outtype_pmodel
    real :: gammastar           ! temperature-dependent photorespiratory compensation point (Pa)
    real :: kmm                 ! Michaelis-Menten coefficient (Pa)
    real :: ci                  ! leaf-internal partial pressure, (Pa)
    real :: chi                 ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
    real :: iwue                ! intrinsic water use efficiency = A / gs = ca - ci = ca ( 1 - chi ) , unitless
    real :: lue                 ! light use efficiency (mol CO2 / mol photon)
    real :: assim               ! leaf-level assimilation rate (mol CO2 m-2 s-1)
    real :: gs                  ! stomatal conductance to CO2, expressed per units absorbed light (mol H2O m-2 m-1 / (mol light m-2))
    real :: gpp                 ! gross primary productivity (g CO2 m-2 d-1)
    real :: vcmax               ! canopy-level maximum carboxylation capacity per unit ground area (mol CO2 m-2 s-1)
    real :: vcmax25             ! canopy-level Vcmax25 (Vcmax normalized to 25 deg C) (mol CO2 m-2 s-1)
    real :: vcmax_unitfapar     ! Vcmax per unit fAPAR (mol CO2 m-2 s-1)
    real :: vcmax_unitiabs      ! Vcmax per unit absorbed light (mol CO2 m-2 s-1 mol-1)
    real :: ftemp_inst_vcmax    ! Instantaneous temperature response factor of Vcmax (unitless)
    real :: ftemp_inst_rd       ! Instantaneous temperature response factor of Rd (unitless)
    real :: rd                  ! Dark respiration (mol CO2 m-2 s-1)
    real :: rd_unitfapar        ! Dark respiration per unit fAPAR (mol CO2 m-2 s-1)
    real :: rd_unitiabs         ! Dark respiration per unit absorbed light (mol CO2 m-2 s-1)
    real :: actnv               ! Canopy-level total metabolic leaf N per unit ground area (g N m-2)
    real :: actnv_unitfapar     ! Metabolic leaf N per unit fAPAR (g N m-2)
    real :: actnv_unitiabs      ! Metabolic leaf N per unit absorbed light (g N m-2 mol-1)
    real :: transp              ! Canopy-level total transpiration rate (g H2O (mol photons)-1)
  end type outtype_pmodel

  type outtype_chi
    real :: chi                 ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
    real :: mj
    real :: mc
    real :: mjoc
  end type outtype_chi

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, allocatable, dimension(:,:) :: outdgpp    ! daily gross primary production [gC/m2/d]
  real, allocatable, dimension(:,:) :: outdrd     ! daily dark respiration [gC/m2/d]
  real, allocatable, dimension(:,:) :: outdtransp ! daily transpiration [mm]

  ! annual
  real, dimension(:,:), allocatable :: outagpp
  real, dimension(:,:), allocatable :: outavcmax        ! canopy-level caboxylation capacity at annual maximum [mol CO2 m-2 s-1]
  real, dimension(:,:), allocatable :: outavcmax_25     ! canopy-level normalised caboxylation capacity at annual maximum [mol CO2 m-2 s-1]
  real, dimension(:,:), allocatable :: outavcmax_leaf   ! leaf-level maximum caboxylation capacity, annual mean of daily values, weighted by daily assimilation rate [mol CO2 m-2 s-1]
  real, dimension(:,:), allocatable :: outalue          ! light use efficiency, mean across growing season, weighted by daily GPP
  real, dimension(:,:), allocatable :: outachi          ! ratio leaf-internal to ambient CO2 partial pressure, mean across growing season, weighted by daily GPP
  real, dimension(:,:), allocatable :: outaci           ! leaf-internal CO2 partial pressure, mean across growing season, weighted by daily GPP (ppm)
  real, dimension(:,:), allocatable :: outags           ! stomatal conductance to H2O, mean across growing season, weighted by daily GPP (mol H2O m-2 s-1)
  real, dimension(:,:), allocatable :: outaiwue         ! intrinsic water use efficiency, weighted by daily GPP [micro-mol CO2 / mol H2O]

  ! These are stored as dayly variables for annual output
  ! at day of year when LAI is at its maximum.
  real, dimension(npft,ndayyear) :: outdvcmax
  real, dimension(npft,ndayyear) :: outdvcmax25

  !----------------------------------------------------------------
  ! Module-specific NetCDF output file and variable names
  !----------------------------------------------------------------
  character(len=256) :: ncoutfilnam_dgpp
  character(len=256) :: ncoutfilnam_agpp

  character(len=*), parameter :: GPP_NAME="gpp"

contains

  subroutine gpp( plant, plant_fluxes, out_pmodel, dppfd, dayl, meanmppfd, wscal, rlmalpha, doy, moy, dtemp, do_soilmstress, do_tempstress, fapar )
    !//////////////////////////////////////////////////////////////////
    ! Wrapper subroutine for invoking function calls to calculate daily
    ! rates (gross primary productivity, dark respiration), and other
    ! corollary predictions of daily means (canopy and leaf-level Vcmax 
    ! and assimilation rate, and leaf-level stomatal conductance).
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !------------------------------------------------------------------
    use md_plant, only: params_pft_plant, plant_type, plant_fluxes_type

    ! arguments
    type(plant_type), dimension(npft), intent(in) :: plant  ! derived type containing plant-related state variables (memory from year to year)
    type(plant_fluxes_type), dimension(npft), intent(inout) :: plant_fluxes  ! derived type containing plant-related flux variables  (overwritten at each time step)
    type(outtype_pmodel), dimension(npft), intent(in)  :: out_pmodel  ! derived type containing outputs of function pmodel()
    real, intent(in) :: dppfd            ! daily total photon flux density, mol m-2
    real, intent(in) :: dayl             ! day length (h)
    real, intent(in) :: meanmppfd        ! monthly mean PPFD (mol m-2 s-1)
    real, dimension(nlu), intent(in) :: wscal            ! relative soil water content (unitless)
    real, dimension(nlu), intent(in) :: rlmalpha         ! rolling mean alpha (mean annual AET/PET, unitless)
    integer, intent(in) :: doy             ! day of year and month of year
    integer, intent(in) :: moy             ! month of year and month of year
    real,    intent(in) :: dtemp           ! this day's air temperature (deg C)
    logical, intent(in) :: do_soilmstress  ! whether empirical soil miosture stress function is applied to GPP
    logical, intent(in) :: do_tempstress   ! whether empirical temperature stress function is applied to GPP
    real, intent(in)    :: fapar  ! fraction of absorbed photosynthetically active radiation (unitless)

    ! local variables
    integer :: pft
    integer :: lu
    real    :: soilmstress
    real    :: ftemp_kphio

    !----------------------------------------------------------------
    ! CALCULATE PREDICTED GPP FROM P-model
    ! using instantaneous (daily) LAI, PPFD, Cramer-Prentice-alpha
    !----------------------------------------------------------------
    do pft=1,npft

      ! land use category (gridcell tile)
      lu = params_pft_plant(pft)%lu_category

      ! Calculate soil moisture stress as a function of soil moisture, mean alpha and vegetation type (grass or not)
      if (do_soilmstress) then
        soilmstress = calc_soilmstress( wscal(lu), rlmalpha(lu), params_pft_plant(pft)%grass )
      else
        soilmstress = 1.0
      end if
      ! print*,'soilmstress ', soilmstress

      ! Include instantaneous temperature effect on quantum yield efficiency
      if (do_tempstress) then
        ftemp_kphio = calc_ftemp_kphio( dtemp )
      else
        ftemp_kphio = 1.0
      end if

      if ( plant(pft)%fpc_grid > 0.0 .and. dayl > 0.0 .and. dtemp > -5.0 ) then

        ! GPP
        ! print*,'md_gpp: ppfd = ', dppfd
        plant_fluxes(pft)%dgpp = calc_dgpp( plant(pft)%fapar_ind, plant(pft)%fpc_grid, dppfd, out_pmodel(pft)%lue, ftemp_kphio, soilmstress )

        ! ! transpiration
        ! ! dtransp(pft) = calc_dtransp( plant(pft)%fapar_ind, plant(pft)%acrown, dppfd, out_pmodel(pft)%transp_unitiabs, ftemp_kphio, soilmstress )
        ! dtransp(pft) = calc_dtransp( plant(pft)%fapar_ind, plant(pft)%acrown, dppfd, out_pmodel(pft)%transp_unitiabs, dtemp )

        ! Dark respiration
        plant_fluxes(pft)%drd = calc_drd( plant(pft)%fapar_ind, plant(pft)%fpc_grid, meanmppfd, out_pmodel(pft)%rd_unitiabs, ftemp_kphio, soilmstress )

        ! Leaf-level assimilation rate
        dassim(pft) = calc_dassim( dppfd, out_pmodel(pft)%lue, dayl, ftemp_kphio, soilmstress )

        ! ! stomatal conductance
        ! dgs(pft) = calc_dgs( dppfd, out_pmodel(pft)%gs_unitiabs, dayl, ftemp_kphio, soilmstress )

        ! Canopy-level Vcmax (actually changes only monthly)
        dvcmax_canop(pft) = calc_vcmax_canop( plant(pft)%fapar_ind, out_pmodel(pft)%vcmax_unitiabs, meanmppfd )

        ! Leaf-level Vcmax
        dvcmax_leaf(pft) = out_pmodel(pft)%vcmax_unitiabs * meanmppfd

      else  

        plant_fluxes(pft)%dgpp    = 0.0
        plant_fluxes(pft)%drd     = 0.0
        plant_fluxes(pft)%dtransp = 0.0
        dvcmax_canop(pft)         = 0.0
        dvcmax_leaf(pft)          = 0.0

      end if 

      ! ! xxx debug
      ! plant_fluxes(pft)%dgpp = fapar

    end do

  end subroutine gpp


  function getlue( co2, dtemp, dvpd, elv ) result( out_pmodel )
    !//////////////////////////////////////////////////////////////////
    ! Wrapper function for invoking function calls of pmodel() to 
    ! calculate light use efficiency and other quantities normalised to
    ! absorbed light as a function of ambient conditions averaged over
    ! a given period (month used within SOFUN).
    !
    ! This is designed for use within SOFUN. For other applications, 
    ! adapt the function call to pmodel() appropriately.
    !------------------------------------------------------------------
    use md_plant, only: params_pft_plant
    use md_sofunutils, only: daily2monthly

    ! arguments
    real, intent(in)                      :: co2      ! atmospheric CO2 (ppm)
    real, dimension(ndayyear), intent(in) :: dtemp    ! daily air temperature (deg C)
    real, dimension(ndayyear), intent(in) :: dvpd     ! daily vapour pressure deficit (Pa)
    real, intent(in)                      :: elv      ! elevation above sea level (m)

    ! function return variable
    type(outtype_pmodel), dimension(npft,nmonth) :: out_pmodel ! P-model output variables for each month and PFT determined beforehand (per unit fAPAR and PPFD only)

    ! local variables
    real, dimension(ndayyear) :: mydtemp
    real, dimension(nmonth)   :: mtemp      ! monthly air temperature (deg C)
    real, dimension(nmonth)   :: mvpd       ! monthly vapour pressure deficit (Pa)
    integer                   :: moy, pft
    integer                   :: doy

    ! ! xxx test
    ! real, dimension(nmonth)   :: mppfd
    ! real :: myco2
    ! real :: myelv

    ! locally used daily temperature (may be different from globally used one)
    mydtemp(:) = dtemp(:)

    ! Get monthly averages
    mtemp(:) = daily2monthly( mydtemp(:), "mean" )
    mvpd(:)  = daily2monthly( dvpd(:), "mean" )

    ! ! xxx try out: -- THIS WORKS PERFECTLY -- 
    ! print*, 'WARNING: TEST INPUT FOR COMPARISON WITH OPTI7.R'
    ! myco2   = 376.0
    ! myelv   = 450.0
    ! mtemp = (/0.4879904, 6.1999985, 7.4999870, 9.6999003, 13.1999913, 19.6999227, 18.6000030, 18.0999577, 13.8999807, 10.7000307, 7.2999217, 4.4999644/)
    ! mvpd  = (/113.0432, 338.4469, 327.1185, 313.8799, 247.9747, 925.9489, 633.8551, 497.6772, 168.7784, 227.1889, 213.0142, 172.6035/)
    ! mppfd = (/223.8286, 315.2295, 547.4822, 807.4035, 945.9020, 1194.1227, 1040.5228, 1058.4161, 814.2580, 408.5199, 268.9183, 191.4482/)

    !------------------------------------------------------------------
    ! Run P-model for monthly averages and store monthly variables 
    ! per unit absorbed light (not corrected for soil moisture)
    !------------------------------------------------------------------
    do pft=1,npft

      do moy=1,nmonth

        out_pmodel(pft,moy) = pmodel( params_pft_gpp(pft)%kphio, fapar = dummy, ppfd = dummy, co2 = co2, tc = mtemp(moy), vpd = mvpd(moy), elv = elv, c4 = params_pft_plant(pft)%c4, method_optci = "prentice14", method_jmaxlim = "wang17" )

      end do
    end do

  end function getlue


  function calc_dgpp( fapar, fpc_grid, dppfd, lue, ftemp_kphio, soilmstress ) result( my_dgpp )
    !//////////////////////////////////////////////////////////////////
    ! Calculates daily GPP given mean daily light use efficiency following
    ! a simple light use efficie model approach.
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: fapar       ! fraction of absorbed photosynthetically active radiation (unitless)
    real, intent(in) :: fpc_grid    ! foliar projective cover, used for dividing grid cell area (unitless)
    real, intent(in) :: dppfd       ! daily total photon flux density (mol m-2)
    real, intent(in) :: lue         ! light use efficiency (g CO2 mol-1)
    real, intent(in) :: ftemp_kphio  ! air temperature (deg C)
    real, intent(in) :: soilmstress ! soil moisture stress factor (unitless)

    ! function return variable
    real :: my_dgpp                 ! Daily total gross primary productivity (gC m-2 d-1)

    ! GPP is light use efficiency multiplied by absorbed light and soil moisture stress function
    my_dgpp = fapar * fpc_grid * dppfd * soilmstress * lue * ftemp_kphio

  end function calc_dgpp


  function calc_dassim( dppfd, lue, daylength, ftemp_kphio, soilmstress ) result( my_dassim )
    !//////////////////////////////////////////////////////////////////
    ! Calculates leaf-level assimilation rate, mean over daylight hours.
    ! Not described in Stocker et al., XXX.
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: dppfd           ! daily total photon flux density (mol m-2)
    real, intent(in) :: lue             ! light use efficiency, (mol CO2 / mol photon)
    real, intent(in) :: daylength       ! day length (h)
    real, intent(in) :: ftemp_kphio      ! this day's air temperature (deg C)
    real, intent(in) :: soilmstress     ! soil moisture stress factor (unitless)

    ! function return variable
    real :: my_dassim                   ! leaf-level assimilation rate, mean over daylight hours ( mol CO2 m-2 s-1 )

    ! Leaf-level assimilation rate, average over daylight hours
    if (daylength>0.0) then
      my_dassim = dppfd * soilmstress * lue * ftemp_kphio / ( 60.0 * 60.0 * daylength )
    else
      my_dassim = 0.0
    end if

  end function calc_dassim


  function calc_dgs( dppfd, dgs_unitiabs, daylength, ftemp_kphio, soilmstress ) result( dgs )
    !//////////////////////////////////////////////////////////////////
    ! Calculates leaf-level stomatal conductance to H2O, mean over daylight hours
    ! Not described in Stocker et al., XXX.
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: dppfd           ! daily total photon flux density, mol m-2
    real, intent(in) :: dgs_unitiabs    ! stomatal conductance per unit absorbed light (mol H2O m-2 s-1 / mol light)
    real, intent(in) :: daylength       ! day length (h)
    real, intent(in) :: ftemp_kphio      ! this day's air temperature, deg C
    real, intent(in) :: soilmstress     ! soil moisture stress factor

    ! function return variable
    real :: dgs                         ! leaf-level stomatal conductance to H2O, mean over daylight hours ( mol H2O m-2 s-1 )

    ! Leaf-level assimilation rate, average over daylight hours
    if (daylength>0.0) then
      dgs = dppfd * soilmstress * dgs_unitiabs * ftemp_kphio / ( 60.0 * 60.0 * daylength )
    else
      dgs = 0.0
    end if
    
  end function calc_dgs


  function calc_drd( fapar, fpc_grid, meanmppfd, rd_unitiabs, ftemp_kphio, soilmstress ) result( my_drd )
    !//////////////////////////////////////////////////////////////////
    ! Calculates daily total dark respiration (Rd) based on monthly mean 
    ! PPFD (assumes acclimation on a monthly time scale).
    ! Not described in Stocker et al., XXX.
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: fapar           ! fraction of absorbed photosynthetically active radiation
    real, intent(in) :: fpc_grid        ! foliar projective cover
    real, intent(in) :: meanmppfd       ! monthly mean PPFD (mol m-2 s-1)
    real, intent(in) :: rd_unitiabs
    real, intent(in) :: ftemp_kphio      ! this day's air temperature, deg C
    real, intent(in) :: soilmstress     ! soil moisture stress factor

    ! function return variable
    real :: my_drd

    ! Dark respiration takes place during night and day (24 hours)
    my_drd = fapar * fpc_grid * meanmppfd * soilmstress * rd_unitiabs * ftemp_kphio * 60.0 * 60.0 * 24.0 * c_molmass

  end function calc_drd


  function calc_dtransp( fapar, acrown, dppfd, transp_unitiabs, ftemp_kphio, soilmstress ) result( my_dtransp )
    !//////////////////////////////////////////////////////////////////
    ! Calculates daily transpiration. 
    ! Exploratory only.
    ! Not described in Stocker et al., XXX.
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: fapar
    real, intent(in) :: acrown
    real, intent(in) :: dppfd              ! daily total photon flux density, mol m-2
    real, intent(in) :: transp_unitiabs
    real, intent(in) :: ftemp_kphio              ! this day's air temperature
    real, intent(in) :: soilmstress        ! soil moisture stress factor

    ! function return variable
    real :: my_dtransp

    ! GPP is light use efficiency multiplied by absorbed light and C-P-alpha
    my_dtransp = fapar * acrown * dppfd * soilmstress * transp_unitiabs * ftemp_kphio * h2o_molmass

  end function calc_dtransp


  function calc_vcmax_canop( fapar, vcmax_unitiabs, meanmppfd ) result( my_vcmax )
    !//////////////////////////////////////////////////////////////////
    ! Calculates canopy-level summed carboxylation capacity (Vcmax). To get
    ! value per unit leaf area, divide by LAI.
    ! Not described in Stocker et al., XXX.
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: fapar
    real, intent(in) :: vcmax_unitiabs
    real, intent(in) :: meanmppfd

    ! function return variable
    real :: my_vcmax    ! canopy-level Vcmax [gCO2/m2-ground/s]

    ! Calculate leafy-scale Rubisco-N as a function of LAI and current LUE
    my_vcmax = fapar * meanmppfd * vcmax_unitiabs

  end function calc_vcmax_canop


  function pmodel( kphio, fapar, ppfd, co2, tc, vpd, elv, c4, method_optci, method_jmaxlim ) result( out_pmodel )
    !//////////////////////////////////////////////////////////////////
    ! Implements the P-model, providing predictions for ci, Vcmax, and 
    ! light use efficiency, etc. 
    ! If fapar and ppfd are provided, calculates GPP, replacing separate
    ! function calc_dgpp().
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: kphio        ! apparent quantum yield efficiency       
    real, intent(in) :: fapar        ! fraction of absorbed photosynthetically active radiation (unitless) 
    real, intent(in) :: ppfd         ! photon flux density (mol/m2)
    real, intent(in) :: co2          ! atmospheric CO2 concentration (ppm)
    real, intent(in) :: tc           ! air temperature (deg C)
    real, intent(in) :: vpd          ! vapor pressure (Pa)
    real, intent(in) :: elv          ! elevation above sea-level (m)
    logical, intent(in) :: c4        ! whether or not C4 photosynthesis pathway is followed. If .false., it's C3.
    character(len=*), intent(in) :: method_optci    ! Method used for deriving optimal ci:ca
    character(len=*), intent(in) :: method_jmaxlim  ! Method used for accounting for Jmax limitation

    ! function return value
    type(outtype_pmodel) :: out_pmodel

    ! local variables
    real :: iabs                ! absorbed photosynthetically active radiation (mol/m2)
    real :: patm                ! atmospheric pressure as a function of elevation (Pa)
    real :: kmm                 ! Michaelis-Menten coefficient (Pa)
    real :: gammastar           ! photorespiratory compensation point - Gamma-star (Pa)
    real :: ca                  ! ambient CO2 partial pressure, (Pa)
    real :: gs                  ! stomatal conductance to H2O, expressed per units absorbed light (mol H2O m-2 m-1 / (mol light m-2))
    real :: ci                  ! leaf-internal partial pressure, (Pa)
    real :: chi                 ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
    real :: ns                  ! viscosity of H2O at ambient temperatures (Pa s)
    real :: ns25                ! viscosity of H2O at 25 deg C (Pa s)
    real :: ns_star             ! viscosity correction factor (unitless)
    real :: mprime              ! factor in light use model with Jmax limitation
    real :: iwue                ! intrinsic water use efficiency = A / gs = ca - ci = ca ( 1 - chi ) , unitless
    real :: lue                 ! light use efficiency (mol CO2 / mol photon)
    real :: gpp                 ! gross primary productivity (g CO2 m-2 d-1)
    real :: vcmax               ! canopy-level maximum carboxylation capacity per unit ground area (mol CO2 m-2 s-1)
    real :: vcmax25             ! canopy-level Vcmax25 (Vcmax normalized to 25 deg C) (mol CO2 m-2 s-1)
    real :: vcmax_unitfapar     ! Vcmax per unit fAPAR (mol CO2 m-2 s-1)
    real :: vcmax25_unitfapar   ! Vcmax25 per unit fAPAR (mol CO2 m-2 s-1)
    real :: vcmax_unitiabs      ! Vcmax per unit absorbed light (mol CO2 m-2 s-1 mol-1)
    real :: vcmax25_unitiabs    ! Vcmax25 per unit absorbed light (mol CO2 m-2 s-1 mol-1)
    real :: ftemp_inst_vcmax    ! Instantaneous temperature response factor of Vcmax (unitless)
    real :: ftemp_inst_rd       ! Instantaneous temperature response factor of Rd (unitless)
    real :: rd                  ! Dark respiration (mol CO2 m-2 s-1)
    real :: rd_unitfapar        ! Dark respiration per unit fAPAR (mol CO2 m-2 s-1)
    real :: rd_unitiabs         ! Dark respiration per unit absorbed light (mol CO2 m-2 s-1)
    real :: actnv               ! Canopy-level total metabolic leaf N per unit ground area (g N m-2)
    real :: actnv_unitfapar     ! Metabolic leaf N per unit fAPAR (g N m-2)
    real :: actnv_unitiabs      ! Metabolic leaf N per unit absorbed light (g N m-2 mol-1)
    real :: transp              ! Canopy-level total transpiration rate (g H2O (mol photons)-1)

    ! local variables for Jmax limitation following Nick Smith's method
    real :: omega, omega_star, vcmax_unitiabs_star, tcref, jmax_over_vcmax, jmax_prime, jvrat
    real, parameter :: theta = 0.85
    real, parameter :: c_cost = 0.05336251

    type(outtype_chi) :: out_optchi

    ! Prevent floating point exception for extremely low temperatures
    if (tc > -20.0) then
      !-----------------------------------------------------------------------
      ! Calculate photosynthesis model parameters depending on temperature, pressure, and CO2.
      !-----------------------------------------------------------------------
      ! atmospheric pressure as a function of elevation (Pa)
      patm = calc_patm( elv )

      ! ambient CO2 partial pression (Pa)
      ca = co2_to_ca( co2, patm )

      ! photorespiratory compensation point - Gamma-star (Pa)
      gammastar = calc_gammastar( tc, patm )

      ! ! XXX PMODEL_TEST: ok
      ! print*,'tc     ', tc
      ! print*,'patm   ', patm
      ! print*,'elv    ', elv
      ! print*,'gstar  ', gammastar

      ! Michaelis-Menten coef. (Pa)
      kmm  = calc_kmm( tc, patm )
      
      ! XXX PMODEL_TEST: ok
      ! print*, 'kmm ', kmm

      ! viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
      ns      = calc_viscosity_h2o( tc, patm )  ! Pa s 
      ns25    = calc_viscosity_h2o( kTo, kPo )  ! Pa s 
      ns_star = ns / ns25                       ! (unitless)

      ! ! XXX PMODEL_TEST: ok
      ! print*, 'ns_star ', ns_star

      !-----------------------------------------------------------------------
      ! Optimal ci
      ! The heart of the P-model: calculate ci:ca ratio (chi) and additional terms
      !-----------------------------------------------------------------------
      if (c4) then

        out_optchi = calc_chi_c4()

      else if (method_optci=="prentice14") then

        !-----------------------------------------------------------------------
        ! B.2 FULL FORMULATION
        !-----------------------------------------------------------------------
        out_optchi = calc_optimal_chi( kmm, gammastar, ns_star, ca, vpd )
        
        ! ! xxx test
        ! print*,'kmm       : ', kmm
        ! print*,'gammastar : ', gammastar 
        ! print*,'ns_star   : ', ns_star
        ! print*,'ca        : ', ca
        ! print*,'vpd       : ', vpd
        ! stop 

      else

        stop 'PMODEL: select valid method'

      end if 

      ! select case (method_optci)

      !   case ("approx")
      !     !-----------------------------------------------------------------------
      !     ! A. APPROXIMATIVE METHOD
      !     !-----------------------------------------------------------------------
      !     out_optchi = lue_approx( tc, vpd, elv, ca, gammastar, ns, kmm )
                    
      !   case ("C3_simpl")
      !     !-----------------------------------------------------------------------
      !     ! B.1 SIMPLIFIED FORMULATION 
      !     !-----------------------------------------------------------------------
      !     out_optchi = lue_vpd_c3_simpl( kmm, gammastar, ns, ca, vpd )

      !   case ("C3_full")
      !     !-----------------------------------------------------------------------
      !     ! B.2 FULL FORMULATION
      !     !-----------------------------------------------------------------------
      !     out_optchi = calc_optimal_chi( kmm, gammastar, ns_star, ca, vpd )

      !   case ("C4")
      !     !-----------------------------------------------------------------------
      !     ! B.2 FULL FORMULATION
      !     !-----------------------------------------------------------------------
      !     out_optchi = lue_c4()

      !   case default

      !     stop 'PMODEL: select valid method'

      ! end select

      ! ratio of leaf internal to ambient CO2
      chi = out_optchi%chi

      ! ! XXX PMODEL_TEST: ok
      ! print*, 'chi ', chi

      ! leaf-internal CO2 partial pressure (Pa)
      ci = out_optchi%chi * ca  

      !-----------------------------------------------------------------------
      ! Corrolary preditions
      !-----------------------------------------------------------------------
      ! ! stomatal conductance
      ! gs = gpp  / ( ca - ci )

      ! intrinsic water use efficiency 
      iwue = ( ca - ci ) / ( 1.6 * patm )

      !-----------------------------------------------------------------------
      ! Vcmax and light use efficiency
      !-----------------------------------------------------------------------
      if (c4) then
        ! Identical to method_jmaxlim = "wang17"

        ! Include effect of Jmax limitation.
        ! In this case, out_optchi%mj = 1, and mprime = 0.669
        mprime = calc_mprime( out_optchi%mj )

        ! Light use efficiency (gpp per unit absorbed light)
        lue = kphio * mprime * c_molmass  ! in g CO2 m-2 s-1 / (mol light m-2 s-1)

        ! Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
        vcmax_unitiabs = kphio * out_optchi%mjoc * mprime / out_optchi%mj

      else if (method_jmaxlim=="wang17") then

        ! Include effect of Jmax limitation
        mprime = calc_mprime( out_optchi%mj )

        ! Light use efficiency (gpp per unit absorbed light)
        lue = kphio * mprime * c_molmass  ! in g CO2 m-2 s-1 / (mol light m-2 s-1)

        ! XXX PMODEL_TEST: ok
        ! print*, 'lue ', lue / c_molmass
        ! stop

        ! Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
        vcmax_unitiabs = kphio * out_optchi%mjoc * mprime / out_optchi%mj

        ! ! xxx test
        ! print*,'kphio           : ', kphio
        ! print*,'out_optchi%mjoc : ', out_optchi%mjoc 
        ! print*,'mprime          : ', mprime
        ! print*,'out_optchi%mj   : ', out_optchi%mj
        ! stop 

      else if (method_jmaxlim=="smith19") then

        ! mc = (ci - gammastar) / (ci + kmm)                       ! Eq. 6
        ! print(paste("mc should be equal: ", mc, out_optchi%mc ) )

        ! mj = (ci - gammastar) / (ci + 2.0 * gammastar)           ! Eq. 8
        ! print(paste("mj should be equal: ", mj, out_optchi%mj ) )

        ! mjoc = (ci + kmm) / (ci + 2.0 * gammastar)               ! mj/mc, used in several instances below
        ! print(paste("mjoc should be equal: ", mjoc, out_optchi%mjoc ) )

        omega = calc_omega( theta = theta, c_cost = c_cost, m = out_optchi%mj )             ! Eq. S4
        omega_star = 1.0 + omega - sqrt( (1.0 + omega)**2 - (4.0 * theta * omega) )       ! Eq. 18
        
        ! calculate Vcmax-star, which corresponds to Vcmax at a reference temperature 'tcref'
        vcmax_unitiabs_star  = kphio * out_optchi%mjoc * omega_star / (8.0 * theta)               ! Eq. 19
        
        ! tcref is the optimum temperature in K, assumed to be the temperature at which Vcmax* is operating. 
        ! tcref is estimated based on its relationship to growth temperature following Kattge & Knorr 2007
        tcref = 0.44 * tc + 24.92

        ! calculated acclimated Vcmax at prevailing growth temperatures
        ftemp_inst_vcmax = calc_ftemp_inst_vcmax( tc, tc, tcref = tcref )
        vcmax_unitiabs = vcmax_unitiabs_star * ftemp_inst_vcmax   ! Eq. 20
        
        ! calculate Jmax
        jmax_over_vcmax = (8.0 * theta * omega) / (out_optchi%mjoc * omega_star)             ! Eq. 15 / Eq. 19
        jmax_prime = jmax_over_vcmax * vcmax_unitiabs 

        ! light use efficiency
        lue = c_molmass * kphio * out_optchi%mj * omega_star / (8.0 * theta) ! * calc_ftemp_inst_vcmax( tc, tc, tcref = tcref )     ! treat theta as a calibratable parameter


      else if (method_jmaxlim=="none") then

        ! Light use efficiency (gpp per unit absorbed light)
        lue = kphio * out_optchi%mj * c_molmass

        ! Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
        vcmax_unitiabs = kphio * out_optchi%mjoc

      else

        stop 'PMODEL: select valid method'

      end if

      ! ! XXX PMODEL_TEST: ok
      ! print*, 'm ', out_optchi%m

      ! ! XXX PMODEL_TEST: ok
      ! print*, 'chi ', chi


      ! ! XXX PMODEL_TEST: ok
      ! print*, 'mprime ', mprime

      ! ! XXX PMODEL_TEST: ok
      ! print*, 'lue ', lue

      !-----------------------------------------------------------------------
      ! Corrolary preditions (This is prelimirary!)
      !-----------------------------------------------------------------------
      ! Vcmax25 (vcmax normalized to 25 deg C)
      ftemp_inst_vcmax  = calc_ftemp_inst_vcmax( tc, tc, tcref = 25.0 )
      vcmax25_unitiabs  = vcmax_unitiabs  / ftemp_inst_vcmax

      ! Dark respiration at growth temperature
      ftemp_inst_rd = calc_ftemp_inst_rd( tc )
      rd_unitiabs  = params_gpp%rd_to_vcmax * (ftemp_inst_rd / ftemp_inst_vcmax) * vcmax_unitiabs 

      ! active metabolic leaf N (canopy-level), mol N/m2-ground (same equations as for nitrogen content per unit leaf area, gN/m2-leaf)
      actnv_unitiabs  = vcmax25_unitiabs  * n_v

      !   ! stomatal conductance to H2O, expressed per unit absorbed light
      !   gs_unitiabs = 1.6 * lue * patm / ( ca - ci )

      if (ppfd /= dummy) then
        !-----------------------------------------------------------------------
        ! Calculate quantities scaling with light assuming fAPAR = 1
        ! representing leaf-level at the top of the canopy.
        !-----------------------------------------------------------------------
        ! ! leaf-level assimilation rate
        ! assim = lue * ppfd

        ! ! stomatal conductance to CO2
        ! gs = assim  / ( ca - ci )

        ! ! Transpiration (E)
        ! transp = 1.6 * gs * vpd

        ! Vcmax normalised per unit fAPAR (assuming fAPAR=1)
        vcmax_unitfapar = ppfd * vcmax_unitiabs

        ! Vcmax25 (vcmax normalized to 25 deg C)
        vcmax25_unitfapar = ppfd * vcmax25_unitiabs

        ! Dark respiration per unit fAPAR (assuming fAPAR=1)
        rd_unitfapar = ppfd * rd_unitiabs

        ! active metabolic leaf N (canopy-level), mol N/m2-ground (same equations as for nitrogen content per unit leaf area, gN/m2-leaf)
        actnv_unitfapar = ppfd * actnv_unitiabs

        if (fapar /= dummy) then
          !-----------------------------------------------------------------------
          ! Calculate quantities scaling with absorbed light
          !-----------------------------------------------------------------------
          ! absorbed photosynthetically active radiation (mol/m2)
          iabs = fapar * ppfd 

          ! XXX PMODEL_TEST: ok
          ! print*, 'iabs ', iabs

          ! Canopy-level quantities 
          ! Defined per unit ground level -> scaling with aborbed light (iabs)
          !-----------------------------------------------------------------------
          ! Gross primary productivity
          gpp = iabs * lue ! in g C m-2 s-1

          ! XXX PMODEL_TEST: ok
          ! print*, 'gpp ', gpp

          ! Vcmax per unit ground area is the product of the intrinsic quantum 
          ! efficiency, the absorbed PAR, and 'n'
          vcmax = iabs * vcmax_unitiabs  ! = iabs * kphio * n 

          ! XXX PMODEL_TEST: ok
          ! print*, 'vcmax ', vcmax

          ! (vcmax normalized to 25 deg C)
          vcmax25 = iabs * vcmax25_unitiabs  ! = factor25_vcmax * vcmax

          ! XXX PMODEL_TEST: ok
          ! print*, 'vcmax25 ', vcmax25

          ! Dark respiration
          rd = iabs * rd_unitiabs ! = rd_to_vcmax * vcmax

          ! XXX PMODEL_TEST: ok
          ! print*, 'rd ', rd

          ! active metabolic leaf N (canopy-level), mol N/m2-ground (same equations as for nitrogen content per unit leaf area, gN/m2-leaf)
          actnv = iabs * actnv_unitiabs ! = vcmax25 * n_v

        else

          gpp     = dummy
          vcmax   = dummy
          vcmax25 = dummy
          rd      = dummy
          actnv   = dummy

        end if

      else

        vcmax_unitfapar   = dummy
        vcmax25_unitfapar = dummy
        rd_unitfapar      = dummy
        actnv_unitfapar   = dummy
        
        gpp               = dummy
        vcmax             = dummy
        vcmax25           = dummy
        rd                = dummy
        actnv             = dummy

      end if
    else
      actnv_unitiabs = 0.0
      actnv_unitfapar = 0.0
      actnv = 0.0
      rd_unitiabs = 0.0
      rd_unitfapar = 0.0
      rd = 0.0
      ftemp_inst_rd = 0.0
      ftemp_inst_vcmax = 0.0
      vcmax_unitiabs = 0.0
      vcmax_unitfapar = 0.0
      vcmax25 = 0.0
      vcmax = 0.0
      gpp = 0.0
      lue = 0.0
      iwue = 0.0
      chi = 0.0
      ci = 0.0
      kmm = 0.0
      gammastar = 0.0

    end if

    ! construct list for output
    out_pmodel%gammastar        = gammastar
    out_pmodel%kmm              = kmm
    out_pmodel%ci               = ci
    out_pmodel%chi              = chi
    out_pmodel%iwue             = iwue
    out_pmodel%lue              = lue
    out_pmodel%gpp              = gpp
    out_pmodel%vcmax            = vcmax
    out_pmodel%vcmax25          = vcmax25
    out_pmodel%vcmax_unitfapar  = vcmax_unitfapar
    out_pmodel%vcmax_unitiabs   = vcmax_unitiabs
    out_pmodel%ftemp_inst_vcmax = ftemp_inst_vcmax
    out_pmodel%ftemp_inst_rd    = ftemp_inst_rd
    out_pmodel%rd               = rd
    out_pmodel%rd_unitfapar     = rd_unitfapar
    out_pmodel%rd_unitiabs      = rd_unitiabs
    out_pmodel%actnv            = actnv
    out_pmodel%actnv_unitfapar  = actnv_unitfapar
    out_pmodel%actnv_unitiabs   = actnv_unitiabs


  end function pmodel


  function calc_optimal_chi( kmm, gammastar, ns_star, ca, vpd ) result( out_optchi )
    !//////////////////////////////////////////////////////////////////
    ! Output:   float, ratio of ci/ca (chi)
    ! Features: Returns an estimate of leaf internal to ambient CO2
    !           partial pressure following the "simple formulation".
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: kmm       ! Pa, Michaelis-Menten coeff.
    real, intent(in) :: gammastar        ! Pa, photores. comp. point (Gamma-star)
    real, intent(in) :: ns_star   ! (unitless) viscosity correction factor for water
    real, intent(in) :: ca        ! Pa, ambient CO2 partial pressure
    real, intent(in) :: vpd       ! Pa, vapor pressure deficit

    ! function return value
    type(outtype_chi) :: out_optchi

    ! local variables
    real :: chi                   ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio (unitless)
    real :: xi
    real :: gamma
    real :: kappa
    real :: mc, mj, mjoc

    ! variable substitutes
    real :: vdcg, vacg, vbkg, vsr

    ! beta = 1.6 * ns_star * vpd * (chi * ca - gammastar) ** 2.0 / ( (kmm + gammastar) * (ca ** 2.0) * (chi - 1.0) ** 2.0 )   ! see Estimation_of_beta.pdf

    ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    xi  = sqrt( ( params_gpp%beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )     ! see Eq. 2 in 'Estimation_of_beta.pdf'
    chi = gammastar / ca + ( 1.0 - gammastar / ca ) * xi / ( xi + sqrt(vpd) )           ! see Eq. 1 in 'Estimation_of_beta.pdf'

    ! consistent with this, directly return light-use-efficiency (m)
    ! see Eq. 13 in 'Simplifying_LUE.pdf'

    ! light use efficiency (m)
    ! m = (ca - gammastar)/(ca + 2.0 * gammastar + 3.0 * gammastar * sqrt( (1.6 * vpd) / (beta * (K + gammastar) / ns_star ) ) )

    ! Define variable substitutes:
    vdcg = ca - gammastar
    vacg = ca + 2.0 * gammastar
    vbkg = params_gpp%beta * (kmm + gammastar)

    ! Check for negatives:
    if (vbkg > 0) then
      vsr = sqrt( 1.6 * ns_star * vpd / vbkg )

      ! Based on the m' formulation (see Regressing_LUE.pdf)
      mj = vdcg / ( vacg + 3.0 * gammastar * vsr )
    end if

    gamma = gammastar / ca
    kappa = kmm / ca

    ! mc
    mc = (chi - gamma) / (chi + kappa)

    ! mj:mv
    mjoc  = (chi + kappa) / (chi + 2 * gamma)

    ! return derived type
    out_optchi%chi  = chi
    out_optchi%mj   = mj
    out_optchi%mc   = mc
    out_optchi%mjoc = mjoc
  
  end function calc_optimal_chi


  function calc_chi_c4() result( out_chi )
    !//////////////////////////////////////////////////////////////////
    ! Output:   float, ratio of ci/ca (chi)
    ! Features: Returns an estimate of leaf internal to ambient CO2
    !           partial pressure following the "simple formulation".
    !-----------------------------------------------------------------------
    ! function return value
    type(outtype_chi) :: out_chi

    ! return derived type
    out_chi%chi  = 9999.9
    out_chi%mj   = 1.0
    out_chi%mc   = 1.0
    out_chi%mjoc = 1.0
  
  end function calc_chi_c4


  function calc_mprime( m ) result( mprime )
    !-----------------------------------------------------------------------
    ! Input:  m   (unitless): factor determining LUE
    ! Output: mpi (unitless): modiefied m accounting for the co-limitation
    !                         hypothesis after Prentice et al. (2014)
    !-----------------------------------------------------------------------
    ! argument
    real, intent(in) :: m

    ! local variables
    real, parameter :: kc = 0.41          ! Jmax cost coefficient

    ! function return variable
    real :: mprime

    ! square of m-prime (mpi)
    mprime = m**2 - kc**(2.0/3.0) * (m**(4.0/3.0))

    ! Check for negatives and take root of square
    if (mprime > 0) then
      mprime = sqrt(mprime)
    else
      print*,'negative mprime (', mprime, '). Setting to zero.'
      mprime = 0.0
    end if 
    
  end function calc_mprime

  function calc_omega( theta, c_cost, m ) result( omega )
    !-----------------------------------------------------------------------
    ! Adopted from Nick Smith's code:
    ! Calculate omega, see Smith et al., 2019 Ecology Letters
    !-----------------------------------------------------------------------
    use md_sofunutils, only: findroot_quadratic

    ! arguments
    real, intent(in) :: theta
    real, intent(in) :: c_cost
    real, intent(in) :: m

    ! function return variable
    real :: omega

    ! local variables
    real :: cm, v, capP, aquad, bquad, cquad, m_star
    real, dimension(2) :: root

    cm = 4.0 * c_cost / m                        ! simplification term for omega calculation
    v  = 1.0/(cm * (1.0 - theta * cm)) - 4.0 * theta ! simplification term for omega calculation
    
    ! account for non-linearities at low m values
    capP  = (((1.0/1.4) - 0.7)**2 / (1.0-theta)) + 3.4
    aquad = -1.0
    bquad = capP
    cquad = -(capP * theta)
    root  = findroot_quadratic( aquad, bquad, cquad )
    m_star = (4.0 * c_cost) / root(1)
    
    if (m < m_star) then
      omega = -( 1.0 - (2 * theta) ) - sqrt( (1.0 - theta) * v)
    else
      omega = -( 1.0 - (2 * theta))  + sqrt( (1.0 - theta) * v)
    end if
    
  end function calc_omega


  ! function findroot_quadratic( aquad, bquad, cquad, return_smallroot ) result( root )
  !   !-----------------------------------------------------------------------
  !   ! Returns the solution for a quadratic function:
  !   ! a + bx + cx^2 = 0
  !   ! Per default returns root2 
  !   !-----------------------------------------------------------------------
  !   ! arguments
  !   real, intent(in) :: aquad, bquad, cquad

  !   ! function return variable
  !   real :: root

  !   ! local variables
  !   real :: d, root1, root2

  !   d = b*b - 4.0*a*c
  !   if (d >= 0.0) then              ! is it solvable?
  !     d     = sqrt(d)
  !     root1 = (-b + d)/(2.0*a)     ! first root
  !     root2 = (-b - d)/(2.0*a)     ! second root
  !   else                            ! complex roots
  !     stop 'findroot_quadratic(): There is no real root.'
  !   end if

  ! end function findroot_quadratic


  function co2_to_ca( co2, patm ) result( ca )
    !-----------------------------------------------------------------------
    ! Output:   - ca in units of Pa
    ! Features: Converts ca (ambient CO2) from ppm to Pa.
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: co2     ! ambient CO2 in units of ppm
    real, intent(in) :: patm    ! monthly atm. pressure, Pa

    ! function return variable
    real :: ca ! ambient CO2 in units of Pa

    ca = ( 1.e-6 ) * co2 * patm         ! Pa, atms. CO2
      
  end function co2_to_ca


  function ca_to_co2( ca, patm ) result( co2 )
    !-----------------------------------------------------------------------
    ! Output:   - co2 in units of Pa
    ! Features: Converts ca (ambient CO2) from Pa to ppm.
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: ca        ! ambient CO2 in units of Pa
    real, intent(in) :: patm      ! monthly atm. pressure, Pa

    ! function return variable
    real :: co2

    co2   = ca * ( 1.e6 ) / patm
    
  end function ca_to_co2


  function calc_kmm( tc, patm ) result( kmm )
    !-----------------------------------------------------------------------
    ! Features: Returns the temperature & pressure dependent Michaelis-Menten
    !           coefficient, K (Pa).
    ! Ref:      Bernacchi et al. (2001), Improved temperature response 
    !           functions for models of Rubisco-limited photosynthesis, 
    !           Plant, Cell and Environment, 24, 253--259.
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc               ! air temperature, deg C 
    real, intent(in) :: patm             ! atmospheric pressure, Pa

    ! local variables
    real, parameter :: dhac = 79430      ! J/mol, Activation energy, Bernacchi et al. (2001)
    real, parameter :: dhao = 36380      ! J/mol, Activation energy, Bernacchi et al. (2001)
    real, parameter :: kc25 = 39.97      ! Pa, assuming 25 deg C & assuming elevation of 227.076 m.a.s.l.
    real, parameter :: ko25 = 27480      ! Pa, assuming 25 deg C & assuming elevation of 227.076 m.a.s.l.
    real, parameter :: kco  = 2.09476d5  ! ppm, US Standard Atmosphere
    real :: kc, ko, po, rat, tk

    ! function return variable
    real :: kmm                           ! temperature & pressure dependent Michaelis-Menten coefficient, K (Pa).

    ! convert to Kelvin
    tk = tc + 273.15

    kc = kc25 * calc_ftemp_arrhenius( tk, dhac )
    ko = ko25 * calc_ftemp_arrhenius( tk, dhao )

    po  = kco * (1.0d-6) * patm ! O2 partial pressure
    kmm = kc * (1.0 + po/ko)

  end function calc_kmm


  function calc_gammastar( tc, patm ) result( gammastar )
    !-----------------------------------------------------------------------
    ! Features: Returns the temperature-dependent photorespiratory 
    !           compensation point, Gamma star (Pascals), based on constants 
    !           derived from Bernacchi et al. (2001) study. Corresponds
    !           to 'calc_gammastar_colin' in pmodel.R.
    ! Ref:      Colin's document
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc                 ! air temperature (degrees C)
    real, intent(in) :: patm               ! air pressure (Pa)

    ! local variables
    real, parameter :: dha    = 37830      ! J/mol, activation energy, Bernacchi et al. (2001)
    real, parameter :: gs25_0 = 4.332      ! Pa, assuming 25 deg C and sea level (1013.25 mbar)

    real :: tk           ! air temperature (Kelvin)
    real :: gammastar25  ! photorespiratory compensation point at 25 deg C and corrected for atmospheric pressure

    ! function return variable
    real :: gammastar   ! gamma-star (Pa)

    gammastar25 = gs25_0 * patm / calc_patm(0.0) 

    ! conversion to temperature in Kelvin
    tk = tc + 273.15
    gammastar = gammastar25 * calc_ftemp_arrhenius( tk, dha )

  end function calc_gammastar


  function calc_ftemp_inst_vcmax( tcleaf, tcgrowth, tcref ) result( fv )
    !-----------------------------------------------------------------------
    ! arguments
    ! tcleaf: temperature (degrees C)
    ! tref: is 'to' in Nick's set it to 25 C (=298.15 K in other cals)
    !
    ! function return variable
    ! fv: temperature response factor, relative to 25 deg C.
    !
    ! Output:   Factor fv to correct for instantaneous temperature response
    !           of Vcmax for:
    !
    !               Vcmax(temp) = fv * Vcmax(25 deg C) 
    !
    ! Ref:      Wang Han et al. (in prep.)
    !-----------------------------------------------------------------------
    use md_params_core, only: kR           ! Universal gas constant, J/mol/K

    ! arguments
    real, intent(in) :: tcleaf
    real, intent(in) :: tcgrowth
    real, intent(in), optional :: tcref

    ! function return variable
    real :: fv

    ! loal parameters
    real, parameter :: Ha    = 71513  ! activation energy (J/mol)
    real, parameter :: Hd    = 200000 ! deactivation energy (J/mol)
    real, parameter :: a_ent = 668.39 ! offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
    real, parameter :: b_ent = 1.07   ! slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
    
    ! local variables
    real :: tkref, tkleaf, dent, fva, fvb, mytcref

    if (present(tcref)) then
      mytcref = tcref
    else
      mytcref = 298.15
    end if

    tkref = mytcref + 273.15  ! to Kelvin

    ! conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C. 
    tkleaf = tcleaf + 273.15

    ! calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
    dent = a_ent - b_ent * tcgrowth   ! 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
    fva = calc_ftemp_arrhenius( tkleaf, Ha, tkref )
    fvb = (1.0 + exp( (tkref * dent - Hd)/(kR * tkref) ) ) / (1.0 + exp( (tkleaf * dent - Hd)/(kR * tkleaf) ) )
    fv  = fva * fvb

  end function calc_ftemp_inst_vcmax


  function calc_ftemp_inst_rd( tc ) result( fr )
    !-----------------------------------------------------------------------
    ! Output:   Factor fr to correct for instantaneous temperature response
    !           of Rd (dark respiration) for:
    !
    !               Rd(temp) = fr * Rd(25 deg C) 
    !
    ! Ref:      Heskel et al. (2016) used by Wang Han et al. (in prep.)
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc      ! temperature (degrees C)

    ! function return variable
    real :: fr                  ! temperature response factor, relative to 25 deg C.

    ! loal parameters
    real, parameter :: apar = 0.1012
    real, parameter :: bpar = 0.0005
    real, parameter :: tk25 = 298.15 ! 25 deg C in Kelvin

    ! local variables
    real :: tk                  ! temperature (Kelvin)

    ! conversion of temperature to Kelvin
    tk = tc + 273.15

    fr = exp( apar * (tc - 25.0) - bpar * (tc**2 - 25.0**2) )
    
  end function calc_ftemp_inst_rd


  function calc_ftemp_arrhenius( tk, dha, tkref ) result( ftemp )
    !-----------------------------------------------------------------------
    ! Calculates the factor to account for the temperature response following 
    ! Arrhenius: 
    !
    !               var(T) = ftemp * var(T=T_ref)
    !
    ! T_ref is 25 deg C (=298.13 K) per default.
    !-----------------------------------------------------------------------
    use md_params_core, only: kR           ! Universal gas constant, J/mol/K

    ! arguments
    real, intent(in) :: tk                 ! temperature (Kelvin)
    real, intent(in) :: dha                ! activation energy (J/mol)
    real, intent(in), optional :: tkref    ! reference temperature 

    ! local variables
    real :: mytkref                        ! reference temperature 

    ! function return variable
    real :: ftemp

    if (present(tkref)) then
      mytkref = tkref
    else
      mytkref = 298.15
    end if

    ftemp = exp( dha * (tk - mytkref) / (mytkref * kR * tk) )

  end function calc_ftemp_arrhenius

  ! XXX REMOVED BECAUSE IT'S NOW IN SOFUNUTILS
  ! function calc_patm( elv ) result( patm )
  !   !-----------------------------------------------------------------------
  !   ! Features: Returns the atmospheric pressure as a function of elevation
  !   !           and standard atmosphere (1013.25 hPa)
  !   ! Depends:  - connect_sql
  !   !           - flux_to_grid
  !   !           - get_data_point
  !   !           - get_msvidx
  !   ! Ref:      Allen et al. (1998)
  !   !-----------------------------------------------------------------------
  !   ! argument
  !   real, intent(in) :: elv           ! elevation above sea level, m

  !   ! local variables
  !   real, parameter :: kPo = 101325   ! standard atmosphere, Pa (Allen, 1973)
  !   real, parameter :: kTo = 298.15   ! base temperature, K (Prentice, unpublished)
  !   real, parameter :: kL  = 0.0065   ! temperature lapse rate, K/m (Allen, 1973)
  !   real, parameter :: kG  = 9.80665  ! gravitational acceleration, m/s**2 (Allen, 1973)
  !   real, parameter :: kR  = 8.3145   ! universal gas constant, J/mol/K (Allen, 1973)
  !   real, parameter :: kMa = 0.028963 ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)

  !   ! function return variable
  !   real :: patm    ! atmospheric pressure at elevation 'elv', Pa 

  !   ! Convert elevation to pressure, Pa:
  !   patm = kPo*(1.0 - kL*elv/kTo)**(kG*kMa/(kR*kL))
    
  ! end function calc_patm


  function calc_density_h2o( tc, patm ) result( density_h2o )
    !-----------------------------------------------------------------------
    ! Features: Calculates density of water at a given temperature and 
    !           pressure using the Tumlirz Equation
    ! Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of 
    !           pure water and sea water, Tech. Rept., Marine Physical 
    !           Laboratory, San Diego, CA.
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc      ! air temperature (tc), degrees C
    real, intent(in) :: patm    ! atmospheric pressure (patm), Pa

    ! local variables
    real :: my_lambda, po, vinf, pbar, vau

    ! function return variable
    real :: density_h2o  ! density of water, kg/m**3

    ! Calculate lambda, (bar cm**3)/g:
    my_lambda = 1788.316 + &
                21.55053*tc + &
            (-0.4695911)*tc*tc + &
           (3.096363e-3)*tc*tc*tc + &
    (-1.0)*(7.341182e-6)*tc*tc*tc*tc

    ! Calculate po, bar
    po = 5918.499 + & 
                58.05267*tc + & 
            (-1.1253317)*tc*tc + & 
          (6.6123869e-3)*tc*tc*tc + & 
    (-1.0)*(1.4661625e-5)*tc*tc*tc*tc

    ! Calculate vinf, cm**3/g
    vinf = 0.6980547 + &
    (-1.0)*(7.435626e-4)*tc + &
           (3.704258e-5)*tc*tc + &
    (-1.0)*(6.315724e-7)*tc*tc*tc + &
           (9.829576e-9)*tc*tc*tc*tc + &
    (-1.0)*(1.197269e-10)*tc*tc*tc*tc*tc + &
          (1.005461e-12)*tc*tc*tc*tc*tc*tc + &
    (-1.0)*(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc + &
           (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc + &
    (-1.0)*(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc

    ! Convert pressure to bars (1 bar = 100000 Pa)
    pbar = (1e-5)*patm
    
    ! Calculate the specific volume (cm**3 g**-1):
    vau = vinf + my_lambda/(po + pbar)

    ! Convert to density (g cm**-3) -> 1000 g/kg; 1000000 cm**3/m**3 -> kg/m**3:
    density_h2o = (1.0e3/vau)

  end function calc_density_h2o


  function calc_viscosity_h2o( tc, patm ) result( viscosity_h2o )
    !-----------------------------------------------------------------------
    ! Features: Calculates viscosity of water at a given temperature and 
    !           pressure.
    ! Depends:  density_h2o
    ! Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. 
    !           Sengers, M. J. Assael, ..., K. Miyagawa (2009) New 
    !           international formulation for the viscosity of H2O, J. Phys. 
    !           Chem. Ref. Data, Vol. 38(2), pp. 101-125.
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc      ! air temperature (tc), degrees C
    real, intent(in) :: patm    ! atmospheric pressure (patm), Pa

    ! local variables
    real, parameter :: tk_ast  = 647.096    ! Kelvin
    real, parameter :: rho_ast = 322.0      ! kg/m**3
    real, parameter :: mu_ast  = 1e-6       ! Pa s

    real, dimension(7,6) :: h_array
    real :: rho                             ! density of water kg/m**3
    real :: tbar, tbarx, tbar2, tbar3, rbar, mu0, mu1, ctbar, mu_bar, &
      coef1, coef2
    integer :: i, j                         ! counter variables

    ! function return variable
    real :: viscosity_h2o

    ! print*,'----- in calc_viscosity_h2o() ------'
    ! print*,'tc ', tc
    ! print*,'patm ', patm

    ! Get the density of water, kg/m**3
    rho = calc_density_h2o(tc, patm)
    ! print*,'rho ', rho

    ! Calculate dimensionless parameters:
    tbar = (tc + 273.15)/tk_ast
    tbarx = tbar**(0.5)
    tbar2 = tbar**2
    tbar3 = tbar**3
    rbar = rho/rho_ast
    ! print*,'rbar ', rbar

    ! Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
    mu0 = 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
    mu0 = 1e2*tbarx/mu0
    ! print*,'mu0 ', mu0

    ! Create Table 3, Huber et al. (2009):
    h_array(1,:) = (/0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0/)  ! hj0
    h_array(2,:) = (/0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573/) ! hj1
    h_array(3,:) = (/-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0/) ! hj2
    h_array(4,:) = (/0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0/) ! hj3
    h_array(5,:) = (/-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0/) ! hj4
    h_array(6,:) = (/0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0/) ! hj5
    h_array(7,:) = (/0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264/) ! hj6

    ! Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
    mu1 = 0.0
    ctbar = (1.0/tbar) - 1.0
    do i=1,6
      coef1 = ctbar**(i-1)
      coef2 = 0.0
      do j=1,7
        ! print*,i, j, ' h_array(j,i): ',h_array(j,i)
        coef2 = coef2 + h_array(j,i) * (rbar - 1.0)**(j-1)
        ! print*,i, j, ' coef2: ',coef2
      end do
      mu1 = mu1 + coef1 * coef2    
    end do
    mu1 = exp( rbar * mu1 )
    ! print*, 'mu1 ', mu1

    ! Calculate mu_bar (Eq. 2, Huber et al., 2009)
    !   assumes mu2 = 1
    mu_bar = mu0 * mu1

    ! Calculate mu (Eq. 1, Huber et al., 2009)
    viscosity_h2o = mu_bar * mu_ast    ! Pa s

    ! print*,'----- END calc_viscosity_h2o() ------'

  end function calc_viscosity_h2o


  function calc_soilmstress( soilm, meanalpha, isgrass ) result( outstress )
    !//////////////////////////////////////////////////////////////////
    ! Calculates empirically-derived stress (fractional reduction in light 
    ! use efficiency) as a function of soil moisture
    ! Input:  soilm (unitless, within [0,1]): daily varying soil moisture
    ! Output: outstress (unitless, within [0,1]): function of alpha to reduce GPP 
    !         in strongly water-stressed months
    !-----------------------------------------------------------------------
    ! argument
    real, intent(in) :: soilm                 ! soil water content (fraction)
    real, intent(in) :: meanalpha             ! mean annual AET/PET, average over multiple years (fraction)
    logical, intent(in), optional :: isgrass  ! vegetation cover information to distinguish sensitivity to low soil moisture

    real, parameter :: x0 = 0.0
    real, parameter :: x1 = 0.6

    ! ! Parameters for approach I (simulation s1a)
    ! real, parameter :: apar = 0.2617121
    ! real, parameter :: bpar = 0.5668587
    ! real, parameter :: apar_grass = 0.2617121
    ! real, parameter :: bpar_grass = 0.5668587
    ! real, parameter :: x0 = 0.125
    ! real, parameter :: x1 = 0.75

    ! ! Parameters for approach IV (simulation s1b)
    ! real, parameter :: apar = 0.1785247
    ! real, parameter :: bpar = 0.4501654
    ! real, parameter :: apar_grass = 0.1007749 
    ! real, parameter :: bpar_grass = 0.0063069
    ! real, parameter :: x0 = 0.0
    ! real, parameter :: x1 = 0.9

    ! ! Parameters for approach II (no simulations done with this)
    ! real, parameter :: apar = 0.0606651 
    ! real, parameter :: bpar = 0.5090085 
    ! real, parameter :: apar_grass = 0.0606651
    ! real, parameter :: bpar_grass = 0.5090085
    ! real, parameter :: x0 = 0.0
    ! real, parameter :: x1 = 0.9

    ! ! Parameters for approach III (simulation s1c)
    ! real, parameter :: apar = 0.0515108 
    ! real, parameter :: bpar = 0.1920844
    ! real, parameter :: apar_grass = 0.0515108
    ! real, parameter :: bpar_grass = 0.1920844
    ! real, parameter :: x0 = 0.0
    ! real, parameter :: x1 = 0.9

    real :: y0, beta

    ! function return variable
    real :: outstress

    ! print*,'soilm: ', soilm

    if (soilm > x1) then
      outstress = 1.0
    else
      ! print*,'soilm_par_a, soilm_par_b, meanalpha', params_gpp%soilm_par_a, params_gpp%soilm_par_b, meanalpha

      y0 = (params_gpp%soilm_par_a + params_gpp%soilm_par_b * meanalpha)

      ! if (present(isgrass)) then
      !   if (isgrass) then
      !     y0 = apar_grass + bpar_grass * meanalpha
      !   else
      !     y0 = apar + bpar * meanalpha
      !   end if
      ! else
      !   y0 = apar + bpar * meanalpha
      ! end if

      beta = (1.0 - y0) / (x0 - x1)**2
      outstress = 1.0 - beta * ( soilm - x1 )**2
      outstress = max( 0.0, min( 1.0, outstress ) )
    end if

  end function calc_soilmstress


  function calc_ftemp_kphio( dtemp ) result( ftemp )
    !////////////////////////////////////////////////////////////////
    ! Calculates the instantaneous temperature response of the quantum
    ! yield efficiency based on Bernacchi et al., 2003 PCE (Equation
    ! and parameter values taken from Appendix B)
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: dtemp

    ! function return variable
    real :: ftemp

    ftemp = 0.352 + 0.022 * dtemp - 3.4d-4 * dtemp**2

  end function calc_ftemp_kphio


  subroutine getpar_modl_gpp()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads module-specific parameters from input file.
    !----------------------------------------------------------------
    use md_interface, only: interface
    use md_sofunutils, only: getparreal
    use md_plant, only: params_pft_plant

    ! local variables
    integer :: pft

    !----------------------------------------------------------------
    ! PFT-independent parameters
    !----------------------------------------------------------------
    ! unit cost of carboxylation
    params_gpp%beta  = getparreal( 'params/params_gpp_pmodel.dat', 'beta' )

    ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
    params_gpp%rd_to_vcmax  = getparreal( 'params/params_gpp_pmodel.dat', 'rd_to_vcmax' )

    ! Apply identical temperature ramp parameter for all PFTs
    if (interface%params_siml%is_calib) then
      params_gpp%soilm_par_a    = interface%params_calib%soilm_par_a     ! is provided through standard input
      params_gpp%soilm_par_b    = interface%params_calib%soilm_par_b     ! is provided through standard input
    else
      params_gpp%soilm_par_a    = getparreal( 'params/params_gpp_pmodel.dat', 'soilm_par_a' )
      params_gpp%soilm_par_b    = getparreal( 'params/params_gpp_pmodel.dat', 'soilm_par_b' )
    end if

    ! PFT-dependent parameter(s)
    do pft=1,npft

      if (interface%params_siml%is_calib) then
        params_pft_gpp(pft)%kphio = interface%params_calib%kphio  ! is provided through standard input
      else
        params_pft_gpp(pft)%kphio = getparreal( 'params/params_gpp_pmodel.dat', 'kphio_'//params_pft_plant(pft)%pftname )
      end if

    end do

    return
 
    999  format (I2.2)

  end subroutine getpar_modl_gpp


  subroutine initio_nc_gpp()
    !////////////////////////////////////////////////////////////////
    ! Initialises module-specific NetCDF output files.
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: init_nc_3D_time, check
    use md_interface, only: interface

    ! local variables
    character(len=256) :: prefix

    character(len=*), parameter :: TITLE = "SOFUN GP-model output, module md_gpp"
    character(len=4)  :: year_char
    character(len=12) :: beta_char
    character(len=12) :: rd_to_vcmax_char
    character(len=12) :: kphio_char
    character(len=12) :: soilm_par_a_char
    character(len=12) :: soilm_par_b_char
    character(len=12) :: tempstress_char
    character(len=12) :: soilmstress_char

    integer :: jpngr, doy

    write(year_char,999) interface%steering%outyear

    ! convert parameter values to charaters
    write(beta_char,888)           params_gpp%beta
    write(rd_to_vcmax_char,888)    params_gpp%rd_to_vcmax
    write(kphio_char,888)          params_pft_gpp(1)%kphio
    write(soilm_par_a_char,888)    params_gpp%soilm_par_a
    write(soilm_par_b_char,888)    params_gpp%soilm_par_b
    if (interface%params_siml%soilmstress) then
      soilmstress_char = ".true."
    else
      soilmstress_char = ".false."
    end if
    if (interface%params_siml%tempstress) then
      tempstress_char = ".true."
    else
      tempstress_char = ".false."
    end if

    prefix = "./output_nc/"//trim(interface%params_siml%runname)

    if ( .not. interface%steering%spinup ) then
      !----------------------------------------------------------------
      ! Annual GPP output file 
      !----------------------------------------------------------------
      if (interface%params_siml%loutgpp) then
        ncoutfilnam_agpp = trim(prefix)//'.'//year_char//".a.gpp.nc"
        print*,'initialising ', trim(ncoutfilnam_agpp), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_agpp), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = GPP_NAME, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual gross primary productivivty", &
                          title    = TITLE, &
                          globatt1_nam = "fapar_source",         globatt1_val = interface%params_siml%fapar_forcing_source, &
                          globatt2_nam = "param_beta",           globatt2_val = beta_char, &
                          globatt3_nam = "param_rd_to_vcmax",    globatt3_val = rd_to_vcmax_char, &
                          globatt4_nam = "param_kphio_GrC3",     globatt4_val = kphio_char,  &
                          globatt5_nam = "param_soilm_par_a",    globatt5_val = soilm_par_a_char,  &
                          globatt6_nam = "param_soilm_par_b",    globatt6_val = soilm_par_b_char,  &
                          globatt7_nam = "soilmstress",          globatt7_val = soilmstress_char,  &
                          globatt8_nam = "tempstress",           globatt8_val = tempstress_char    &
                          )
      end if

      if ( interface%steering%outyear>=interface%params_siml%daily_out_startyr .and. &
        interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

        !----------------------------------------------------------------
        ! Daily GPP output file 
        !----------------------------------------------------------------
        if (interface%params_siml%loutdgpp) then
          ncoutfilnam_dgpp = trim(prefix)//'.'//year_char//".d.gpp.nc"
          print*,'initialising ', trim(ncoutfilnam_dgpp), '...'
          call init_nc_3D_time(  filnam  = trim(ncoutfilnam_dgpp), &
                            nlon     = interface%domaininfo%nlon, &
                            nlat     = interface%domaininfo%nlat, &
                            lon      = interface%domaininfo%lon, &
                            lat      = interface%domaininfo%lat, &
                            outyear  = interface%steering%outyear, &
                            outdt    = interface%params_siml%outdt, &
                            outnt    = interface%params_siml%outnt, &
                            varnam   = GPP_NAME, &
                            varunits = "gC m-2 d-1", &
                            longnam  = "daily gross primary productivivty", &
                            title    = TITLE, &
                            globatt1_nam = "fapar_source",      globatt1_val = interface%params_siml%fapar_forcing_source, &
                            globatt2_nam = "param_beta",        globatt2_val = beta_char, &
                            globatt3_nam = "param_rd_to_vcmax", globatt3_val = rd_to_vcmax_char, &
                            globatt4_nam = "param_kphio_GrC3",  globatt4_val = kphio_char,  &
                            globatt5_nam = "param_soilm_par_a", globatt5_val = soilm_par_a_char,  &
                            globatt6_nam = "param_soilm_par_b", globatt6_val = soilm_par_b_char,  &
                            globatt7_nam = "soilmstress",       globatt7_val = soilmstress_char,  &
                            globatt8_nam = "tempstress",        globatt8_val = tempstress_char    &
                            )
        end if

      end if

    end if

    888  format (F12.6)
    999  format (I4.4)
    
  end subroutine initio_nc_gpp


  subroutine initoutput_gpp( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises module-specific output variables
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: ngridcells

    ! daily
    if (interface%steering%init.and.interface%params_siml%loutdgpp    ) allocate( outdgpp(interface%params_siml%outnt,ngridcells) )
    if (interface%steering%init.and.interface%params_siml%loutdrd    ) allocate( outdrd(interface%params_siml%outnt,ngridcells) )
    if (interface%steering%init.and.interface%params_siml%loutdtransp) allocate( outdtransp(interface%params_siml%outnt,ngridcells) )

    if (interface%params_siml%loutdgpp )   outdgpp(:,:)    = 0.0
    if (interface%params_siml%loutdrd    ) outdrd(:,:)     = 0.0
    if (interface%params_siml%loutdtransp) outdtransp(:,:) = 0.0

    ! annual
    if (interface%params_siml%loutgpp) then

      if (interface%steering%init) then
        allocate( outagpp       (npft,ngridcells) )
        allocate( outavcmax     (npft,ngridcells) )
        allocate( outavcmax_25  (npft,ngridcells) )
        allocate( outavcmax_leaf(npft,ngridcells) )
        allocate( outalue       (npft,ngridcells) )
        allocate( outachi       (npft,ngridcells) )
        allocate( outaci        (npft,ngridcells) )
        allocate( outags        (npft,ngridcells) )
        allocate( outaiwue      (npft,ngridcells) )
      end if

      outagpp(:,:)        = 0.0
      outavcmax(:,:)      = 0.0
      outavcmax_25(:,:)   = 0.0
      outavcmax_leaf(:,:) = 0.0
      outachi(:,:)        = 0.0
      outaiwue(:,:)       = 0.0
      outalue(:,:)        = 0.0
      outaci(:,:)         = 0.0
      outags(:,:)         = 0.0
    
    end if

  end subroutine initoutput_gpp


  subroutine getout_daily_gpp( out_pmodel, plant_fluxes, jpngr, doy )
    !////////////////////////////////////////////////////////////////
    ! Called daily to gather daily output variables.
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use md_interface
    use md_plant, only: plant_fluxes_type

    ! argument
    type(outtype_pmodel), dimension(npft), intent(in)    :: out_pmodel
    type(plant_fluxes_type), dimension(npft), intent(in) :: plant_fluxes
    integer, intent(in)                                    :: jpngr
    integer, intent(in)                                    :: doy

    ! local
    integer :: it

    !----------------------------------------------------------------
    ! DAILY FOR HIGH FREQUENCY OUTPUT
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    it = floor( real( doy - 1 ) / real( interface%params_siml%outdt ) ) + 1

    ! sum over PFT
    if (interface%params_siml%loutdgpp   ) outdgpp(it,jpngr)    = outdgpp(it,jpngr)    + sum(plant_fluxes(:)%dgpp)    / real( interface%params_siml%outdt )
    if (interface%params_siml%loutdrd    ) outdrd(it,jpngr)     = outdrd(it,jpngr)     + sum(plant_fluxes(:)%drd)     / real( interface%params_siml%outdt )
    if (interface%params_siml%loutdtransp) outdtransp(it,jpngr) = outdtransp(it,jpngr) + sum(plant_fluxes(:)%dtransp) / real( interface%params_siml%outdt )

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables
    !----------------------------------------------------------------
    ! store all daily values for outputting annual maximum
    ! if (npft>1) stop 'getout_daily_gpp not implemented for npft>1'

    outdvcmax(1,doy)      = dvcmax_canop(1)
    outdvcmax25(1,doy)    = out_pmodel(1)%ftemp_inst_vcmax * dvcmax_canop(1)

    ! weighted by daily GPP
    if (interface%params_siml%loutgpp) then

      outagpp(:,jpngr)        = outagpp(:,jpngr)        + plant_fluxes(:)%dgpp

      outachi       (:,jpngr) = outachi       (:,jpngr) + out_pmodel(1)%chi  * plant_fluxes(:)%dgpp
      outaci        (:,jpngr) = outaci        (:,jpngr) + out_pmodel(1)%ci   * plant_fluxes(:)%dgpp
      outags        (:,jpngr) = outags        (:,jpngr) + dgs(:)             * plant_fluxes(:)%dgpp
      outavcmax_leaf(:,jpngr) = outavcmax_leaf(:,jpngr) + dvcmax_leaf(1)     * plant_fluxes(:)%dgpp
      outaiwue      (:,jpngr) = outaiwue      (:,jpngr) + out_pmodel(1)%iwue * plant_fluxes(:)%dgpp

      if (doy==ndayyear) then
        if (sum(outagpp(:,jpngr))==0.0) then
          outachi       (:,jpngr) = dummy
          outaiwue      (:,jpngr) = dummy
          outaci        (:,jpngr) = dummy
          outags        (:,jpngr) = dummy
          outavcmax_leaf(:,jpngr) = dummy
        else
          outachi       (:,jpngr) = outachi       (:,jpngr) / outagpp(:,jpngr)
          outaiwue      (:,jpngr) = outaiwue      (:,jpngr) / outagpp(:,jpngr)
          outaci        (:,jpngr) = outaci        (:,jpngr) / outagpp(:,jpngr)
          outags        (:,jpngr) = outags        (:,jpngr) / outagpp(:,jpngr)
          outavcmax_leaf(:,jpngr) = outavcmax_leaf(:,jpngr) / outagpp(:,jpngr)
        end if
      end if

    end if


  end subroutine getout_daily_gpp


  subroutine getout_annual_gpp( jpngr )
    !////////////////////////////////////////////////////////////////
    ! Called once a year to gather annual output variables.
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr

    ! local variables
    integer :: pft

    ! outanrlarea(jpngr) = anrlarea
    if (interface%params_siml%loutgpp) then
      ! xxx to do: get vcmax at annual maximum (of monthly values)
      do pft=1,npft
        outavcmax(pft,jpngr)    = maxval(outdvcmax(pft,:))
        outavcmax_25(pft,jpngr) = maxval(outdvcmax25(pft,:))
      end do
    end if

  end subroutine getout_annual_gpp


  ! subroutine writeout_ascii_gpp()
  !   !/////////////////////////////////////////////////////////////////////////
  !   ! Write module-specific ASCII output.
  !   !
  !   ! This is designed for use within SOFUN and requires arguments as
  !   ! derived-types, defined elsewhere. For other applications, implement 
  !   ! the function calls (e.g., calc_dgpp()) differently and 
  !   ! comment/delete this subroutine.
  !   !-------------------------------------------------------------------------
  !   use md_interface

  !   ! local variables
  !   real :: itime
  !   integer :: it, jpngr

  !   ! xxx implement this: sum over gridcells? single output per gridcell?
  !   if (maxgrid>1) stop 'writeout_ascii_gpp: think of something ...'
  !   jpngr = 1

  !   !-------------------------------------------------------------------------
  !   ! DAILY OUTPUT
  !   !-------------------------------------------------------------------------
  !   if ( .not. interface%steering%spinup &
  !        .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
  !        .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

  !     ! Write daily output only during transient simulation
  !     do it=1,interface%params_siml%outnt

  !       ! Define 'itime' as a decimal number corresponding to day in the year + year
  !       itime = real(interface%steering%outyear) + real( it - 1 ) * interface%params_siml%outdt / real( ndayyear )
        
  !       if (interface%params_siml%loutdgpp  )  write(101,999) itime, outdgpp(it,jpngr)
  !       if (interface%params_siml%loutdrd    ) write(135,999) itime, outdrd(it,jpngr)
  !       if (interface%params_siml%loutdtransp) write(114,999) itime, outdtransp(it,jpngr)

  !     end do

  !   end if

  !   !-------------------------------------------------------------------------
  !   ! ANNUAL OUTPUT
  !   ! Write annual value, summed over all PFTs / LUs
  !   ! xxx implement taking sum over PFTs (and gridcells) in this land use category
  !   !-------------------------------------------------------------------------
  !   if (interface%params_siml%loutgpp) then
  
  !     itime = real(interface%steering%outyear)

  !     write(310,999) itime, sum(outagpp(:,jpngr))
  !     write(323,999) itime, sum(outavcmax(:,jpngr))
  !     write(654,999) itime, sum(outavcmax_25(:,jpngr))
  !     write(653,999) itime, sum(outalue(:,jpngr))
  !     write(652,999) itime, sum(outachi(:,jpngr))
  !     write(658,999) itime, sum(outaiwue(:,jpngr)) * 1e6 / 1.6 ! converting from unitless to micro-mol CO2 / mol H2O
  !     write(655,999) itime, sum(outaci(:,jpngr))
  !     write(656,999) itime, sum(outags(:,jpngr))
  !     write(657,999) itime, sum(outavcmax_leaf(:,jpngr))

  !   end if

  !   return
    
  !   999 format (F20.8,F20.8)

  ! end subroutine writeout_ascii_gpp


  subroutine writeout_nc_gpp()
    !/////////////////////////////////////////////////////////////////////////
    ! Writes module-specific NetCDF output
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !-------------------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: write_nc_2D, write_nc_3D_time, check
    use md_interface, only: interface

    if ( .not. interface%steering%spinup ) then
      !-------------------------------------------------------------------------
      ! Annual GPP
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutgpp) print*,'writing ', trim(ncoutfilnam_agpp), '...'
      if (interface%params_siml%loutgpp) call write_nc_2D( trim(ncoutfilnam_agpp), &
                                                              GPP_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              sum( outagpp(:,:), dim=1 ) &
                                                              )

      if (       interface%steering%outyear>=interface%params_siml%daily_out_startyr &
           .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then
        !-------------------------------------------------------------------------
        ! Daily GPP
        !-------------------------------------------------------------------------
        if (interface%params_siml%loutdgpp) print*,'writing ', trim(ncoutfilnam_dgpp), '...'
        if (interface%params_siml%loutdgpp) call write_nc_3D_time( trim(ncoutfilnam_dgpp), &
                                                                GPP_NAME, &
                                                                interface%domaininfo%maxgrid, &
                                                                interface%domaininfo%nlon, &
                                                                interface%domaininfo%nlat, &
                                                                interface%grid(:)%ilon, &
                                                                interface%grid(:)%ilat, &
                                                                interface%params_siml%outnt, &
                                                                interface%grid(:)%dogridcell, &
                                                                outdgpp(:,:) &
                                                                )
      end if
    end if

  end subroutine writeout_nc_gpp


end module md_gpp
