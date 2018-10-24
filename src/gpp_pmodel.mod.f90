module md_gpp
  !////////////////////////////////////////////////////////////////
  ! P-MODEL GPP MODULE
  ! Contains P model functions adopted from GePiSaT
  ! This module ontains the "main" subroutine 'pmodel' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'pmodel' must contain this list 
  ! of subroutines (names that way).
  !   - getpar_modl_gpp
  !   - initio_gpp
  !   - initoutput_gpp
  !   - getout_daily_gpp
  !   - writeout_ascii_gpp
  !   - pmodel
  ! Required module-independent model state variables (necessarily 
  ! updated by 'pmodel') are:
  !   - xxx
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: nmonth, npft, nlu, c_molmass, h2o_molmass, maxgrid, ndayyear

  implicit none

  private
  public params_pft_gpp, getpar_modl_gpp, initio_gpp, initoutput_gpp, &
    gpp, getlue, pmodel, getout_daily_gpp, getout_annual_gpp, &
    writeout_ascii_gpp, outtype_pmodel, calc_tempstress, calc_dgpp, calc_drd, &
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
  real, parameter :: temp1 = 5.0           ! temperature above which ramp function is 1.0 (linear between temp0 and temp1) (deg C)

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
    real :: temp_ramp_edge
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
 
  !-----------------------------------------------------------------------
  ! Email from Tyler (10.3.2015):
  ! I was estimating values of β based on the Wang Han approximation equation 
  ! of χ using both the simplified and "more precise" expressions for χ and ξ 
  ! (Prentice et al., 2014, Ecology Letters).  After examination, Colin and I 
  ! noticed that the value of β is not significantly influenced by the 
  ! expressions for χ and ξ. Since then, Colin has theorised the use of a 
  ! "ground state" universal value of β, which is derived from the Wang Han 
  ! equation at sea level (i.e., z = 0 m and Patm = 101325 Pa), standard temp-
  ! erature (i.e., Tair = 25 deg C) and a non-influencial VPD (i.e., 
  ! D = 1000 Pa). Based on these climatological values, the following were 
  ! calculated:
  !   a. Γ* = 4.220 Pa
  !   b. K = 70.842 Pa
  !   c. η* = 1.0
  !   d. χ = 0.767
  !   e. β = 244.033
  ! Results from modelled versus "observed" monthly GPP based on the universal 
  ! value of β are promising. Colin and I are currently in the works on the next 
  ! set of improvements, which, as I far as I know, will be based on this uni-
  ! versal value of β.
  !-----------------------------------------------------------------------

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC, PRIVATE VARIABLES
  !----------------------------------------------------------------
  ! Function return variables as derived types
  type outtype_pmodel
    real :: gpp
    real :: gstar                 ! photorespiratory compensation point - Gamma-star (Pa)
    real :: chi                   ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
    real :: ci                    ! leaf-internal partial pressure, (Pa)
    real :: ca                    ! ambient partial pressure, (Pa)
    real :: iwue                  ! intrinsic water use efficiency (unitless)
    real :: gs_unitiabs           ! stomatal conductance to H2O, expressed per units absorbed light (mol H2O m-2 m-1 / (mol light m-2))
    real :: vcmax                 ! maximum carboxylation capacity per unit ground area (mol CO2 m-2 s-1)
    real :: vcmax25               ! Vcmax25 (vcmax normalized to 25 deg C) (mol CO2 m-2 s-1)
    real :: vcmax_unitfapar       ! Vcmax per fAPAR (mol CO2 m-2 s-1)
    real :: vcmax_unitiabs        ! Vcmax per unit absorbed light (xxx units)
    real :: ftemp_inst_vcmax      ! Instantaneous temperature response factor of Vcmax (unitless)
    real :: ftemp_inst_rd         ! Instantaneous temperature response factor of Rd (unitless)
    real :: rd                    ! Dark respiration (mol CO2 m-2 s-1)
    real :: rd_unitfapar          ! Dark respiration per fAPAR (mol CO2 m-2 s-1)
    real :: rd_unitiabs           ! Dark respiration per unit absorbed light (mol CO2 m-2 s-1)
    real :: actnv 
    real :: actnv_unitfapar 
    real :: actnv_unitiabs 
    real :: lue                   ! light use efficiency (mol CO2 / mol photon)
    ! real :: transp                ! transpiration [g H2O (mol photons)-1]
    ! real :: transp_unitfapar      ! transpiration per unit fAPAR [g H2O (mol photons)-1]
    ! real :: transp_unitiabs       ! transpiration per unit light absorbed light [g H2O (mol photons)-1]
  end type outtype_pmodel

  type outtype_lue
    real :: chi                   ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
    real :: m
    real :: n
  end type outtype_lue

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, allocatable, dimension(:,:) :: outdgpp    ! daily gross primary production [gC/m2/d]
  real, allocatable, dimension(:,:) :: outdrd     ! daily dark respiration [gC/m2/d]
  real, allocatable, dimension(:,:) :: outdtransp ! daily transpiration [mm]

  ! ! monthly
  ! real, allocatable, dimension(:,:,:) :: outmgpp    ! monthly gross primary production [gC/m2/mo.]
  ! real, allocatable, dimension(:,:,:) :: outmrd     ! monthly dark respiration [gC/m2/mo.]
  ! real, allocatable, dimension(:,:,:) :: outmtransp ! monthly transpiration [mm]

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

  subroutine gpp( plant, plant_fluxes, out_pmodel, dppfd, dayl, meanmppfd, wscal, rlmalpha, doy, moy, dtemp, do_soilmstress, do_tempstress )
    !//////////////////////////////////////////////////////////////////
    ! Calculates daily GPP (gC/m2/d) from monthly acclimated photosynth-
    ! etic parameters (P-model output) and actual daily PPFD and soil
    ! moisture stress (Cramer-Prentice Alpha).
    ! Alternatively ready from data (prescribed GPP array).
    !
    ! Output:   
    ! - gpp (gC/m2/d)   : gross primary production
    !
    !------------------------------------------------------------------
    use md_plant, only: params_pft_plant, plant_type, plant_fluxes_type

    ! arguments
    type(plant_type), dimension(npft), intent(in) :: plant
    type(plant_fluxes_type), dimension(npft), intent(inout) :: plant_fluxes
    type(outtype_pmodel), dimension(npft), intent(in)  :: out_pmodel
    real, intent(in) :: dppfd            ! daily total photon flux density, mol m-2
    real, intent(in) :: dayl             ! day length (h)
    real, intent(in) :: meanmppfd        ! monthly mean PPFD (mol m-2 s-1)
    real, dimension(nlu), intent(in) :: wscal            ! relative soil water content (unitless)
    real, dimension(nlu), intent(in) :: rlmalpha         ! rolling mean alpha (mean annual AET/PET, unitless)
    integer, intent(in) :: doy             ! day of year and month of year
    integer, intent(in) :: moy             ! month of year and month of year
    real,    intent(in) :: dtemp           ! this day's air temperature
    logical, intent(in) :: do_soilmstress  ! whether empirical soil miosture stress function is applied to GPP
    logical, intent(in) :: do_tempstress   ! whether empirical temperature stress function is applied to GPP

    ! local variables
    integer :: pft
    integer :: lu
    real    :: soilmstress
    real    :: tempstress

    !----------------------------------------------------------------
    ! CALCULATE PREDICTED GPP FROM P-model
    ! using instantaneous (daily) LAI, PPFD, Cramer-Prentice-alpha
    !----------------------------------------------------------------
    do pft=1,npft

      lu = params_pft_plant(pft)%lu_category

      ! Calculate soil moisture stress as a function of soil moisture, mean alpha and vegetation type (grass or not)
      if (do_soilmstress) then
        soilmstress = calc_soilmstress( wscal(lu), rlmalpha(lu), params_pft_plant(pft)%grass )
      else
        soilmstress = 1.0
      end if

      ! Calculate low temperature stress
      if (do_tempstress) then
        tempstress = calc_tempstress( dtemp )
      else
        tempstress = 1.0
      end if

      if ( plant(pft)%fpc_grid>0.0 .and. doy>0.0) then

        ! GPP
        plant_fluxes(pft)%dgpp = calc_dgpp( plant(pft)%fapar_ind, plant(pft)%fpc_grid, dppfd, out_pmodel(pft)%lue, tempstress, soilmstress )

        ! ! transpiration
        ! ! dtransp(pft) = calc_dtransp( plant(pft)%fapar_ind, plant(pft)%acrown, dppfd, out_pmodel(pft)%transp_unitiabs, tempstress, soilmstress )
        ! dtransp(pft) = calc_dtransp( plant(pft)%fapar_ind, plant(pft)%acrown, dppfd, out_pmodel(pft)%transp_unitiabs, dtemp )

        ! Dark respiration
        plant_fluxes(pft)%drd = calc_drd( plant(pft)%fapar_ind, plant(pft)%fpc_grid, meanmppfd, out_pmodel(pft)%rd_unitiabs, tempstress, soilmstress )

        ! Leaf-level assimilation rate
        dassim(pft) = calc_dassim( dppfd, out_pmodel(pft)%lue, dayl, tempstress, soilmstress )

        ! stomatal conductance
        dgs(pft) = calc_dgs( dppfd, out_pmodel(pft)%gs_unitiabs, dayl, tempstress, soilmstress )

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

    end do

  end subroutine gpp


  function getlue( co2, dtemp, dvpd, elv ) result( out_pmodel )
    !//////////////////////////////////////////////////////////////////
    ! Averages ambient conditions across given time scale and 
    ! calls pmodel() with averaged values.
    ! xxx todo: abandon this and do the averaging and loop directly in biosphere_()
    !------------------------------------------------------------------
    use md_params_core, only: ndayyear, nlu
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

    ! xxx test
    real, dimension(nmonth)   :: mppfd
    real :: myco2
    real :: myelv

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

        if ( params_pft_plant(pft)%c4 ) then
          ! C4: use infinite CO2 for ci (note lower quantum efficiency 'kphio' parameter for C4)
          out_pmodel(pft,moy) = pmodel( params_pft_gpp(pft)%kphio, -9999.0, -9999.0, 3.0 * co2, mtemp(moy), mvpd(moy), elv, "C4" )
        else
          ! C3
          out_pmodel(pft,moy) = pmodel( params_pft_gpp(pft)%kphio, -9999.0, -9999.0, co2, mtemp(moy), mvpd(moy), elv, "C3_full" )

          ! ! XXX PMODEL_TEST:
          ! out_pmodel(pft,moy) = pmodel( pft, 1.0, mppfd(moy), myco2, mtemp(moy), mvpd(moy), myelv, "C3_full" )
        end if

      end do
    end do

  end function getlue


  function calc_dgpp( fapar, fpc_grid, dppfd, lue, tempstress, soilmstress ) result( my_dgpp )
    !//////////////////////////////////////////////////////////////////
    ! Calculates daily GPP
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: fapar       ! fraction of absorbed photosynthetically active radiation
    real, intent(in) :: fpc_grid    ! foliar projective cover
    real, intent(in) :: dppfd       ! daily total photon flux density, mol m-2
    real, intent(in) :: lue         ! light use efficiency
    real, intent(in) :: tempstress  ! this day's air temperature, deg C 
    real, intent(in) :: soilmstress ! soil moisture stress factor

    ! function return variable
    real :: my_dgpp                 ! Daily total gross primary productivity (gC m-2 d-1)

    ! GPP is light use efficiency multiplied by absorbed light and soil moisture stress function
    my_dgpp = fapar * fpc_grid * dppfd * soilmstress * lue * tempstress * c_molmass

  end function calc_dgpp


  function calc_dassim( dppfd, lue, daylength, tempstress, soilmstress ) result( my_dassim )
    !//////////////////////////////////////////////////////////////////
    ! Calculates leaf-level assimilation rate, mean over daylight hours ( mol CO2 m-2 s-1 )
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: dppfd           ! daily total photon flux density, mol m-2
    real, intent(in) :: lue             ! light use efficiency, mol CO2 / mol photon
    real, intent(in) :: daylength       ! day length (h)
    real, intent(in) :: tempstress      ! this day's air temperature, deg C
    real, intent(in) :: soilmstress     ! soil moisture stress factor

    ! function return variable
    real :: my_dassim                   ! leaf-level assimilation rate, mean over daylight hours ( mol CO2 m-2 s-1 )

    ! Leaf-level assimilation rate, average over daylight hours
    my_dassim = dppfd * soilmstress * lue * tempstress / ( 60.0 * 60.0 * daylength )

  end function calc_dassim


  function calc_dgs( dppfd, dgs_unitiabs, daylength, tempstress, soilmstress ) result( dgs )
    !//////////////////////////////////////////////////////////////////
    ! Calculates leaf-level stomatal conductance to H2O, mean over daylight hours ( mol H2O m-2 s-1 )
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: dppfd           ! daily total photon flux density, mol m-2
    real, intent(in) :: dgs_unitiabs    ! stomatal conductance per unit absorbed light (mol H2O m-2 s-1 / mol light)
    real, intent(in) :: daylength       ! day length (h)
    real, intent(in) :: tempstress      ! this day's air temperature, deg C
    real, intent(in) :: soilmstress     ! soil moisture stress factor

    ! function return variable
    real :: dgs                         ! leaf-level stomatal conductance to H2O, mean over daylight hours ( mol H2O m-2 s-1 )

    ! Leaf-level assimilation rate, average over daylight hours
    dgs = dppfd * soilmstress * dgs_unitiabs * tempstress / ( 60.0 * 60.0 * daylength )

  end function calc_dgs


  function calc_drd( fapar, fpc_grid, meanmppfd, rd_unitiabs, tempstress, soilmstress ) result( my_drd )
    !//////////////////////////////////////////////////////////////////
    ! Calculates daily dark respiration (Rd) based on monthly mean 
    ! PPFD (assumes acclimation on a monthly time scale).
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: fapar           ! fraction of absorbed photosynthetically active radiation
    real, intent(in) :: fpc_grid        ! foliar projective cover
    real, intent(in) :: meanmppfd       ! monthly mean PPFD (mol m-2 s-1)
    real, intent(in) :: rd_unitiabs
    real, intent(in) :: tempstress      ! this day's air temperature, deg C
    real, intent(in) :: soilmstress     ! soil moisture stress factor

    ! function return variable
    real :: my_drd

    ! Dark respiration takes place during night and day (24 hours)
    my_drd = fapar * fpc_grid * meanmppfd * soilmstress * rd_unitiabs * tempstress * 60.0 * 60.0 * 24.0 * c_molmass

  end function calc_drd


  function calc_dtransp( fapar, acrown, dppfd, transp_unitiabs, tempstress, soilmstress ) result( my_dtransp )
    !//////////////////////////////////////////////////////////////////
    ! Calculates daily GPP
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: fapar
    real, intent(in) :: acrown
    real, intent(in) :: dppfd              ! daily total photon flux density, mol m-2
    real, intent(in) :: transp_unitiabs
    real, intent(in) :: tempstress              ! this day's air temperature
    real, intent(in) :: soilmstress        ! soil moisture stress factor

    ! function return variable
    real :: my_dtransp

    ! GPP is light use efficiency multiplied by absorbed light and C-P-alpha
    my_dtransp = fapar * acrown * dppfd * soilmstress * transp_unitiabs * tempstress * h2o_molmass

  end function calc_dtransp


  function calc_vcmax_canop( fapar, vcmax_unitiabs, meanmppfd ) result( my_vcmax )
    !//////////////////////////////////////////////////////////////////
    ! Calculates canopy-level carboxylation capacity (Vcmax). To get
    ! value per unit leaf area, divide by LAI.
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


  function pmodel( kphio, fpar, ppfd, co2, tc, vpd, elv, method ) result( out_pmodel )
    !//////////////////////////////////////////////////////////////////
    ! Implements the P-model, providing predictions for ci, Vcmax, and 
    ! light use efficiency, etc. If fpar and ppfd are provided, calculates GPP.
    !------------------------------------------------------------------
    ! arguments
    real, intent(in) :: kphio        ! apparent quantum yield efficiency       
    real, intent(in) :: fpar         ! monthly fraction of absorbed photosynthetically active radiation (unitless) 
    real, intent(in) :: ppfd         ! monthly photon flux density (mol/m2)
    real, intent(in) :: co2          ! atmospheric CO2 concentration (ppm)
    real, intent(in) :: tc           ! monthly air temperature (deg C)
    real, intent(in) :: vpd          ! mean monthly vapor pressure (Pa) -- CRU data is in hPa
    real, intent(in) :: elv          ! elevation above sea-level (m)
    character(len=*), intent(in) :: method

    ! function return value
    type(outtype_pmodel) :: out_pmodel

    ! local variables
    real :: ppfdabs                  ! absorbed photosynthetically active radiation (mol/m2)
    real :: patm                     ! atmospheric pressure as a function of elevation (Pa)
    real :: chi                      ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
    real :: ci                       ! leaf-internal partial pressure, (Pa)
    real :: iwue                     ! intrinsic water use efficiency = A / gs = ca - ci = ca ( 1 - chi ) , unitless
    real :: gs_unitiabs              ! stomatal conductance to H2O, expressed per unit absorbed light (mol H2O m-2 s-1 / (mol light m-2 -1))
    real :: ca                       ! ambient CO2 partial pressure (Pa)
    real :: gstar                    ! photorespiratory compensation point - Gamma-star (Pa)
    real :: fa                       ! function of alpha to reduce GPP in strongly water-stressed months (unitless)
    real :: kmm                      ! Michaelis-Menten coefficient (Pa)
    real :: ns                       ! viscosity of H2O at ambient temperatures (Pa s)
    real :: ns25                     ! viscosity of H2O at 25 deg C (Pa s)
    real :: ns_star                  ! viscosity correction factor (unitless)
    real :: mprime                   ! factor in light use model with Jmax limitation
    real :: assim                    ! assimilation rate per unit ground area, ecosystem scale (mol m-2 s-1)
    real :: assim_unitfapar          ! assimilation rate per unit leaf area, leaf scale, representative for top-canopy leaves (mol m-2 s-1)
    real :: lue                      ! Light use efficiency = assimilation rate per unit aborbed light
    real :: vcmax                    ! Vcmax per unit ground area (mol CO2 m-2 s-1)
    real :: vcmax_unitfapar          ! Vcmax per fAPAR (mol CO2 m-2 s-1)
    real :: vcmax_unitiabs           ! Vcmax per unit absorbed light (mol CO2 m-2 s-1 [mol PPFD]-1)
    real :: vcmax25                  ! Vcmax25 (vcmax normalized to 25 deg C) (mol CO2 m-2 s-1)
    real :: vcmax25_unitfapar        ! Vcmax25 per fAPAR (mol CO2 m-2 s-1)
    real :: vcmax25_unitiabs         ! Vcmax25 per unit absorbed light
    real :: rd                       ! Dark respiration (mol CO2 m-2 s-1)
    real :: rd_unitfapar             ! Dark respiration per fAPAR (mol CO2 m-2 s-1)
    real :: rd_unitiabs              ! Dark respiration per unit absorbed light (mol CO2 m-2 s-1)
    real :: rd25                     ! Dark respiration normalised to 25 deg C (mol CO2 m-2 s-1)
    real :: rd25_unitfapar           ! Dark respiration per fAPAR normalised to 25 deg C per fAPAR (mol CO2 m-2 s-1)
    real :: rd25_unitiabs            ! Dark respiration per absorbed light normalised to 25 deg C per unit absorbed light (mol CO2 m-2 s-1)
    real :: ftemp_inst_vcmax         ! Instantaneous temperature response factor of Vcmax (unitless)
    real :: ftemp_inst_rd            ! Instantaneous temperature response factor of Rd (unitless)
    real :: actnv 
    real :: actnv_unitfapar 
    real :: actnv_unitiabs 
    ! real :: transp                   ! transpiration [g H2O (mol photons)-1]
    ! real :: transp_unitfapar         ! transpiration per unit fAPAR [g H2O (mol photons)-1]
    ! real :: transp_unitiabs          ! transpiration per unit light absorbed light [g H2O (mol photons)-1]

    type(outtype_lue) :: out_lue

    ! xxx discuss this: approapriate?
    if (tc > 0.0) then

      ! absorbed photosynthetically active radiation (mol/m2)
      ppfdabs = fpar * ppfd

      ! atmospheric pressure as a function of elevation (Pa)
      patm = calc_patm( elv )

      ! ambient CO2 partial pression (Pa)
      ca = co2_to_ca( co2, patm )

      ! photorespiratory compensation point - Gamma-star (Pa)
      gstar = calc_gstar( tc )

      ! XXX PMODEL_TEST: ok
      ! print*, 'gstar ', gstar

      ! Michaelis-Menten coef. (Pa)
      kmm  = calc_k( tc, patm )

      ! XXX PMODEL_TEST: ok
      ! print*, 'kmm ', kmm

      ! viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
      ns      = calc_viscosity_h2o( tc, patm )  ! Pa s 
      ns25    = calc_viscosity_h2o( kTo, kPo )  ! Pa s 
      ns_star = ns / ns25                       ! (unitless)

      ! XXX PMODEL_TEST: ok
      ! print*, 'ns_star ', ns_star

      select case (method)

        case ("approx")
          !-----------------------------------------------------------------------
          ! A. APPROXIMATIVE METHOD
          !-----------------------------------------------------------------------
          out_lue = lue_approx( tc, vpd, elv, ca, gstar, ns, kmm )
                    
        case ("C3_simpl")
          !-----------------------------------------------------------------------
          ! B.1 SIMPLIFIED FORMULATION 
          !-----------------------------------------------------------------------
          out_lue = lue_vpd_c3_simpl( kmm, gstar, ns, ca, vpd )

        case ("C3_full")
          !-----------------------------------------------------------------------
          ! B.2 FULL FORMULATION
          !-----------------------------------------------------------------------
          out_lue = lue_vpd_c3_full( kmm, gstar, ns_star, ca, vpd )

        case ("C4")
          !-----------------------------------------------------------------------
          ! B.2 FULL FORMULATION
          !-----------------------------------------------------------------------
          out_lue = lue_c4()

        case default

          stop 'PMODEL: select valid method'

      end select

      ! LUE-functions return m, n, and chi
      chi = out_lue%chi

      ! XXX PMODEL_TEST: ok
      ! print*, 'm ', out_lue%m

      ! XXX PMODEL_TEST: ok
      ! print*, 'chi ', chi

      !-----------------------------------------------------------------------
      ! Calculate function return variables
      !-----------------------------------------------------------------------

      ! GPP per unit ground area is the product of the intrinsic quantum 
      ! efficiency, the absorbed PAR, the function of alpha (drought-reduction),
      ! and 'm'
      mprime = calc_mprime( out_lue%m )

      ! Light use efficiency (assimilation rate per unit absorbed light)
      lue = kphio * mprime  ! in mol CO2 m-2 s-1 / (mol light m-2 s-1)

      ! Gross primary productivity = ecosystem-level assimilation rate (per unit ground area)
      assim = ppfdabs * lue ! in mol CO2 m-2 s-1

      ! Leaf-level assimilation rate (per unit leaf area), representative for top-canopy leaves
      assim_unitfapar = ppfd * lue  ! in mol m-2 s-1

      ! ! XXX PMODEL_TEST: ok
      ! print*, 'lue ', lue
      ! print*, 'chi ', chi

      ! leaf-internal CO2 partial pressure (Pa)
      ci = chi * ca

      ! stomatal conductance to H2O, expressed per unit absorbed light
      gs_unitiabs = 1.6 * lue * patm / ( ca - ci )

      ! Vcmax per unit ground area is the product of the intrinsic quantum 
      ! efficiency, the absorbed PAR, and 'n'
      vcmax = ppfdabs * kphio * out_lue%n

      ! Vcmax normalised per unit fAPAR (assuming fAPAR=1)
      vcmax_unitfapar = ppfd * kphio * out_lue%n 

      ! Vcmax normalised per unit absorbed PPFD (assuming ppfdabs=1)
      vcmax_unitiabs = kphio * out_lue%n 

      ! Vcmax25 (vcmax normalized to 25 deg C)
      ftemp_inst_vcmax  = calc_ftemp_inst_vcmax( tc )
      vcmax25           = vcmax           / ftemp_inst_vcmax
      vcmax25_unitfapar = vcmax_unitfapar / ftemp_inst_vcmax
      vcmax25_unitiabs  = vcmax_unitiabs  / ftemp_inst_vcmax

      ! ! Dark respiration, normalised to 25 deg C proportional to Vcmax normalised to 25 deg C
      ! rd25           = params_gpp%rd_to_vcmax * vcmax25
      ! rd25_unitfapar = params_gpp%rd_to_vcmax * vcmax25_unitfapar
      ! rd25_unitiabs  = params_gpp%rd_to_vcmax * vcmax25_unitiabs 

      ! Dark respiration at growth temperature
      ftemp_inst_rd = calc_ftemp_inst_rd( tc )
      ! rd            = ftemp_inst_rd * rd25
      ! rd_unitfapar  = ftemp_inst_rd * rd25_unitfapar
      ! rd_unitiabs   = ftemp_inst_rd * rd25_unitiabs 

      ! print*,'fr/fv: ', ftemp_inst_rd / ftemp_inst_vcmax
      rd           = params_gpp%rd_to_vcmax * (ftemp_inst_rd / ftemp_inst_vcmax) * vcmax
      rd_unitfapar = params_gpp%rd_to_vcmax * (ftemp_inst_rd / ftemp_inst_vcmax) * vcmax_unitfapar
      rd_unitiabs  = params_gpp%rd_to_vcmax * (ftemp_inst_rd / ftemp_inst_vcmax) * vcmax_unitiabs 

      ! active metabolic leaf N (canopy-level), mol N/m2-ground (same equations as for nitrogen content per unit leaf area, gN/m2-leaf)
      actnv           = vcmax25           * n_v
      actnv_unitfapar = vcmax25_unitfapar * n_v
      actnv_unitiabs  = vcmax25_unitiabs  * n_v

      ! ! Transpiration (E)
      ! ! Using 
      ! ! - E = 1.6 gs D
      ! ! - gs = A / (ca (1-chi))
      ! ! (- chi = ci / ca)
      ! ! => E = f
      ! transp           = (1.6 * ppfdabs * kphio * fa * mprime * vpd) / (ca - ci)   ! gpp = ppfdabs * kphio * fa * m
      ! transp_unitfapar = (1.6 * ppfd * kphio * fa * mprime * vpd) / (ca - ci)
      ! transp_unitiabs  = (1.6 * 1.0  * kphio * fa * mprime * vpd) / (ca - ci)

      ! Construct derived type for output
      out_pmodel%gpp              = assim * c_molmass
      out_pmodel%gstar            = gstar
      out_pmodel%chi              = chi
      out_pmodel%ci               = co2 * chi  ! return value 'out_pmodel%ci' is used for output in units of ppm. 
      out_pmodel%ca               = ca
      out_pmodel%iwue             = ( ca - ci ) / ( 1.6 * patm )
      out_pmodel%gs_unitiabs      = gs_unitiabs
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
      out_pmodel%lue              = lue
      ! out_pmodel%transp           = transp          
      ! out_pmodel%transp_unitfapar = transp_unitfapar
      ! out_pmodel%transp_unitiabs  = transp_unitiabs 

    else

      ! for monthly mean temperatures below 0 deg C:
      out_pmodel%gpp              = 0.0
      out_pmodel%gstar            = 0.0
      out_pmodel%chi              = 0.0
      out_pmodel%ci               = 0.0
      out_pmodel%iwue             = 0.0
      out_pmodel%gs_unitiabs      = 0.0
      out_pmodel%vcmax            = 0.0
      out_pmodel%vcmax25          = 0.0
      out_pmodel%vcmax_unitfapar  = 0.0
      out_pmodel%vcmax_unitiabs   = 0.0
      out_pmodel%ftemp_inst_vcmax = 0.0
      out_pmodel%ftemp_inst_rd    = 0.0
      out_pmodel%rd               = 0.0
      out_pmodel%rd_unitfapar     = 0.0
      out_pmodel%rd_unitiabs      = 0.0
      out_pmodel%actnv            = 0.0
      out_pmodel%actnv_unitfapar  = 0.0
      out_pmodel%actnv_unitiabs   = 0.0
      out_pmodel%lue              = 0.0
      ! out_pmodel%transp           = 0.0
      ! out_pmodel%transp_unitfapar = 0.0
      ! out_pmodel%transp_unitiabs  = 0.0

    end if

  end function pmodel


  function calc_soilmstress( soilm, meanalpha, isgrass ) result( outstress )
    !-----------------------------------------------------------------------
    ! Calculates empirically-derived stress (fractional reduction in light 
    ! use efficiency) as a function of soil moisture
    ! Input:  soilm (unitless, within [0,1]): daily varying soil moisture
    ! Output: outstress (unitless, within [0,1]): function of alpha to reduce GPP 
    !         in strongly water-stressed months
    !-----------------------------------------------------------------------
    ! argument
    real, intent(in) :: soilm       ! soil water content (fraction)
    real, intent(in) :: meanalpha   ! mean annual AET/PET (fraction)
    logical, intent(in), optional :: isgrass 

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


  subroutine getpar_modl_gpp()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads waterbalance module-specific parameters 
    ! from input file
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
      params_gpp%temp_ramp_edge = interface%params_calib%temp_ramp_edge  ! is provided through standard input
      params_gpp%soilm_par_a    = interface%params_calib%soilm_par_a     ! is provided through standard input
      params_gpp%soilm_par_b    = interface%params_calib%soilm_par_b     ! is provided through standard input
    else
      params_gpp%temp_ramp_edge = getparreal( 'params/params_gpp_pmodel.dat', 'temp_ramp_edge' )
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


  function lue_approx( temp, vpd, elv, ca, gstar, ns_star, kmm ) result( out_lue )
    !//////////////////////////////////////////////////////////////////
    ! Output:   list: 'm' (unitless), 'chi' (unitless)
    ! Returns list containing light use efficiency (m) and ci/ci ratio 
    ! (chi) based on the approximation of the theoretical relationships
    ! of chi with temp, vpd, and elevation. Is now based on SI units as 
    ! inputs.
    !------------------------------------------------------------------

    ! arguments
    real, intent(in) :: temp      ! deg C, air temperature
    real, intent(in) :: vpd       ! Pa, vapour pressure deficit
    real, intent(in) :: elv       ! m, elevation above sea level
    real, intent(in) :: ca        ! Pa, ambient CO2 partial pressure
    real, intent(in) :: gstar ! Pa, photores. comp. point (Gamma-star)
    real, intent(in) :: ns_star   ! (unitless) viscosity correction factor for water
    real, intent(in) :: kmm       ! Pa, Michaelis-Menten coeff.

    ! function return value
    type(outtype_lue) :: out_lue

    ! local variables
    ! real :: beta_wh
    real :: whe                ! value of "Wang-Han Equation"
    real :: gamma
    real :: chi                ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio (unitless)
    real :: m

    ! ! variable substitutes
    ! real :: vdcg, vacg, vbkg, vsr

    ! Wang-Han Equation
    whe = exp( &
      1.19 &
      + 0.0545 * ( temp - 25.0 ) &
      - 0.5 * log( vpd ) &   ! convert vpd from Pa to kPa 
      - 8.15e-5 * elv &      ! convert elv from m to km
      )

    ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    chi = whe / ( 1.0 + whe )

    !  m
    gamma = gstar / ca
    m = (chi - gamma) / (chi + 2 * gamma)

    ! xxx try
    ! ! beta derived from chi and xi (see Estimation_of_beta.pdf). Uses empirical chi
    ! beta_wh = ( 1.6 * ns_star * vpd * ( chi * ca - gstar )**2 ) / ( ca**2 * ( chi - 1.0 )**2 * ( kmm + gstar ) )

    ! ! Define variable substitutes:
    ! vdcg = ca - gstar
    ! vacg = ca + 2.0 * gstar
    ! vbkg = beta_wh * (kmm + gstar)

    ! ! Check for negatives:
    ! if (vbkg > 0) then
    !   vsr = sqrt( 1.6 * ns_star * vpd / vbkg )

    !   ! Based on the m' formulation (see Regressing_LUE.pdf)
    !   m = vdcg / ( vacg + 3.0 * gstar * vsr )
    ! end if
    ! xxx

    ! return derived type
    out_lue%chi = chi
    out_lue%m = m
    out_lue%n = -9999
  
  end function lue_approx


  function lue_vpd_c3_simpl( kmm, gstar, ns_star, ca, vpd ) result( out_lue )
    !//////////////////////////////////////////////////////////////////
    ! Output:   float, ratio of ci/ca (chi)
    ! Returns an estimate of leaf internal to ambient CO2
    ! partial pressure following the "simple formulation".
    !-----------------------------------------------------------------------

    ! arguments
    real, intent(in) :: kmm       ! Pa, Michaelis-Menten coeff.
    real, intent(in) :: gstar     ! Pa, photores. comp. point (Gamma-star)
    real, intent(in) :: ns_star   ! (unitless) viscosity correction factor for water
    real, intent(in) :: ca        ! Pa, ambient CO2 partial pressure
    real, intent(in) :: vpd       ! Pa, vapor pressure deficit

    ! function return value
    type(outtype_lue) :: out_lue

    ! local variables
    real :: xi
    real :: gamma
    real :: kappa
    real :: chi                   ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio (unitless)
    real :: m
    real :: n

    ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    xi  = sqrt( params_gpp%beta * kmm / (1.6 * ns_star))
    chi = xi / (xi + sqrt(vpd))

    ! light use efficiency (m)
    ! consistent with this, directly return light-use-efficiency (m)
    m = ( xi * (ca - gstar) - gstar * sqrt( vpd ) ) / ( xi * (ca + 2.0 * gstar) + 2.0 * gstar * sqrt( vpd ) )

    ! n 
    gamma = gstar / ca
    kappa = kmm / ca
    n = (chi + kappa) / (chi + 2 * gamma)

    ! return derived type
    out_lue%chi=chi
    out_lue%m=m
    out_lue%n=n
      
  end function lue_vpd_c3_simpl


  function lue_vpd_c3_full( kmm, gstar, ns_star, ca, vpd ) result( out_lue )
    !//////////////////////////////////////////////////////////////////
    ! Output:   float, ratio of ci/ca (chi)
    ! Features: Returns an estimate of leaf internal to ambient CO2
    !           partial pressure following the "simple formulation".
    !-----------------------------------------------------------------------

    ! arguments
    real, intent(in) :: kmm       ! Pa, Michaelis-Menten coeff.
    real, intent(in) :: gstar        ! Pa, photores. comp. point (Gamma-star)
    real, intent(in) :: ns_star   ! (unitless) viscosity correction factor for water
    real, intent(in) :: ca        ! Pa, ambient CO2 partial pressure
    real, intent(in) :: vpd       ! Pa, vapor pressure deficit

    ! function return value
    type(outtype_lue) :: out_lue

    ! local variables
    real :: chi                   ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio (unitless)
    real :: xi
    real :: gamma
    real :: kappa
    real :: m
    real :: n

    ! variable substitutes
    real :: vdcg, vacg, vbkg, vsr

    ! beta = 1.6 * ns_star * vpd * (chi * ca - gstar) ** 2.0 / ( (kmm + gstar) * (ca ** 2.0) * (chi - 1.0) ** 2.0 )   ! see Estimation_of_beta.pdf

    ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    xi  = sqrt( ( params_gpp%beta * ( kmm + gstar ) ) / ( 1.6 * ns_star ) )     ! see Eq. 2 in 'Estimation_of_beta.pdf'
    chi = gstar / ca + ( 1.0 - gstar / ca ) * xi / ( xi + sqrt(vpd) )           ! see Eq. 1 in 'Estimation_of_beta.pdf'

    ! consistent with this, directly return light-use-efficiency (m)
    ! see Eq. 13 in 'Simplifying_LUE.pdf'

    ! light use efficiency (m)
    ! m = (ca - gstar)/(ca + 2.0 * gstar + 3.0 * gstar * sqrt( (1.6 * vpd) / (beta * (K + gstar) / ns_star ) ) )

    ! Define variable substitutes:
    vdcg = ca - gstar
    vacg = ca + 2.0 * gstar
    vbkg = params_gpp%beta * (kmm + gstar)

    ! Check for negatives:
    if (vbkg > 0) then
      vsr = sqrt( 1.6 * ns_star * vpd / vbkg )

      ! Based on the m' formulation (see Regressing_LUE.pdf)
      m = vdcg / ( vacg + 3.0 * gstar * vsr )
    end if

    ! n 
    gamma = gstar / ca
    kappa = kmm / ca
    n = (chi + kappa) / (chi + 2 * gamma)

    ! return derived type
    out_lue%chi=chi
    out_lue%m=m
    out_lue%n=n
  
  end function lue_vpd_c3_full


  function lue_c4() result( out_lue )
    !//////////////////////////////////////////////////////////////////
    ! Output:   float, ratio of ci/ca (chi)
    ! Features: Returns an estimate of leaf internal to ambient CO2
    !           partial pressure following the "simple formulation".
    !-----------------------------------------------------------------------
    use md_params_core, only: dummy

    ! function return value
    type(outtype_lue) :: out_lue

    ! return derived type
    out_lue%chi=dummy
    out_lue%m=1
    out_lue%n=1
  
  end function lue_c4


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


  function calc_k( tc, patm ) result( k )
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
    real, parameter :: kc25 = 39.97      ! Pa, assuming 25 deg C & 98.716 kPa
    real, parameter :: ko25 = 2.748e4    ! Pa, assuming 25 deg C & 98.716 kPa
    real, parameter :: dhac = 79430      ! J/mol
    real, parameter :: dhao = 36380      ! J/mol
    real, parameter :: kR   = 8.3145     ! J/mol/K
    real, parameter :: kco  = 2.09476e5  ! ppm, US Standard Atmosphere

    real :: kc, ko, po

    ! function return variable
    real :: k               ! temperature & pressure dependent Michaelis-Menten coefficient, K (Pa).

    kc = kc25 * exp( dhac * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) ) 
    ko = ko25 * exp( dhao * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) ) 

    po     = kco * (1e-6) * patm ! O2 partial pressure
    k      = kc * (1.0 + po/ko)

  end function calc_k


  function calc_gstar( tc ) result( gstar )
    !-----------------------------------------------------------------------
    ! Features: Returns the temperature-dependent photorespiratory 
    !           compensation point, Gamma star (Pascals), based on constants 
    !           derived from Bernacchi et al. (2001) study. Corresponds
    !           to 'calc_gstar_colin' in pmodel.R.
    ! Ref:      Colin's document
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc             ! air temperature (degrees C)

    ! local variables
    real, parameter :: gs25 = 4.220    ! Pa, assuming 25 deg C & 98.716 kPa)
    real, parameter :: kR   = 8.3145   ! J/mol/K
    real, parameter :: dha  = 37830    ! J/mol

    real :: tk                         ! air temperature (Kelvin)

    ! function return variable
    real :: gstar   ! gamma-star (Pa)

    !! conversion to temperature in Kelvin
    tk = tc + 273.15

    gstar = gs25 * exp( ( dha / kR ) * ( 1.0/298.15 - 1.0/tk ) )
    
  end function calc_gstar


  function calc_ftemp_inst_vcmax( tc ) result( fv )
    !-----------------------------------------------------------------------
    ! Output:   Factor fv to correct for instantaneous temperature response
    !           of Vcmax for:
    !
    !               Vcmax(temp) = fv * Vcmax(25 deg C) 
    !
    ! Ref:      Wang Han et al. (in prep.)
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc            ! temperature (degrees C)

    ! function return variable
    real :: fv                        ! temperature response factor, relative to 25 deg C.

    ! loal parameters
    real, parameter :: Ha    = 71513  ! activation energy (J/mol)
    real, parameter :: Hd    = 200000 ! deactivation energy (J/mol)
    real, parameter :: Rgas  = 8.3145 ! universal gas constant (J/mol/K)
    real, parameter :: a_ent = 668.39 ! offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
    real, parameter :: b_ent = 1.07   ! slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
    real, parameter :: tk25  = 298.15 ! 25 deg C in Kelvin

    ! local variables
    real :: tk                        ! temperature (Kelvin)
    real :: dent                      ! entropy change (J/mol)

    ! conversion of temperature to Kelvin
    tk = tc + 273.15

    ! calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
    dent = a_ent - b_ent * tc

    fv = exp( (Ha * (tk - tk25))/(tk * tk25 * Rgas) ) * (1 + exp( (tk25 * dent - Hd)/(Rgas * tk25) ) )/(1 + exp( (tk * dent - Hd)/(Rgas * tk) ) )
    
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
    real, parameter :: tk25  = 298.15 ! 25 deg C in Kelvin

    ! local variables
    real :: tk                  ! temperature (Kelvin)

    ! conversion of temperature to Kelvin
    tk = tc + 273.15

    fr = exp( apar * (tc - 25.0) - bpar * (tc**2 - 25.0**2) )
    
  end function calc_ftemp_inst_rd


  function calc_patm( elv ) result( patm )
    !-----------------------------------------------------------------------
    ! Features: Returns the atmospheric pressure as a function of elevation
    !           and standard atmosphere (1013.25 hPa)
    ! Depends:  - connect_sql
    !           - flux_to_grid
    !           - get_data_point
    !           - get_msvidx
    ! Ref:      Allen et al. (1998)
    !-----------------------------------------------------------------------
    ! argument
    real, intent(in) :: elv           ! elevation above sea level, m

    ! local variables
    real, parameter :: kPo = 101325   ! standard atmosphere, Pa (Allen, 1973)
    real, parameter :: kTo = 298.15   ! base temperature, K (Prentice, unpublished)
    real, parameter :: kL = 0.0065    ! temperature lapse rate, K/m (Allen, 1973)
    real, parameter :: kG = 9.80665   ! gravitational acceleration, m/s**2 (Allen, 1973)
    real, parameter :: kR = 8.3145    ! universal gas constant, J/mol/K (Allen, 1973)
    real, parameter :: kMa = 0.028963 ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)

    ! function return variable
    real :: patm    ! atmospheric pressure at elevation 'elv', Pa 

    ! Convert elevation to pressure, Pa:
    patm = kPo*(1.0 - kL*elv/kTo)**(kG*kMa/(kR*kL))
    
  end function calc_patm


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
    density_h2o = (1e3/vau)

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


  function calc_tempstress( dtemp ) result( ftemp )
    !////////////////////////////////////////////////////////////////
    ! Simple temperature inihibtion function for photosynthesis
    ! 0 at 0 deg C, 1 at 10 deg C
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: dtemp

    ! function return variable
    real :: ftemp

    real, parameter :: temp0 = 0.0

    ! ftemp is a linear ramp down from 1.0 at 12 deg C to 0.0 at 0 deg C
    if (params_gpp%temp_ramp_edge<=0.0) then
      ! no temperature ramp. GPP set to 0 if temp below zero in function pmodel()
      ftemp = 1.0
    else
      ! linear temperature ramp from 0 at 0 deg C to 1.0 at <temp_ramp_edge>
      ftemp = max( 0.0, min( 1.0, (dtemp - temp0) / params_gpp%temp_ramp_edge ) )
    end if

  end function calc_tempstress


  function sigm_gpp_lotemp( dtemp ) result( ftemp )
    !////////////////////////////////////////////////////////////////
    ! Simple temperature inihibtion function for photosynthesis
    ! 0 at 0 deg C, 1 at 10 deg C
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: dtemp

    ! function return variable
    real :: ftemp


    ftemp = max(  0.0, &
                  min(  1.0, &
                        ( ( ( 1.0 / ( 1.0 + exp( -0.7 * dtemp ) ) ) - 0.5 ) * 2.0 ) &
                      ) &
                ) 

  end function sigm_gpp_lotemp


  subroutine initio_gpp()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use md_interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !----------------------------------------------------------------
    ! DAILY OUTPUT
    !----------------------------------------------------------------
    ! GPP
    if (interface%params_siml%loutdgpp) then
      filnam=trim(prefix)//'.d.gpp.out'
      print*,'filnam ', filnam
      open(101,file=filnam,err=888,status='unknown')
    end if 

    ! RD
    if (interface%params_siml%loutdrd) then
      filnam=trim(prefix)//'.d.rd.out'
      open(135,file=filnam,err=888,status='unknown')
    end if 

    ! TRANSPIRATION
    if (interface%params_siml%loutdtransp) then
      filnam=trim(prefix)//'.d.transp.out'
      open(114,file=filnam,err=888,status='unknown')
    end if

    ! xxx: not yet included in writeout_ascii_gpp()
    ! !----------------------------------------------------------------
    ! ! MONTHLY OUTPUT
    ! !----------------------------------------------------------------
    ! ! GPP
    ! if (interface%params_siml%loutdgpp) then     ! monthly and daily output switch are identical
    !   filnam=trim(prefix)//'.m.gpp.out'
    !   open(151,file=filnam,err=888,status='unknown')
    ! end if 

    ! ! RD
    ! if (interface%params_siml%loutdrd) then     ! monthly and daily output switch are identical
    !   filnam=trim(prefix)//'.m.rd.out'
    !   open(152,file=filnam,err=888,status='unknown')
    ! end if 

    ! ! TRANSP
    ! if (interface%params_siml%loutdtransp) then     ! monthly and daily output switch are identical
    !   filnam=trim(prefix)//'.m.transp.out'
    !   open(153,file=filnam,err=888,status='unknown')
    ! end if 

    !----------------------------------------------------------------
    ! ANNUAL OUTPUT
    !----------------------------------------------------------------
    if (interface%params_siml%loutgpp) then

      ! GPP 
      filnam=trim(prefix)//'.a.gpp.out'
      open(310,file=filnam,err=888,status='unknown')

      ! VCMAX (canopy-level, annual maximum) (mol m-2 s-1)
      filnam=trim(prefix)//'.a.vcmax.out'
      open(323,file=filnam,err=888,status='unknown')

      ! 25degC-normalised VCMAX (annual maximum) (mol m-2 s-1)
      filnam=trim(prefix)//'.a.vcmax25.out'
      open(654,file=filnam,err=888,status='unknown')

      ! chi = ci:ca (annual mean, weighted by monthly PPFD) (unitless)
      filnam=trim(prefix)//'.a.chi.out'
      open(652,file=filnam,err=888,status='unknown')

      ! LUE (annual  mean, weighted by monthly PPFD) (unitless)
      filnam=trim(prefix)//'.a.lue.out'
      open(653,file=filnam,err=888,status='unknown')

      ! ci: leaf-internal CO2 partial pressure (Pa)
      filnam=trim(prefix)//'.a.ci.out'
      open(655,file=filnam,err=888,status='unknown')

      ! gs: stomatal conductance
      filnam=trim(prefix)//'.a.gs.out'
      open(656,file=filnam,err=888,status='unknown')

      ! VCMAX (leaf-level, annual mean) (mol m-2 s-1)
      filnam=trim(prefix)//'.a.vcmax_mean.out'
      open(657,file=filnam,err=888,status='unknown')

      ! intrinsic water use efficiency
      filnam=trim(prefix)//'.a.iwue.out'
      open(658,file=filnam,err=888,status='unknown')

    end if

    return

    888  stop 'INITIO_GPP: error opening output files'

  end subroutine initio_gpp


  subroutine initio_nc_gpp()
    !////////////////////////////////////////////////////////////////
    ! Opens NetCDF output files.
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
    character(len=12) :: temp_ramp_edge_char
    character(len=12) :: soilm_par_a_char
    character(len=12) :: soilm_par_b_char

    integer :: jpngr, doy

    write(year_char,999) interface%steering%outyear

    ! convert parameter values to charaters
    write(beta_char,888)           params_gpp%beta
    write(rd_to_vcmax_char,888)    params_gpp%rd_to_vcmax
    write(kphio_char,888)          params_pft_gpp(1)%kphio
    write(temp_ramp_edge_char,888) params_gpp%temp_ramp_edge
    write(soilm_par_a_char,888)    params_gpp%soilm_par_a
    write(soilm_par_b_char,888)    params_gpp%soilm_par_b

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
                          globatt5_nam = "param_temp_ramp_edge", globatt5_val = temp_ramp_edge_char,  &
                          globatt6_nam = "param_soilm_par_a",    globatt6_val = soilm_par_a_char,  &
                          globatt7_nam = "param_soilm_par_a",    globatt7_val = soilm_par_b_char  &
                          )
      end if

      if ( interface%steering%outyear>=interface%params_siml%daily_out_startyr .and. &
        interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then
        !----------------------------------------------------------------
        ! Daily GPP output file 
        !----------------------------------------------------------------
        if (interface%params_siml%lncoutdgpp) then
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
                            globatt5_nam = "param_temp_ramp_edge", globatt5_val = temp_ramp_edge_char,  &
                            globatt6_nam = "param_soilm_par_a",    globatt6_val = soilm_par_a_char,  &
                            globatt7_nam = "param_soilm_par_a",    globatt7_val = soilm_par_b_char  &
                            )
        end if

      end if

    end if

    888  format (F12.6)
    999  format (I4.4)
    
  end subroutine initio_nc_gpp


  subroutine initoutput_gpp( ngridcells )
    !////////////////////////////////////////////////////////////////
    !  Initialises waterbalance-specific output variables
    !----------------------------------------------------------------
    use md_params_core, only: npft, ndayyear, nmonth
    use md_interface

    ! arguments
    integer, intent(in) :: ngridcells

    ! daily
    if ( interface%steering%init.and.interface%params_siml%loutdgpp )   allocate( outdgpp(interface%params_siml%outnt,ngridcells) )
    if (interface%steering%init.and. interface%params_siml%loutdrd    ) allocate( outdrd(interface%params_siml%outnt,ngridcells) )
    if (interface%steering%init.and. interface%params_siml%loutdtransp) allocate( outdtransp(interface%params_siml%outnt,ngridcells) )

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
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
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
          stop 'annual GPP is zero'
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
    !  SR called once a year to gather annual output variables.
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


  subroutine writeout_ascii_gpp()
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use md_interface

    ! local variables
    real :: itime
    integer :: it, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii_gpp: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    !-------------------------------------------------------------------------
    if ( .not. interface%steering%spinup &
         .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
         .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      ! Write daily output only during transient simulation
      do it=1,interface%params_siml%outnt

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(interface%steering%outyear) + real( it - 1 ) * interface%params_siml%outdt / real( ndayyear )
        
        if (interface%params_siml%loutdgpp  )  write(101,999) itime, outdgpp(it,jpngr)
        if (interface%params_siml%loutdrd    ) write(135,999) itime, outdrd(it,jpngr)
        if (interface%params_siml%loutdtransp) write(114,999) itime, outdtransp(it,jpngr)

      end do

    end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutgpp) then
  
      itime = real(interface%steering%outyear)

      write(310,999) itime, sum(outagpp(:,jpngr))
      write(323,999) itime, sum(outavcmax(:,jpngr))
      write(654,999) itime, sum(outavcmax_25(:,jpngr))
      write(653,999) itime, sum(outalue(:,jpngr))
      write(652,999) itime, sum(outachi(:,jpngr))
      write(658,999) itime, sum(outaiwue(:,jpngr)) * 1e6 / 1.6 ! converting from unitless to micro-mol CO2 / mol H2O
      write(655,999) itime, sum(outaci(:,jpngr))
      write(656,999) itime, sum(outags(:,jpngr))
      write(657,999) itime, sum(outavcmax_leaf(:,jpngr))

    end if

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_gpp


  subroutine writeout_nc_gpp()
    !/////////////////////////////////////////////////////////////////////////
    ! Write NetCDF output
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
        if (interface%params_siml%lncoutdgpp) print*,'writing ', trim(ncoutfilnam_dgpp), '...'
        if (interface%params_siml%lncoutdgpp) call write_nc_3D_time( trim(ncoutfilnam_dgpp), &
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
