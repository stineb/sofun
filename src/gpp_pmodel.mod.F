module gpp_pmodel
  !////////////////////////////////////////////////////////////////
  ! P-MODEL GPP MODULE
  ! Contains P model functions adopted from GePiSaT
  !----------------------------------------------------------------
  !
  ! gpp_pmodel
  ! __version__ 1.0
  !
  ! written by Benjamin Stocker, adopted from Python code written by Tyler Davis
  ! Imperial College London
  !
  ! 2015-03-21 -- created
  ! 2015-03-21 -- last updated
  !
  ! ------------
  ! description:
  ! ------------
  ! This collection of functions implements the Global ecosystem
  ! in Space and Time (GePiSaT) model of the terrestrial biosphere. 
  !
  ! GePiSaT takes a simplistic approach to modelling terrestrial gross primary 
  ! production (GPP) by making  use of the optimality principle of vegetation 
  ! minimizing the summed costs associated with maintaining carbon fixation and 
  ! water transport capabilities (Prentice et al., 2014).
  !
  ! ----------
  ! changelog:
  ! ----------
  ! VERSION 1.0.0
  !  - file created [2015-03-21]
  !
  ! VERSION 1.0.1
  !  - 'lue_approx' function now based on SI units.
  !
  !
  ! This module ontains the "main" subroutine 'pmodel' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'pmodel' must contain this list 
  ! of subroutines (names that way).
  !   - getpar_modl_pmodel
  !   - initio_pmodel
  !   - initoutput_pmodel
  !   - getout_daily_pmodel
  !   - getout_monthly_pmodel
  !   - writeout_ascii_pmodel
  !   - pmodel
  ! Required module-independent model state variables (necessarily 
  ! updated by 'pmodel') are:
  !   - xxx
  !----------------------------------------------------------------
  implicit none
  ! private
  ! public gpp, outtype_pmodel

  ! pmodel PARAMETERS
  ! hard-wired here because well known (xxx some are not well known => read these from file xxx)
  real, parameter :: kphio = 0.093      ! quantum efficiency (Long et al., 1993)
  real, parameter :: kPo   = 101325.0   ! standard atmosphere, Pa (Allen, 1973)
  real, parameter :: kTo   = 25.0       ! base temperature, deg C (Prentice, unpublished)

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
  real, parameter :: beta = 244.033
  
  !-----------------------------------------------------------------------
  ! Metabolic N ratio (N per unit Vcmax)
  ! Reference: Harrison et al., 2009, Plant, Cell and Environment; Eq. 3
  !-----------------------------------------------------------------------
  real, parameter :: mol_weight_rubisco = 5.5e5    ! molecular weight of Rubisco, (g R)(mol R)-1
  real, parameter :: n_conc_rubisco     = 1.14e-2  ! N concentration in rubisco, (mol N)(g R)-1
  real, parameter :: cat_turnover_per_site = 2.33  ! catalytic turnover rate per site at 25 deg C, (mol CO2)(mol R sites)-1; use 2.33 instead of (3.5) as not all Rubisco is active (see Harrison et al., 2009)  
  real, parameter :: cat_sites_per_mol_R   = 8.0   ! number of catalytic sites per mol R, (mol R sites)(mol R)-1

  ! Metabolic N ratio (mol N s (mol CO2)-1 )
  real, parameter :: n_v = mol_weight_rubisco * n_conc_rubisco / ( cat_turnover_per_site * cat_sites_per_mol_R )

  ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  real, parameter :: rd_to_vcmax = 0.015
  
  ! parameters for Narea -- under construction
  ! sla = 0.0014       ! specific leaf area (m2/gC)

  ! N in cell walls: Slope of WN~LMA is 0.0002 mol N / g leaf mass (Hikosaka&Shigeno, 2009)
  ! With 0.5 g C / g leaf mass and 14 g N / mol N: n_cw = 0.0056 g N / g C

  ! real, parameter :: ncw = 0.0056          ! N:C ratio in cell walls, working hypothesis: leaf N is solely determined by Vcmax25
  ! n_v  = 1.0/40.96    ! gN µmol-1 s-1. Value 40.96 is 'sv' in Table 2 in Kattge et al., 2009, GCB, C3 herbaceous
  ! -- under construction


  ! MODULE-SPECIFIC VARIABLES
  ! These are not required outside module 'pmodel', but are used in different SRs of this module
  ! real :: xxx                    ! xxx

  ! DAILY OUTPUT VARIABLES
  ! xxx add jpngr dim


  ! MONTHLY OUTPUT VARIABLES
  ! xxx add jpngr dim

  ! Function return variables as derived types
  type outtype_pmodel
    real :: gpp
    real :: vcmax
    real :: vcmax25
    real :: vcmax_unitfapar
    real :: vcmax_unitiabs
    real :: factor25_vcmax
    real :: rd
    real :: rd_unitfapar 
    real :: rd_unitiabs 
    real :: actnv 
    real :: actnv_unitfapar 
    real :: actnv_unitiabs 
    real :: lue
  end type outtype_pmodel

  type outtype_lue
    real :: chi
    real :: m
    real :: n
  end type outtype_lue

contains

  subroutine gpp( jpngr, fpar, ppfd, co2, tc, vpd, elv, method )
    !//////////////////////////////////////////////////////////////////
    ! Output:   gpp (gC/m2/d)   : gross primary production (xxx month -> day)
    !------------------------------------------------------------------
    use fluxes, only: dgpp

    ! arguments
    integer, intent(in) :: jpngr    ! gridcell number
    real, intent(in)    :: fpar     ! monthly fraction of absorbed photosynthetically active radiation (unitless) 
    real, intent(in)    :: ppfd     ! photon flux density (xxx)
    real, intent(in)    :: co2      ! atmospheric CO2 (ppm)
    real, intent(in)    :: tc       ! air temperature (deg C)
    real, intent(in)    :: vpd      ! vapour pressure deficit (Pa)
    real, intent(in)    :: elv      ! elevation above sea level (m)
    character(len=*), intent(in), optional :: method

    ! local variables
    type(outtype_pmodel) :: out_pmodel
    character(len=10) :: mymethod

    if (present(method)) then
      mymethod = method
    else  
      mymethod = "full"
    end if

    if (trim(mymethod)=="prescr") then
      !----------------------------------------------------------------
      ! PRESCRIBED GPP
      !----------------------------------------------------------------
      dgpp(:) = dgpp_data(:)

    else
      !----------------------------------------------------------------
      ! PREDICTED GPP FROM CO2, FAPAR
      !----------------------------------------------------------------
      out_pmodel   = pmodel( fpar, ppfd, co2, tc, cpalpha, vpd, elv, trim(mymethod) )


      out_pmodel   = pmodel( fpar, mppfd(moy), co2, mtemp(moy), cpalpha, mvpd(moy), elv, "full" )
      out_gpp(moy) = out_pmodel%gpp
      out_rd(moy)  = out_pmodel%rd



    end if

  end subroutine gpp


  function pmodel( fpar, ppfd, co2, tc, cpalpha, vpd, elv, method ) result( out_pmodel )
    !//////////////////////////////////////////////////////////////////
    ! Output:   gpp (mol/m2/month)   : gross primary production
    !------------------------------------------------------------------

    ! arguments
    real, intent(in) :: fpar         ! monthly fraction of absorbed photosynthetically active radiation (unitless) 
    real, intent(in) :: ppfd         ! monthly photon flux density (mol/m2)
    real, intent(in) :: co2          ! atmospheric CO2 concentration (ppm)
    real, intent(in) :: tc           ! monthly air temperature (deg C)
    real, intent(in) :: cpalpha      ! monthly Cramer-Prentice-alpha (unitless, within [0,1.26]) 
    real, intent(in) :: vpd          ! mean monthly vapor pressure (Pa) -- CRU data is in hPa
    real, intent(in) :: elv          ! elevation above sea-level (m)
    character(len=*), intent(in) :: method

    ! local variables
    real :: iabs                     ! absorbed photosynthetically active radiation (mol/m2)
    real :: patm                     ! atmospheric pressure as a function of elevation (Pa)
    real :: ca                       ! ambient CO2 partial pression (Pa)
    real :: gs                       ! photorespiratory compensation point - Gamma-star (Pa)
    real :: fa                       ! function of alpha to reduce GPP in strongly water-stressed months (unitless)
    real :: kmm                      ! Michaelis-Menten coefficient (Pa)
    real :: ns                       ! viscosity of H2O at ambient temperatures (Pa s)
    real :: ns25                     ! viscosity of H2O at 25 deg C (Pa s)
    real :: ns_star                  ! viscosity correction factor (unitless)
    real :: m, n

    real :: gpp                      ! assimilation (mol m-2 s-1)
    real :: lue                      ! Light use efficiency
    real :: vcmax                    ! Vcmax per unit ground area (mol m-2 s-1)
    real :: vcmax_unitfapar
    real :: vcmax_unitiabs
    real :: vcmax25                  ! Vcmax25 (vcmax normalized to 25 deg C)
    real :: vcmax25_unitfapar
    real :: vcmax25_unitiabs
    real :: rd                       ! Dark respiration (mol m-2 s-1)
    real :: rd_unitfapar 
    real :: rd_unitiabs 
    real :: factor25_vcmax           ! correction factor to normalise Vcmax to 25 deg C
    real :: actnv 
    real :: actnv_unitfapar 
    real :: actnv_unitiabs 


    ! function return value
    type(outtype_pmodel) :: out_pmodel

    ! local variables
    type(outtype_lue) :: out_lue


    ! absorbed photosynthetically active radiation (mol/m2)
    iabs = fpar * ppfd

    ! atmospheric pressure as a function of elevation (Pa)
    patm = calc_patm( elv )

    ! ambient CO2 partial pression (Pa)
    ca   = co2_to_ca( co2, patm )

    ! photorespiratory compensation point - Gamma-star (Pa)
    gs   = calc_gstar( tc )

    ! function of alpha to reduce GPP in strongly water-stressed months (unitless)
    fa   = calc_fa( cpalpha )

    ! Michaelis-Menten coef. (Pa)
    kmm  = calc_k( tc, patm )

    ! viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
    ns      = calc_viscosity_h2o( tc, patm )  ! Pa s 
    ns25    = calc_viscosity_h2o( kTo, kPo )  ! Pa s 
    ns_star = ns / ns25                       ! (unitless)

    select case (method)

      case ("approx")
        !-----------------------------------------------------------------------
        ! A. APPROXIMATIVE METHOD
        !-----------------------------------------------------------------------
        out_lue = lue_approx( tc, vpd, elv, ca, gs )

      case ("simpl")
        !-----------------------------------------------------------------------
        ! B.1 SIMPLIFIED FORMULATION 
        !-----------------------------------------------------------------------
        out_lue = lue_vpd_simpl( kmm, gs, ns, ca, vpd, beta  )

      case ("full")
        !-----------------------------------------------------------------------
        ! B.2 FULL FORMULATION
        !-----------------------------------------------------------------------
        out_lue = lue_vpd_full( kmm, gs, ns_star, ca, vpd, beta  )

      case default

        stop 'PMODEL: select valid method'

    end select

    ! LUE-functions return m, n, and chi
    m = out_lue%m
    n = out_lue%n

    !-----------------------------------------------------------------------
    ! Calculate function return variables
    !-----------------------------------------------------------------------

    ! GPP per unit ground area is the product of the intrinsic quantum 
    ! efficiency, the absorbed PAR, the function of alpha (drought-reduction),
    ! and 'm'
    m   = calc_mprime( m )
    gpp = iabs * kphio * fa * m

    ! Light use efficiency
    lue = kphio * fa * m 

    ! Vcmax per unit ground area is the product of the intrinsic quantum 
    ! efficiency, the absorbed PAR, and 'n'
    vcmax = iabs * kphio * n

    ! Vcmax normalised per unit fAPAR (assuming fAPAR=1)
    vcmax_unitfapar = ppfd * kphio * n 

    ! Vcmax normalised per unit absorbed PPFD (assuming iabs=1)
    vcmax_unitiabs = kphio * n 

    ! Vcmax25 (vcmax normalized to 25 deg C)
    factor25_vcmax    = calc_vcmax25( 1.0, tc )
    vcmax25           = factor25_vcmax * vcmax
    vcmax25_unitfapar = factor25_vcmax * vcmax_unitfapar
    vcmax25_unitiabs  = factor25_vcmax * vcmax_unitiabs

    ! Dark respiration
    rd = rd_to_vcmax * vcmax

    ! Dark respiration per unit fAPAR (assuming fAPAR=1)
    rd_unitfapar = rd_to_vcmax * vcmax_unitfapar

    ! Dark respiration per unit absorbed PPFD (assuming iabs=1)
    rd_unitiabs = rd_to_vcmax * vcmax_unitiabs

    ! active metabolic leaf N (canopy-level), mol N/m2-ground (same equations as for nitrogen content per unit leaf area, gN/m2-leaf)
    actnv = vcmax25 * n_v
    actnv_unitfapar = vcmax25_unitfapar * n_v
    actnv_unitiabs  = vcmax25_unitiabs  * n_v


    ! Construct derived type for output
    out_pmodel%gpp             = gpp
    out_pmodel%vcmax           = vcmax
    out_pmodel%vcmax25         = vcmax25
    out_pmodel%vcmax_unitfapar = vcmax_unitfapar
    out_pmodel%vcmax_unitiabs  = vcmax_unitiabs
    out_pmodel%factor25_vcmax  = factor25_vcmax
    out_pmodel%rd              = rd
    out_pmodel%rd_unitfapar    = rd_unitfapar 
    out_pmodel%rd_unitiabs     = rd_unitiabs 
    out_pmodel%actnv           = actnv 
    out_pmodel%actnv_unitfapar = actnv_unitfapar 
    out_pmodel%actnv_unitiabs  = actnv_unitiabs 
    out_pmodel%lue             = lue
  
  end function pmodel


  function calc_n_rubisco_area( vcmax25 ) result( n_rubisco_area )
    !-----------------------------------------------------------------------
    ! Returns Rubisco N content per unit leaf area for a given Vcmax.
    ! Reference: Harrison et al., 2009, Plant, Cell and Environment; Eq. 3
    !-----------------------------------------------------------------------

    ! argument
    real, intent(in) :: vcmax25                      ! leaf level Vcmax  at 25 deg C, (mol CO2) m-2 s-1

    ! function retrurn value
    real, intent(out) :: n_rubisco_area              ! Rubisco N content per unit leaf area, (g N)(m-2 leaf)

    ! local variables
    real :: n_v                                      ! Rubisco N per unit Vcmax (xxx units xxx)

    real, parameter :: mol_weight_rubisco = 5.5e5    ! molecular weight of Rubisco, (g R)(mol R)-1
    real, parameter :: n_conc_rubisco     = 1.14e-2  ! N concentration in rubisco, (mol N)(g R)-1
    real, parameter :: mol_weight_n       = 14.0067  ! molecular weight of N, (g N)(mol N)-1
    real, parameter :: cat_turnover_per_site = 3.5   ! catalytic turnover rate per site at 25 deg C, (mol CO2)(mol R sites)-1
    real, parameter :: cat_sites_per_mol_R   = 8.0   ! number of catalytic sites per mol R, (mol R sites)(mol R)-1

    ! Metabolic N ratio
    n_v = mol_weight_rubisco * n_conc_rubisco * mol_weight_n / ( cat_turnover_per_site * cat_sites_per_mol_R )

    n_rubisco_area = vcmax25 * n_v

  end function calc_n_rubisco_area


  function lue_approx( temp, vpd, elv, ca, gs ) result( out_lue )
    !//////////////////////////////////////////////////////////////////
    ! Output:   list: 'm' (unitless), 'chi' (unitless)
    ! Returns list containing light use efficiency (m) and ci/ci ratio 
    ! (chi) based on the approximation of the theoretical relationships
    ! of chi with temp, vpd, and elevation. Is now based on SI units as 
    ! inputs.
    !------------------------------------------------------------------

    ! arguments
    real, intent(in) :: temp   ! deg C, air temperature
    real, intent(in) :: vpd    ! Pa, vapour pressure deficit
    real, intent(in) :: elv    ! m, elevation above sea level
    real, intent(in) :: ca     ! Pa, ambient CO2 partial pressure
    real, intent(in) :: gs     ! Pa, photores. comp. point (Gamma-star)

    ! function return value
    type(outtype_lue) :: out_lue

    ! local variables
    real :: whe                ! value of "Wang-Han Equation"
    real :: gamma
    real :: chi                ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio (unitless)
    real :: m

    ! Wang-Han Equation
    whe = exp( &
      4.644 &
      + 0.0545 * ( temp - 25.0 ) &
      - 0.5 * log( vpd ) &   ! convert vpd from Pa to kPa 
      - 8.15e-5 * elv &      ! convert elv from m to km
      )

    ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    chi = whe / ( 1.0 + whe )

    !  m
    gamma = gs / ca
    m = (chi - gamma) / (chi + 2 * gamma)

    ! return derived type
    out_lue%chi=chi
    out_lue%m=m
    out_lue%n=-9999
  
  end function lue_approx


  function lue_vpd_simpl( kmm, gs, ns_star, ca, vpd, beta ) result( out_lue )
    !//////////////////////////////////////////////////////////////////
    ! Output:   float, ratio of ci/ca (chi)
    ! Returns an estimate of leaf internal to ambient CO2
    ! partial pressure following the "simple formulation".
    !-----------------------------------------------------------------------

    ! arguments
    real, intent(in) :: kmm       ! Pa, Michaelis-Menten coeff.
    real, intent(in) :: gs        ! Pa, photores. comp. point (Gamma-star)
    real, intent(in) :: ns_star   ! (unitless) viscosity correction factor for water
    real, intent(in) :: ca        ! Pa, ambient CO2 partial pressure
    real, intent(in) :: vpd       ! Pa, vapor pressure deficit
    real, intent(in) :: beta      ! unit cost xxx

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
    xi  = sqrt( beta * kmm / (1.6 * ns_star))
    chi = xi / (xi + sqrt(vpd))

    ! light use efficiency (m)
    ! consistent with this, directly return light-use-efficiency (m)
    m = ( xi * (ca - gs) - gs * sqrt( vpd ) ) / ( xi * (ca + 2.0 * gs) + 2.0 * gs * sqrt( vpd ) )

    ! n 
    gamma = gs / ca
    kappa = kmm / ca
    n = (chi + kappa) / (chi + 2 * gamma)

    ! return derived type
    out_lue%chi=chi
    out_lue%m=m
    out_lue%n=n
      
  end function lue_vpd_simpl


  function lue_vpd_full( kmm, gs, ns_star, ca, vpd, beta ) result( out_lue )
    !//////////////////////////////////////////////////////////////////
    ! Output:   float, ratio of ci/ca (chi)
    ! Features: Returns an estimate of leaf internal to ambient CO2
    !           partial pressure following the "simple formulation".
    !-----------------------------------------------------------------------

    ! arguments
    real, intent(in) :: kmm       ! Pa, Michaelis-Menten coeff.
    real, intent(in) :: gs        ! Pa, photores. comp. point (Gamma-star)
    real, intent(in) :: ns_star   ! (unitless) viscosity correction factor for water
    real, intent(in) :: ca        ! Pa, ambient CO2 partial pressure
    real, intent(in) :: vpd       ! Pa, vapor pressure deficit
    real, intent(in) :: beta      ! unit cost xxx

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

    ! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    xi  = sqrt( (beta * ( kmm + gs ) ) / ( 1.6 * ns_star ) )
    chi = gs / ca + ( 1.0 - gs / ca ) * xi / ( xi + sqrt(vpd) )

    ! consistent with this, directly return light-use-efficiency (m)
    ! see Eq. 13 in 'Simplifying_LUE.pdf'

    ! light use efficiency (m)
    ! m = (ca - gs)/(ca + 2.0 * gs + 3.0 * gs * sqrt( (1.6 * vpd) / (beta * (K + gs) / ns_star ) ) )

    ! Define variable substitutes:
    vdcg = ca - gs
    vacg = ca + 2.0 * gs
    vbkg = beta * (kmm + gs)

    ! Check for negatives:
    if (vbkg > 0) then
      vsr = sqrt( 1.6 * ns_star * vpd / vbkg )

      ! Based on the m' formulation (see Regressing_LUE.pdf)
      m = vdcg / ( vacg + 3.0 * gs * vsr )
    end if

    ! n 
    gamma = gs / ca
    kappa = kmm / ca
    n = (chi + kappa) / (chi + 2 * gamma)

    ! return derived type
    out_lue%chi=chi
    out_lue%m=m
    out_lue%n=n
  
  end function lue_vpd_full


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
    real, intent(out) :: mprime

    ! square of m-prime (mpi)
    mprime = m**2 - kc**(2.0/3.0) * (m**(4.0/3.0))

    ! Check for negatives and take root of square
    if (mprime > 0) mprime = sqrt(mprime) 
    
  end function calc_mprime


  function calc_fa( cpalpha ) result( fa )
    !-----------------------------------------------------------------------
    ! Input:  cpalpha (unitless, within [0,1.26]): monthly Cramer-Prentice-alpha
    ! Output: fa (unitless, within [0,1]): function of alpha to reduce GPP 
    !                                      in strongly water-stressed months
    !-----------------------------------------------------------------------
    ! argument
    real, intent(in) :: cpalpha

    ! function return variable
    real, intent(out) :: fa

    fa = ( cpalpha / 1.26 )**(0.25)
    
  end function calc_fa


  function co2_to_ca( co2, patm ) result( ca )
    !-----------------------------------------------------------------------
    ! Output:   - ca in units of Pa
    ! Features: Converts ca (ambient CO2) from ppm to Pa.
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: co2     ! ambient CO2 in units of ppm
    real, intent(in) :: patm    ! monthly atm. pressure, Pa

    ! function return variable
    real, intent(out) :: ca ! ambient CO2 in units of Pa

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
    real, intent(out) :: co2

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
    real, intent(out) :: k               ! temperature & pressure dependent Michaelis-Menten coefficient, K (Pa).

    kc = kc25 * exp( dhac * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) ) 
    ko = ko25 * exp( dhao * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) ) 

    po     = kco * (1e-6) * patm ! O2 partial pressure
    k      = kc * (1.0 + po/ko)

  end function calc_k


  function calc_gstar( tc ) result( gstar )
    !-----------------------------------------------------------------------
    ! Features: Returns the temperature-dependent photorespiratory 
    !           compensation point, Gamma star (Pascals), based on constants 
    !           derived from Bernacchi et al. (2001) study.
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
    real, intent(out) :: gstar   ! gamma-star (Pa)

    !! conversion to temperature in Kelvin
    tk = tc + 273.15

    gstar = gs25 * exp( ( dha / kR ) * ( 1.0/298.15 - 1.0/tk ) )
    
  end function calc_gstar


  function calc_vcmax25( vcmax, tc ) result( vcmax25 )
    !-----------------------------------------------------------------------
    ! Output:   vcmax25  : Vcmax at 25 deg C
    ! Features: Returns the temperature-corrected Vcmax at 25 deg C
    ! Ref:      Analogue function like 'calc_gstar_gepisat'
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: vcmax   ! Vcmax at a given temperature tc 
    real, intent(in) :: tc      ! air temperature (degrees C)

    ! loal variables
    real, parameter :: dhav = 65330    ! J/mol
    real, parameter :: kR   = 8.3145   ! J/mol/K

    real :: tk                         ! air temperature (Kelvin)

    ! function return variable
    real, intent(out) :: vcmax25  ! Vcmax at 25 deg C 

    !! conversion to temperature in Kelvin
    tk = tc + 273.15

    vcmax25 = vcmax * exp( -dhav/kR * (1/298.15 - 1/tk) )
    
  end function calc_vcmax25


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
    real, parameter :: kR = 8.3143    ! universal gas constant, J/mol/K (Allen, 1973)
    real, parameter :: kMa = 0.028963 ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)

    ! function return variable
    real, intent(out) :: patm    ! atmospheric pressure at elevation 'elv', Pa 

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
    real, intent(out) :: density_h2o  ! density of water, kg/m**3

    ! Calculate lambda, (bar cm**3)/g:
    my_lambda = 1788.316 + &
            21.55053*tc + &
          -0.4695911*tc*tc + &
       (3.096363e-3)*tc*tc*tc + &
      -(7.341182e-6)*tc*tc*tc*tc

    ! Calculate po, bar
    po = 5918.499 + & 
             58.05267*tc + & 
           -1.1253317*tc*tc + & 
       (6.6123869e-3)*tc*tc*tc + & 
      -(1.4661625e-5)*tc*tc*tc*tc

    ! Calculate vinf, cm**3/g
    vinf = 0.6980547 + &
      -(7.435626e-4)*tc + &
       (3.704258e-5)*tc*tc + &
      -(6.315724e-7)*tc*tc*tc + &
       (9.829576e-9)*tc*tc*tc*tc + &
     -(1.197269e-10)*tc*tc*tc*tc*tc + &
      (1.005461e-12)*tc*tc*tc*tc*tc*tc + &
     -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc + &
       (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc + &
     -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc

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
    real, tbar, tbarx, tbar2, tbar3, rbar, mu0, mu1, ctbar, mu_bar, &
      coef1, coef2
    integer :: i, j                         ! counter variables

    ! function return variable
    real :: viscosity_h2o

    ! Get the density of water, kg/m**3
    rho = calc_density_h2o(tc, patm)

    ! Calculate dimensionless parameters:
    tbar = (tc + 273.15)/tk_ast
    tbarx = tbar**(0.5)
    tbar2 = tbar**2
    tbar3 = tbar**3
    rbar = rho/rho_ast

    ! Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
    mu0 = 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
    mu0 = 1e2*tbarx/mu0

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
        coef2 = coef2 + h_array(j,i) * (rbar - 1.0)**(j-1)
      end do
      mu1 = mu1 + coef1 * coef2    
    end do
    mu1 = exp( rbar * mu1 )

    ! Calculate mu_bar (Eq. 2, Huber et al., 2009)
    !   assumes mu2 = 1
    mu_bar = mu0 * mu1

    ! Calculate mu (Eq. 1, Huber et al., 2009)
    viscosity_h2o = mu_bar * mu_ast    ! Pa s

  end function calc_viscosity_h2o

end module gpp_pmodel
