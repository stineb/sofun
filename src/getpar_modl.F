subroutine getpar_modl
  !////////////////////////////////////////////////////////////////
  !  Subroutine reads model parameters from input file.
  !  It was necessary to separate this SR from module params_modl
  !  because this SR uses module waterbal_stash, which also uses
  !  params_modl.
  !----------------------------------------------------------------
  use params_modl
  use waterbal_stash, only: getpar_modl_waterbal
  use nuptake_fun, only: getpar_modl_nuptake
  use treegeometry, only: get_geometry_sapl
  use phenology_lpx, only: getpar_pft_summerphenology

  implicit none

  ! local variables
  integer :: pft
  integer, parameter :: npar = 23
  real, dimension(npar) :: params_array

  integer, parameter :: npar_pft = 21
  real, dimension(npar_pft) :: pft_params_array

  real, dimension(npft) :: lu_category_prov    ! land use category associated with PFT (provisional)
  real, dimension(3) :: out_get_geometry_sapl

  !----------------------------------------------------------------
  ! NON-PFT DEPENDENT PARAMETERS
  ! read parameter input file and store values in single array
  !----------------------------------------------------------------
  open(unit=01,file='params_modl.dat',status='OLD')      
  read (01,*) params_array
  close (01)

  ! above-ground fast (leaf) litter decay rate [1/d]
  klitt_af10 = params_array(1)/ndayyear

  ! above-ground slow (woody) litter decay rate [1/d] 
  klitt_as10 = params_array(2)/ndayyear

  ! below-ground (root) litter decay rate [1/d] 
  klitt_bg10 = params_array(3)/ndayyear

  ! exudates decay rate [1/d]
  kexu10 = params_array(4)/ndayyear

  ! fast soil pool decay rate [1/d]
  ksoil_fs10 = params_array(5)/ndayyear

  ! slow soil pool decay rate [1/d]
  ksoil_sl10 = params_array(6)/ndayyear

  ! factor for "Manzoni Equation" (XPXXX) [1]
  ntoc_crit1 = params_array(7)

  ! exponent for "Manzoni Equation" (XPXXX) [1]
  ntoc_crit2 = params_array(8)

  ! C:N ratio of microbial biomass [1]
  cton_microb = params_array(9)

  ! fraction of litter input to fast soil pool [1]
  fastfrac = params_array(10)

  ! constant ratio between leaf and sapwood C/N ratios (Friend et al.,2008 in XP08)
  rcton_LS = params_array(11)

  ! constant ratio between leaf and root C/N ratios (Friend et al., 2008 in XP08)
  rcton_LR = params_array(12)

  ! maximum nitrification rate
  MAXNITR = params_array(13)

  ! maximum NO from nitrification (day-1)
  NON = params_array(14)

  ! maximum N2O from nitrification (day-1)
  N2ON = params_array(15)

  ! Michaelis-Menten coefficient [gN/m2]. Use this value if soil represents top 100 cm 
  KN = params_array(16)

  ! Michaelis-Menten coefficient [gC/m2]. Use this value if soil represents top 100 cm 
  KDOC = params_array(17)

  ! DOCMAX
  DOCMAX = params_array(18)

  ! Fraction of denitrification lost as N2O. Range of possible values: 0.002 - 0.047 (Xu-Ri and Prentice, 2008)
  DNITR2N2O = params_array(19)

  ! canopy light extinction coefficient
  kbeer = params_array(20)

  ! additional respiration coefficient, unknown origin ("proposed by Annet Wolf 2008")
  kresp = params_array(21)

  ! fraction of N retained at leaf abscission 
  F_NRETAIN = params_array(22)        

  ! maximum fractional plant coverage of trees (sum of all tree PFTs)
  fpc_tree_max = params_array(23)

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC PARAMETERS
  !----------------------------------------------------------------
  call getpar_modl_waterbal('params_waterbal_stash.dat')
  call getpar_modl_nuptake('params_nuptake_fun.dat')


  !----------------------------------------------------------------
  ! PFT DEPENDENT PARAMETERS
  ! read parameter input file and store values in single array
  !----------------------------------------------------------------
  open(unit=02,file='params_pft.dat',status='OLD')      
  read (02,*) pft_params_array
  close (02)

  do pft=1,npft 

    ! C:N ratio of soil organic matter [1]
    cton_soil(pft) = pft_params_array(1) 

    ! C:N ratio of new production
    cton_pro(pft) = pft_params_array(6) 

    ! land use category associated with PFT (provisional)
    lu_category_prov(pft) = pft_params_array(2)
    if (lu_category_prov(pft)==1.0) then
      lu_category(pft) = lunat
      islu(pft,lunat) = .true.
    else
      islu(pft,lunat) = .false.
    end if

    ! leaf type: broadleaved (1), needleleaved (2), grass (3) or moss (4)
    if (pft_params_array(3)==3.0) then
      grass(pft) = .true.
      tree(pft)  = .false.
    else
      tree(pft) = .true.
      grass(pft) = .false.
    end if 

    ! sapling (or grass on initialisation) LAI
    lai0(pft) = pft_params_array(4) 

    ! leaf longevity [days], read in as [years]
    long_leaf(pft) = pft_params_array(5) * ndayyear 

    ! specific leaf area [m2/gC]
    ! defined as a function of 'long_leaf'
    ! Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
    ! Equation based on Reich et al 1997, Fig 1f:
    ! This is a mathematically identical expression as in Sitch et al., 2003, Eq. 6
    sla(pft) = 2.0e-4 * exp( 6.15 - 0.46 * log( long_leaf(pft) * 12.0 ) ) 

    ! PFT identification code (same as in LPX)
    ! 3 = Temperate Needle Evergreen   TeNE
    ! 8 = C3 grass                     TeH
    pftcode(pft) = int(pft_params_array(7)) 

    ! sapwood turnover period [days], read in as [years]
    long_sapw(pft) = pft_params_array(8) * ndayyear 

    ! root turnover period [days], read in as [years]
    long_root(pft) = pft_params_array(9) * ndayyear 

    ! maintenance respiration coefficient
    rresp_maint(pft) = pft_params_array(10) 

    ! N-fixing plant? (0=false, 1=true)
    if (pft_params_array(11)==0.) then
      nfixer(pft) = .false.
    else
      nfixer(pft) = .true.
    endif 

    ! growth respiration coefficient
    rresp_growth(pft) = pft_params_array(12) 

    ! Reinicker-p for geometry
    reinickerp = pft_params_array(13) ! put into pft-loop to have all geometry parameters in one place 

    ! wood density (gC * m-3)
    wooddens(pft) = pft_params_array(14) 

    ! allometry parameters
    allom1(pft) = pft_params_array(15)
    allom2(pft) = pft_params_array(16)
    allom3 = pft_params_array(17) ! put into pft-loop to have all geometry parameters in one place 

    ! maximum crown area
    crownarea_max = pft_params_array(18) ! put into pft-loop to have all geometry parameters in one place 

    ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
    latosa(pft) = pft_params_array(19)

    ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    woodtosapw(pft) = pft_params_array(20)

    ! leaf to root ratio under non-water stressed conditionss
    lmtorm0(pft) = pft_params_array(21)

    ! leaf (, sapwood, and heartwood) mass of sapling
    if (grass(pft)) then
      lm_sapl(pft)%c%c12 = lai0(pft) / sla(pft)
    else
      ! call function returning (lm_sapl, sm_sapl, hm_sapl)
      out_get_geometry_sapl = get_geometry_sapl( &
          lai0(pft), sla(pft), latosa(pft), reinickerp, allom1(pft) &
        , allom2(pft), allom3, wooddens(pft), woodtosapw(pft), lmtorm0(pft) &
        )            
      lm_sapl(pft)%c%c12 = out_get_geometry_sapl(1)
      sm_sapl(pft)%c%c12 = out_get_geometry_sapl(2)
      hm_sapl(pft)%c%c12 = out_get_geometry_sapl(3)
    end if
    
    ! root mass of sapling is a function of leaf mass            
    rm_sapl(pft)%c%c12 = ( 1.0 / lmtorm0(pft) ) * lm_sapl(pft)%c%c12

    ! N content of sapling
    lm_sapl(pft)%n%n14 = lm_sapl(pft)%c%c12 / cton_pro(pft)
    rm_sapl(pft)%n%n14 = rm_sapl(pft)%c%c12 / cton_pro(pft) 
    if (tree(pft)) then
      sm_sapl(pft)%n%n14 = sm_sapl(pft)%c%c12 / cton_pro(pft)
      hm_sapl(pft)%n%n14 = hm_sapl(pft)%c%c12 / cton_pro(pft)
    endif

  end do

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC PARAMETERS
  !----------------------------------------------------------------
  call getpar_pft_summerphenology('params_pft_summerphenology.dat')


end subroutine getpar_modl
