module params_modl
  !////////////////////////////////////////////////////////////////
  !  Module contains (constrainable) model parameters.
  !  Model parameters adopted here are from LPX C3 grass PFT
  !  Litter and soil turnover parameters are divided by 365 to 
  !  convert from [1/yr] to [1/d].
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  use params_core
  use classdefs

  implicit none

  ! NON PFT-DEPENDENT PARAMETERS
  real :: klitt_af10        !  1 above-ground fast (leaf) litter decay rate [1/d]
  real :: klitt_as10        !  2 above-ground slow (woody) litter decay rate [1/d] 
  real :: klitt_bg10        !  3 below-ground (root) litter decay rate [1/d] 
  real :: kexu10            !  4 exudates decay rate [1/d]
  real :: ksoil_fs10        !  5 fast soil pool decay rate [1/d]
  real :: ksoil_sl10        !  6 slow soil pool decay rate [1/d]
  real :: ntoc_crit1        !  7 factor for "Manzoni Equation" (XPXXX) [1]
  real :: ntoc_crit2        !  8 exponent for "Manzoni Equation" (XPXXX) [1]
  real :: cton_microb       !  9 C:N ratio of microbial biomass [1]
  real :: fastfrac          ! 10 fraction of litter input to fast soil pool [1]
  real :: rcton_LS          ! 11 constant ratio between leaf and sapwood C/N ratios (Friend et al.,2008 in XP08)
  real :: rcton_LR          ! 12 constant ratio between leaf and root C/N ratios (Friend et al., 2008 in XP08)
  real :: MAXNITR           ! 13 maximum nitrification rate
  real :: NON               ! 14 maximum NO from nitrification (day-1)
  real :: N2ON              ! 15 maximum N2O from nitrification (day-1)
  real :: KN                ! 16 Michaelis-Menten coefficient [gN/m2]. Use this value if soil represents top 100 cm 
  real :: KDOC              ! 17 Michaelis-Menten coefficient [gC/m2]. Use this value if soil represents top 100 cm 
  real :: DOCMAX            ! 18 DOCMAX
  real :: DNITR2N2O         ! 19 Fraction of denitrification lost as N2O. Range of possible values: 0.002 - 0.047 (Xu-Ri and Prentice, 2008)
  real :: kbeer             ! 20 canopy light extinction coefficient
  real :: kresp             ! 21 additional respiration coefficient, unknown origin ("proposed by Annet Wolf 2008")
  real :: F_NRETAIN         ! 29 fraction of N retained at leaf abscission 
  real :: fpc_tree_max      ! 30 maximum fractional plant coverage of trees

  ! PFT-DEPENDENT PARAMETERS
  real, dimension(npft) :: cton_soil           ! C:N ratio of soil organic matter [1]
  integer, dimension(npft) :: lu_category      ! land use category associated with PFT
  logical, dimension(npft,nlu) :: islu         ! islu(ipft,ilu) is true if ipft belongs to ilu
  logical, dimension(npft) :: grass            ! boolean for growth form 'grass'
  logical, dimension(npft) :: tree             ! boolean for growth form 'tree'
  real, dimension(npft) :: lai0                ! sapling (or grass on initialisation) LAI (=lai_sapl)
  type(orgpool), dimension(npft) :: lm_sapl    ! sapling leaf mass
  type(orgpool), dimension(npft) :: sm_sapl    ! sapling leaf mass
  type(orgpool), dimension(npft) :: hm_sapl    ! sapling leaf mass
  type(orgpool), dimension(npft) :: rm_sapl    ! sapling root mass
  real, dimension(npft) :: long_leaf           ! leaf longeivity [days]
  real, dimension(npft) :: long_sapw           ! sapwood longeivity [days]
  real, dimension(npft) :: long_root           ! root longeivity [days]
  real, dimension(npft) :: sla                 ! specific leaf area [m2/gC]
  real, dimension(npft) :: cton_pro            ! C:N ratio of new production [1]
  integer, dimension(npft) :: pftcode          ! code for identification of PFT (1=C3 grass)
  real, dimension(npft) :: rresp_maint         ! maintenance respiration coefficient
  real, dimension(npft) :: rresp_growth        ! growth respiration coefficient
  logical, dimension(npft) :: nfixer           ! whether plant is capable of symbiotically fixing N
  real :: reinickerp                           ! Reinicker-p for geometry
  real, dimension(npft) :: wooddens            ! wood density (gC * m-3)
  real, dimension(npft) :: allom1, allom2      ! allometry parameters
  real :: allom3                               ! allometry parameter
  real :: crownarea_max                        ! maximum crown area
  real, dimension(npft) :: latosa              ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
  real, dimension(npft) :: woodtosapw          ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
  real, dimension(npft) :: lmtorm0             ! leaf to root ratio under non-water stressed conditions


!----------------------------------------------------------------
! Parameter values from Makela et al., 2008, Table 1
!----------------------------------------------------------------
#if _gpp_makela
  ! Amount of foliage capturing 50% of maximum C gain [g/m2] (given as 2500 kg-leaf-DW/ha)
  real, parameter :: K_f = 500.0

  ! non-photosynthtic or structural N [gN/gC], assuming 0.5 gC/gDW (given as 0.009 kgN/kgDW)
  real, parameter :: n_struct = 0.018

  ! reference photosynthetic N [gN/gC], assuming 0.5 gC/gDW (given as 0.002 kgN/kgDW)
  real, parameter :: n_ref = 0.004

  ! N-saturated specific rate of photosynthesis [gC/g-foliage-C/yr] (given as 8.0 kgC/kgDW/yr)
  ! xxx try: This is per year. Adjust it to get daily value: 16.0/365.0
  real, parameter :: sigma_fm0 = 0.04383562

#endif

#if _allocation_makela
   ! 'Steady-state' pipe length coefficient [m/(gC/gN)], given as 2800 m/(kg-N/kg-DW), assuming 0.5 gC/g-DW
   real, parameter :: c_H = 5600.0

   ! Sapwood weight per unit foliage and pipe length [m-1]
   real, parameter :: alpha_w = 0.8

#endif

#if _nuptake_makela
  ! Amount of roots capturing 50% of available N [g/m2] (given as 2000 kg-leaf-DW/ha)
  real, parameter :: K_r = 400.0
#endif


end module params_modl
