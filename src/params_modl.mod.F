module _params_modl
  !////////////////////////////////////////////////////////////////
  !  Module contains (constrainable) model parameters.
  !  Model parameters adopted here are from LPX C3 grass PFT
  !  Litter and soil turnover parameters are divided by 365 to 
  !  convert from [1/yr] to [1/d].
  ! Copyright (C) 2015  Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core

  implicit none

  ! NON PFT-DEPENDENT PARAMETERS
  real :: kbeer             ! canopy light extinction coefficient
  real :: F_NRETAIN         ! fraction of N retained at leaf abscission 
  real :: fpc_tree_max      ! maximum fractional plant coverage of trees
  real :: growtheff         ! growth respiration coefficient = yield factor [unitless]

  ! PFT-DEPENDENT PARAMETERS
  real, dimension(npft)        :: k_decay_leaf        ! leaf decay constant [days-1]
  real, dimension(npft)        :: k_decay_sapw        ! sapwood decay constant [days-1]
  real, dimension(npft)        :: k_decay_root        ! root decay constant [days-1]
  real, dimension(npft)        :: cton_soil           ! C:N ratio of soil organic matter [gC/gN]

  integer, dimension(npft)     :: lu_category      ! land use category associated with PFT
  integer, dimension(npft)     :: pftcode          ! code for identification of PFT (1=C3 grass)
  logical, dimension(npft,nlu) :: islu             ! islu(ipft,ilu) is true if ipft belongs to ilu
  logical, dimension(npft)     :: grass            ! boolean for growth form 'grass'
  logical, dimension(npft)     :: tree             ! boolean for growth form 'tree'
  logical, dimension(npft)     :: nfixer              ! whether plant is capable of symbiotically fixing N

  real                         :: reinickerp          ! Reinicker-p for geometry
  real, dimension(npft)        :: wooddens            ! wood density (gC * m-3)
  real, dimension(npft)        :: allom1, allom2      ! allometry parameters
  real                         :: allom3              ! allometry parameter
  real, dimension(npft)        :: crownarea_max       ! maximum crown area
  real, dimension(npft)        :: latosa              ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
  real, dimension(npft)        :: woodtosapw          ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
  real, dimension(npft)        :: lmtorm0             ! leaf to root ratio under non-water stressed conditions

  real, dimension(npft)        :: r_cton_root         ! C:N ratio in roots
  real, dimension(npft)        :: r_ntoc_root         ! N:C ratio in roots (inverse of 'r_cton_root')

end module _params_modl
