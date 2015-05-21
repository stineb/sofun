module _vars_core_two
  !////////////////////////////////////////////////////////////////
  !  Module contains all pool variables. These are handled as state
  !  variables of the model and contain a dummy spatial dimension
  !  (to be extended). Daily output for pool variables are done by
  !  copying value of pool to a daily variable. Always use 'p' as
  !  a prefix of variable name.
  !----------------------------------------------------------------
  use _params_core
  use _classdefs

  implicit none

  type(orgpool), dimension(npft,maxgrid)  :: pleaf           ! leaf biomass [gC/ind.] (=lm_ind)
  type(orgpool), dimension(npft,maxgrid)  :: proot           ! root biomass [gC/ind.] (=rm_ind)
  type(orgpool), dimension(npft,maxgrid)  :: psapw           ! sapwood biomass [gC/ind.] (=sm_ind)
  type(orgpool), dimension(npft,maxgrid)  :: pwood           ! heartwood (non-living) biomass [gC/ind.] (=hm_ind)
  type(orgpool), dimension(npft,maxgrid)  :: plabl           ! labile pool, temporary storage of N and C [gC/ind.] (=bm_inc but contains also N) 

  type(carbon),  dimension(npft,maxgrid)  :: pexud           ! exudates pool (very short turnover) [gC/m2]

  type(orgpool), dimension(npft,maxgrid)  :: plitt_af        ! above-ground litter, fast turnover [gC/m2]
  type(orgpool), dimension(npft,maxgrid)  :: plitt_as        ! above-ground litter, slow turnover [gC/m2]
  type(orgpool), dimension(npft,maxgrid)  :: plitt_bg        ! below-ground litter [gC/m2]

  type(orgpool), dimension(nlu,maxgrid)   :: psoil_sl        ! soil organic matter, fast turnover [gC/m2]
  type(orgpool), dimension(nlu,maxgrid)   :: psoil_fs        ! soil organic matter, slow turnover [gC/m2]

  type(nitrogen), dimension(nlu,maxgrid)  :: pninorg         ! total inorganic N pool (sum of NO3 and NH4) [gC/m2]

  real, dimension(nlu,maxgrid) :: dwtot                      ! soil water mass [mm = kg/m2]

  real, dimension(npft,maxgrid) :: r_cton_leaf, r_ntoc_leaf  ! leaf C:N [gC/gN] and N:C ratio [gN/gC]
  real, dimension(npft,maxgrid) :: lma, sla                  ! leaf mass per area [g-dry mass/m2], specific leaf area [m2/g-dry mass]

  real, dimension(npft,maxgrid) :: lai_ind
  real, dimension(npft,maxgrid) :: fpc_ind


end module _vars_core_two

