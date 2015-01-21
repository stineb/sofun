module pools
  !////////////////////////////////////////////////////////////////
  !  Module contains all pool variables. These are handled as state
  !  variables of the model and contain a dummy spatial dimension
  !  (to be extended). Daily output for pool variables are done by
  !  copying value of pool to a daily variable. Always use 'p' as
  !  a prefix of variable name.
  !----------------------------------------------------------------
  use params_core
  use classdefs

  implicit none

  type(orgpool), dimension(npft,maxgrid)  :: pleaf    ! leaf biomass [gC/ind.] (=lm_ind)
  type(orgpool), dimension(npft,maxgrid)  :: proot    ! root biomass [gC/ind.] (=rm_ind)
  type(orgpool), dimension(npft,maxgrid)  :: psapw    ! sapwood biomass [gC/ind.] (=sm_ind)
  type(orgpool), dimension(npft,maxgrid)  :: pwood    ! heartwood (non-living) biomass [gC/ind.] (=hm_ind)
  type(orgpool), dimension(npft,maxgrid)  :: plabl    ! labile pool, temporary storage of N and C [gC/ind.] (=bm_inc but contains also N) 

  type(carbon),  dimension(npft,maxgrid)  :: pexud    ! exudates pool (very short turnover) [gC/m2]

  type(orgpool), dimension(npft,maxgrid)  :: plitt_af ! above-ground litter, fast turnover [gC/m2]
  type(orgpool), dimension(npft,maxgrid)  :: plitt_as ! above-ground litter, slow turnover [gC/m2]
  type(orgpool), dimension(npft,maxgrid)  :: plitt_bg ! below-ground litter [gC/m2]

  type(orgpool), dimension(nlu,maxgrid)   :: psoil_sl ! soil organic matter, fast turnover [gC/m2]
  type(orgpool), dimension(nlu,maxgrid)   :: psoil_fs ! soil organic matter, slow turnover [gC/m2]

  type(nitrogen), dimension(nlu,maxgrid)  :: nh4      ! ammonium [gN/m2]
  type(nitrogen), dimension(nlu,maxgrid)  :: no3      ! nitrate [gN/m2]
  type(nitrogen), dimension(nlu,maxgrid)  :: no2      ! NO2 [gN/m2]

  real, dimension(nlu,maxgrid)            :: no_w,  no_d    ! NO in wet and dry microsites (split done in ntransform) [gN/m2]
  real, dimension(nlu,maxgrid)            :: n2o_w, n2o_d   ! N2O in wet and dry microsites (split done in ntransform) [gN/m2]
  real, dimension(nlu,maxgrid)            :: n2_w           ! N2 in wet microsites (split done in ntransform) [gN/m2]

end module pools

