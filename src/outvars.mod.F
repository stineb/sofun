module outvars
!////////////////////////////////////////////////////////////////
!  Module contains all daily variables. These variables have no 
!  spatial dimension. At the end of the daily loop, values may be 
!  copied to spatial arrays.
!----------------------------------------------------------------
  use params_core
  use classdefs

  implicit none

  !----------------------------------------------------------------
  ! Daily output variables
  !----------------------------------------------------------------
  real, dimension(npft,ndayyear,maxgrid) :: outdgpp
  real, dimension(npft,ndayyear,maxgrid) :: outdnpp
  real, dimension(npft,ndayyear,maxgrid) :: outdnup
  real, dimension(npft,ndayyear,maxgrid) :: outdCleaf
  real, dimension(npft,ndayyear,maxgrid) :: outdcex

  real, dimension(nlu,ndayyear,maxgrid)  :: outdninorg
  real, dimension(nlu,ndayyear,maxgrid)  :: outdnetmin
  real, dimension(nlu,ndayyear,maxgrid)  :: outdnfixfree  ! N fixation by free-living organisms [gN/m2/d]
  real, dimension(nlu,ndayyear,maxgrid)  :: outdtemp_soil
  
  !----------------------------------------------------------------
  ! Monthly output variables
  !----------------------------------------------------------------


  !----------------------------------------------------------------
  ! Annual output variables
  !----------------------------------------------------------------
  real, dimension(npft,maxgrid) :: outagpp
  real, dimension(npft,maxgrid) :: outanpp
  real, dimension(npft,maxgrid) :: outanup
  real, dimension(npft,maxgrid) :: outaCveg
  real, dimension(npft,maxgrid) :: outaClitt
  real, dimension(nlu,maxgrid)  :: outaCsoil
  real, dimension(npft,maxgrid) :: outaCalloc
  real, dimension(npft,maxgrid) :: outaNalloc
  real, dimension(nlu,maxgrid)  :: outanreq      ! N required from litter -> soil transfer [gN/m2/d]
  real, dimension(npft,maxgrid) :: outaClit2soil
  real, dimension(npft,maxgrid) :: outaNlit2soil
  real, dimension(npft,maxgrid) :: outaCveg2lit
  real, dimension(npft,maxgrid) :: outaNveg2lit

end module outvars
