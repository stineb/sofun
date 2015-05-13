module _statevars
  !////////////////////////////////////////////////////////////////
  ! Module contains required state variables apart from those 
  ! contained in '_pools' and '_fluxes'. 
  !----------------------------------------------------------------
  use _params_core, only: nlu, maxgrid

  implicit none

  real, dimension(nlu,maxgrid) :: dtemp_soil          ! soil temperature [deg C] xxx really maxgrid necessary? xxx
  real, dimension(nlu) :: aet_d                       ! daily AET (mm)
  real, dimension(nlu) :: ncpa_d                      ! daily normalised Cramer-Prentice-Alpha (AET/EET) {0,1}
  real, dimension(nlu) :: dro                         ! daily runoff [mm]

  real, dimension(nlu) :: ddoc                         ! surrogate for dissolved organic carbon


end module _statevars

