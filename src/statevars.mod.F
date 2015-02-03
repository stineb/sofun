module statevars
  !////////////////////////////////////////////////////////////////
  ! Module contains required state variables apart from those 
  ! contained in 'pools' and 'fluxes'. 
  !----------------------------------------------------------------
  use params_core, only: nlu, maxgrid

  implicit none

  real, dimension(nlu,maxgrid) :: dtemp_soil          ! soil temperature [deg C] xxx really maxgrid necessary? xxx
  real, dimension(nlu) :: aet_d                       ! daily AET (mm)
  real, dimension(nlu) :: cpa_d                       ! daily Cramer-Prentice-Alpha (AET/EET)
  real, dimension(nlu) :: dro                         ! daily runoff [mm]

  real, dimension(nlu) :: ddoc                         ! surrogate for dissolved organic carbon


end module statevars

