subroutine getout_daily( jpngr, moy, doy )
  !////////////////////////////////////////////////////////////////
  ! SR called daily to sum up daily output variables.
  ! Note that output variables are collected only for those variables
  ! that are global anyway (e.g., outdcex). Others are not made 
  ! global just for this, but are collected inside the subroutine 
  ! where they are defined.
  !----------------------------------------------------------------
  use _params_core, only: ndayyear
  use _params_modl, only: npft, lu_category
  use _vars_core
  use _outvars
  use _pools
  use _statevars
  use _waterbal, only: getout_daily_waterbal
  use _nuptake, only: getout_daily_nuptake
  use _littersom, only: getout_daily_littersom

  ! ARGUMENTS
  integer, intent(in) :: jpngr
  integer, intent(in) :: moy
  integer, intent(in) :: doy

  ! LOCAL VARIABLES
  integer :: pft

  !print*,'netmin_litt ',netmin_litt(1)
  !print*,'netmin_soil ',netmin_soil(1)

  !----------------------------------------------------------------
  ! DAILY
  ! Collect daily output variables
  ! so far not implemented for isotopes
  !----------------------------------------------------------------
  outdgpp(:,doy,jpngr)   = dgpp(:)
  outdnpp(:,doy,jpngr)   = dnpp(:)%c12
  outdnup(:,doy,jpngr)   = dnup(:)%n14
  outdCleaf(:,doy,jpngr) = pleaf(:,jpngr)%c%c12

  outdninorg(:,doy,jpngr)   = pninorg(:,jpngr)%n14
  outdtemp_soil(:,doy,jpngr)= dtemp_soil(:,jpngr)

  !----------------------------------------------------------------
  ! ANNUAL SUM OVER DAILY VALUES
  ! Collect annual output variables as sum of daily values
  !----------------------------------------------------------------
  outagpp(:,jpngr)    = outagpp(:,jpngr) + dgpp(:)
  outanpp(:,jpngr)    = outanpp(:,jpngr) + dnpp(:)%c12
  outanup(:,jpngr)    = outanup(:,jpngr) + dnup(:)%n14
  outaninorg(:,jpngr) = outaNinorg(:,jpngr) + pninorg(:,jpngr)%n14 / ndayyear

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC VARIABLES
  !----------------------------------------------------------------
  call getout_daily_waterbal( jpngr, moy, doy )
  call getout_daily_nuptake( jpngr, moy, doy )
  call getout_daily_littersom( jpngr, moy, doy )

end subroutine getout_daily


subroutine getout_monthly( jpngr, moy )
  !////////////////////////////////////////////////////////////////
  !  SR called once a year to gather monthly output variables.
  !----------------------------------------------------------------
  use _waterbal, only: getout_monthly_waterbal
  implicit none

  ! ARGUMENTS
  integer, intent(in) :: jpngr
  integer, intent(in) :: moy

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC VARIABLES
  !----------------------------------------------------------------
  call getout_monthly_waterbal( jpngr, moy )

  return

end subroutine getout_monthly


subroutine getout_annual( jpngr )
  !////////////////////////////////////////////////////////////////
  !  SR called once a year to gather annual output variables.
  !----------------------------------------------------------------
  use _params_modl, only: npft, lu_category
  use _pools
  use _outvars
  use _littersom, only: getout_annual_littersom

  implicit none

  ! ARGUMENTS
  integer, intent(in) :: jpngr

  outaCveg(:,jpngr)  = pleaf(:,jpngr)%c%c12 + proot(:,jpngr)%c%c12 !+ psapw(:,jpngr)%c%c12 + pwood(:,jpngr)%c%c12

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC VARIABLES
  !----------------------------------------------------------------
  call getout_annual_littersom( jpngr )


  return

end subroutine getout_annual

