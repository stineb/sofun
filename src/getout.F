subroutine getout_daily( jpngr, moy, doy )
  !////////////////////////////////////////////////////////////////
  !  SR called daily to sum up daily output variables.
  !----------------------------------------------------------------
  use params_modl, only: npft, lu_category
  use fluxes
  use outvars
  use pools
  use waterbal_stash, only: getout_daily_waterbal
  use nuptake_fun, only: getout_daily_nuptake

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
  outdcex(:,doy,jpngr)   = dcex(:)

  do pft=1,npft
    lu=lu_category(pft)
    outdnetmin(lu,doy,jpngr) = outdnetmin(lu,doy,jpngr) + netmin_litt(pft)
  end do

  outdnetmin(:,doy,jpngr)   = outdnetmin(:,doy,jpngr) + netmin_soil(:)
  outdninorg(:,doy,jpngr)   = nh4(:,jpngr)%n14 + no3(:,jpngr)%n14
  outdnfixfree(:,doy,jpngr) = dnfix_free(:) 

  !----------------------------------------------------------------
  ! ANNUAL SUM OVER DAILY VALUES
  ! Collect annual output variables as sum of daily values
  !----------------------------------------------------------------
  outagpp(:,jpngr) = outagpp(:,jpngr) + dgpp(:)
  outanpp(:,jpngr) = outanpp(:,jpngr) + dnpp(:)%c12
  outanup(:,jpngr) = outanup(:,jpngr) + dnup(:)%n14

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC VARIABLES
  !----------------------------------------------------------------
  call getout_daily_waterbal( jpngr, moy, doy )
  call getout_daily_nuptake( jpngr, moy, doy )

end subroutine getout_daily


subroutine getout_monthly( jpngr, moy )
  !////////////////////////////////////////////////////////////////
  !  SR called once a year to gather monthly output variables.
  !----------------------------------------------------------------
  use waterbal_stash, only: getout_monthly_waterbal
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
  use params_modl, only: npft, lu_category
  use pools
  use outvars

  implicit none

  ! ARGUMENTS
  integer, intent(in) :: jpngr

  outaClitt(:,jpngr) = plitt_af(:,jpngr)%c%c12 + plitt_as(:,jpngr)%c%c12 + plitt_bg(:,jpngr)%c%c12
  outaCsoil(:,jpngr) = psoil_sl(:,jpngr)%c%c12 + psoil_fs(:,jpngr)%c%c12


end subroutine getout_annual

