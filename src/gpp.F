subroutine gpp( jpngr, co2, dgpp_data )
  !////////////////////////////////////////////////////////////////
  !  SR reads GPP from input file (for the time being)
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  use fluxes
  use params_siml
  use params_modl
  use pools, only: pleaf
  use treegeometry, only: nind
  implicit none

  ! optional set of arguments if 'dgpp' is to be predicted
  integer, intent(in), optional :: jpngr                  ! gridcell number
  real, intent(in), optional    :: co2                    ! atmospheric CO2 concentration [ppm]

  ! optional set of arguments if 'dgpp' is to be prescribed
  real, dimension(npft), intent(in), optional :: dgpp_data

  ! local variables 
  integer :: pft                     ! pft counter

#if _gpp_makela
  real :: n_photo
  real :: sigma_fm 
#endif


  if (present(dgpp_data)) then
    !----------------------------------------------------------------
    ! PRESCRIBED GPP
    !----------------------------------------------------------------
    dgpp(:) = dgpp_data(:)

  else
    !----------------------------------------------------------------
    ! PREDICTED GPP FROM CO2, FAPAR
    !----------------------------------------------------------------
   ! necessary because of keyword arguments used in function call
!       interface
!         function getvalreal( filename, realyear, day, dm, mo )
!           character(len=*), intent(in) :: filename
!           integer, intent(in) :: realyear
!           integer, intent(in), optional :: day ! day in year (1:365)
!           integer, intent(in), optional :: dm  ! day in month (1:31)
!           integer, intent(in), optional :: mo  ! month in year (1:12)
!         end function getvalreal
!       end interface
#if _gpp_makela
    !----------------------------------------------------------------
    ! parameter values for Pinus Sylvestris (their table 1)
    !----------------------------------------------------------------
    ! Actual photosynthetic N concentration in foliage, reduced by non-photosynthtic/structural N
    do pft=1,npft
      n_photo = max( pleaf(pft,jpngr)%n%n14 * nind(pft,jpngr) - n_struct, 0.0 )

      ! light-saturated specific rate of photosynthesis
      sigma_fm = ( sigma_fm0 * n_photo ) / ( n_photo + n_ref )

      ! Eq.6 in Makela et al., 2008
      dgpp(pft) = ( sigma_fm * pleaf(pft,jpngr)%c%c12 * nind(pft,jpngr) * K_f ) &
            / ( pleaf(pft,jpngr)%c%c12 * nind(pft,jpngr) + K_f )

    enddo
#else

  stop 'GPP: no GPP parametrisation selected. no default defined. abort.'

#endif

  end if

end subroutine gpp

