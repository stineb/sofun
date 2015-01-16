module input_site
  !////////////////////////////////////////////////////////////////
  ! This module contains subroutines and variables needed when pre-
  ! scribing certain variables instead of simulating them online.
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  use params_core
  use params_siml, only: nyeartrend, sitename
  use forcing_siterun, only: read1year_daily

  implicit none

#if _prescr_gpp_site
  !real, dimension(nyeartrend,ndayyear,npft) :: dgpp_data
  real, dimension(:,:,:), allocatable :: dgpp_data
#endif

contains

  subroutine getinput_site
    !////////////////////////////////////////////////////////////////
    !  Subroutine reads to prescribe variables that are otherwise 
    !  (per default) calculated online
    ! This is not appropriate here as its only for one year
    !----------------------------------------------------------------

    ! local variables
    integer :: day, mo, dm, pft, yr

#if _prescr_gpp_site
    ! PRESCRIBED DAILY GPP FOR ONE YEAR
    write(0,*) 'prescribe daily GPP ...'
    allocate(dgpp_data(nyeartrend,ndayyear,npft))
    dgpp_data(:,:,:) = 0.0
    do yr=1,nyeartrend
      dgpp_data(yr,:,1) = read1year_daily(sitename//'_daily_gpp_med_STANDARD.txt')
    end do
    write(0,*) '... done'

    !allocate(dgpp_data(nyeartrend,ndayyear,npft))
    !dgpp_data(:,:,:) = 0.0
    !do yr=1,nyeartrend
    !  day=0
    !  do mo=1,nmonth
    !    do dm=1,ndaymonth(mo)
    !      day=day+1
    !      do pft=1,npft
    !        dgpp_data(yr,day,pft) = getvalreal_STANDARD( &
    !          sitename//'_daily_gpp_med_STANDARD.txt', mo=mo, dm=dm &
    !          )
    !      end do
    !    end do
    !  enddo
    !enddo
    !write(0,*) '... done'
#endif

    return

  end subroutine getinput_site

end module input_site
