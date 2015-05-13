module input_site
  !////////////////////////////////////////////////////////////////
  ! This module contains subroutines and variables needed when pre-
  ! scribing certain variables instead of simulating them online.
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  use _params_core, only: npft, maxgrid

  implicit none

#ifdef _prescr_gpp_site
  !real, dimension(nyeartrend,ndayyear,npft) :: dgpp_data
  real, dimension(:,:,:), allocatable :: dgpp_data
#endif

#ifdef _fix_veg
  real, dimension(npft,maxgrid) :: fpc_grid_data
#endif

contains

  subroutine getinput_site
    !////////////////////////////////////////////////////////////////
    !  Subroutine reads to prescribe variables that are otherwise 
    !  (per default) calculated online
    ! This is not appropriate here as its only for one year
    !----------------------------------------------------------------
    use _params_siml, only: nyeartrend, sitename
    use _forcing_siterun, only: read1year_daily  ! xxx put all reading functions/subroutines into a single module
    use _sofunutils, only: getreal  ! xxx put all reading functions/subroutines into a single module
    use _params_core, only: ndayyear

    implicit none

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

#ifdef _fix_veg
    ! Get prescribed PFT selection 
    ! xxx try: make npft a definable variable depending on this selection 
    ! Get prescribed fractional plant cover (FPC) for each PFT
    call getreal(sitename//".parameter",'in_fpc_grid_1',fpc_grid_data(1,1))
    call getreal(sitename//".parameter",'in_fpc_grid_2',fpc_grid_data(2,1))
    call getreal(sitename//".parameter",'in_fpc_grid_3',fpc_grid_data(3,1))
    call getreal(sitename//".parameter",'in_fpc_grid_4',fpc_grid_data(4,1))
    call getreal(sitename//".parameter",'in_fpc_grid_5',fpc_grid_data(5,1))
    call getreal(sitename//".parameter",'in_fpc_grid_6',fpc_grid_data(6,1))
    call getreal(sitename//".parameter",'in_fpc_grid_7',fpc_grid_data(7,1))
    call getreal(sitename//".parameter",'in_fpc_grid_8',fpc_grid_data(8,1))
    call getreal(sitename//".parameter",'in_fpc_grid_9',fpc_grid_data(9,1))
#endif

    return

  end subroutine getinput_site

end module input_site
