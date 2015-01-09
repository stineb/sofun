      module input_site
        !////////////////////////////////////////////////////////////////
        ! This module contains subroutines and variables needed when pre-
        ! scribing certain variables instead of simulating them online.
        !----------------------------------------------------------------
#include "sofun_module_control.inc"
        use params_core
        use params_siml, only: nyeartrend, sitename
        use forcing

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
        !----------------------------------------------------------------

        ! local variables
        integer :: day, mo, dm, pft, yr

#if _prescr_gpp_site
        ! Prescribe daily GPP
        write(0,*) 'prescribe daily GPP ...'
        allocate(dgpp_data(nyeartrend,ndayyear,npft))
        dgpp_data(:,:,:) = 0.
        !print*,'length(dgpp_data)[1]', size(dgpp_data(:,1,1))
        !print*,'length(dgpp_data)[2]', size(dgpp_data(1,:,1))
        !print*,'length(dgpp_data)[3]', size(dgpp_data(1,1,:))

        do yr=1,nyeartrend
          day=0
          do mo=1,nmonth
            do dm=1,ndaymonth(mo)
              day=day+1
              do pft=1,npft
                  dgpp_data(yr,day,pft) = getvalreal_STANDARD( &
                    sitename//'_daily_gpp_med_STANDARD.txt', mo=mo, dm=dm &
                    )
              end do
            end do
          enddo
        enddo
        write(0,*) '... done'
#endif

      end subroutine getinput_site

      end module input_site
