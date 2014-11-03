      subroutine gpp( co2, year, realyear, mo, dm, day )
      !////////////////////////////////////////////////////////////////
      !  SR reads GPP from input file (for the time being)
      !----------------------------------------------------------------
#include "sofun_module_control.inc"
      use forcing
      use dailyvars
      use params_siml
#if _prescr_gpp_site
      use input_site, only: dgpp_data
#endif
      implicit none

      ! arguments
      real,    intent(in) :: co2                    ! atmospheric CO2 concentration [ppm]
      integer, intent(in), optional :: year         ! simulation year
      integer, intent(in), optional :: realyear     ! year AD
      integer, intent(in), optional :: mo           ! month of year (1:12)
      integer, intent(in), optional :: dm           ! day of month (1:31)
      integer, intent(in), optional :: day          ! day of year (1:365)

      ! local variables 
      integer :: pft                     ! pft counter

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


#if _prescr_gpp_site
      ! Prescribe GPP for each day for site-scale simulation
      do pft=1,npft
        dgpp(pft) = carbon( dgpp_data( transyear, day, pft) )
      enddo
#else
      ! Calculate GPP online
      ! xxx try: this is just a placeholder for an explicit GPP parameterisation
      do pft=1,npft
          ! read daily GPP value from input file
          dgpp(pft) = carbon( getvalreal_STANDARD( sitename//'_daily_gpp_med_STANDARD.txt', mo=mo, dm=dm ) )
      end do  
#endif

      end subroutine gpp

