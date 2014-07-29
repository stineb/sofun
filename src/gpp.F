      subroutine gpp( co2, realyear, mo, dm )
      !////////////////////////////////////////////////////////////////
      !  SR reads GPP from input file (for the time being)
      !----------------------------------------------------------------
      use forcing
      use dailyvars
      use params_siml

      implicit none

      ! arguments
      real, intent(in) :: co2            ! atmospheric CO2 concentration [ppm]
      integer, intent(in) :: realyear     ! year AD
      integer, intent(in) :: mo           ! month of year (1:12)
      integer, intent(in) :: dm           ! day of month (1:31)

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

      do pft=1,npft
    
          ! read daily GPP value from input file
          dgpp(pft) = carbon( getvalreal_STANDARD( sitename//'_daily_gpp_STANDARD.txt', realyear, mo=mo, dm=dm ) )
          
          !print*,'realyear, mo, dm, dgpp ',realyear,mo,dm,dgpp

      end do  

      end subroutine gpp