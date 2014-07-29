       module forcing
       !////////////////////////////////////////////////////////////////
       !  Module contains forcing variables (climate, co2, ...), and
       !  subroutines used to read forcing input files
       !----------------------------------------------------------------
         implicit none
         !real, allocatable :: in_co2(:)
         real :: pco2, ndep
!         real, dimension(:), allocatable :: dtemp, dprec
         real, dimension(365) :: dtemp, dprec

         contains

         subroutine getco2(realyear)
         !////////////////////////////////////////////////////////////////
         !  SR reads this year's atmospheric CO2 from input
         !----------------------------------------------------------------
           implicit none
           integer, intent(in) :: realyear

           pco2 = 399. 

         end subroutine getco2

         subroutine getclimate(realyear)
         !////////////////////////////////////////////////////////////////
         !  SR reads this year's daily temperature and precipitation.
         !----------------------------------------------------------------
           use params_precompiling
           implicit none
           integer, intent(in) :: realyear
           
!           allocate(dtemp(ndayyear))
!           allocate(dprec(ndayyear))

           dtemp(:) = 12.
           dprec(:) = 2.7  ! ~1000 mm / 365 d

         end subroutine getclimate

         subroutine getndep(realyear)
         !////////////////////////////////////////////////////////////////
         !  SR reads this year's annual ndeposition and distributes it
         !  over days according to daily precipitation.
         !----------------------------------------------------------------
           implicit none
           integer, intent(in) :: realyear
           
           ndep    = 1.

         end subroutine getndep

       end module forcing
