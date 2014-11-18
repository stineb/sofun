      subroutine gpp( jpngr, co2, year, realyear, mo, dm, day )
      !////////////////////////////////////////////////////////////////
      !  SR reads GPP from input file (for the time being)
      !----------------------------------------------------------------
#include "sofun_module_control.inc"
      use forcing
      use statevars
      use params_siml
      use pools, only: pleaf
#if _prescr_gpp_site
      use input_site, only: dgpp_data
#endif
      implicit none

      ! arguments
      integer, intent(in) :: jpngr                  ! gridcell number
      real,    intent(in) :: co2                    ! atmospheric CO2 concentration [ppm]
      integer, intent(in), optional :: year         ! simulation year
      integer, intent(in), optional :: realyear     ! year AD
      integer, intent(in), optional :: mo           ! month of year (1:12)
      integer, intent(in), optional :: dm           ! day of month (1:31)
      integer, intent(in), optional :: day          ! day of year (1:365)

      ! local variables 
      integer :: pft                     ! pft counter

#if _gpp_makela
      real :: n_photo
      real :: sigma_fm 

      !----------------------------------------------------------------
      ! Parameter values from Makela et al., 2008, Table 1
      !----------------------------------------------------------------
      ! Amount of foliage capturing 50% of maximum C gain [g/m2] (given as 2500 kg-leaf-DW/ha)
      real, parameter :: K_f = 500.0

      ! non-photosynthtic or structural N [gN/gC], assuming 0.5 gC/gDW (given as 0.009 kgN/kgDW)
      real, parameter :: n_struct = 0.018

      ! reference photosynthetic N [gN/gC], assuming 0.5 gC/gDW (given as 0.002 kgN/kgDW)
      real, parameter :: n_ref = 0.004

      ! N-saturated specific rate of photosynthesis [gC/g-foliage-C/yr] (given as 8.0 kgC/kgDW/yr)
      ! xxx try: This is per year. Adjust it to get daily value: 16.0/365.0
      real, parameter :: sigma_fm0 = 0.04383562

#endif
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
      ! Prescribe GPP for each day for site-scale simulation as stored in array 'gpp_data'
      do pft=1,npft
        dgpp(pft) = dgpp_data( transyear, day, pft )
      enddo

#elif _gpp_makela
      !----------------------------------------------------------------
      ! Calculate GPP online
      ! parameter values for Pinus Sylvestris (their table 1)
      !----------------------------------------------------------------
      ! Actual photosynthetic N concentration in foliage, reduced by non-photosynthtic/structural N
      n_photo = max( pleaf(pft,jpngr)%n%n14 * nind(pft,jpngr) - n_struct, 0.0 )

      ! light-saturated specific rate of photosynthesis
      sigma_fm = ( sigma_fm0 * n_photo ) / ( n_photo + n_ref )

      ! Eq.6 in Makela et al., 2008
      dgpp(pft) = ( sigma_fm * pleaf(pft,jpngr)%c%c12 * nind(pft,jpngr) * K_f ) / ( pleaf(pft,jpngr)%c%c12 * nind(pft,jpngr) + K_f )

#else
      ! xxx try: this is just a placeholder for an explicit GPP parameterisation
      do pft=1,npft
          ! read daily GPP value from input file
          dgpp(pft) = getvalreal_STANDARD( sitename//'_daily_gpp_med_STANDARD.txt', mo=mo, dm=dm )
      end do  

#endif

      end subroutine gpp

