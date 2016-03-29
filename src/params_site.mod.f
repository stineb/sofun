# 1 "params_site.mod.f90"
module md_params_site
!////////////////////////////////////////////////////////////////
!  Module contains site-specific parameters read by getpar_site
! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
! contact: b.stocker@imperial.ac.uk
!----------------------------------------------------------------
  use md_params_siml, only: sitename
  use md_params_core, only: maxgrid

  implicit none

  real    :: lon_site, lat_site, elv_site
  logical :: lTeBS, lGrC3, lGrC4
! logical :: lcultfile
  integer, dimension(maxgrid) :: soilcode_field

contains

  subroutine getpar_site()
!////////////////////////////////////////////////////////////////
!  SR reads site-specific parameters from <sitename>.parameter
!----------------------------------------------------------------
    use md_sofunutils, only: getparreal, getparlogical, getparint

! local variables
  
    write(0,*) 'reading site parameter file ', 'site_paramfils/'//trim(sitename)//'.parameter ...'

    lon_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'longitude' )
    lat_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'latitude' )
    elv_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'altitude' )

! order is important and must be equal to order in site parameter file
! as well as in getpar_modl()
    lTeBS = getparlogical( 'site_paramfils/'//trim(sitename)//'.parameter', 'lTeBS' )
    lGrC3 = getparlogical( 'site_paramfils/'//trim(sitename)//'.parameter', 'lGrC3' )
    lGrC4 = getparlogical( 'site_paramfils/'//trim(sitename)//'.parameter', 'lGrC4' )

! dimension length (maxgrid) is 1 for site-scale simulations
    soilcode_field(1) = getparint( 'site_paramfils/'//trim(sitename)//'.parameter', 'soilcode' )

! ! For C4 grasses, check if cultivar file is available for this site
! if (lGrC4) then
!   inquire( file='site_paramfils/site_cultivars/site_cultivars_'//trim(sitename)//'.txt', exist=lcultfile )
! end if

  end subroutine getpar_site

end module md_params_site
