module params_site
  !////////////////////////////////////////////////////////////////
  !  Module contains site-specific parameters read by getpar_site
  !----------------------------------------------------------------
  use params_siml
  implicit none
  real :: lon_site, lat_site, elv_site

contains

  subroutine getpar_site
    !////////////////////////////////////////////////////////////////
    !  SR reads site-specific parameters from <sitename>.parameter
    !----------------------------------------------------------------
    implicit none

    write(0,*) 'reading site parameter file ', sitename//".parameter ..."
    call getreal(sitename//".parameter",'longitude',lon_site)
    call getreal(sitename//".parameter",'latitude',lat_site)
    call getreal(sitename//".parameter",'altitude',elv_site)

  end subroutine getpar_site

end module params_site
