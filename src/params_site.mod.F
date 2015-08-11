module _params_site
  !////////////////////////////////////////////////////////////////
  !  Module contains site-specific parameters read by getpar_site
  ! Copyright (C) 2015  Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  implicit none

  real :: lon_site, lat_site, elv_site

contains

  subroutine getpar_site()
    !////////////////////////////////////////////////////////////////
    !  SR reads site-specific parameters from <sitename>.parameter
    !----------------------------------------------------------------
    use _sofunutils, only: getreal
    use _params_siml, only: sitename

    write(0,*) 'reading site parameter file ', sitename//".parameter ..."
    call getreal(sitename//".parameter",'longitude',lon_site)
    call getreal(sitename//".parameter",'latitude',lat_site)
    call getreal(sitename//".parameter",'altitude',elv_site)

  end subroutine getpar_site

end module _params_site
