module _params_site
  !////////////////////////////////////////////////////////////////
  !  Module contains site-specific parameters read by getpar_site
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _sofunutils, only: getparreal
  use _params_siml, only: sitename
  implicit none

  real :: lon_site, lat_site, elv_site

contains

  subroutine getpar_site()
    !////////////////////////////////////////////////////////////////
    !  SR reads site-specific parameters from <sitename>.parameter
    !----------------------------------------------------------------
    write(0,*) 'reading site parameter file ', sitename//".parameter ..."
    lon_site = getparreal( sitename//".parameter", 'longitude' )
    lat_site = getparreal( sitename//".parameter", 'latitude' )
    elv_site = getparreal( sitename//".parameter", 'altitude' )

  end subroutine getpar_site

end module _params_site
