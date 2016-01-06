module _params_site
  !////////////////////////////////////////////////////////////////
  !  Module contains site-specific parameters read by getpar_site
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _sofunutils, only: getparreal, getparlogical
  use _params_siml, only: sitename

  implicit none

  real    :: lon_site, lat_site, elv_site
  logical :: lTeBS, lGrC3

contains

  subroutine getpar_site()
    !////////////////////////////////////////////////////////////////
    !  SR reads site-specific parameters from <sitename>.parameter
    !----------------------------------------------------------------
    write(0,*) 'reading site parameter file ', sitename//".parameter ..."

    lon_site = getparreal( sitename//".parameter", 'longitude' )
    lat_site = getparreal( sitename//".parameter", 'latitude' )
    elv_site = getparreal( sitename//".parameter", 'altitude' )

    ! order is important and must be equal to order in site parameter file
    ! as well as in getpar_modl()
    lTeBS = getparlogical( sitename//".parameter", 'lTeBS' )
    lGrC3 = getparlogical( sitename//".parameter", 'lGrC3' )

  end subroutine getpar_site

end module _params_site
