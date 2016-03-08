module _gridvars
  !////////////////////////////////////////////////////////////////
  ! Module contains global variables defining the model grid.
  ! A module 'gridvars_*' must contain the following subroutines:
  ! - getgrid
  !
  ! ... and define the following variables that are global within
  ! 'sofun' (but passed on to 'biosphere' as arguments).
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core

  implicit none

  real, dimension(maxgrid)  :: lon        ! longitude vector/field (degrees E)              
  real, dimension(maxgrid)  :: lat        ! latitude vector/field (degrees N)             
  real, dimension(maxgrid)  :: elv        ! elevation vector/field (m above sea level)                  

contains

  subroutine getgrid
    !////////////////////////////////////////////////////////////////
    ! Defines grid variables
    !----------------------------------------------------------------
    use _params_site, only: lon_site, lat_site, elv_site

    implicit none

    ! xxx try. otherwise loop over sites and allocate values for each 
    ! site into vectors containing all sites
    lon(1) = lon_site
    lat(1) = lat_site
    elv(1) = elv_site

    return

  end subroutine getgrid

end module _gridvars

