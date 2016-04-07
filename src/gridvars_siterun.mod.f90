module md_grid
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
  use md_params_core

  implicit none

  private
  public gridtype, getgrid

  type gridtype
    real :: lon
    real :: lat
    real :: elv
  end type gridtype

contains

  function getgrid( sitename ) result( out_grid )
    !////////////////////////////////////////////////////////////////
    ! Defines grid variables
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    character(len=*) :: sitename

    type( gridtype ), dimension(maxgrid) :: out_grid

    ! xxx try. otherwise loop over sites and allocate values for each 
    ! site into vectors containing all sites
    out_grid(:)%lon = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'longitude' )
    out_grid(:)%lat = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'latitude' )
    out_grid(:)%elv = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'altitude' )

  end function getgrid

end module md_grid

