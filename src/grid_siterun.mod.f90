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
  use md_params_domain, only: type_params_domain

  implicit none

  private
  public gridtype, getgrid, domaininfo_type

  type domaininfo_type
    integer :: nlon
    integer :: nlat
    integer :: maxgrid
    real :: lon_start
    real :: lat_start
    real :: dlon
    real :: dlat
    integer, dimension(:,:), allocatable :: gridarray
  end type domaininfo_type


  type gridtype
    real :: lon
    real :: lat
    real :: elv
    integer :: soilcode
  end type gridtype

contains

  function get_domaininfo( params_domain ) result( domaininfo )
    !////////////////////////////////////////////////////////////////
    ! Gets domain information needed to allocate size of arrays
    !----------------------------------------------------------------
    use md_params_domain, only: type_params_domain

    ! arguments
    type( type_params_domain ) :: params_domain

    ! function return variable
    type( domaininfo_type ) :: domaininfo

    domaininfo%nlon = 1
    domaininfo%nlat = 1

    domaininfo%lon_start = params_domain%lon_site
    domaininfo%lat_start = params_domain%lat_site

    domaininfo%maxgrid = 1

  end function get_domaininfo


  function getgrid( domaininfo, params_domain ) result( out_grid )
    !////////////////////////////////////////////////////////////////
    ! Defines grid variables
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal
    use md_params_domain, only: type_params_domain

    ! arguments
    type( domaininfo_type ) :: domaininfo
    type( type_params_domain ) :: params_domain

    ! function return variable
    type( gridtype ), allocatable, dimension(:) :: out_grid

    ! local variables
    integer :: jpngr, ilon, ilat
    integer :: maxgrid_test

    print*,'reading grid ...'

    ! xxx try. otherwise loop over sites and allocate values for each 
    ! site into vectors containing all sites
    out_grid(:)%lon = domaininfo%lon_start
    out_grid(:)%lat = domaininfo%lat_start
    out_grid(:)%elv = params_domain%elv_site

    print*,'...done'

  end function getgrid

end module md_grid

