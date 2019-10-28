module md_grid
  !////////////////////////////////////////////////////////////////
  ! Module for handling variables defining the model grid. For site-
  ! scale simulation this is trivial but still necesssary.
  !
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core
  use md_params_domain, only: type_params_domain

  implicit none

  private
  public gridtype, getgrid, domaininfo_type, get_domaininfo

  ! note that for site scale simulations 'domaininfo_type' contains also 'elv', 'soilcode', and 'whc'
  type domaininfo_type
    integer :: nlon
    integer :: nlat
    real    :: dlon
    real    :: dlat
    integer :: maxgrid
    integer, dimension(:,:), allocatable :: gridarray
    real, dimension(:), allocatable :: lon
    real, dimension(:), allocatable :: lat
    real :: elv
    real :: whc
    integer :: soilcode
    real :: landarea
    character(len=256) :: domain_name  ! This is the site name for site-scale simulations or the character identifyier defining the resolution for global simulations
  end type domaininfo_type

  type gridtype
    integer :: ilon
    integer :: ilat
    real :: lon
    real :: lat
    real :: elv
    real :: landfrac
    real :: area
    logical :: dogridcell
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

    domaininfo%dlon = 0
    domaininfo%dlat = 0

    domaininfo%maxgrid = 1

    allocate( domaininfo%gridarray(1,1) )
    allocate( domaininfo%lon(1) )
    allocate( domaininfo%lat(1) )

    domaininfo%gridarray(1,1) = 0.0

    ! Copy domain parameters
    domaininfo%lon(1)      = params_domain%lon_site
    domaininfo%lat(1)      = params_domain%lat_site
    domaininfo%elv         = params_domain%elv_site
    domaininfo%whc         = params_domain%whc_site
    domaininfo%soilcode    = params_domain%soilcode_site

    domaininfo%domain_name = params_domain%domain_name

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

    allocate( out_grid(domaininfo%maxgrid) )

    out_grid(1)%ilon = 1
    out_grid(1)%ilat = 1

    out_grid(1)%landfrac = 1.0
    out_grid(1)%area     = 1.0

    out_grid(1)%dogridcell = .true.
    
    out_grid(1)%lon = domaininfo%lon(1)
    out_grid(1)%lat = domaininfo%lat(1)
    out_grid(1)%elv = params_domain%elv_site

  end function getgrid

end module md_grid

