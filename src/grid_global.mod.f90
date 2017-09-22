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
  public gridtype, getgrid, domaininfo_type, get_domaininfo

  type domaininfo_type
    integer :: nlon
    integer :: nlat
    integer :: maxgrid
    real :: lon_start
    real :: lat_start
    real :: dlon
    real :: dlat
    integer, dimension(:,:), allocatable :: gridarray
    real, dimension(:), allocatable :: lon
    real, dimension(:), allocatable :: lat
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
    type( type_params_domain ), intent(in) :: params_domain

    ! function return variable
    type( domaininfo_type ) :: domaininfo

    ! local variables
    integer, dimension(:,:), allocatable :: myarr
    integer :: ilon, ilat, jpngr

    if ( trim(params_domain%filnam_landmask)=="landmaskfile_global_halfdeg.nc" ) then

      print*,'this is a halfdeg simulation'

      allocate( domaininfo%gridarray(4,4) )
      allocate( domaininfo%lon(4) )
      allocate( domaininfo%lat(4) )

      domaininfo%nlon = 4
      domaininfo%nlat = 4

      domaininfo%lon_start = 0.25
      domaininfo%lat_start = 30.25

      domaininfo%dlat = 0.5
      domaininfo%dlon = 0.5

      ! data domaininfo%gridarray(1,:)/0,0,0,0/ 
      ! data domaininfo%gridarray(2,:)/0,0,0,0/ 
      ! data domaininfo%gridarray(3,:)/1,1,1,1/
      ! data domaininfo%gridarray(4,:)/1,1,1,1/
      do ilon=1,domaininfo%nlon
        do ilat=1,domaininfo%nlat
          if (ilat<3) then
            domaininfo%gridarray(ilon,ilat) = 0
          else
            domaininfo%gridarray(ilon,ilat) = 1
          end if
        end do
      end do

      do ilon=1,domaininfo%nlon
        domaininfo%lon(ilon) = domaininfo%lon_start + (ilon-1) * domaininfo%dlon
      end do

      do ilat=1,domaininfo%nlat
        domaininfo%lat(ilat) = domaininfo%lat_start + (ilat-1) * domaininfo%dlat
      end do

    else if ( trim(params_domain%filnam_landmask)=="landmaskfile_global_1x1deg.nc" ) then

      print*,'this is a 1x1deg simulation'

      allocate( domaininfo%gridarray(2,2) )
      allocate( domaininfo%lon(2) )
      allocate( domaininfo%lat(2) )

      domaininfo%nlon = 2
      domaininfo%nlat = 2

      domaininfo%lon_start = 0.5
      domaininfo%lat_start = 30.5

      domaininfo%dlat = 1.0
      domaininfo%dlon = 1.0

      ! data domaininfo%gridarray(1,:)/0,0/ 
      ! data domaininfo%gridarray(2,:)/1,1/
      do ilon=1,domaininfo%nlon
        do ilat=1,domaininfo%nlat
          if (ilat>1) then
            domaininfo%gridarray(ilon,ilat) = 0
          else
            domaininfo%gridarray(ilon,ilat) = 1
          end if
        end do
      end do
      
      do ilon=1,domaininfo%nlon
        domaininfo%lon(ilon) = domaininfo%lon_start + (ilon-1) * domaininfo%dlon
      end do

      do ilat=1,domaininfo%nlat
        domaininfo%lat(ilat) = domaininfo%lat_start + (ilat-1) * domaininfo%dlat
      end do

    else

      print*,'Error: landmask file name unknown'
      stop

    end if

    ! Get total number of land gridcells (value 1 in landmaskfile)
    jpngr = 0
    do ilon=1,domaininfo%nlon
      do ilat=1,domaininfo%nlat
        if (domaininfo%gridarray(ilon,ilat)==1) then
          jpngr = jpngr + 1
        end if
      end do
    end do
    domaininfo%maxgrid = jpngr

  end function get_domaininfo


  function getgrid( domaininfo ) result( out_grid )
    !////////////////////////////////////////////////////////////////
    ! Defines grid variables
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal
    use md_params_domain, only: type_params_domain

    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo

    ! function return variable
    type( gridtype ), allocatable, dimension(:) :: out_grid

    ! local variables
    integer :: jpngr, ilon, ilat
    integer :: maxgrid_test

    allocate( out_grid(domaininfo%maxgrid ) )

    ! xxx test
    jpngr = 0
    do ilon=1,domaininfo%nlon
      do ilat=1,domaininfo%nlat
        if (domaininfo%gridarray(ilon,ilat)==1) then
          jpngr=jpngr+1
          out_grid(jpngr)%lon = domaininfo%lon(ilon)
          out_grid(jpngr)%lat = domaininfo%lat(ilat)
        end if
      end do
    end do

    out_grid(:)%elv = 100.0
    out_grid(:)%soilcode = 1

  end function getgrid

end module md_grid

