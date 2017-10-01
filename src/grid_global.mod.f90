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
  use netcdf

  implicit none

  private
  public gridtype, getgrid, domaininfo_type, get_domaininfo

  type domaininfo_type
    integer :: nlon
    integer :: nlat
    real    :: dlon
    real    :: dlat
    integer :: maxgrid
    integer, dimension(:,:), allocatable :: gridarray
    real, dimension(:), allocatable :: lon
    real, dimension(:), allocatable :: lat
    real :: landarea
  end type domaininfo_type

  type gridtype
    integer :: ilon
    integer :: ilat
    real :: lon
    real :: lat
    real :: elv
    real :: landfrac
    real :: area
    integer :: soilcode
    logical :: dogridcell
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

    ! This will be the netCDF ID for the file and data variable.
    integer :: ncid, varid, status, latdimid, londimid

    print*,'landmask file: ', './input/global/grid/'//trim(params_domain%filnam_landmask)
    call check( nf90_open( './input/global/grid/'//trim(params_domain%filnam_landmask), NF90_NOWRITE, ncid ) )

    ! if ( trim(params_domain%filnam_landmask)=="landmaskfile_global_halfdeg.nc" ) then

    !   print*,'this is a halfdeg simulation'

    !   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
    !   call check( nf90_open( './input/global/grid/gicew_halfdeg.cdf', NF90_NOWRITE, ncid ) )

    ! else if ( trim(params_domain%filnam_landmask)=="landmaskfile_global_1x1deg.nc" ) then

    !   print*,'this is a 1x1deg simulation'

    !   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
    !   call check( nf90_open( './input/global/grid/gicew_1x1deg.cdf', NF90_NOWRITE, ncid ) )

    ! else

    !   print*,'Error: landmask file name unknown'
    !   stop

    ! end if

    ! get dimension ID for latitude
    status = nf90_inq_dimid( ncid, "lat", latdimid )
    if ( status /= nf90_noerr ) then
      status = nf90_inq_dimid( ncid, "latitude", latdimid )
      if ( status /= nf90_noerr ) then
        status = nf90_inq_dimid( ncid, "LAT", latdimid )
        if ( status /= nf90_noerr ) then
          status = nf90_inq_dimid( ncid, "LATITUDE", latdimid )
          if ( status /= nf90_noerr ) then
            print*,'Error: Unknown latitude name.'
            stop
          end if
        end if
      end if
    end if

    ! Get latitude information: nlat
    call check( nf90_inquire_dimension( ncid, latdimid, len = domaininfo%nlat ) )

    ! get dimension ID for longitude
    status = nf90_inq_dimid( ncid, "lon", londimid )
    if ( status /= nf90_noerr ) then
      status = nf90_inq_dimid( ncid, "longitude", londimid )
      if ( status /= nf90_noerr ) then
        status = nf90_inq_dimid( ncid, "LON", londimid )
        if ( status /= nf90_noerr ) then
          status = nf90_inq_dimid( ncid, "LONGITUDE", londimid )
          if ( status /= nf90_noerr ) then
            print*,'Error: Unknown latitude name.'
            stop
          end if
        end if
      end if
    end if

    ! Get latitude information: nlon
    call check( nf90_inquire_dimension( ncid, londimid, len = domaininfo%nlon ) )

    ! Allocate array sizes now knowing nlon and nlat 
    allocate( domaininfo%gridarray(domaininfo%nlon,domaininfo%nlat) )
    allocate( domaininfo%lon(domaininfo%nlon) )
    allocate( domaininfo%lat(domaininfo%nlat) )

    ! Get longitude and latitude values
    call check( nf90_get_var( ncid, londimid, domaininfo%lon ) )
    call check( nf90_get_var( ncid, latdimid, domaininfo%lat ) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid( ncid, "GICEW", varid ) )

    ! Read the grid data
    call check( nf90_get_var( ncid, varid, domaininfo%gridarray ) )

    ! Get total number of land gridcells (value < 1.0 in gicew landmaskfile)
    jpngr = 0
    do ilon=1,domaininfo%nlon
      do ilat=1,domaininfo%nlat
        if (domaininfo%gridarray(ilon,ilat)<1) then
          jpngr = jpngr + 1
        end if
      end do
    end do
    domaininfo%maxgrid = jpngr

    ! get resolution
    domaininfo%dlon = domaininfo%lon(2) - domaininfo%lon(1)
    domaininfo%dlat = domaininfo%lat(2) - domaininfo%lat(1)

  end function get_domaininfo


  function getgrid( domaininfo ) result( out_grid )
    !////////////////////////////////////////////////////////////////
    ! Defines grid variables
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal, area
    use md_params_domain, only: type_params_domain

    ! arguments
    type( domaininfo_type ), intent(inout) :: domaininfo

    ! function return variable
    type( gridtype ), allocatable, dimension(:) :: out_grid

    ! local variables
    integer :: jpngr, ilon, ilat
    integer :: maxgrid_test
    real :: landarea = 0.0

    allocate( out_grid(domaininfo%maxgrid ) )

    jpngr = 0
    do ilon=1,domaininfo%nlon
      do ilat=1,domaininfo%nlat
        if (domaininfo%gridarray(ilon,ilat)<1) then
          jpngr=jpngr+1
          out_grid(jpngr)%ilon = ilon
          out_grid(jpngr)%ilat = ilat
          out_grid(jpngr)%lon = domaininfo%lon(ilon)
          out_grid(jpngr)%lat = domaininfo%lat(ilat)
          out_grid(jpngr)%landfrac = 1.0 - domaininfo%gridarray(ilon,ilat)
          out_grid(jpngr)%area = area( domaininfo%lat(ilat), dx=domaininfo%dlon, dy=domaininfo%dlat )
          landarea = landarea + out_grid(jpngr)%landfrac * out_grid(jpngr)%area
        end if
      end do
    end do

    out_grid(:)%dogridcell = .true.
    out_grid(:)%elv = 100.0
    out_grid(:)%soilcode = 1

    ! complement domaininfo here in order to avoid passing too many large arrays around with domaininfo
    domaininfo%landarea = landarea

  end function getgrid


  subroutine check( status )
    !/////////////////////////////////////////////////////////////////////////
    ! Auxiliary subroutine handling NetCDF 
    !-------------------------------------------------------------------------
    use netcdf
    integer, intent (in) :: status
    if ( status /= nf90_noerr ) then 
      print *, trim( nf90_strerror(status) )
      stop "Stopped"
    end if
  end subroutine check  

end module md_grid

