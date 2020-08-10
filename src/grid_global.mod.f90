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
    real, dimension(:,:), allocatable :: gridarray
    real, dimension(:), allocatable :: lon
    real, dimension(:), allocatable :: lat
    real :: landarea
    character(len=256) :: domain_name  ! This is the site name for site-scale simulations or the character identifyier defining the resolution for global simulations
  end type domaininfo_type

  type gridtype
    logical :: dogridcell
    integer :: ilon
    integer :: ilat
    real :: lon
    real :: lat
    real :: elv
    real :: landfrac
    real :: area
    real :: nu               ! true anomaly (orbital parameter), recalculated each year for each gridcell in solar()
    real :: lambda           ! true longitude (orbital parameter), recalculated each year for each gridcell in solar()
    real :: decl_angle       ! declination angle (degrees)
    real :: dayl             ! day length (s), is updated daily in waterbal SR
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

    print*,'getting grid from landmask file: ', './input/global/grid/'//trim(params_domain%filnam_landmask)
    call check( nf90_open( './input/global/grid/'//trim(params_domain%filnam_landmask), NF90_NOWRITE, ncid ) )

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

    ! Read the grid data (warning: this here produces some IEEE exception)
    call check( nf90_get_var( ncid, varid, domaininfo%gridarray ) )

    ! Get total number of land gridcells (value < 1.0 in gicew landmaskfile)
    jpngr = 0
    do ilon=1,domaininfo%nlon
      do ilat=1,domaininfo%nlat
        if (domaininfo%gridarray(ilon,ilat)<1.0) then
          jpngr = jpngr + 1

          ! ! search particular
          ! if (ilon==551 .and. ilat==290) then
          !   print*,'jpngr', jpngr
          !   stop
          ! end if

        end if
      end do
    end do
    domaininfo%maxgrid = jpngr

    ! get resolution
    domaininfo%dlon = domaininfo%lon(2) - domaininfo%lon(1)
    domaininfo%dlat = domaininfo%lat(2) - domaininfo%lat(1)

    print*,'... done'

  end function get_domaininfo


  function getgrid( domaininfo, params_domain ) result( out_grid )
    !////////////////////////////////////////////////////////////////
    ! Defines grid variables
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal, area
    use md_params_domain, only: type_params_domain

    ! arguments
    type( domaininfo_type ), intent(inout) :: domaininfo
    type( type_params_domain ), intent(in) :: params_domain

    ! function return variable
    type( gridtype ), allocatable, dimension(:) :: out_grid

    ! local variables
    integer :: jpngr, ilon, ilat
    integer :: maxgrid_test
    real :: landarea = 0.0
    character(len=256) :: filnam
    integer :: ncid, varid, latdimid, londimid
    integer:: nlon_arr, nlat_arr, ilat_arr, ilon_arr, nrec_arr
    real :: dlon_elv, dlat_elv                           ! resolution in longitude and latitude in climate input files
    integer, dimension(100000), save :: ilon_tmp, ilat_tmp
    real, dimension(:), allocatable :: lon_arr, lat_arr  ! longitude and latitude vectors from climate NetCDF files
    real, dimension(:,:), allocatable :: elv_arr         ! elevation, array read from NetCDF file in m
    real :: ncfillvalue                                  ! _FillValue attribute in NetCDF file
    integer :: nmissing                                  ! number of land cells where climate data is not available

    if (domaininfo%maxgrid>100000) stop 'problem for ilon and ilat length'

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

    ! complement domaininfo here in order to avoid passing too many large arrays around with domaininfo
    domaininfo%landarea = landarea

    !----------------------------------------------------------------    
    ! Get elevation data from WATCH-WFDEI
    !----------------------------------------------------------------    
    ! Get associations of elevation-array gridcells to jpngr (ilon, ilat)
    !----------------------------------------------------------------    
    print*,'getting grid ...'

    call check( nf90_open( "./input/global/grid/"//trim(params_domain%filnam_topography), NF90_NOWRITE, ncid ) )

    ! get dimension ID for latitude
    call check( nf90_inq_dimid( ncid, "lat", latdimid ) )

    ! Get latitude information: nlat
    call check( nf90_inquire_dimension( ncid, latdimid, len = nlat_arr ) )

    ! get dimension ID for longitude
    call check( nf90_inq_dimid( ncid, "lon", londimid ) )

    ! Get latitude information: nlon
    call check( nf90_inquire_dimension( ncid, londimid, len = nlon_arr ) )

    ! for index association, get ilon and ilat vectors
    ! Allocate array sizes now knowing nlon and nlat 
    allocate( lon_arr(nlon_arr) )
    allocate( lat_arr(nlat_arr) )

    ! Get longitude and latitude values
    call check( nf90_get_var( ncid, londimid, lon_arr ) )
    call check( nf90_get_var( ncid, latdimid, lat_arr ) )

    ! Check if the resolution of the climate input files is identical to the model grid resolution
    dlon_elv = lon_arr(2) - lon_arr(1)
    dlat_elv = lat_arr(2) - lat_arr(1)
    
    if (dlon_elv/=domaininfo%dlon) stop 'Longitude resolution of elevation input file not identical with model grid.'
    if (dlat_elv/=domaininfo%dlat) stop 'latitude resolution of elevation input file not identical with model grid.'

    do jpngr=1,domaininfo%maxgrid

      ilon_arr = 1

      do while (out_grid(jpngr)%lon/=lon_arr(ilon_arr))
        ilon_arr = ilon_arr + 1
      end do
      ilon_tmp(jpngr) = ilon_arr

      ilat_arr = 1
      do while (out_grid(jpngr)%lat/=lat_arr(ilat_arr))
        ilat_arr = ilat_arr + 1
      end do
      ilat_tmp(jpngr) = ilat_arr

    end do

    ! allocate size of output array
    allocate( elv_arr(nlon_arr,nlat_arr) )

    ! Get the varid of the data variable, based on its name
    call check( nf90_inq_varid( ncid, "elevation", varid ) )

    ! Read the full array data
    call check( nf90_get_var( ncid, varid, elv_arr ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid, varid, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    call check( nf90_close( ncid ) )

    ! read from array to define grid type 
    nmissing = 0
    do jpngr=1,domaininfo%maxgrid

      if ( elv_arr(ilon_tmp(jpngr),ilat_tmp(jpngr))/=ncfillvalue ) then
        out_grid(jpngr)%elv = elv_arr(ilon_tmp(jpngr),ilat_tmp(jpngr))
      else
        nmissing = nmissing + 1
        out_grid(jpngr)%elv = dummy
        out_grid(jpngr)%dogridcell = .false.
      end if

    end do

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( elv_arr )

    print*,'... done.'

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

