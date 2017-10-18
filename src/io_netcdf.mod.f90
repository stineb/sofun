module md_io_netcdf
  !////////////////////////////////////////////////////////////////
  ! Some additional wrapper functions
  !----------------------------------------------------------------    
  use netcdf
  use md_params_core, only: dummy

  implicit none
  private
  ! public arrsize_2D_type, get_arrsize_2D, 
  public check, init_nc_3D, write_nc_3D

  ! type arrsize_2D_type
  !   integer :: nlon
  !   integer :: nlat
  ! end type arrsize_2D_type

contains

  subroutine init_nc_3D( filnam, nlon, nlat, lon, lat, outyear, varnam, varunits, longnam, title )
    !////////////////////////////////////////////////////////////////
    ! Subroutine to initialise a NetCDF file with one variable
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear

    ! arguments
    character(len=*), intent(in) :: filnam

    integer, intent(in) :: nlon
    integer, intent(in) :: nlat
    
    real,    dimension(nlon), intent(in) :: lon
    real,    dimension(nlat), intent(in) :: lat
    integer, intent(in) :: outyear

    character(len=*), intent(in) :: varnam
    character(len=*), intent(in) :: varunits
    character(len=*), intent(in) :: longnam
    character(len=*), intent(in) :: title

    ! local variables
    integer :: ncid
    integer, parameter :: nz = 1
    integer, parameter :: zvals = 1
    integer, parameter :: ndims = 4
    integer :: dimids(ndims)

    integer :: londimid, latdimid, zdimid, tdimid
    integer :: varid_lat, varid_lon, varid_z, varid_t
    integer :: varid_var
    integer :: doy

    character(len=*), parameter :: LAT_NAME  = "lat"
    character(len=*), parameter :: LON_NAME  = "lon"
    character(len=*), parameter :: Z_NAME    = "z"
    character(len=*), parameter :: T_NAME    = "time"
    character(len=*), parameter :: UNITS     = "units"
    character(len=*), parameter :: LAT_UNITS = "degrees_north"
    character(len=*), parameter :: LON_UNITS = "degrees_east"
    character(len=*), parameter :: T_UNITS   = "days since 2001-1-1 0:0:0"

    integer, dimension(ndayyear) :: tvals

    ! create time values as integers counting days since 1 Jan 2001 (assuming no leap years)
    tvals = (/ ( ((outyear-2001)*ndayyear+doy-1), doy = 1, ndayyear) /)

    call check( nf90_create( trim(filnam), NF90_CLOBBER, ncid ) )

    ! Define the dimensions. NetCDF will hand back an ID for each. 
    call check( nf90_def_dim( ncid, LON_NAME, nlon,     londimid  ) )
    call check( nf90_def_dim( ncid, LAT_NAME, nlat,     latdimid  ) )
    call check( nf90_def_dim( ncid, Z_NAME,   nz,       zdimid    ) )
    call check( nf90_def_dim( ncid, T_NAME,   ndayyear, tdimid  ) )

    ! Define the coordinate variables. They will hold the coordinate
    ! information, that is, the latitudes and longitudes. A varid is
    ! returned for each.
    call check( nf90_def_var( ncid, LAT_NAME, NF90_REAL, latdimid,  varid_lat ) )
    call check( nf90_def_var( ncid, LON_NAME, NF90_REAL, londimid,  varid_lon ) )
    call check( nf90_def_var( ncid, Z_NAME,   NF90_INT,  zdimid,    varid_z   ) )
    call check( nf90_def_var( ncid, T_NAME,   NF90_INT , tdimid,    varid_t   ) )

    ! Assign units attributes to coordinate var data. This attaches a
    ! text attribute to each of the coordinate variables, containing the
    ! units.
    call check( nf90_put_att( ncid, varid_lat, UNITS, LAT_UNITS ) )
    call check( nf90_put_att( ncid, varid_lon, UNITS, LON_UNITS ) )
    call check( nf90_put_att( ncid, varid_t  , UNITS, T_UNITS   ) )
    call check( nf90_put_att( ncid, varid_t  , "calendar", "noleap"   ) )

    ! The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids = (/ londimid, latdimid, zdimid, tdimid /)

    ! Define the variable. The type of the variable in this case is
    ! NF90_DOUBLE.
    call check( nf90_def_var( ncid, varnam, NF90_REAL, dimids, varid_var ) )

    ! Define some attributes
    ! variable-specific
    call check( nf90_put_att( ncid, varid_var, UNITS, varunits ) )
    call check( nf90_put_att( ncid, varid_var, "_FillValue", dummy ) )
    call check( nf90_put_att( ncid, varid_var, "long_name", longnam ) )

    ! global
    call check( nf90_put_att( ncid, NF90_GLOBAL, "title", title ) )

    ! End define mode. This tells netCDF we are done defining metadata.
    call check( nf90_enddef( ncid ) )

    ! Write the coordinate variable data. This will put the latitudes
    ! and longitudes of our data grid into the netCDF file.
    call check( nf90_put_var( ncid, varid_lat, lat   ) )
    call check( nf90_put_var( ncid, varid_lon, lon   ) )
    call check( nf90_put_var( ncid, varid_z,   zvals ) )
    call check( nf90_put_var( ncid, varid_t,   tvals ) )

    ! Close the file. This frees up any internal netCDF resources
    ! associated with the file, and flushes any buffers.
    call check( nf90_close( ncid ) )

  end subroutine init_nc_3D


  subroutine write_nc_3D( filnam, varnam, maxgrid, nlon, nlat, ilon, ilat, ndays, out )
    !////////////////////////////////////////////////////////////////
    ! Subroutine to put data into an already opened NetCDF file
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filnam
    character(len=*), intent(in) :: varnam
    integer, intent(in) :: maxgrid, nlon, nlat, ndays
    integer, dimension(maxgrid), intent(in) :: ilon, ilat
    real, dimension(ndays,maxgrid), intent(in) :: out

    ! local variables
    real, dimension(:,:,:,:), allocatable :: outarr
    integer :: jpngr
    integer :: ncid, varid

    allocate( outarr(nlon,nlat,1,ndays) )
    outarr(:,:,:,:) = dummy        

    ! Populate output array
    do jpngr=1,maxgrid
        outarr(ilon(jpngr),ilat(jpngr),1,:) = out(:,jpngr)
    end do

    ! open NetCDF output file
    call check( nf90_open( trim(filnam), NF90_WRITE, ncid ) )

    ! Get the varid of the data variable, based on its name
    call check( nf90_inq_varid( ncid, varnam, varid ) )

    ! write the data into the file
    call check( nf90_put_var( ncid, varid, outarr(:,:,:,:) ) )

    ! close NetCDF output file
    call check( nf90_close( ncid ) )

    ! deallocate memory
    deallocate( outarr )

  end subroutine write_nc_3D

  ! function get_arrsize_2D( filnam ) result( out_arrsize )
  !   !----------------------------------------------------------------    
  !   ! returns the longitude and latitude dimension lengths of a file
  !   !----------------------------------------------------------------    
  !   ! arguments
  !   character(len=*) :: filnam

  !   ! function return variable
  !   type( arrsize_2D_type ) :: out_arrsize

  !   ! local variables
  !   integer :: ncid

  !   call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid ) )

  !   ! get dimension ID for latitude
  !   status = nf90_inq_dimid( ncid, "lat", latdimid )
  !   if ( status /= nf90_noerr ) then
  !     status = nf90_inq_dimid( ncid, "latitude", latdimid )
  !     if ( status /= nf90_noerr ) then
  !       status = nf90_inq_dimid( ncid, "LAT", latdimid )
  !       if ( status /= nf90_noerr ) then
  !         status = nf90_inq_dimid( ncid, "LATITUDE", latdimid )
  !         if ( status /= nf90_noerr ) then
  !           print*,'Error: Unknown latitude name.'
  !           stop
  !         end if
  !       end if
  !     end if
  !   end if

  !   ! Get latitude information: nlat
  !   call check( nf90_inquire_dimension( ncid, latdimid, len = out_arrsize%nlat ) )

  !   ! get dimension ID for longitude
  !   status = nf90_inq_dimid( ncid, "lon", londimid )
  !   if ( status /= nf90_noerr ) then
  !     status = nf90_inq_dimid( ncid, "longitude", londimid )
  !     if ( status /= nf90_noerr ) then
  !       status = nf90_inq_dimid( ncid, "LON", londimid )
  !       if ( status /= nf90_noerr ) then
  !         status = nf90_inq_dimid( ncid, "LONGITUDE", londimid )
  !         if ( status /= nf90_noerr ) then
  !           print*,'Error: Unknown latitude name.'
  !           stop
  !         end if
  !       end if
  !     end if
  !   end if

  !   ! Get latitude information: nlon
  !   call check( nf90_inquire_dimension( ncid, londimid, len = out_arrsize%nlon ) )

  !   call check( nf90_close( ncid ) )

  ! end function get_arrsize_2D


  ! function my_nf90_get_var_unknownsize_2D( lonname, latname, varname ) result( outarr )
  !   !////////////////////////////////////////////////////////////////
  !   ! Function to return a 2D array of unknown size
  !   !----------------------------------------------------------------    
  !   ! arguments
  !   character(len=*), intent(in) :: lonname
  !   character(len=*), intent(in) :: latname
  !   character(len=*), intent(in) :: varname

  !   ! function return variable
  !   real, dimension(:,:), allocatable :: outarr

  !   ! local variables
  !   integer :: londimid, latdimid, varid
  !   integer :: nlon, nlat

  !   ! get dimension IDs
  !   call check( nf90_inq_dimid( ncid, lonname, londimid ) )   
  !   call check( nf90_inq_dimid( ncid, latname, latdimid ) )   

  !   ! get dimension lengths
  !   call check( nf90_get_var( ncid, londimid, nlon ) )
  !   call check( nf90_get_var( ncid, londimid, nlat ) )

  !   ! allocate size of output array
  !   allocate( outarr(nlon,nlat) )

  !   ! get data
  !   call check( nf90_inq_varid( ncid, varname, varid ) )
  !   call check( nf90_get_var( ncid, varid, outarr ) )

  ! end function my_nf90_get_var_unknownsize_2D


  ! function my_nf90_get_var_unknownsize_3D( lonname, latname, recname, varname ) result( outarr )
  !   !////////////////////////////////////////////////////////////////
  !   ! Function to return a 3D array of unknown size
  !   !----------------------------------------------------------------    
  !   ! arguments
  !   character(len=*), intent(in) :: lonname
  !   character(len=*), intent(in) :: latname
  !   character(len=*), intent(in) :: recname
  !   character(len=*), intent(in) :: varname

  !   ! function return variable
  !   real, dimension(:,:,:), allocatable :: outarr

  !   ! local variables
  !   integer :: londimid, latdimid, tdimid, varid
  !   integer :: nlon, nlat, nrec

  !   ! get dimension IDs
  !   call check( nf90_inq_dimid( ncid, lonname, londimid ) )   
  !   call check( nf90_inq_dimid( ncid, latname, latdimid ) )   
  !   call check( nf90_inq_dimid( ncid, recname, tdimid ) )   

  !   ! get dimension lengths
  !   call check( nf90_get_var( ncid, londimid, nlon ) )
  !   call check( nf90_get_var( ncid, londimid, nlat ) )
  !   call check( nf90_get_var( ncid, londimid, nrec ) )

  !   ! allocate size of output array
  !   allocate( outarr(nlon,nlat,nrec) )

  !   ! get data
  !   call check( nf90_inq_varid( ncid, varname, varid ) )
  !   call check( nf90_get_var( ncid, varid, outarr ) )

  ! end function my_nf90_get_var_unknownsize_3D


  subroutine check( status )
    
    integer, intent (in) :: status
    
    if ( status /= nf90_noerr ) then 
      print *, trim( nf90_strerror(status) )
      stop "Stopped"
    end if

  end subroutine check  

end module md_io_netcdf