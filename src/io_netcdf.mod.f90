module md_io_netcdf
  !////////////////////////////////////////////////////////////////
  ! Some additional wrapper functions
  !----------------------------------------------------------------    
  use netcdf
  use md_params_core, only: dummy

  implicit none
  private
  public check, init_nc_2D, init_nc_3D_time, write_nc_2D, write_nc_3D_time, init_nc_3D_pft, write_nc_3D_pft

contains

  subroutine init_nc_3D_time( filnam, nlon, nlat, lon, lat, outyear, outdt, outnt, varnam, varunits, longnam, title, &
    globatt1_nam, globatt1_val, &
    globatt2_nam, globatt2_val, &
    globatt3_nam, globatt3_val, &
    globatt4_nam, globatt4_val, &
    globatt5_nam, globatt5_val  &
    )
    !////////////////////////////////////////////////////////////////
    ! Subroutine to initialise a NetCDF file with one variable and lon/lat/time
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear

    ! arguments
    character(len=*), intent(in) :: filnam
    integer, intent(in) :: nlon
    integer, intent(in) :: nlat
    real,    dimension(nlon), intent(in) :: lon
    real,    dimension(nlat), intent(in) :: lat
    integer, intent(in) :: outyear, outdt, outnt
    character(len=*), intent(in) :: varnam
    character(len=*), intent(in) :: varunits
    character(len=*), intent(in) :: longnam
    character(len=*), intent(in) :: title
    character(len=*), intent(in), optional :: globatt1_nam, globatt2_nam, globatt3_nam, globatt4_nam, globatt5_nam
    character(len=*), intent(in), optional :: globatt1_val, globatt2_val, globatt3_val, globatt4_val, globatt5_val

    ! local variables
    integer :: ncid
    integer, parameter :: nz = 1
    integer, parameter :: zvals = 1
    integer, parameter :: ndims = 4
    integer :: dimids(ndims)

    integer :: londimid, latdimid, zdimid, tdimid
    integer :: varid_lat, varid_lon, varid_z, varid_t
    integer :: varid_var
    integer :: it

    character(len=*), parameter :: LAT_NAME  = "lat"
    character(len=*), parameter :: LON_NAME  = "lon"
    character(len=*), parameter :: Z_NAME    = "z"
    character(len=*), parameter :: T_NAME    = "time"
    character(len=*), parameter :: UNITS     = "units"
    character(len=*), parameter :: LAT_UNITS = "degrees_north"
    character(len=*), parameter :: LON_UNITS = "degrees_east"
    character(len=*), parameter :: T_UNITS   = "days since 2001-1-1 0:0:0"

    integer, dimension(outnt) :: tvals

    ! create time values as integers counting days since 1 Jan 2001 (assuming no leap years)
    tvals = (/ ( ((outyear-2001)*ndayyear+(it-1)*outdt), it = 1, outnt) /)

    call check( nf90_create( trim(filnam), NF90_CLOBBER, ncid ) )

    ! Define the dimensions. NetCDF will hand back an ID for each. 
    call check( nf90_def_dim( ncid, LON_NAME, nlon,     londimid  ) )
    call check( nf90_def_dim( ncid, LAT_NAME, nlat,     latdimid  ) )
    call check( nf90_def_dim( ncid, Z_NAME,   nz,       zdimid    ) )
    call check( nf90_def_dim( ncid, T_NAME,   outnt,    tdimid  ) )

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
    ! variable-specific attributes
    call check( nf90_put_att( ncid, varid_var, UNITS, varunits ) )
    call check( nf90_put_att( ncid, varid_var, "_FillValue", dummy ) )
    call check( nf90_put_att( ncid, varid_var, "long_name", longnam ) )
    call check( nf90_put_att( ncid, varid_var, "output_periodicity_d", outdt ) )

    ! global attributes
    call check( nf90_put_att( ncid, NF90_GLOBAL, "title", title ) )
    if (present(globatt1_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt1_nam, globatt1_val ) )
    if (present(globatt2_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt2_nam, globatt2_val ) )
    if (present(globatt3_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt3_nam, globatt3_val ) )
    if (present(globatt4_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt4_nam, globatt4_val ) )
    if (present(globatt5_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt5_nam, globatt5_val ) )

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

  end subroutine init_nc_3D_time


  subroutine init_nc_3D_pft( filnam, nlon, nlat, lon, lat, outnz, varnam, varunits, longnam, title, &
    globatt1_nam, globatt1_val, &
    globatt2_nam, globatt2_val, &
    globatt3_nam, globatt3_val, &
    globatt4_nam, globatt4_val, &
    globatt5_nam, globatt5_val  &
    )
    !////////////////////////////////////////////////////////////////
    ! Subroutine to initialise a NetCDF file with one variable and lon/lat/time
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear

    ! arguments
    character(len=*), intent(in) :: filnam
    integer, intent(in) :: nlon
    integer, intent(in) :: nlat
    real,    dimension(nlon), intent(in) :: lon
    real,    dimension(nlat), intent(in) :: lat
    integer, intent(in) :: outnz
    character(len=*), intent(in) :: varnam
    character(len=*), intent(in) :: varunits
    character(len=*), intent(in) :: longnam
    character(len=*), intent(in) :: title
    character(len=*), intent(in), optional :: globatt1_nam, globatt2_nam, globatt3_nam, globatt4_nam, globatt5_nam
    character(len=*), intent(in), optional :: globatt1_val, globatt2_val, globatt3_val, globatt4_val, globatt5_val

    ! local variables
    integer :: ncid
    integer, parameter :: nt = 1
    integer, parameter :: tvals = 1
    integer, parameter :: ndims = 4
    integer :: dimids(ndims)

    integer :: londimid, latdimid, zdimid, tdimid
    integer :: varid_lat, varid_lon, varid_z, varid_t
    integer :: varid_var
    integer :: it

    character(len=*), parameter :: LAT_NAME  = "lat"
    character(len=*), parameter :: LON_NAME  = "lon"
    character(len=*), parameter :: Z_NAME    = "pft"
    character(len=*), parameter :: T_NAME    = "t"
    character(len=*), parameter :: UNITS     = "units"
    character(len=*), parameter :: LAT_UNITS = "degrees_north"
    character(len=*), parameter :: LON_UNITS = "degrees_east"
    character(len=*), parameter :: Z_UNITS   = "PFT number"

    integer, dimension(outnz) :: zvals

    ! create time values as integers counting days since 1 Jan 2001 (assuming no leap years)
    zvals = (/ ( it, it = 1, outnz) /)

    call check( nf90_create( trim(filnam), NF90_CLOBBER, ncid ) )

    ! Define the dimensions. NetCDF will hand back an ID for each. 
    call check( nf90_def_dim( ncid, LON_NAME, nlon,     londimid  ) )
    call check( nf90_def_dim( ncid, LAT_NAME, nlat,     latdimid  ) )
    call check( nf90_def_dim( ncid, Z_NAME,   outnz,    zdimid    ) )
    call check( nf90_def_dim( ncid, T_NAME,   nt,       tdimid    ) )

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
    call check( nf90_put_att( ncid, varid_z  , UNITS, Z_UNITS   ) )

    ! The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids = (/ londimid, latdimid, zdimid, tdimid /)

    ! Define the variable. The type of the variable in this case is
    ! NF90_DOUBLE.
    call check( nf90_def_var( ncid, varnam, NF90_REAL, dimids, varid_var ) )

    ! Define some attributes
    ! variable-specific attributes
    call check( nf90_put_att( ncid, varid_var, UNITS, varunits ) )
    call check( nf90_put_att( ncid, varid_var, "_FillValue", dummy ) )
    call check( nf90_put_att( ncid, varid_var, "long_name", longnam ) )

    ! global attributes
    call check( nf90_put_att( ncid, NF90_GLOBAL, "title", title ) )
    if (present(globatt1_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt1_nam, globatt1_val ) )
    if (present(globatt2_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt2_nam, globatt2_val ) )
    if (present(globatt3_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt3_nam, globatt3_val ) )
    if (present(globatt4_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt4_nam, globatt4_val ) )
    if (present(globatt5_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt5_nam, globatt5_val ) )

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

  end subroutine init_nc_3D_pft


  subroutine init_nc_2D( filnam, nlon, nlat, lon, lat, varnam, varunits, longnam, title, &
    globatt1_nam, globatt1_val, &
    globatt2_nam, globatt2_val, &
    globatt3_nam, globatt3_val, &
    globatt4_nam, globatt4_val, &
    globatt5_nam, globatt5_val  &
    )
    !////////////////////////////////////////////////////////////////
    ! Subroutine to initialise a NetCDF file with one variable and lon/lat (no time)
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear

    ! arguments
    character(len=*), intent(in) :: filnam
    integer, intent(in) :: nlon
    integer, intent(in) :: nlat
    real,    dimension(nlon), intent(in) :: lon
    real,    dimension(nlat), intent(in) :: lat
    character(len=*), intent(in) :: varnam
    character(len=*), intent(in) :: varunits
    character(len=*), intent(in) :: longnam
    character(len=*), intent(in) :: title
    character(len=*), intent(in), optional :: globatt1_nam, globatt2_nam, globatt3_nam, globatt4_nam, globatt5_nam
    character(len=*), intent(in), optional :: globatt1_val, globatt2_val, globatt3_val, globatt4_val, globatt5_val

    ! local variables
    integer :: ncid
    integer, parameter :: nz = 1
    integer, parameter :: zvals = 1
    integer, parameter :: nt = 1
    integer, parameter :: tvals = 1
    integer, parameter :: ndims = 4
    integer :: dimids(ndims)

    integer :: londimid, latdimid, zdimid, tdimid
    integer :: varid_lat, varid_lon, varid_z, varid_t
    integer :: varid_var

    character(len=*), parameter :: LAT_NAME  = "lat"
    character(len=*), parameter :: LON_NAME  = "lon"
    character(len=*), parameter :: Z_NAME    = "z"
    character(len=*), parameter :: T_NAME    = "t"
    character(len=*), parameter :: UNITS     = "units"
    character(len=*), parameter :: LAT_UNITS = "degrees_north"
    character(len=*), parameter :: LON_UNITS = "degrees_east"

    call check( nf90_create( trim(filnam), NF90_CLOBBER, ncid ) )

    ! Define the dimensions. NetCDF will hand back an ID for each. 
    call check( nf90_def_dim( ncid, LON_NAME, nlon, londimid  ) )
    call check( nf90_def_dim( ncid, LAT_NAME, nlat, latdimid  ) )
    call check( nf90_def_dim( ncid, Z_NAME,   nz,   zdimid    ) )
    call check( nf90_def_dim( ncid, T_NAME,   nt,   tdimid  ) )

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

    ! global attributes
    call check( nf90_put_att( ncid, NF90_GLOBAL, "title", title ) )
    if (present(globatt1_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt1_nam, globatt1_val ) )
    if (present(globatt2_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt2_nam, globatt2_val ) )
    if (present(globatt3_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt3_nam, globatt3_val ) )
    if (present(globatt4_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt4_nam, globatt4_val ) )
    if (present(globatt5_nam)) call check( nf90_put_att( ncid, NF90_GLOBAL, globatt5_nam, globatt5_val ) )

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

  end subroutine init_nc_2D


  subroutine write_nc_3D_time( filnam, varnam, maxgrid, nlon, nlat, ilon, ilat, outnt, dogridcell, out )
    !////////////////////////////////////////////////////////////////
    ! Subroutine to put data into an already opened NetCDF file
    ! Output array has dimensions (lon, lat, time).
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filnam
    character(len=*), intent(in) :: varnam
    integer, intent(in) :: maxgrid, nlon, nlat, outnt
    integer, dimension(maxgrid), intent(in) :: ilon, ilat
    logical, dimension(maxgrid), intent(in) :: dogridcell
    real, dimension(outnt,maxgrid), intent(in) :: out

    ! local variables
    real, dimension(:,:,:,:), allocatable :: outarr
    integer :: jpngr
    integer :: ncid, varid

    allocate( outarr(nlon,nlat,1,outnt) )
    outarr(:,:,:,:) = dummy        

    ! Populate output array
    do jpngr=1,maxgrid
      if (dogridcell(jpngr)) outarr(ilon(jpngr),ilat(jpngr),1,:) = out(:,jpngr)
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

  end subroutine write_nc_3D_time


  subroutine write_nc_3D_pft( filnam, varnam, maxgrid, nlon, nlat, ilon, ilat, outnz, dogridcell, out )
    !////////////////////////////////////////////////////////////////
    ! Subroutine to put data into an already opened NetCDF file
    ! Output array has dimensions (lon, lat, time).
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filnam
    character(len=*), intent(in) :: varnam
    integer, intent(in) :: maxgrid, nlon, nlat, outnz
    integer, dimension(maxgrid), intent(in) :: ilon, ilat
    logical, dimension(maxgrid), intent(in) :: dogridcell
    real, dimension(outnz,maxgrid), intent(in) :: out

    ! local variables
    real, dimension(:,:,:,:), allocatable :: outarr
    integer :: jpngr
    integer :: ncid, varid

    allocate( outarr(nlon,nlat,outnz,1) )
    outarr(:,:,:,:) = dummy        

    ! Populate output array
    do jpngr=1,maxgrid
      if (dogridcell(jpngr)) outarr(ilon(jpngr),ilat(jpngr),:,1) = out(:,jpngr)
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

  end subroutine write_nc_3D_pft


  subroutine write_nc_2D( filnam, varnam, maxgrid, nlon, nlat, ilon, ilat, dogridcell, out )
    !////////////////////////////////////////////////////////////////
    ! Subroutine to put data into an already opened NetCDF file
    ! Output array has dimensions (lon, lat).
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filnam
    character(len=*), intent(in) :: varnam
    integer, intent(in) :: maxgrid, nlon, nlat
    integer, dimension(maxgrid), intent(in) :: ilon, ilat
    logical, dimension(maxgrid), intent(in) :: dogridcell
    real, dimension(maxgrid), intent(in) :: out

    ! local variables
    real, dimension(:,:,:,:), allocatable :: outarr
    integer :: jpngr
    integer :: ncid, varid

    allocate( outarr(nlon,nlat,1,1) )
    outarr(:,:,:,:) = dummy        

    ! Populate output array
    do jpngr=1,maxgrid
      if (dogridcell(jpngr)) outarr(ilon(jpngr),ilat(jpngr),1,1) = out(jpngr)
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

  end subroutine write_nc_2D


  subroutine check( status )
    
    integer, intent (in) :: status
    
    if ( status /= nf90_noerr ) then 
      print *, trim( nf90_strerror(status) )
      stop "Stopped"
    end if

  end subroutine check  

end module md_io_netcdf
