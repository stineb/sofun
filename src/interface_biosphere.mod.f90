module md_interface

  use md_params_core, only: maxgrid, nlu, ndayyear, dummy
  use md_grid, only: gridtype, domaininfo_type
  use md_forcing, only: landuse_type, climate_type, ninput_type
  use md_params_domain, only: type_params_domain
  use md_params_soil, only: paramtype_soil
  use md_params_siml, only: outtype_steering, paramstype_siml

  implicit none

  private
  public interfacetype_biosphere, interface, initoutput_forcing, initio_forcing, &
    initio_nc_forcing, getout_daily_forcing, writeout_ascii_forcing, writeout_nc_forcing

  type interfacetype_biosphere
    integer                                             :: year
    real                                                :: pco2
    type( gridtype )      , dimension(:),   allocatable :: grid
    type( paramtype_soil ), dimension(:),   allocatable :: soilparams
    type( landuse_type)   , dimension(:),   allocatable :: landuse
    type( climate_type )  , dimension(:),   allocatable :: climate
    type( ninput_type)    , dimension(:),   allocatable :: ninput_field
    real                  , dimension(:,:), allocatable :: dfapar_field
    type( domaininfo_type )                             :: domaininfo
    type( outtype_steering )                            :: steering
    type( paramstype_siml )                             :: params_siml
  end type interfacetype_biosphere

  !----------------------------------------------------------------
  ! Interface instance is created here 
  ! (instead of locally defined and passed on as argument. Both are 
  ! ok but this has the advantage that unknown-size arguments are
  ! avoided).
  !----------------------------------------------------------------
  type( interfacetype_biosphere ) :: interface

  !----------------------------------------------------------------
  ! Module-specific daily output variables
  !----------------------------------------------------------------
  real, allocatable, dimension(:,:) :: outdtemp

  !----------------------------------------------------------------
  ! Module-specific NetCDF output file and variable names
  !----------------------------------------------------------------
  character(len=256) :: ncoutfilnam_temp
  character (len = *), parameter :: TEMP_NAME="temperature"

  ! !----------------------------------------------------------------
  ! ! Module-specific annual output variables
  ! !----------------------------------------------------------------
  ! real, dimension(maxgrid)     :: outatemp
  ! real, dimension(nlu,maxgrid) :: outanin


contains

  subroutine initoutput_forcing( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear

    ! arguments
    integer, intent(in) :: ngridcells

    ! Allocate memory for daily output variables
    if ( interface%steering%init .and. interface%params_siml%loutdtemp ) allocate( outdtemp(ndayyear,ngridcells) )
    if ( interface%params_siml%loutdtemp ) outdtemp(:,:) = 0.0

    ! ! annual output variables
    ! if (interface%params_siml%loutforcing) then
    !   outatemp(:) = 0.0
    !   outanin (:,:) = 0.0
    ! end if

  end subroutine initoutput_forcing


  subroutine initio_forcing()
    !////////////////////////////////////////////////////////////////
    ! Opens ascii output files.
    !----------------------------------------------------------------
    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !////////////////////////////////////////////////////////////////
    ! DAILY OUTPUT: OPEN ASCII OUTPUT FILES 
    !----------------------------------------------------------------
    ! DAILY MEAN TEMPERATURE (DEG C)
    if (interface%params_siml%loutdtemp) then
      filnam=trim(prefix)//'.d.temp.out'
      open(950,file=filnam,err=999,status='unknown')
    end if 

    return

    999  stop 'INITIO: error opening output files'

  end subroutine initio_forcing


  subroutine initio_nc_forcing()
    !////////////////////////////////////////////////////////////////
    ! Opens NetCDF output files.
    !----------------------------------------------------------------
    use netcdf

    ! local variables
    character(len=256) :: prefix

    integer :: ncid
    integer, parameter :: ndims = 4

    integer :: londimid, latdimid, doydimid, yeardimid

    character (len = *), parameter :: LAT_NAME = "latitude"
    character (len = *), parameter :: LON_NAME = "longitude"
    character (len = *), parameter :: DOY_NAME = "doy"
    character (len = *), parameter :: YEAR_NAME = "year"

    ! In addition to the latitude and longitude dimensions, we will also
    ! create latitude and longitude netCDF variables which will hold the
    ! actual latitudes and longitudes. Since they hold data about the
    ! coordinate system, the netCDF term for these is: "coordinate
    ! variables."
    integer :: varid_lat, varid_lon, varid_doy, varid_year

    ! It's good practice for each variable to carry a "units" attribute.
    character (len = *), parameter :: UNITS = "units"
    character (len = *), parameter :: TEMP_UNITS = "degrees Celsius"
    character (len = *), parameter :: LAT_UNITS = "degrees_north"
    character (len = *), parameter :: LON_UNITS = "degrees_east"

    integer :: varid_temp
    integer :: dimids(ndims)
    integer :: jpngr, doy
    integer, dimension(ndayyear) :: doy_vals

    doy_vals = (/ (doy, doy = 1, ndayyear) /)

    if (interface%params_siml%lncoutdtemp) then

      prefix = "./output_nc/"//trim(interface%params_siml%runname)

      ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
      ! overwrite this file, if it already exists.
      ncoutfilnam_temp = trim(prefix)//'.d.temp.nc'
      call check( nf90_create( trim(ncoutfilnam_temp), NF90_CLOBBER, ncid ) )

      ! Define the dimensions. NetCDF will hand back an ID for each. 
      call check( nf90_def_dim( ncid, LON_NAME,  interface%domaininfo%nlon, londimid  ) )
      call check( nf90_def_dim( ncid, LAT_NAME,  interface%domaininfo%nlat, latdimid  ) )
      call check( nf90_def_dim( ncid, DOY_NAME,  ndayyear,                  doydimid  ) )
      call check( nf90_def_dim( ncid, YEAR_NAME, 1,                         yeardimid ) )

      ! Define the coordinate variables. They will hold the coordinate
      ! information, that is, the latitudes and longitudes. A varid is
      ! returned for each.
      call check( nf90_def_var( ncid, LAT_NAME,  NF90_REAL, latdimid,  varid_lat ) )
      call check( nf90_def_var( ncid, LON_NAME,  NF90_REAL, londimid,  varid_lon ) )
      call check( nf90_def_var( ncid, DOY_NAME,  NF90_INT,  doydimid,  varid_doy ) )
      call check( nf90_def_var( ncid, YEAR_NAME, NF90_INT,  yeardimid, varid_year ) )

      ! Assign units attributes to coordinate var data. This attaches a
      ! text attribute to each of the coordinate variables, containing the
      ! units.
      call check( nf90_put_att( ncid, varid_lat, UNITS, LAT_UNITS ) )
      call check( nf90_put_att( ncid, varid_lon, UNITS, LON_UNITS ) )

      ! The dimids array is used to pass the IDs of the dimensions of
      ! the variables. Note that in fortran arrays are stored in
      ! column-major format.
      ! Option A
      ! dimids =  (/ latdimid, londimid /)
      ! Option B
      dimids =  (/ londimid, latdimid, doydimid, yeardimid /)

      ! Define the variable. The type of the variable in this case is
      ! NF90_DOUBLE.
      call check( nf90_def_var( ncid, TEMP_NAME, NF90_REAL, dimids, varid_temp ) )

      ! Define some attributes
      ! variable-specific
      call check( nf90_put_att( ncid, varid_temp, UNITS, TEMP_UNITS ) )
      call check( nf90_put_att( ncid, varid_temp, "_FillValue", dummy ) )
      call check( nf90_put_att( ncid, varid_temp, "long_name", "daily average 2 m temperature" ) )

      ! global
      call check( nf90_put_att( ncid, NF90_GLOBAL, "title", "SOFUN GP-model output, module md_interface" ) )

      ! End define mode. This tells netCDF we are done defining metadata.
      call check( nf90_enddef( ncid ) )

      ! Write the coordinate variable data. This will put the latitudes
      ! and longitudes of our data grid into the netCDF file.
      call check( nf90_put_var( ncid, varid_lat,  interface%domaininfo%lat ) )
      call check( nf90_put_var( ncid, varid_lon,  interface%domaininfo%lon ) )
      call check( nf90_put_var( ncid, varid_doy,  doy_vals ) )
      call check( nf90_put_var( ncid, varid_year, interface%steering%outyear ) )

      ! Close the file. This frees up any internal netCDF resources
      ! associated with the file, and flushes any buffers.
      call check( nf90_close( ncid ) )

    end if

  end subroutine initio_nc_forcing


  subroutine getout_daily_forcing( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    if (interface%params_siml%loutdtemp) outdtemp(doy,jpngr) = interface%climate(jpngr)%dtemp(doy)

    ! !----------------------------------------------------------------
    ! ! ANNUAL SUM OVER DAILY VALUES
    ! ! Collect annual output variables as sum of daily values
    ! !----------------------------------------------------------------
    ! if (interface%params_siml%loutforcing) then
    !   outatemp(jpngr)  = outatemp(jpngr)  + interface%climate(jpngr)%dtemp(doy) / ndayyear
    !   outanin(:,jpngr) = outanin(:,jpngr) + interface%ninput_field(jpngr)%dtot(doy)
    ! end if

  end subroutine getout_daily_forcing


  subroutine writeout_ascii_forcing()
    !/////////////////////////////////////////////////////////////////////////
    ! Write daily ASCII output
    !-------------------------------------------------------------------------
    ! use md_params_siml, only: spinup, interface%params_siml%daily_out_startyr, &
    use md_params_core, only: ndayyear

    ! local variables
    real :: itime
    integer :: doy, moy, jpngr
    real, dimension(ndayyear) :: outdtemp_tot

    outdtemp_tot(:) = 0.0

    if (nlu>1) stop 'Output only for one LU category implemented.'

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutdtemp) then
      ! if ( .not. interface%steering%spinup &
      !   .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
      !   .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

        ! Write daily output only during transient simulation
        do doy=1,ndayyear

          ! Get weighted average
          do jpngr=1,size(interface%grid)
            outdtemp_tot(doy) = outdtemp_tot(doy) + outdtemp(doy,jpngr) * interface%grid(jpngr)%landfrac * interface%grid(jpngr)%area
          end do
          outdtemp_tot(doy) = outdtemp_tot(doy) / interface%domaininfo%landarea

          ! Define 'itime' as a decimal number corresponding to day in the year + year
          itime = real( interface%steering%outyear ) + real( doy - 1 ) / real( ndayyear )
          
          write(950,999) itime, outdtemp_tot(doy)

        end do
      ! end if
    end if

    return

    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_forcing


  subroutine writeout_nc_forcing()
    !/////////////////////////////////////////////////////////////////////////
    ! Write NetCDF output
    !-------------------------------------------------------------------------
    use netcdf

    ! local variables
    integer :: doy, jpngr
    integer :: ncid
    integer :: varid_temp

    real, dimension(:,:,:,:), allocatable :: outarr

    ! if ( .not. interface%steering%spinup &
    !       .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
    !       .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      if (interface%params_siml%lncoutdtemp) then

        allocate( outarr(interface%domaininfo%nlon,interface%domaininfo%nlat,ndayyear,1) )
        outarr(:,:,:,:) = dummy        

        ! open NetCDF output file
        call check( nf90_open( trim(ncoutfilnam_temp), NF90_WRITE, ncid ) )

        ! Get the varid of the data variable, based on its name
        call check( nf90_inq_varid( ncid, TEMP_NAME, varid_temp ) )

        ! Write the data, gridcell by gridcell
        do jpngr=1,size(interface%grid)
          if (interface%grid(jpngr)%dogridcell) then

            ! do doy=1,ndayyear
              ! option A
              ! call check( nf90_put_var( ncid, varid_temp, interface%climate(jpngr)%dtemp(1), start = (/ interface%grid(jpngr)%ilat, interface%grid(jpngr)%ilon /) ) )

              ! option B
              ! call check( nf90_put_var( ncid, varid_temp, outdtemp(doy,jpngr), start = (/ doy, interface%grid(jpngr)%ilon, interface%grid(jpngr)%ilat /) ) )            
            ! end do
            ! call check( nf90_put_var( ncid, varid_temp, outdtemp(:,jpngr), start = (/ interface%grid(jpngr)%ilon, interface%grid(jpngr)%ilat, 1 /), count = (/ 1, 1, ndayyear /) ) )

            ! populate array
            outarr(interface%grid(jpngr)%ilon,interface%grid(jpngr)%ilat,:,1) = outdtemp(:,jpngr)

          end if
        end do

        ! write the data into the file
        call check( nf90_put_var( ncid, varid_temp, outarr(:,:,:,:) ) )

        ! close NetCDF output file
        call check( nf90_close( ncid ) )

        ! deallocate memory
        deallocate( outarr )

      end if

      stop 'ok, written NetCDF file now.'

    ! end if

  end subroutine writeout_nc_forcing


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

end module md_interface
