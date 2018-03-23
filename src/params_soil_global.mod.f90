module md_params_soil
  !////////////////////////////////////////////////////////////////
  ! Module containing soil parameters and functions to read them
  !----------------------------------------------------------------
  use md_grid, only: gridtype, domaininfo_type
  use netcdf
  use md_io_netcdf, only: check

  implicit none

  private
  public paramtype_soil, getsoil

  type paramtype_soil
    real :: whc            ! water holding capacity of the soil, from SoilGrids as (FC-PWP) * fgravel * min( 2m, soildepth )
    real :: perc_k1
    real :: thdiff_wp
    real :: thdiff_whc15
    real :: thdiff_fc
    real :: forg
    real :: por
    real :: fsand
    real :: fclay
    real :: fsilt
  end type

contains

  function getsoil( domaininfo, grid ) result( params_soil_field )
    !////////////////////////////////////////////////////////////////
    ! Function returns the field of soil parameters
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(in) :: grid

    ! function return variable
    type(paramtype_soil), dimension(domaininfo%maxgrid) :: params_soil_field

    ! local variables
    type(paramtype_soil), dimension(9) :: soilparams_per_code

    integer :: ncid, varid
    integer :: latdimid, londimid
    integer :: nlat_arr, nlon_arr
    real, allocatable, dimension(:)     :: lon_arr
    real, allocatable, dimension(:)     :: lat_arr
    real, allocatable, dimension(:,:) :: soilparams_arr

    integer :: isoilcode, jpngr, ilon_arr, ilat_arr, moy, dom, doy
    integer, dimension(domaininfo%maxgrid) :: ilon
    integer, dimension(domaininfo%maxgrid) :: ilat
    integer :: fileyear, read_idx
    real :: tmp
    real :: ncfillvalue
    real :: dlat, dlon
    character(len=3), parameter :: lonname = "LON"
    character(len=3), parameter :: latname = "LAT" 
    character(len=100), parameter :: varname_whc = "WHC_FILLED"
    character(len=100), parameter :: varname_soilcode = "SOILTYPE"
    character(len=100), parameter :: filnam_whc = "./input/global/soil/whc_soilgrids_halfdeg_FILLED.nc"
    character(len=100), parameter :: filnam_soiltype = "./input/global/soil/soil_type_hwsd_halfdeg.cdf"

    !----------------------------------------------------------------  
    ! Get soil parameters per code
    !----------------------------------------------------------------  
    do isoilcode=1,1
      soilparams_per_code(isoilcode) = get_soilparams_per_code( isoilcode )
    end do

    !----------------------------------------------------------------  
    ! Get soil type (1-9) from HWSD (does not define WHC!)
    !----------------------------------------------------------------
    print*,'getting soil type from ', trim(filnam_soiltype), '...'

    ! Read arrays of all months of current year from file  
    call check( nf90_open( trim(filnam_soiltype), NF90_NOWRITE, ncid ) )

    ! get dimension ID for latitude
    call check( nf90_inq_dimid( ncid, trim(latname), latdimid ) )

    ! Get latitude information: nlat
    call check( nf90_inquire_dimension( ncid, latdimid, len = nlat_arr ) )

    ! get dimension ID for longitude
    call check( nf90_inq_dimid( ncid, trim(lonname), londimid ) )

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
    dlon = lon_arr(2) - lon_arr(1)
    dlat = lat_arr(2) - lat_arr(1)

    if (dlon/=domaininfo%dlon) stop 'Longitude resolution of soil input file is not identical with model grid.'
    if (dlat/=domaininfo%dlat) stop 'latitude resolution of soil input file is not identical with model grid.'

    ! get index associations
    do jpngr=1,domaininfo%maxgrid
      ilon_arr = 1
      do while (grid(jpngr)%lon/=lon_arr(ilon_arr))
        ilon_arr = ilon_arr + 1
      end do
      ilon(jpngr) = ilon_arr

      ilat_arr = 1
      do while (grid(jpngr)%lat/=lat_arr(ilat_arr))
        ilat_arr = ilat_arr + 1
      end do
      ilat(jpngr) = ilat_arr
    end do

    ! allocate size of output array
    allocate( soilparams_arr(nlon_arr,nlat_arr) )

    ! Get the varid of the data variable, based on its name
    print*,trim(varname_soilcode)
    call check( nf90_inq_varid( ncid, trim(varname_soilcode), varid ) )

    ! Read the array
    call check( nf90_get_var( ncid, varid, soilparams_arr, start=(/1, 1/), count=(/nlon_arr, nlat_arr/) ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid, varid, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    call check( nf90_close( ncid ) )

    ! read from array to define grid type 
    do jpngr=1,domaininfo%maxgrid
      tmp = soilparams_arr(ilon(jpngr),ilat(jpngr))
      if ( tmp/=ncfillvalue ) then
        params_soil_field(jpngr) = soilparams_per_code( int(tmp) )
      end if
    end do

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( soilparams_arr )
    deallocate( lon_arr )
    deallocate( lat_arr )

    !----------------------------------------------------------------  
    ! Get WHC from SoilGrids 
    !----------------------------------------------------------------  
    print*,'getting water holding capacity from ', trim(filnam_whc), '...'

    ! Read arrays of all months of current year from file  
    call check( nf90_open( trim(filnam_whc), NF90_NOWRITE, ncid ) )

    ! get dimension ID for latitude
    call check( nf90_inq_dimid( ncid, trim(latname), latdimid ) )

    ! Get latitude information: nlat
    call check( nf90_inquire_dimension( ncid, latdimid, len = nlat_arr ) )

    ! get dimension ID for longitude
    call check( nf90_inq_dimid( ncid, trim(lonname), londimid ) )

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
    dlon = lon_arr(2) - lon_arr(1)
    dlat = lat_arr(2) - lat_arr(1)

    if (dlon/=domaininfo%dlon) stop 'Longitude resolution of soil input file is not identical with model grid.'
    if (dlat/=domaininfo%dlat) stop 'latitude resolution of soil input file is not identical with model grid.'

    ! get index associations
    do jpngr=1,domaininfo%maxgrid
      ilon_arr = 1
      do while (grid(jpngr)%lon/=lon_arr(ilon_arr))
        ilon_arr = ilon_arr + 1
      end do
      ilon(jpngr) = ilon_arr

      ilat_arr = 1
      do while (grid(jpngr)%lat/=lat_arr(ilat_arr))
        ilat_arr = ilat_arr + 1
      end do
      ilat(jpngr) = ilat_arr
    end do

    ! allocate size of output array
    allocate( soilparams_arr(nlon_arr,nlat_arr) )

    ! Get the varid of the data variable, based on its name
    call check( nf90_inq_varid( ncid, trim(varname_whc), varid ) )

    ! Read the array
    call check( nf90_get_var( ncid, varid, soilparams_arr, start=(/1, 1/), count=(/nlon_arr, nlat_arr/) ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid, varid, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    call check( nf90_close( ncid ) )

    ! read from array to define grid type 
    do jpngr=1,domaininfo%maxgrid
      tmp = soilparams_arr(ilon(jpngr),ilat(jpngr))
      if ( tmp/=ncfillvalue ) then
        params_soil_field(jpngr)%whc = tmp
      else
        print*, 'WARNING: No WHC data for (ilon, ilat)=', ilon(jpngr), ',', ilat(jpngr)
        print*, 'Assuming WHC = 250 mm'
        params_soil_field(jpngr)%whc = 250.0
      end if
    end do

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( soilparams_arr )
    deallocate( lon_arr )
    deallocate( lat_arr )

    return
    999  format (I2.2)

  end function getsoil


  function get_soilparams_per_code( soilcode ) result( params_soil )
    !////////////////////////////////////////////////////////////////
    ! Function returns soil parameter values given soil code (except WHC)
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    integer, intent(in) :: soilcode

    ! local variables
    character(len=2) :: soilcode_char

    ! function return variable
    type(paramtype_soil) :: params_soil

    write( soilcode_char, 999 ) soilcode

    params_soil%whc          = -9999.9

    params_soil%perc_k1      = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'perc_k1' )
    params_soil%thdiff_wp    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_wp' )
    params_soil%thdiff_whc15 = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_whc15' )
    params_soil%thdiff_fc    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_fc' )
    params_soil%forg         = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'forg' )
    params_soil%por          = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'por' )
    params_soil%fsand        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsand' )
    params_soil%fclay        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fclay' )
    params_soil%fsilt        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsilt' )

    return
    999  format (I2.2)

  end function get_soilparams_per_code

end module md_params_soil
