module md_params_soil
  !////////////////////////////////////////////////////////////////
  ! Module containing soil parameters and functions to read them
  !----------------------------------------------------------------
  use netcdf
  use md_io_netcdf, only: check

  implicit none

  private
  public paramtype_soil, getsoil_field

  type paramtype_soil
    real :: whc
  end type

contains

  function getsoil_field( grid ) result( params_soil_field )
    !////////////////////////////////////////////////////////////////
    ! Function returns array containing all soil parameter values
    !----------------------------------------------------------------
    use md_params_core, only: maxgrid
    use md_grid, only: gridtype

    ! arguments
    type( gridtype ), intent(in), dimension(maxgrid) :: grid

    ! local variables
    integer :: jpngr

    ! function return variable
    type(paramtype_soil), dimension(maxgrid) :: params_soil_field

    if (maxgrid>1) stop 'in getsoil_field: think of something'

    do jpngr=1,maxgrid
      params_soil_field(jpngr) = getsoil( 1 )
    end do

  end function getsoil_field


  function getsoil( soilcode ) result( params_soil )
    !////////////////////////////////////////////////////////////////
    ! Function returns soil parameter values given soil code
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    integer, intent(in) :: soilcode

    ! local variables
    character(len=2) :: soilcode_char

    ! function return variable
    type(paramtype_soil) :: params_soil

    write( soilcode_char, 999 ) soilcode

    params_soil%perc_k1      = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'perc_k1' )
    params_soil%whc_eff      = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'whc_eff' )
    params_soil%thdiff_wp    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_wp' )
    params_soil%thdiff_whc15 = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_whc15' )
    params_soil%thdiff_fc    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_fc' )
    params_soil%forg         = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'forg' )
    params_soil%whc_blwp     = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'whc_blwp' )
    params_soil%por          = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'por' )
    params_soil%fsand        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsand' )
    params_soil%fclay        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fclay' )
    params_soil%fsilt        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsilt' )

    return
    999  format (I2.2)

  end function getsoil


  function getsoil( domaininfo, grid, year, fapar_forcing_source ) result( fapar_field )
    !////////////////////////////////////////////////////////////////
    ! Reads fAPAR from fapar3g data file.
    ! Assumes fAPAR=0 for cells with missing data
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(in) :: grid
    integer, intent(in) :: year
    character(len=*), intent(in) :: fapar_forcing_source

    ! function return variable
    real, dimension(ndayyear,domaininfo%maxgrid) :: fapar_field

    ! local variables
    integer :: ncid, varid
    integer :: latdimid, londimid
    integer :: nlat_arr, nlon_arr
    real, allocatable, dimension(:)     :: lon_arr
    real, allocatable, dimension(:)     :: lat_arr
    real, allocatable, dimension(:,:,:) :: fapar_arr

    integer :: jpngr, ilon_arr, ilat_arr, moy, dom, doy
    integer, dimension(domaininfo%maxgrid) :: ilon
    integer, dimension(domaininfo%maxgrid) :: ilat
    integer :: fileyear, read_idx
    real :: tmp
    real :: ncfillvalue
    real :: dlat, dlon
    character(len=100) :: lonname, latname, varname
    integer :: firstyr_data, nyrs_data
    character(len=100) :: filnam

    !----------------------------------------------------------------  
    ! Set file-specific variables
    !----------------------------------------------------------------    
    if (fapar_forcing_source=="evi_modis") then

      ! fAPAR data from MODIS EVI
      firstyr_data = 2001
      nyrs_data = 15
      lonname ="LON"
      latname = "LAT"
      varname = "EVI_FILLED"
      filnam = "./input/global/fapar/modis_vegetation__LPDAAC__v5__0.5deg_FILLED.nc"

    else if (fapar_forcing_source=="fapar3g") then
      
      ! fAPAR data from fAPAR3g
      firstyr_data = 1982
      nyrs_data = 35
      lonname ="LON"
      latname = "LAT"
      varname = "FAPAR_FILLED"
      filnam = "./input/global/fapar/fAPAR3g_v2_1982_2016_FILLED.nc"

    else

      stop 'getsoil: argument fapar_forcing_source is invalid'

    end if

    !----------------------------------------------------------------  
    ! Read arrays of all months of current year from file  
    !----------------------------------------------------------------    
    call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid ) )

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

    if (dlon/=domaininfo%dlon) stop 'Longitude resolution of fapar (modis evi) input file not identical with model grid.'
    if (dlat/=domaininfo%dlat) stop 'latitude resolution of fapar (modis evi) input file not identical with model grid.'

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
    allocate( fapar_arr(nlon_arr,nlat_arr,nmonth) )

    ! Get the varid of the data variable, based on its name
    call check( nf90_inq_varid( ncid, trim(varname), varid ) )

    ! Read the array, only current year
    read_idx = ( min( max( year - firstyr_data + 1, 1 ), nyrs_data ) - 1 ) * nmonth + 1
    call check( nf90_get_var( ncid, varid, fapar_arr, start=(/1, 1, 1, read_idx/), count=(/nlon_arr, nlat_arr, 1, nmonth/) ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid, varid, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    call check( nf90_close( ncid ) )

    ! read from array to define grid type 
    do jpngr=1,domaininfo%maxgrid
      doy = 0
      do moy=1,nmonth
        do dom=1,ndaymonth(moy)
          doy = doy + 1
          tmp = fapar_arr(ilon(jpngr),ilat(jpngr),moy)
          if ( tmp/=ncfillvalue ) then
            fapar_field(doy,jpngr) = tmp
          else
            fapar_field(doy,jpngr) = 0.0
          end if
        end do
      end do
    end do

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( fapar_arr )

  end function getsoil

end module md_params_soil
