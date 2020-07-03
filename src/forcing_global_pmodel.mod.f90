module md_forcing
  !////////////////////////////////////////////////////////////////
  ! Module contains forcing variables (climate, co2, ...), and
  ! subroutines used to read forcing input files for a specific year
  ! ('forcingyear'), specifically for site scale simulations.
  ! This module is only used on the level of 'sofun', but not
  ! within 'biosphere', as these variables are passed on to 'biosphere'
  ! as arguments.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: nmonth, ndaymonth, lunat, ndayyear, maxgrid, nlu, dummy
  use md_sofunutils, only: daily2monthly, read1year_daily, read1year_monthly, &
    getvalreal, monthly2daily_weather, monthly2daily, calc_patm
  use md_grid, only: gridtype, domaininfo_type
  use netcdf
  use md_io_netcdf, only: check

  implicit none

  private
  public getco2, getninput, ninput_type, gettot_ninput, getfapar, &
    getclimate, getlanduse, landuse_type, climate_type, get_fpc_grid, &
    vegcover_type

  type climate_type
    real :: dtemp  ! deg C
    real :: dprec  ! mm d-1
    real :: dsnow  ! mm d-1 water equivalents
    real :: dfsun  ! unitless
    real :: dvpd   ! Pa
    real :: dtmin  ! deg C
    real :: dtmax  ! deg C
    real :: dppfd  ! mol m-2 d-1
    real :: dpatm  ! Pa
    real :: dnetrad! W m-2
  end type climate_type

  type vegcover_type
    real :: dfapar ! fraction of absorbed photosynthetically active radiation
  end type vegcover_type

  type landuse_type
    real, dimension(nlu)         :: lu_area
    logical, dimension(ndayyear) :: do_grharvest
  end type landuse_type

  type ninput_type
    real, dimension(ndayyear) :: dnoy
    real, dimension(ndayyear) :: dnhx
    real, dimension(ndayyear) :: dtot
  end type ninput_type

contains

  function getco2( runname, domaininfo, forcingyear, const_co2_year, firstyeartrend, co2_forcing_file ) result( pco2 )
    !////////////////////////////////////////////////////////////////
    !  Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    type( domaininfo_type ), intent(in) :: domaininfo
    integer, intent(in) :: forcingyear
    integer, intent(in) :: const_co2_year
    integer, intent(in) :: firstyeartrend
    character(len=*), intent(in) :: co2_forcing_file

    ! function return variable
    real :: pco2

    ! local variables 
    integer :: readyear

    if (const_co2_year/=int(dummy)) then
      readyear = const_co2_year
    else  
      readyear = forcingyear
    end if
    ! write(0,*) 'GETCO2: use CO2 data of year ', readyear
    print*,'CO2: reading file:', 'input/global/co2/'//trim(co2_forcing_file)
    pco2 = getvalreal( 'global/co2/'//trim(co2_forcing_file), readyear )

  end function getco2


  function getninput( ntype, runname, domaininfo, forcingyear, firstyeartrend, const_ninput_year, ninput_noy_forcing_file, ninput_nhx_forcing_file, climate ) result( out_getninput )
    !////////////////////////////////////////////////////////////////
    ! Dummy function
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: ntype   ! either "nfert" or "ndep"
    character(len=*), intent(in) :: runname
    type( domaininfo_type ), intent(in) :: domaininfo
    integer, intent(in)          :: forcingyear
    integer, intent(in) :: firstyeartrend
    integer, intent(in) :: const_ninput_year
    character(len=*), intent(in) :: ninput_noy_forcing_file
    character(len=*), intent(in) :: ninput_nhx_forcing_file
    type( climate_type ), dimension(maxgrid), intent(in) :: climate

    ! function return variable
    type( ninput_type ), dimension(maxgrid) :: out_getninput 

    ! local variables
    integer :: jpngr
    
    do jpngr=1,maxgrid
      out_getninput(jpngr)%dnoy(:) = dummy
      out_getninput(jpngr)%dnhx(:) = dummy
      out_getninput(jpngr)%dtot(:) = dummy
    end do

  end function getninput


  function gettot_ninput( ninput1, ninput2 ) result( out_gettot_ninput )
    !////////////////////////////////////////////////////////////////
    ! Function returns totals of two ninput type variables with 
    ! dimension maxgrid
    !----------------------------------------------------------------
    ! arguments
    type( ninput_type ), dimension(maxgrid), intent(in) :: ninput1, ninput2 

    ! local variables
    integer :: jpngr

    ! function return variable
    type( ninput_type ), dimension(maxgrid) :: out_gettot_ninput 

    do jpngr=1,maxgrid
      out_gettot_ninput(jpngr)%dnoy(:) = dummy
      out_gettot_ninput(jpngr)%dnhx(:) = dummy
      out_gettot_ninput(jpngr)%dtot(:) = dummy
    end do

  end function gettot_ninput
  

  function getfapar( domaininfo, grid, year, fapar_forcing_source ) result( out_vegcover )
    !////////////////////////////////////////////////////////////////
    ! Reads fAPAR from fapar3g data file.
    ! Assumes fAPAR=0 for cells with missing data
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid
    integer, intent(in) :: year
    character(len=*), intent(in) :: fapar_forcing_source

    ! function return variable
    type( vegcover_type ), dimension(ndayyear,domaininfo%maxgrid) :: out_vegcover

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
    character(len=256) :: filnam

    !----------------------------------------------------------------  
    ! Set file-specific variables
    !----------------------------------------------------------------    
    print*,'fapar_forcing_source: *', trim(fapar_forcing_source), '*' 
    if (trim(fapar_forcing_source)=="evi_modis") then

      ! fAPAR data from MODIS EVI ZMAW data file
      print*,'Using MODIS EVI from ./input/global/fapar/file modis_vegetation__LPDAAC__v5__0.5deg_FILLED.nc ...'
      firstyr_data = 2001
      nyrs_data = 15
      lonname ="LON"
      latname = "LAT"
      varname = "EVI_FILLED"
      filnam = "./input/global/fapar/modis_vegetation__LPDAAC__v5__0.5deg_FILLED.nc"

    else if (trim(fapar_forcing_source)=="fapar3g" .or. trim(fapar_forcing_source)=="fAPAR3g") then
      
      ! fAPAR data from fAPAR3g
      firstyr_data = 1982
      nyrs_data = 35
      lonname ="LON"
      latname = "LAT"
      varname = "FAPAR_FILLED"
      filnam = "./input/global/fapar/fAPAR3g_v2_1982_2016_FILLED.nc"

    else if (trim(fapar_forcing_source)=="fpar_modis") then

      ! fAPAR data from MODIS FPAR ZMAW data file
      firstyr_data = 2000
      nyrs_data = 19
      lonname ="lon"
      latname = "lat"
      varname = "fpar"
      filnam = "./input/global/fapar/MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__2000_2018__MON__fv0.02.nc"

    else

      print*,'fapar_forcing_source: ', fapar_forcing_source
      stop 'getfapar: argument fapar_forcing_source is invalid'

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

    ! file has no data for january 2000
    if (trim(fapar_forcing_source)=="fpar_modis") then
      read_idx = max(1, read_idx - 1)
    end if

    if (trim(fapar_forcing_source)=="evi_modis") then

      call check( nf90_get_var( ncid, varid, fapar_arr, start=(/1, 1, read_idx/), count=(/nlon_arr, nlat_arr, nmonth/) ) )

    else if (trim(fapar_forcing_source)=="fapar3g" .or. trim(fapar_forcing_source)=="fAPAR3g") then

      call check( nf90_get_var( ncid, varid, fapar_arr, start=(/1, 1, 1, read_idx/), count=(/nlon_arr, nlat_arr, 1, nmonth/) ) )

    else if (trim(fapar_forcing_source)=="fpar_modis") then

      call check( nf90_get_var( ncid, varid, fapar_arr, start=(/1, 1, read_idx/), count=(/nlon_arr, nlat_arr, nmonth/) ) )

    else

      print*,'fapar_forcing_source: ', fapar_forcing_source
      stop 'getfapar: argument fapar_forcing_source is invalid'

    end if

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
            out_vegcover(doy,jpngr)%dfapar = tmp
          else
            out_vegcover(doy,jpngr)%dfapar = 0.0
          end if
        end do
      end do
    end do

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( fapar_arr )

  end function getfapar


  function get_fpc_grid( domaininfo, grid, params_siml ) result( fpc_grid_field )
    !////////////////////////////////////////////////////////////////
    ! Function returns the fractional land cover by vegetation types 
    ! based on the 10 IGBP types in the input file (MODIS Landcover)
    ! 1: ENF: type2 = "evergreen needleleaf forest" ;
    ! 2: EBF: type3 = "evergreen broadleaf forest" ;
    ! 3: DNF: type4 = "deciduous needleleaf forest" ;
    ! 4: DBF: type5 = "deciduous broadleaf forest" ;
    ! 5: MF:  type6 = "mixed forest" ;
    ! 6: SHR: type7+type8 = "closed shrublands" + "open shrublands";
    ! 7: SAV: type9+type10 = "savannas" plus "woody savannas"
    ! 8: GRA: type11 = "grasslands" ;
    ! 9: WET: type12 = "permanent wetlands" ;
    ! 10:CRO: type13 + type15 = "croplands" + "cropland (natural vegetation mosaic)";
    !----------------------------------------------------------------
    use md_params_siml, only: paramstype_siml
    use md_params_core, only: npft, eps

    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(in) :: grid
    type( paramstype_siml ), intent(in) :: params_siml

    ! function return variable
    real, dimension(npft,domaininfo%maxgrid) :: fpc_grid_field

    ! local variables
    integer :: ncid, varid
    integer :: latdimid, londimid, pftdimid
    integer :: nlat_arr, nlon_arr, npft_in
    real, allocatable, dimension(:)     :: lon_arr
    real, allocatable, dimension(:)     :: lat_arr
    real, allocatable, dimension(:,:,:) :: vegtype_arr

    integer :: i, pft, jpngr, ilon_arr, ilat_arr, n_noinfo
    integer, dimension(domaininfo%maxgrid) :: ilon
    integer, dimension(domaininfo%maxgrid) :: ilat
    integer :: fileyear, read_idx
    real, allocatable, dimension(:) :: tmp
    real :: ncfillvalue
    real :: dlat, dlon
    character(len=3), parameter :: lonname = "lon"
    character(len=3), parameter :: latname = "lat"
    character(len=100), parameter :: dimname_pft = "z"
    character(len=100), parameter :: varname = "pftcover"
    character(len=100), parameter :: filnam = "./input/global/landcover/modis_landcover_halfdeg_2010_FILLED.nc"

    !----------------------------------------------------------------  
    ! Get vegetation cover information from file
    !----------------------------------------------------------------
    print*,'getting vegetation cover from ', trim(filnam), ' ...'

    ! Read arrays of all months of current year from file  
    call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid ) )

    ! get dimension ID for latitude
    call check( nf90_inq_dimid( ncid, trim(latname), latdimid ) )

    ! Get latitude information: nlat
    call check( nf90_inquire_dimension( ncid, latdimid, len = nlat_arr ) )

    ! get dimension ID for longitude
    call check( nf90_inq_dimid( ncid, trim(lonname), londimid ) )

    ! Get latitude information: nlon
    call check( nf90_inquire_dimension( ncid, londimid, len = nlon_arr ) )

    ! get dimension ID for PFT
    call check( nf90_inq_dimid( ncid, trim(dimname_pft), pftdimid ) )

    ! Get PFT information: number of PFTs
    call check( nf90_inquire_dimension( ncid, pftdimid, len = npft_in ) )

    ! for index association, get ilon and ilat vectors
    ! Allocate array sizes now knowing nlon and nlat 
    allocate( lon_arr(nlon_arr) )
    allocate( lat_arr(nlat_arr) )

    ! print*,'size(lon_arr)', size(lon_arr)
    ! print*,'size(lat_arr)', size(lat_arr)

    ! ! Get longitude and latitude values
    ! print*,'1'
    ! call check( nf90_get_var( ncid, londimid, lon_arr ) )
    ! print*,'2'
    ! call check( nf90_get_var( ncid, latdimid, lat_arr ) )
    ! print*,'3'

    ! xxx try:
    lon_arr = (/ (i, i = 1,nlon_arr) /)
    lon_arr = (lon_arr - 1) * 0.5 - 180.0 + 0.25 
    lat_arr = (/ (i, i = 1,nlat_arr) /)
    lat_arr = (lat_arr - 1) * 0.5 - 90.0 + 0.25 

    ! Check if the resolution of the climate input files is identical to the model grid resolution
    dlon = lon_arr(2) - lon_arr(1)
    dlat = lat_arr(2) - lat_arr(1)

    if (dlon/=domaininfo%dlon) stop 'Longitude resolution of FPC input file is not identical with model grid.'
    if (dlat/=domaininfo%dlat) stop 'latitude resolution of FPC input file is not identical with model grid.'

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
    allocate( vegtype_arr(nlon_arr,nlat_arr,npft_in) )
    allocate( tmp(npft_in))

    ! Get the varid of the data variable, based on its name
    call check( nf90_inq_varid( ncid, trim(varname), varid ) )

    ! Read the array
    call check( nf90_get_var( ncid, varid, vegtype_arr, start=(/1, 1, 1/), count=(/nlon_arr, nlat_arr, npft_in/) ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid, varid, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    call check( nf90_close( ncid ) )

    ! read from array to define land cover 'fpc_grid_field'
    n_noinfo = 0
    do jpngr=1,domaininfo%maxgrid
      
      tmp = vegtype_arr(ilon(jpngr),ilat(jpngr),:)
      
      if ( (tmp(1)==ncfillvalue .or. sum(tmp(:))<eps) .and. grid(jpngr)%dogridcell ) then

        n_noinfo = n_noinfo + 1
        fpc_grid_field(:,jpngr) = 1.0 / real( npft )

      else        

        fpc_grid_field(:,jpngr) = 0.0

        if (npft==2) then
          if (jpngr==1) print*,'GET_FPC_GRID: npft=2 ==> assuming distinction between grasslands/croplands and others'

          pft = 0
          if ( params_siml%lTrE ) then
            ! xxx dirty: call all non-grass vegetation types 'TrE', see indeces above
            pft = pft + 1
  
            ! TrE defined as: 1: ENF, 2: EBF, 3: DNF, 4: DBF, 5: MF:, 6: SHR, 7: SAV, 9: WET     
            fpc_grid_field(pft,jpngr) = sum( tmp(1:7) ) + tmp(9)
          end if

          if ( params_siml%lGr3 ) then
            ! xxx dirty: call all grass vegetation types 'Gr3'
            pft = pft + 1

            ! Gr3 defined as: 8: GRA, 10:CRO: 
            fpc_grid_field(pft,jpngr) = tmp(8) + tmp(10)
          end if

        else if (npft==1) then

          if (jpngr==1) print*,'GET_FPC_GRID: npft=1 ==> assuming no distinction between vegetation types'
          ! 1: ENF, 2: EBF, 3: DNF, 4: DBF, 5: MF:, 6: SHR, 7: SAV, 9: WET     
          pft = 1
          fpc_grid_field(pft,jpngr) = sum( tmp(:) )

        else

          stop 'GET_FPC_GRID: only implemented for npft = 1 or 2.'

        end if


      end if

      if ( abs(sum(fpc_grid_field(:,jpngr)) - 1.0)>eps .and. sum(fpc_grid_field(:,jpngr)) > 0.0 ) &
        fpc_grid_field(:,jpngr) = fpc_grid_field(:,jpngr) / sum( fpc_grid_field(:,jpngr) )

    end do

    print*,'GET_FPC_GRID: number of gridcells with no fpc_grid info:', n_noinfo

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( vegtype_arr )
    deallocate( lon_arr )
    deallocate( lat_arr )

    return
    999  format (I2.2)

  end function get_fpc_grid


  function getclimate( domaininfo, grid, init, climateyear, in_ppfd, in_netrad ) result ( out_climate )
    !////////////////////////////////////////////////////////////////
    ! This function invokes file format specific "sub-functions/routines"
    ! to read from NetCDF. This nesting is necessary because this 
    ! cannot be done file-specific in SR sofun, but can be done here
    ! as this module is compilation-specific (only for global simulations)
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid
    logical, intent(in) :: init
    integer, intent(in) :: climateyear
    logical, intent(in) :: in_ppfd
    logical, intent(in) :: in_netrad

    ! function return variable
    type( climate_type ), dimension(ndayyear,domaininfo%maxgrid) :: out_climate

    ! first get ccov, tmin, and tmax
    out_climate(:,:) =  getclimate_cru( &
                                      domaininfo, &
                                      grid, &
                                      init, &
                                      climateyear &
                                      )

    ! then complement climate derived-type with remaining variables
    call getclimate_wfdei( &
                          domaininfo, &
                          grid, &
                          init, &
                          climateyear, &
                          in_ppfd,  &
                          in_netrad, &
                          out_climate(:,:) &
                          )

    out_climate(:,:)%dpatm = calc_patm( grid(1)%elv )

  end function getclimate


  subroutine getclimate_wfdei( domaininfo, grid, init, climateyear, in_ppfd, in_netrad, inout_climate )
    !////////////////////////////////////////////////////////////////
    ! SR reads this year's daily temperature and precipitation.
    ! Read year-2013 data after 2013
    !----------------------------------------------------------------
    use md_params_core, only: kfFEC

    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid
    logical, intent(in) :: init
    integer, intent(in) :: climateyear
    logical, intent(in) :: in_ppfd
    logical, intent(in) :: in_netrad
    type( climate_type ), dimension(ndayyear,domaininfo%maxgrid), intent(inout) :: inout_climate

    ! local variables
    integer :: doy, dom, moy
    integer :: jpngr = 1
    character(len=4) :: climateyear_char
    character(len=256) :: filnam
    character(len=2) :: moy_char
    integer :: ncid_temp, ncid_prec, ncid_snow, ncid_humd, ncid_fsun, ncid_nrad, ncid_ppfd
    integer :: varid_temp, varid_prec, varid_snow, varid_humd, varid_fsun, varid_nrad, varid_ppfd
    integer :: latdimid, londimid, recdimid, status
    integer, dimension(100000), save :: ilon, ilat
    integer, save :: nlon_arr, nlat_arr, ilat_arr, ilon_arr, nrec_arr
    real, dimension(:,:,:), allocatable :: temp_arr      ! temperature, array read from NetCDF file in K
    real, dimension(:,:,:), allocatable :: prec_arr      ! precipitation, array read from NetCDF file in kg/m2/s
    real, dimension(:,:,:), allocatable :: snow_arr      ! snow fall, array read from NetCDF file in kg/m2/s
    real, dimension(:,:,:), allocatable :: qair_arr      ! specific humidity, array read from NetCDF file 
    real, dimension(:,:,:), allocatable :: fsun_arr      ! sunshine fraction, array read from NetCDF file 
    real, dimension(:,:,:), allocatable :: nrad_arr      ! net radiation, array read from NetCDF file 
    real, dimension(:,:,:), allocatable :: rswd_arr      ! photosynthetic photon flux density, array read from NetCDF file 
    real, dimension(:), allocatable :: lon_arr, lat_arr  ! longitude and latitude vectors from climate NetCDF files
    real :: dlon_clim, dlat_clim                         ! resolution in longitude and latitude in climate input files
    real :: ncfillvalue                                  ! _FillValue attribute in NetCDF file
    integer :: nmissing                                  ! number of land cells where climate data is not available
    character(len=5) :: recname = "tstep"
    logical, parameter :: verbose = .false.

    ! create 4-digit string for year  
    write(climateyear_char,999) climateyear

    if (verbose) print*,'Start getclimate_wfdei() ...'

    if (domaininfo%maxgrid>100000) stop 'problem for ilon and ilat length'

    !----------------------------------------------------------------    
    ! Get longitude and latitude information from WATCH-WFDEI file
    !----------------------------------------------------------------    
    if (init) then

      write(moy_char,888) moy
      filnam = './input/global/climate/temp/Tair_daily_WFDEI_'//climateyear_char//'01.nc'

      ! out_arrsize_2D = get_arrsize_2D( filnam )

      !print*,'Opening climate file to get lon and lat: ', trim(filnam)
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_temp ) )

      ! get dimension ID for latitude
      status = nf90_inq_dimid( ncid_temp, "lat", latdimid )
      if ( status /= nf90_noerr ) then
        status = nf90_inq_dimid( ncid_temp, "latitude", latdimid )
        if ( status /= nf90_noerr ) then
          status = nf90_inq_dimid( ncid_temp, "LAT", latdimid )
          if ( status /= nf90_noerr ) then
            status = nf90_inq_dimid( ncid_temp, "LATITUDE", latdimid )
            if ( status /= nf90_noerr ) then
              print*,'Error: Unknown latitude name.'
              stop
            end if
          end if
        end if
      end if

      ! Get latitude information: nlat
      call check( nf90_inquire_dimension( ncid_temp, latdimid, len = nlat_arr ) )

      ! get dimension ID for longitude
      status = nf90_inq_dimid( ncid_temp, "lon", londimid )
      if ( status /= nf90_noerr ) then
        status = nf90_inq_dimid( ncid_temp, "longitude", londimid )
        if ( status /= nf90_noerr ) then
          status = nf90_inq_dimid( ncid_temp, "LON", londimid )
          if ( status /= nf90_noerr ) then
            status = nf90_inq_dimid( ncid_temp, "LONGITUDE", londimid )
            if ( status /= nf90_noerr ) then
              print*,'Error: Unknown latitude name.'
              stop
            end if
          end if
        end if
      end if

      ! Get latitude information: nlon
      call check( nf90_inquire_dimension( ncid_temp, londimid, len = nlon_arr ) )

      ! for index association, get ilon and ilat vectors
      ! Allocate array sizes now knowing nlon and nlat 
      allocate( lon_arr(nlon_arr) )
      allocate( lat_arr(nlat_arr) )

      ! Get longitude and latitude values
      call check( nf90_get_var( ncid_temp, londimid, lon_arr ) )
      call check( nf90_get_var( ncid_temp, latdimid, lat_arr ) )

      call check( nf90_close( ncid_temp ) )

      ! Check if the resolution of the climate input files is identical to the model grid resolution
      dlon_clim = lon_arr(2) - lon_arr(1)
      dlat_clim = lat_arr(2) - lat_arr(1)
      
      if (dlon_clim/=domaininfo%dlon) stop 'Longitude resolution of climate input file not identical with model grid.'
      if (dlat_clim/=domaininfo%dlat) stop 'latitude resolution of climate input file not identical with model grid.'

      !----------------------------------------------------------------    
      ! Get associations of climate-array gridcells to jpngr (ilon, ilat)
      !----------------------------------------------------------------    
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
    end if

    !----------------------------------------------------------------    
    ! Read climate fields for each month (and day) this year
    !----------------------------------------------------------------
    doy = 0
    monthloop: do moy=1,nmonth

      write(moy_char,888) moy

      ! open NetCDF files to get ncid_*
      ! temperature
      filnam = './input/global/climate/temp/Tair_daily_WFDEI_'//climateyear_char//moy_char//'.nc'
      if (verbose) print*,'opening ', filnam
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_temp ) )

      ! precipitation (rain)
      filnam = './input/global/climate/prec/Rainf_daily_WFDEI_CRU_'//climateyear_char//moy_char//'.nc'
      if (verbose) print*,'opening ', filnam
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_prec ) )

      ! precipitation (snow)
      filnam = './input/global/climate/prec/Snowf_daily_WFDEI_CRU_'//climateyear_char//moy_char//'.nc'
      if (verbose) print*,'opening ', filnam
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_snow ) )

      ! VPD from Qair
      filnam = './input/global/climate/humd/Qair_daily_WFDEI_'//climateyear_char//moy_char//'.nc'
      if (verbose) print*,'opening ', filnam
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_humd ) )

      ! PPFD from SWdown
      if (in_ppfd) then
        filnam = './input/global/climate/srad/SWdown_daily_WFDEI_'//climateyear_char//moy_char//'.nc'
        if (verbose) print*,'opening ', filnam
        call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_ppfd ) )
      end if

      ! get dimension IDs
      call check( nf90_inq_dimid( ncid_temp, recname, recdimid ) )
      call check( nf90_inquire_dimension( ncid_temp, recdimid, len = nrec_arr ) )

      ! allocate size of output array
      allocate( temp_arr(nlon_arr,nlat_arr,nrec_arr) )
      allocate( prec_arr(nlon_arr,nlat_arr,nrec_arr) )
      allocate( snow_arr(nlon_arr,nlat_arr,nrec_arr) )
      allocate( qair_arr(nlon_arr,nlat_arr,nrec_arr) )
      if (in_ppfd) allocate( rswd_arr(nlon_arr,nlat_arr,nrec_arr) )

      ! Get the varid of the data variable, based on its name
      call check( nf90_inq_varid( ncid_temp, "Tair",  varid_temp ) )
      call check( nf90_inq_varid( ncid_prec, "Rainf", varid_prec ) )
      call check( nf90_inq_varid( ncid_snow, "Snowf", varid_snow ) )
      call check( nf90_inq_varid( ncid_humd, "Qair",  varid_humd ) )
      if (in_ppfd) call check( nf90_inq_varid( ncid_ppfd, "SWdown", varid_ppfd ) )

      ! Read the full array data
      call check( nf90_get_var( ncid_temp, varid_temp, temp_arr ) )
      call check( nf90_get_var( ncid_prec, varid_prec, prec_arr ) )
      call check( nf90_get_var( ncid_snow, varid_snow, snow_arr ) )
      call check( nf90_get_var( ncid_humd, varid_humd, qair_arr ) )
      if (in_ppfd) call check( nf90_get_var( ncid_ppfd, varid_ppfd, rswd_arr ) )

      ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
      call check( nf90_get_att( ncid_temp, varid_temp, "_FillValue", ncfillvalue ) )

      ! close NetCDF files
      call check( nf90_close( ncid_temp ) )
      call check( nf90_close( ncid_prec ) )
      call check( nf90_close( ncid_snow ) )
      call check( nf90_close( ncid_humd ) )
      if (in_ppfd) call check( nf90_close( ncid_ppfd ) )

      ! read from array to define climate type 
      domloop: do dom=1,ndaymonth(moy)
        
        doy = doy + 1

        nmissing = 0
        gridcellloop: do jpngr=1,domaininfo%maxgrid

          if ( temp_arr(ilon(jpngr),ilat(jpngr),dom)/=ncfillvalue ) then
            
            ! required input variables
            inout_climate(doy,jpngr)%dtemp = temp_arr(ilon(jpngr),ilat(jpngr),dom) - 273.15  ! conversion from Kelving to Celsius
            inout_climate(doy,jpngr)%dprec = prec_arr(ilon(jpngr),ilat(jpngr),dom) * 60.0 * 60.0 * 24.0  ! kg/m2/s -> mm/day
            inout_climate(doy,jpngr)%dsnow = snow_arr(ilon(jpngr),ilat(jpngr),dom) * 60.0 * 60.0 * 24.0  ! kg/m2/s -> mm/day
            inout_climate(doy,jpngr)%dvpd  = calc_vpd( qair_arr(ilon(jpngr),ilat(jpngr),dom), inout_climate(doy,jpngr)%dtemp, inout_climate(doy,jpngr)%dtmin, inout_climate(doy,jpngr)%dtmax, grid(jpngr)%elv )

            ! optional input variables
            if (in_ppfd) then
              inout_climate(doy,jpngr)%dppfd = 1.0e-6 * rswd_arr(ilon(jpngr),ilat(jpngr),dom) * 60.0 * 60.0 * 24.0 * kfFEC ! W m-2 -> mol m-2 d-1
            else
              inout_climate(doy,jpngr)%dppfd = dummy
            end if

            ! if ( in_netrad .and. in_ppfd ) then
            !   inout_climate(doy,jpngr)%dfsun = dummy
            ! else
            !   inout_climate(doy,jpngr)%dfsun = 1111
            ! end if

            if (in_netrad) then
              inout_climate(doy,jpngr)%dnetrad = 1111.0
            else
              inout_climate(doy,jpngr)%dnetrad = dummy
            end if

          else
            nmissing = nmissing + 1
            inout_climate(doy,jpngr)%dtemp = dummy
            inout_climate(doy,jpngr)%dprec = dummy
            inout_climate(doy,jpngr)%dsnow = dummy
            inout_climate(doy,jpngr)%dppfd = dummy
            inout_climate(doy,jpngr)%dvpd = dummy
            grid(jpngr)%dogridcell = .false.
          end if

        end do gridcellloop

      end do domloop

      ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
      deallocate( temp_arr )
      deallocate( prec_arr )
      deallocate( snow_arr )
      deallocate( qair_arr )
      if (in_ppfd) deallocate( rswd_arr )

    end do monthloop

    if (init) print*,'number of land cells without climate data: ', nmissing
    if (verbose) print*,'done with getclimate_wfdei().'
    return
    888  format (I2.2)
    999  format (I4.4)

  end subroutine getclimate_wfdei


  function getclimate_cru( domaininfo, grid, init, climateyear ) result ( out_climate )
    !////////////////////////////////////////////////////////////////
    ! SR reads this year's daily temperature and precipitation.
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid
    logical, intent(in) :: init
    integer, intent(in) :: climateyear

    ! function return variable
    type( climate_type ), dimension(ndayyear,domaininfo%maxgrid) :: out_climate

    ! local variables
    integer :: doy, dom, moy, read_idx
    integer :: jpngr = 1
    integer :: ncid_ccov, ncid_tmin, ncid_tmax
    integer :: varid_ccov, varid_tmin, varid_tmax
    integer :: latdimid, londimid
    integer, dimension(100000), save :: ilon, ilat
    integer, save :: nlon_arr, nlat_arr, ilat_arr, ilon_arr, nrec_arr
    real, dimension(:,:,:), allocatable :: ccov_arr, tmin_arr, tmax_arr
    real, dimension(:), allocatable :: lon_arr, lat_arr  ! longitude and latitude vectors from climate NetCDF files
    real :: dlon_clim, dlat_clim                         ! resolution in longitude and latitude in climate input files
    real :: ncfillvalue                                  ! _FillValue attribute in NetCDF file
    integer :: nmissing                                  ! number of land cells where climate data is not available
    real :: ccov, tmin, tmax
    character(len=5) :: recname = "tstep"
    integer, parameter :: firstyr_cru = 1901
    integer, parameter :: nyrs_cru = 116
    character(len=256), parameter :: filnam_ccov = './input/global/climate/ccov/cru_ts4.01.1901.2016.cld.dat.nc'
    character(len=256), parameter :: filnam_tmin = './input/global/climate/ccov/cru_ts4.01.1901.2016.tmn.dat.nc'
    character(len=256), parameter :: filnam_tmax = './input/global/climate/ccov/cru_ts4.01.1901.2016.tmx.dat.nc'
    character(len=3),   parameter :: varnamm_ccov = 'cld'
    character(len=3),   parameter :: varnamm_tmin = 'tmn'
    character(len=3),   parameter :: varnamm_tmax = 'tmx'
    logical, parameter :: verbose = .false.

    if (domaininfo%maxgrid>100000) stop 'problem for ilon and ilat length'
    if (verbose) print*,'start with getclimate_cru()....'

    !----------------------------------------------------------------    
    ! Get longitude and latitude information from CRU file
    ! Only for one of the files because lon and lat are identical in all CRU TS files.
    !----------------------------------------------------------------    
    if (init) then

      ! Open the files
      if (verbose) print*,'opening CRU cloud cover file ', trim(filnam_ccov), '...'
      call check( nf90_open( trim(filnam_ccov), NF90_NOWRITE, ncid_ccov ) )

      ! get dimension ID for latitude
      call check( nf90_inq_dimid( ncid_ccov, "lat", latdimid ) )

      ! Get latitude information: nlat
      call check( nf90_inquire_dimension( ncid_ccov, latdimid, len = nlat_arr ) )

      ! get dimension ID for longitude
      call check( nf90_inq_dimid( ncid_ccov, "lon", londimid ) )

      ! Get latitude information: nlon
      call check( nf90_inquire_dimension( ncid_ccov, londimid, len = nlon_arr ) )

      ! for index association, get ilon and ilat vectors
      ! Allocate array sizes now knowing nlon and nlat 
      allocate( lon_arr(nlon_arr) )
      allocate( lat_arr(nlat_arr) )

      ! Get longitude and latitude values
      call check( nf90_get_var( ncid_ccov, londimid, lon_arr ) )
      call check( nf90_get_var( ncid_ccov, latdimid, lat_arr ) )

      call check( nf90_close( ncid_ccov ) )

      ! Check if the resolution of the climate input files is identical to the model grid resolution
      dlon_clim = lon_arr(2) - lon_arr(1)
      dlat_clim = lat_arr(2) - lat_arr(1)
      
      if (dlon_clim/=domaininfo%dlon) stop 'Longitude resolution of cloud cover (CRU) input file not identical with model grid.'
      if (dlat_clim/=domaininfo%dlat) stop 'latitude resolution of cloud cover (CRU) input file not identical with model grid.'

      !----------------------------------------------------------------    
      ! Get associations of climate-array gridcells to jpngr (ilon, ilat)
      !----------------------------------------------------------------    
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

    end if

    !----------------------------------------------------------------    
    ! Read climate fields for each month (and day) this year
    !----------------------------------------------------------------
    ! open files
    if (verbose) print*,'opening CRU cloud cover file ', trim(filnam_ccov), '...'
    call check( nf90_open( trim(filnam_ccov), NF90_NOWRITE, ncid_ccov ) )

    if (verbose) print*,'opening CRU daily minimum temperature file ', trim(filnam_tmin), '...'
    call check( nf90_open( trim(filnam_tmin), NF90_NOWRITE, ncid_tmin ) )

    if (verbose) print*,'opening CRU daily maximum temperature file ', trim(filnam_tmax), '...'
    call check( nf90_open( trim(filnam_tmax), NF90_NOWRITE, ncid_tmax ) )

    ! allocate size of output array
    allocate( ccov_arr(nlon_arr,nlat_arr,nmonth) )
    allocate( tmin_arr(nlon_arr,nlat_arr,nmonth) )
    allocate( tmax_arr(nlon_arr,nlat_arr,nmonth) )

    ! Get the varid_ccov of the data variable, based on its name
    if (verbose) print*,'inquiring variable IDs ...'
    call check( nf90_inq_varid( ncid_ccov, trim(varnamm_ccov),  varid_ccov ) )
    call check( nf90_inq_varid( ncid_tmin, trim(varnamm_tmin),  varid_tmin ) )
    call check( nf90_inq_varid( ncid_tmax, trim(varnamm_tmax),  varid_tmax ) )

    ! Read the full array data (years before 1901 are set to 1901, years after)
    if (verbose) print*,'reading data ...'
    read_idx = ( min( max( climateyear - firstyr_cru + 1, 1 ), nyrs_cru ) - 1 ) * nmonth + 1
    call check( nf90_get_var( ncid_ccov, varid_ccov, ccov_arr, start = (/1,1,read_idx/), count = (/nlon_arr, nlat_arr, nmonth/) ) )
    call check( nf90_get_var( ncid_tmin, varid_tmin, tmin_arr, start = (/1,1,read_idx/), count = (/nlon_arr, nlat_arr, nmonth/) ) )
    call check( nf90_get_var( ncid_tmax, varid_tmax, tmax_arr, start = (/1,1,read_idx/), count = (/nlon_arr, nlat_arr, nmonth/) ) )

    ! Get _FillValue from file (assuming that all are the same for CRU)
    call check( nf90_get_att( ncid_ccov, varid_ccov, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    if (verbose) print*,'closing files ...'
    call check( nf90_close( ncid_ccov ) )
    call check( nf90_close( ncid_tmin ) )
    call check( nf90_close( ncid_tmax ) )

    ! read from array to define grid type 
    if (verbose) print*,'populating climate derived-type ...'
    gridcellloop: do jpngr=1,domaininfo%maxgrid
      doy = 0
      monthloop: do moy=1,nmonth
        ccov = ccov_arr(ilon(jpngr),ilat(jpngr),moy)
        tmin = tmin_arr(ilon(jpngr),ilat(jpngr),moy)
        tmax = tmax_arr(ilon(jpngr),ilat(jpngr),moy)
        domloop: do dom=1,ndaymonth(moy)
          doy = doy + 1
          if ( ccov/=ncfillvalue ) then
            out_climate(doy,jpngr)%dfsun = ( 100.0 - ccov ) / 100.0
            out_climate(doy,jpngr)%dtmin = tmin
            out_climate(doy,jpngr)%dtmax = tmax
          else
            out_climate(doy,jpngr)%dfsun = dummy
            out_climate(doy,jpngr)%dtmin = dummy
            out_climate(doy,jpngr)%dtmax = dummy
          end if
        end do domloop
      end do monthloop
    end do gridcellloop

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( ccov_arr )
    deallocate( tmin_arr )
    deallocate( tmax_arr )

    if (verbose) print*,'... done with getclimate_cru().'

  end function getclimate_cru


  function getlanduse( runname, domaininfo, forcingyear, do_grharvest_forcing_file, const_lu_year, firstyeartrend ) result( out_landuse )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's annual landuse state and harvesting regime (day of above-ground harvest)
    ! Grass harvest forcing file is read for specific year, if none is available,
    ! use earliest forcing file available. 
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    type( domaininfo_type ), intent(in) :: domaininfo
    integer, intent(in)          :: forcingyear
    character(len=*), intent(in), optional :: do_grharvest_forcing_file
    integer, intent(in) :: const_lu_year
    integer, intent(in) :: firstyeartrend

    ! local variables
    integer :: doy
    integer :: findyear
    real, dimension(ndayyear) :: tmp
    character(len=4) :: landuseyear_char
    character(len=245) :: filnam
    integer :: readyear
    logical :: file_exists

    ! function return variable
    type( landuse_type ) :: out_landuse

    ! xxx dummy
    out_landuse%lu_area(lunat)  = 1.0
    out_landuse%do_grharvest(:) = .false.

  end function getlanduse


  function calc_vpd( qair, temp, tmin, tmax, elv ) result( vpd )
    !////////////////////////////////////////////////////////////////////////
    ! Calculates vapor pressure deficit, given air temperature and assuming
    ! standard atmosphere, corrected for elevation above sea level.
    ! Ref:      Allen et al. (1998)
    !-----------------------------------------------------------------------
    use md_sofunutils, only: calc_patm
    use md_params_core, only: kR, kMv, kMa

    ! arguments
    real, intent(in)    :: qair    ! specific humidity (g g-1)
    real, intent(in)    :: temp    ! daily mean air temperature (deg C), daily varying from WATCH-WFDEI (ACTUALLY NOT USED)
    real, intent(inout) :: tmin    ! daily minimum temperature (deg C), constant by month from CRU
    real, intent(inout) :: tmax    ! daily maximum temperature (deg C), constant by month from CRU
    real, intent(in)    :: elv     ! elevation above sea level (m)

    ! function return variable
    real :: vpd         ! vapor pressure deficit as the mean of vpd_min and vpd_max

    ! local variables
    real :: wair         ! mass mixing ratio of water vapor to dry air (dimensionless)
    real :: patm         ! atmopheric pressure (Pa)
    real :: rv           ! specific gas constant of water vapor (J g-1 K-1)
    real :: rd           ! specific gas constant of dry air (J g-1 K-1)
    real :: eact         ! actual water vapor pressure (Pa)
    real :: esat_min     ! saturation water vapor pressure (Pa) based on daily minimum temperature
    real :: esat_max     ! saturation water vapor pressure (Pa) based on daily maximum temperature
    real :: vpd_min      ! vapor pressure deficit based on daily minimum temperature (Pa)
    real :: vpd_max      ! vapor pressure deficit based on daily maximum temperature (Pa)
    real :: halfamp_temp ! half of the daily temperature amplitude

    ! ! xxx don't do this because it probably overestimates extremes a lot
    ! ! take daily temperature half-amplitude from CRU and add it to temp to re-define tmin and tmax
    ! print*,'BEFORE: tmin, tmax ', tmin, tmax
    ! halfamp_temp = (tmax - tmin) / 2.0
    ! tmin = temp - halfamp_temp
    ! tmax = temp + halfamp_temp
    ! print*,'AFTER : tmin, tmax ', tmin, tmax
    ! print*,'----------------------------------'

    ! calculate the mass mising ratio of water vapor to dry air (dimensionless)
    wair = qair / ( 1 - qair )

    ! calculate atmopheric pressure (Pa) assuming standard conditions at sea level (elv=0)
    patm = calc_patm( elv )
    
    ! calculate water vapor pressure 
    rv = kR / kMv
    rd = kR / kMa
    eact = patm * wair * rv / (rd + wair * rv)

    ! calculate saturation water vapour pressure in Pa
    esat_min = 611.0 * exp( (17.27 * tmin)/(tmin + 237.3) )
    esat_max = 611.0 * exp( (17.27 * tmax)/(tmax + 237.3) )

    ! VPD is the difference between actual and saturation vapor pressure
    vpd_min = esat_min - eact
    vpd_max = esat_max - eact
    vpd = (vpd_min + vpd_max) / 2.0 

    ! Set negative VPD to zero
    vpd = max( 0.0, vpd )

  end function calc_vpd

end module md_forcing

