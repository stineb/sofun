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
    getvalreal, monthly2daily_weather, monthly2daily
  use md_grid, only: gridtype, domaininfo_type
  use netcdf
  use md_io_netcdf, only: check

  implicit none

  private
  public getco2, getninput, ninput_type, gettot_ninput, getfapar, &
    getclimate, getlanduse, landuse_type, climate_type, get_fpc_grid

  type climate_type
    real, dimension(ndayyear) :: dtemp  ! deg C
    real, dimension(ndayyear) :: dprec  ! mm d-1
    real, dimension(ndayyear) :: dfsun  ! unitless
    real, dimension(ndayyear) :: dvpd   ! Pa
    real, dimension(ndayyear) :: dppfd  ! mol m-2 d-1
    real, dimension(ndayyear) :: dnetrad! W m-2
  end type climate_type

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
      out_gettot_ninput(jpngr)%dnoy(:) = ninput1(jpngr)%dnoy(:) + ninput2(jpngr)%dnoy(:)
      out_gettot_ninput(jpngr)%dnhx(:) = ninput1(jpngr)%dnhx(:) + ninput2(jpngr)%dnhx(:)
      out_gettot_ninput(jpngr)%dtot(:) = ninput1(jpngr)%dtot(:) + ninput2(jpngr)%dtot(:)
    end do

  end function gettot_ninput
  

  function getfapar( domaininfo, grid, year, fapar_forcing_source ) result( fapar_field )
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
    type( climate_type ), dimension(domaininfo%maxgrid) :: out_climate

    out_climate(:) = getclimate_wfdei( &
                                      domaininfo, &
                                      grid, &
                                      init, &
                                      climateyear, &
                                      in_ppfd,  &
                                      in_netrad &
                                      )

    call getclimate_cru( &
                        domaininfo, &
                        grid, &
                        init, &
                        climateyear, &
                        out_climate(:) &
                        )

  end function getclimate


  function getclimate_wfdei( domaininfo, grid, init, climateyear, in_ppfd, in_netrad ) result ( out_climate )
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

    ! function return variable
    type( climate_type ), dimension(domaininfo%maxgrid) :: out_climate

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

    ! create 4-digit string for year  
    write(climateyear_char,999) climateyear

    if (domaininfo%maxgrid>100000) stop 'problem for ilon and ilat length'

    !----------------------------------------------------------------    
    ! Get longitude and latitude information from WATCH-WFDEI file
    !----------------------------------------------------------------    
    if (init) then

      write(moy_char,888) moy
      filnam = './input/global/climate/temp/Tair_daily_WFDEI_'//climateyear_char//'01.nc'

      ! out_arrsize_2D = get_arrsize_2D( filnam )

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

      ! Open the file
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_temp ) )

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
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_temp ) )

      ! precipitation (rain)
      filnam = './input/global/climate/prec/Rainf_daily_WFDEI_CRU_'//climateyear_char//moy_char//'.nc'
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_prec ) )

      ! precipitation (snow)
      filnam = './input/global/climate/prec/Snowf_daily_WFDEI_CRU_'//climateyear_char//moy_char//'.nc'
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_snow ) )

      ! VPD from Qair
      filnam = './input/global/climate/humd/Qair_daily_WFDEI_'//climateyear_char//moy_char//'.nc'
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_humd ) )

      ! PPFD from SWdown
      if (in_ppfd) then
        filnam = './input/global/climate/srad/SWdown_daily_WFDEI_'//climateyear_char//moy_char//'.nc'
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
            out_climate(jpngr)%dtemp(doy) = temp_arr(ilon(jpngr),ilat(jpngr),dom) - 273.15  ! conversion from Kelving to Celsius
            out_climate(jpngr)%dprec(doy) = ( prec_arr(ilon(jpngr),ilat(jpngr),dom) + snow_arr(ilon(jpngr),ilat(jpngr),dom) ) * 60.0 * 60.0 * 24.0  ! kg/m2/s -> mm/day
            out_climate(jpngr)%dvpd(doy)  = calc_vpd( qair_arr(ilon(jpngr),ilat(jpngr),dom), out_climate(jpngr)%dtemp(doy), grid(jpngr)%elv )
            
            ! optional input variables
            if (in_ppfd) then
              out_climate(jpngr)%dppfd(doy) = 1.0e-6 * rswd_arr(ilon(jpngr),ilat(jpngr),dom) * 60.0 * 60.0 * 24.0 * kfFEC ! W m-2 -> mol m-2 d-1
            else
              out_climate(jpngr)%dppfd(doy) = dummy
            end if

            ! if ( in_netrad .and. in_ppfd ) then
            !   out_climate(jpngr)%dfsun(doy) = dummy
            ! else
            !   out_climate(jpngr)%dfsun(doy) = 1111
            ! end if

            ! if (in_netrad) then
            !   out_climate(jpngr)%dnetrad(:) = 1111
            ! else
            !   out_climate(jpngr)%dnetrad(:) = dummy
            ! end if

          else
            nmissing = nmissing + 1
            out_climate(jpngr)%dtemp(doy) = dummy
            out_climate(jpngr)%dprec(doy) = dummy
            out_climate(jpngr)%dppfd(doy) = dummy
            out_climate(jpngr)%dvpd (doy) = dummy
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

    return
    888  format (I2.2)
    999  format (I4.4)

  end function getclimate_wfdei


  subroutine getclimate_cru( domaininfo, grid, init, climateyear, inout_climate )
    !////////////////////////////////////////////////////////////////
    ! SR reads this year's daily temperature and precipitation.
    ! Read year-2013 data after 2013
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid
    logical, intent(in) :: init
    integer, intent(in) :: climateyear
    type( climate_type ), dimension(domaininfo%maxgrid), intent(inout) :: inout_climate

    ! local variables
    integer :: doy, dom, moy, read_idx
    integer :: jpngr = 1
    integer :: ncid_ccov
    integer :: varid_ccov
    integer :: latdimid, londimid
    integer, dimension(100000), save :: ilon, ilat
    integer, save :: nlon_arr, nlat_arr, ilat_arr, ilon_arr, nrec_arr
    real, dimension(:,:,:), allocatable :: ccov_arr      ! temperature, array read from NetCDF file in K
    real, dimension(:), allocatable :: lon_arr, lat_arr  ! longitude and latitude vectors from climate NetCDF files
    real :: dlon_clim, dlat_clim                         ! resolution in longitude and latitude in climate input files
    real :: ncfillvalue                                  ! _FillValue attribute in NetCDF file
    integer :: nmissing                                  ! number of land cells where climate data is not available
    real :: tmp
    character(len=5) :: recname = "tstep"
    integer, parameter :: firstyr_cru = 1901
    integer, parameter :: nyrs_cru = 114
    character(len=256), parameter :: filnam = './input/global/climate/ccov/cru_ts3.23.1901.2014.cld.dat.nc'

    if (domaininfo%maxgrid>100000) stop 'problem for ilon and ilat length'

    !----------------------------------------------------------------    
    ! Get longitude and latitude information from CRU file
    !----------------------------------------------------------------    
    if (init) then

      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_ccov ) )

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

      ! Open the file
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_ccov ) )

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
    ! open NetCDF files to get ncid_*
    ! temperature
    call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid_ccov ) )

    ! allocate size of output array
    allocate( ccov_arr(nlon_arr,nlat_arr,nmonth) )

    ! Get the varid of the data variable, based on its name
    call check( nf90_inq_varid( ncid_ccov, "cld",  varid_ccov ) )

    ! Read the full array data (years before 1901 are set to 1901, years after)
    read_idx = ( min( max( climateyear - firstyr_cru + 1, 1 ), nyrs_cru ) - 1 ) * nmonth + 1
    call check( nf90_get_var( ncid_ccov, varid_ccov, ccov_arr, start = (/1,1,read_idx/), count = (/nlon_arr, nlat_arr, nmonth/) ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid_ccov, varid_ccov, "_FillValue", ncfillvalue ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid_ccov, varid_ccov, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    call check( nf90_close( ncid_ccov ) )

    ! read from array to define grid type 
    gridcellloop: do jpngr=1,domaininfo%maxgrid
      doy = 0
      monthloop: do moy=1,nmonth
        domloop: do dom=1,ndaymonth(moy)
          doy = doy + 1
          tmp = ccov_arr(ilon(jpngr),ilat(jpngr),moy)
          if ( tmp/=ncfillvalue ) then
            inout_climate(jpngr)%dfsun(doy) = ( 100.0 - tmp ) / 100.0
          else
            inout_climate(jpngr)%dfsun(doy) = dummy
          end if
        end do domloop
      end do monthloop
    end do gridcellloop

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( ccov_arr )

    return
    888  format (I2.2)
    999  format (I4.4)

  end subroutine getclimate_cru


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


  function calc_vpd( qair, temp, elv ) result( vpd )
    !////////////////////////////////////////////////////////////////////////
    ! Calculates vapor pressure deficit, given air temperature and assuming
    ! standard atmosphere, corrected for elevation above sea level.
    ! Ref:      Allen et al. (1998)
    !-----------------------------------------------------------------------
    use md_sofunutils, only: calc_patm
    use md_params_core, only: kR, kMv, kMa

    ! arguments
    real, intent(in) :: qair    ! specific humidity (g g-1)
    real, intent(in) :: temp    ! temperature (degrees Celsius)
    real, intent(in) :: elv     ! elevation above sea level (m)

    ! function return variable
    real :: vpd                 ! vapor pressure deficit (Pa)

    ! local variables
    real :: wair    ! mass mising ratio of water vapor to dry air (dimensionless)
    real :: patm    ! atmopheric pressure (Pa)
    real :: rv      ! specific gas constant of water vapor (J g-1 K-1)
    real :: rd      ! specific gas constant of dry air (J g-1 K-1)
    real :: eact    ! actual water vapor pressure (Pa)
    real :: esat    ! saturation water vapor pressure (Pa)


    ! calculate the mass mising ratio of water vapor to dry air (dimensionless)
    wair = qair / ( 1 - qair )

    ! calculate atmopheric pressure (Pa) assuming standard conditions at sea level (elv=0)
    patm = calc_patm( elv )

    ! calculate water vapor pressure 
    rv = kR / kMv
    rd = kR / kMa
    eact = patm * wair * rv / (rd + wair * rv)

    ! calculate saturation water vapour pressure in Pa
    esat = 611.0 * exp( (17.27 * temp)/(temp + 237.3) )

    ! VPD is the difference between actual and saturation vapor pressure
    vpd = esat - eact

    ! Set negative VPD to zero
    vpd = max( 0.0, vpd )

  end function calc_vpd

end module md_forcing

