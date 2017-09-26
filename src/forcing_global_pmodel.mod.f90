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

  implicit none

  private
  public getco2, getninput, ninput_type, gettot_ninput, getfapar, getclimate, getlanduse, landuse_type, climate_type

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

  function getco2( runname, sitename, forcingyear, const_co2_year, firstyeartrend, co2_forcing_file ) result( pco2 )
    !////////////////////////////////////////////////////////////////
    !  Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
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
    pco2 = getvalreal( 'sitedata/co2/'//trim(sitename)//'/'//trim(co2_forcing_file), readyear )

  end function getco2


  function getninput( ntype, runname, sitename, forcingyear, firstyeartrend, const_ninput_year, ninput_noy_forcing_file, ninput_nhx_forcing_file, climate ) result( out_getninput )
    !////////////////////////////////////////////////////////////////
    ! Dummy function
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: ntype   ! either "nfert" or "ndep"
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
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


  function getfapar( runname, sitename, ngridcells, grid, forcingyear, fapar_forcing_source ) result( fapar_field )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
    integer, intent(in) :: ngridcells
    type( gridtype ), dimension(ngridcells), intent(in) :: grid
    integer, intent(in)          :: forcingyear
    character(len=*), intent(in) :: fapar_forcing_source

    ! function return variable
    real, dimension(ndayyear,ngridcells) :: fapar_field

    ! local variables 
    integer :: jpngr
    integer :: readyear
    character(len=4) :: faparyear_char

    if (trim(fapar_forcing_source)=='NA') then
      ! If in simulation parameter file 'NA' is specified for 'fapar_forcing_source', then set fapar_field to dummy value
      do jpngr=1,ngridcells
        fapar_field(:,jpngr) = dummy
      end do

    else
      ! Prescribed. Read monthly fAPAR value from file
      do jpngr=1,ngridcells
        ! create 4-digit string for year  
        write(faparyear_char,999) min( max( 2000, forcingyear ), 2014 )
        fapar_field(:,jpngr) = 0.75
      end do

    end if

    return
    999  format (I4.4)

  end function getfapar


  function getclimate( sitename, ngridcells, grid, init, climateyear, in_ppfd, in_netrad ) result ( out_climate )
    !////////////////////////////////////////////////////////////////
    ! SR reads this year's daily temperature and precipitation.
    ! Read year-2013 data after 2013
    !----------------------------------------------------------------    
    ! arguments
    character(len=*), intent(in) :: sitename
    integer, intent(in) :: ngridcells
    type( gridtype ), dimension(ngridcells), intent(in) :: grid
    logical, intent(in) :: init
    integer, intent(in) :: climateyear
    logical, intent(in) :: in_ppfd
    logical, intent(in) :: in_netrad

    ! function return variable
    type( climate_type ), dimension(ngridcells) :: out_climate

    ! local variables
    integer :: doy, dom, moy
    integer :: jpngr = 1
    character(len=4) :: climateyear_char
    character(len=256) :: filnam
    character(len=2) :: moy_char
    integer :: ncid, varid, latdimid, londimid, recdimid, status
    integer, dimension(100000), save :: ilon, ilat
    integer :: ilat_arr, ilon_arr, nlat_arr, nlon_arr, nrec_arr
    real, dimension(:,:,:), allocatable :: temp_arr
    real, dimension(:), allocatable :: lon_arr, lat_arr
    character(len=5) :: recname = "tstep"

    ! create 4-digit string for year  
    write(climateyear_char,999) climateyear

    if (ngridcells>100000) stop 'problem for ilon and ilat length'

    print*,'1'

    if (init) then

      write(moy_char,888) moy
      filnam = './input/global/climate/temp/Tair_daily_WFDEI_'//climateyear_char//'01.nc'
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid ) )

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
      call check( nf90_inquire_dimension( ncid, latdimid, len = nlat_arr ) )

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
      call check( nf90_inquire_dimension( ncid, londimid, len = nlon_arr ) )

      ! for index association, get ilon and ilat vectors
      ! Allocate array sizes now knowing nlon and nlat 
      allocate( lon_arr(nlon_arr) )
      allocate( lat_arr(nlat_arr) )

      ! Get longitude and latitude values
      call check( nf90_get_var( ncid, londimid, lon_arr ) )
      call check( nf90_get_var( ncid, latdimid, lat_arr ) )

      do jpngr=1,ngridcells

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

    doy = 0
    do moy=1,nmonth

      write(moy_char,888) moy

      ! xxx test
      write(moy_char,888) 1

      filnam = './input/global/climate/temp/Tair_daily_WFDEI_'//climateyear_char//moy_char//'.nc'
      call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid ) )

      ! get dimension IDs
      call check( nf90_inq_dimid( ncid, recname, recdimid ) )
      call check( nf90_inquire_dimension( ncid, recdimid, len = nrec_arr ) )

      ! allocate size of output array
      allocate( temp_arr(nlon_arr,nlat_arr,nrec_arr) )

      ! Get the varid of the data variable, based on its name.
      call check( nf90_inq_varid( ncid, "Tair", varid ) )

      ! Read the full array data
      call check( nf90_get_var( ncid, varid, temp_arr ) )

      ! read from array to define climate type 
      do dom=1,ndaymonth(moy)
        
        doy = doy + 1

        do jpngr=1,ngridcells

          out_climate(jpngr)%dtemp(doy) = temp_arr(ilon(jpngr),ilat(jpngr),dom)

        end do

      end do

      stop 'ok'

    end do

    ! xxx test
    do jpngr=1,ngridcells

      out_climate(jpngr)%dprec(:) = 1111
      out_climate(jpngr)%dvpd(:)  = 1111
      if (in_ppfd) then
        out_climate(jpngr)%dppfd(:) = 1111
      else
        out_climate(jpngr)%dppfd(:) = dummy
      end if
      if (in_netrad) then
        out_climate(jpngr)%dnetrad(:) = 1111
      else
        out_climate(jpngr)%dnetrad(:) = dummy
      end if
      if ( in_netrad .and. in_ppfd ) then
        out_climate(jpngr)%dfsun(:) = dummy
      else
        out_climate(jpngr)%dfsun(:) = 1111
      end if

    end do


    return
    888  format (I2.2)
    999  format (I4.4)

  end function getclimate


  function getlanduse( runname, sitename, forcingyear, do_grharvest_forcing_file, const_lu_year, firstyeartrend ) result( out_landuse )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's annual landuse state and harvesting regime (day of above-ground harvest)
    ! Grass harvest forcing file is read for specific year, if none is available,
    ! use earliest forcing file available. 
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
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


  subroutine check( status )
    
    integer, intent (in) :: status
    
    if ( status /= nf90_noerr ) then 
      print *, trim( nf90_strerror(status) )
      stop "Stopped"
    end if

  end subroutine check     

end module md_forcing

