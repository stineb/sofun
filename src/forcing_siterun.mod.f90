module md_forcing_siterun
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

  implicit none

  private
  public getco2, getninput, ninput_type, gettot_ninput, getfapar, getclimate_site, &
    getlanduse, landuse_type, climate_type

  type climate_type
    real, dimension(ndayyear) :: dtemp  ! deg C
    real, dimension(ndayyear) :: dprec  ! mm d-1
    real, dimension(ndayyear) :: dsnow  ! mm d-1 water equivalents
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
    pco2 = getvalreal( 'sitedata/co2/'//trim(domaininfo%domain_name)//'/'//trim(co2_forcing_file), readyear )

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
    real                      :: aninput_noy
    real                      :: aninput_nhx
    real, dimension(ndayyear) :: dprec_rel
    integer                   :: jpngr
    integer                   :: readyear
    real, dimension(ndayyear) :: dninput_noy
    real, dimension(ndayyear) :: dninput_nhx
    
    if (const_ninput_year/=int(dummy)) then
      readyear = const_ninput_year
    else  
      readyear = forcingyear
    end if

    ! xxx try
    if (ntype=="ndep") readyear = min( readyear, 2009 )
    readyear = max( 1850, readyear )
    ! write(0,*) 'GETNINPUT: use '//trim(ntype)//' data of year ', readyear, '...'

    ! aninput = getvalreal( trim(input_dir)//trim(ninput_forcing_file), readyear )
    aninput_noy = getvalreal( 'sitedata/'//trim(ntype)//'/'//trim(sitename)//'/'//trim(ninput_noy_forcing_file), readyear )
    aninput_nhx = getvalreal( 'sitedata/'//trim(ntype)//'/'//trim(sitename)//'/'//trim(ninput_nhx_forcing_file), readyear )

    do jpngr=1,maxgrid
      dprec_rel(:)               = climate(jpngr)%dprec(:)/sum(climate(jpngr)%dprec(:))
      out_getninput(jpngr)%dnoy(:) = aninput_noy * dprec_rel(:)
      out_getninput(jpngr)%dnhx(:) = aninput_nhx * dprec_rel(:)
      ! out_getninput(jpngr)%dnoy(:) = aninput_noy / 365.0
      ! out_getninput(jpngr)%dnhx(:) = aninput_nhx / 365.0
      out_getninput(jpngr)%dtot(:) = out_getninput(jpngr)%dnoy(:) + out_getninput(jpngr)%dnhx(:)
    end do

    ! print*,'out_getninput(jpngr) ', sum(out_getninput(1)%dnoy(:)), sum(out_getninput(1)%dnhx(:)), sum(out_getninput(1)%dtot(:))
    ! write(0,*) 'GETNINPUT: done'

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


  function getfapar( domaininfo, grid, year, fapar_forcing_source ) result( out_vegcover )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(in) :: grid
    integer, intent(in) :: year
    character(len=*), intent(in) :: fapar_forcing_source

    ! function return variable
    type( vegcover_type ), dimension(domaininfo%maxgrid) :: out_vegcover

    ! local variables 
    integer :: readyear
    character(len=4) :: faparyear_char

    if (trim(fapar_forcing_source)=='NA') then
      ! If in simulation parameter file 'NA' is specified for 'fapar_forcing_source'
      out_vegcover(1)%dfapar(:) = dummy

    else
      ! Prescribed. Read monthly fAPAR value from file
      ! create 4-digit string for year  
      write(faparyear_char,999) min( max( 2000, year ), 2014 )
      out_vegcover(1)%dfapar(:) = read1year_daily( 'sitedata/fapar/'//trim(domaininfo%domain_name)//'/'//faparyear_char//'/'//'dfapar_'//trim(domaininfo%domain_name)//'_'//faparyear_char//'.txt' )

      ! ! "Correct" fAPAR
      ! print*,"WARNING: normalising fAPAR to within 0.12 and 1.0."
      ! out_vegcover(1)%dfapar(:) = max((out_vegcover(1)%dfapar(:) - 0.12), 0.0)/(1.0 - 0.12)

    end if

    return
    999  format (I4.4)

  end function getfapar


  function getclimate_site( sitename, climateyear ) result ( out_climate )
    !////////////////////////////////////////////////////////////////
    ! SR reads this year's daily temperature and precipitation.
    ! Read year-2013 data after 2013
    !----------------------------------------------------------------    
    ! arguments
    character(len=*), intent(in) :: sitename
    integer, intent(in) :: climateyear
    logical, intent(in) :: in_ppfd
    logical, intent(in) :: in_netrad

    ! local variables
    character(len=4) :: climateyear_char
    character(len=256) :: snowf_name
    logical :: snowf_exists

    ! function return variable
    type( climate_type ), dimension(domaininfo%maxgrid) :: out_climate

    ! create 4-digit string for year  
    write(climateyear_char,999) climateyear

    out_climate(1)%dtemp(:) = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dtemp_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    out_climate(1)%dprec(:) = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dprec_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    out_climate(1)%dvpd(:)  = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dvpd_' //trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    if (in_ppfd) then
      out_climate(1)%dppfd(:) = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dppfd_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    else
      out_climate(1)%dppfd(:) = dummy
    end if
    if (in_netrad) then
      out_climate(1)%dnetrad(:) = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dnetrad_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    else
      out_climate(1)%dnetrad(:) = dummy
    end if
    if ( in_netrad .and. in_ppfd ) then
      out_climate(1)%dfsun(:) = dummy
    else
      out_climate(1)%dfsun(:) = read1year_daily('sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dfsun_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt')
    end if

    ! read snow file if it exists
    snowf_name = 'sitedata/climate/'//trim(domaininfo%domain_name)//'/'//climateyear_char//'/'//'dsnow_'//trim(domaininfo%domain_name)//'_'//climateyear_char//'.txt'
    INQUIRE(FILE = trim(snowf_name), EXIST = snowf_exists)
    if (snowf_exists) then
      out_climate(1)%dsnow(:) = read1year_daily(snowf_name)
    else
      out_climate(1)%dsnow(:) = 0.0
    end if

    return
    999  format (I4.4)

  end function getclimate


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
    out_landuse%lu_area(lunat) = 1.0

    if (const_lu_year/=int(dummy)) then
      readyear = const_lu_year
    else
      readyear = forcingyear
    end if    

    ! get harvest data for forcing year
    if (present(do_grharvest_forcing_file)) then

      ! create 4-digit string for year  
      write(landuseyear_char,999) readyear

      filnam = 'sitedata/landuse/'//trim(sitename)//'/'//landuseyear_char//'/'//trim(do_grharvest_forcing_file)//'_'//trim(sitename)//'_'//landuseyear_char//'.txt'
      inquire( file='./input/'//trim(filnam), exist=file_exists )

      if ( file_exists ) then
        ! found data file
        ! write(0,*) 'GETLANDUSE: use harvest data for year ', readyear
        tmp(:) = read1year_daily( trim(filnam) )
      else
        
        ! find first year with data available
        findyear = readyear
        do while ( .not. file_exists .and. findyear<2501 )
          findyear = findyear - 1
          write(landuseyear_char,999) findyear
          filnam = 'sitedata/landuse/'//trim(sitename)//'/'//landuseyear_char//'/'//trim(do_grharvest_forcing_file)//'_'//trim(sitename)//'_'//landuseyear_char//'.txt'
          inquire( file='./input/'//trim(filnam), exist=file_exists )
          if ( file_exists ) tmp(:) = read1year_daily( trim(filnam) )
        end do   
        write(0,*) 'GETLANDUSE: found harvest data for year  ', landuseyear_char

      end if

      ! translate zeros and ones to boolean
      do doy=1,ndayyear
        if (tmp(doy)==1.0) then
          out_landuse%do_grharvest(doy) = .true.
        else
          out_landuse%do_grharvest(doy) = .false.
        end if
      end do
    end if

    return
    999  format (I4.4)

  end function getlanduse


  function calc_vpd( tc, vap, tmin, tmax ) result( vpd )
    !-----------------------------------------------------------------------
    ! Output:   mean monthly vapor pressure deficit, Pa (vpd)
    ! Features: Returns mean monthly vapor pressure deficit
    ! Ref:      Eq. 5.1, Abtew and Meleese (2013), Ch. 5 Vapor Pressure 
    !           Calculation Methods, in Evaporation and Evapotranspiration: 
    !           Measurements and Estimations, Springer, London.
    !             vpd = 0.611*exp[ (17.27 tc)/(tc + 237.3) ] - ea
    !             where:
    !                 tc = average daily air temperature, deg C
    !                 vap  = actual vapor pressure, kPa
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc            ! mean monthly temperature, deg C
    real, intent(in) :: vap             ! mean monthly vapor pressure, hPa (because CRU data is in hPa)
    real, intent(in), optional :: tmin  ! (optional) mean monthly min daily air temp, deg C 
    real, intent(in), optional :: tmax  ! (optional) mean monthly max daily air temp, deg C 

    ! local variables
    real :: my_tc

    ! function return variable
    real :: vpd       !  mean monthly vapor pressure deficit, Pa

    if ( present(tmin) .and. present(tmax) ) then
      my_tc = 0.5 * (tmin + tmax)
    else
      my_tc = tc
    end if

    ! calculate VPD in units of kPa
    vpd = ( 0.611 * exp( (17.27 * my_tc)/(my_tc + 237.3) ) - 0.10 * vap )    

    ! this empirical equation may lead to negative values for VPD (happens very rarely). assume positive...
    vpd = max( 0.0, vpd )

    !! convert to Pa
    vpd = vpd * 1.0e3

  end function calc_vpd


end module md_forcing_siterun


