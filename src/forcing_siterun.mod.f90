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
    real, dimension(ndayyear) :: dtemp
    real, dimension(ndayyear) :: dprec
    real, dimension(ndayyear) :: dfsun
    real, dimension(ndayyear) :: dvpd
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
    ! Function reads this year's annual ninputosition and distributes it
    ! over days according to daily precipitation.
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


  function getfapar( runname, sitename, forcingyear, fapar_forcing_source ) result( fapar_field )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
    integer, intent(in) :: forcingyear
    character(len=*), intent(in) :: fapar_forcing_source

    ! function return variable
    real, dimension(nmonth,maxgrid) :: fapar_field

    ! local variables 
    integer :: jpngr
    integer :: readyear
    character(len=4) :: faparyear_char

    if (trim(fapar_forcing_source)=='NA') then
      ! If in simulation parameter file 'NA' is specified for 'fapar_forcing_source', then set fapar_field to dummy value
      do jpngr=1,maxgrid
        fapar_field(:,jpngr) = dummy
      end do

    else
      ! Prescribed. Read monthly fAPAR value from file
      do jpngr=1,maxgrid
        ! create 4-digit string for year  
        write(faparyear_char,999) min( max( 2000, forcingyear ), 2014 )
        fapar_field(:,jpngr) = read1year_monthly( 'sitedata/fapar/'//trim(sitename)//'/'//faparyear_char//'/'//'fapar_'//trim(fapar_forcing_source)//'_'//trim(sitename)//'_'//faparyear_char//'.txt' )
      end do
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

    ! local variables
    integer :: day
    integer :: jpngr = 1
    real, dimension(ndayyear) :: dvapr
    character(len=4) :: climateyear_char

    ! function return variable
    type( climate_type ), dimension(maxgrid) :: out_climate

    if (climateyear>2013) then
      ! write(0,*) 'GETCLIMATE_SITE: held climate fixed at year 2013'
      write(climateyear_char,999) 2013
    else
      ! create 4-digit string for year  
      write(climateyear_char,999) climateyear
    end if
    ! filnam_dtemp = 'sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dtemp_'//trim(sitename)//'_'//climateyear_char//'.txt'
    ! filnam_dprec = 'sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dprec_'//trim(sitename)//'_'//climateyear_char//'.txt'
    ! filnam_dfsun = 'sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dfsun_'//trim(sitename)//'_'//climateyear_char//'.txt'
    ! filnam_dvapr = 'sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dvapr_'//trim(sitename)//'_'//climateyear_char//'.txt'
    
    ! write(0,*) 'GETCLIMATE_SITE: use climate data of year ', climateyear_char

    jpngr = 1

    out_climate(jpngr)%dtemp(:) = read1year_daily('sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dtemp_'//trim(sitename)//'_'//climateyear_char//'.txt')
    out_climate(jpngr)%dprec(:) = read1year_daily('sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dprec_'//trim(sitename)//'_'//climateyear_char//'.txt')
    out_climate(jpngr)%dfsun(:) = read1year_daily('sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dfsun_'//trim(sitename)//'_'//climateyear_char//'.txt')

    dvapr(:) = read1year_daily('sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dvapr_'//trim(sitename)//'_'//climateyear_char//'.txt')

    ! calculate daily VPD based on daily vapour pressure and temperature data
    do day=1,ndayyear
      out_climate(jpngr)%dvpd(day) = calc_vpd( out_climate(jpngr)%dtemp(day), dvapr(day) )
    end do

    return
    
    999  format (I4.4)

  end function getclimate_site


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

