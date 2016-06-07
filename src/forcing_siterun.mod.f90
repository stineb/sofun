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
  use md_params_core, only: nmonth, ndaymonth, lunat, ndayyear, maxgrid, nlu
  use md_sofunutils, only: daily2monthly, read1year_daily, read1year_monthly, &
    getvalreal, monthly2daily_weather, monthly2daily

  implicit none

  private
  public getco2, getndep, ndep_type, getfapar, getclimate_site, getlanduse, landuse_type, climate_type

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

  type ndep_type
    real, dimension(ndayyear) :: dnoy
    real, dimension(ndayyear) :: dnhx
    real, dimension(ndayyear) :: dtot
  end type ndep_type

contains

  function getco2( runname, sitename, forcingyear, const_co2, firstyeartrend, co2_forcing_file ) result( pco2 )
    !////////////////////////////////////////////////////////////////
    !  Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
    integer, intent(in) :: forcingyear
    logical, intent(in) :: const_co2
    integer, intent(in) :: firstyeartrend
    character(len=*), intent(in) :: co2_forcing_file

    ! function return variable
    real, intent(out) :: pco2

    ! local variables 
    integer :: readyear

    if (const_co2) then
      readyear = firstyeartrend
    else  
      readyear = forcingyear
    end if
    write(0,*) 'GETCO2: use CO2 data of year ', readyear
    pco2 = getvalreal( 'sitedata/co2/'//trim(sitename)//'/'//trim(co2_forcing_file), readyear )

  end function getco2


  function getndep( runname, sitename, forcingyear, firstyeartrend, const_ndep, ndep_noy_forcing_file, ndep_nhx_forcing_file, climate ) result( out_getndep )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's annual ndeposition and distributes it
    !  over days according to daily precipitation.
    !----------------------------------------------------------------
    use md_params_core, only: dummy

    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
    integer, intent(in)          :: forcingyear
    integer, intent(in) :: firstyeartrend
    logical, intent(in) :: const_ndep
    character(len=*), intent(in) :: ndep_noy_forcing_file
    character(len=*), intent(in) :: ndep_nhx_forcing_file
    type( climate_type ), dimension(maxgrid), intent(in) :: climate

    ! function return variable
    type( ndep_type ), dimension(maxgrid) :: out_getndep 

    ! local variables
    real                      :: andep_noy
    real                      :: andep_nhx
    real, dimension(ndayyear) :: dprec_rel
    integer                   :: jpngr
    integer                   :: readyear
    real, dimension(ndayyear) :: dndep_noy
    real, dimension(ndayyear) :: dndep_nhx
    
    if (const_ndep) then
      readyear = firstyeartrend
    else  
      readyear = forcingyear
    end if

    ! xxx try
    ! readyear = max( 1992, min( readyear, 2009 ) )
    readyear = max( 1850, readyear )
    write(0,*) 'GETNDEP: use N fertilisation data of year ', readyear

    ! andep = getvalreal( trim(input_dir)//trim(ndep_forcing_file), readyear )
    andep_noy = getvalreal( 'sitedata/ndep/'//trim(sitename)//'/'//trim(ndep_noy_forcing_file), readyear )
    andep_nhx = getvalreal( 'sitedata/ndep/'//trim(sitename)//'/'//trim(ndep_nhx_forcing_file), readyear )

    ! Distribute annual Ndep to days by daily precipitation

    do jpngr=1,maxgrid
      dprec_rel(:)               = climate(jpngr)%dprec(:)/sum(climate(jpngr)%dprec(:))
      out_getndep(jpngr)%dnoy(:) = andep_noy * dprec_rel(:)
      out_getndep(jpngr)%dnhx(:) = andep_nhx * dprec_rel(:)
      ! out_getndep(jpngr)%dnoy(:) = andep_noy / 365.0
      ! out_getndep(jpngr)%dnhx(:) = andep_nhx / 365.0
      out_getndep(jpngr)%dtot(:) = out_getndep(jpngr)%dnoy(:) + out_getndep(jpngr)%dnhx(:)
    end do

    ! print*,'out_getndep(jpngr) ', sum(out_getndep(1)%dnoy(:)), sum(out_getndep(1)%dnhx(:)), sum(out_getndep(1)%dtot(:))

  end function getndep


  function getfapar( runname, sitename, forcingyear ) result( fapar_field )
    !////////////////////////////////////////////////////////////////
    ! Function reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    use md_params_core, only: dummy

    ! arguments
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: sitename
    integer, intent(in) :: forcingyear

    ! function return variable
    real, dimension(nmonth,maxgrid) :: fapar_field

    ! local variables 
    integer :: jpngr
    integer :: readyear
    character(len=4) :: faparyear_char

    do jpngr=1,maxgrid
      ! create 4-digit string for year  
      write(faparyear_char,999) min( max( 2000, forcingyear ), 2014 )
      fapar_field(:,jpngr) = read1year_monthly( 'sitedata/fapar/'//trim(sitename)//'/'//faparyear_char//'/'//'fapar_modis_'//trim(sitename)//'_'//faparyear_char//'.txt' )
    end do

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
      write(0,*) 'GETCLIMATE_SITE: held climate fixed at year 2013'
      write(climateyear_char,999) 2013
    ! else if (climateyear<1993) then
    !   write(0,*) 'GETCLIMATE_SITE: held climate fixed at year 1993'
    !   write(climateyear_char,999) 1993
    else
      ! create 4-digit string for year  
      write(climateyear_char,999) climateyear
    end if
    ! filnam_dtemp = 'sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dtemp_'//trim(sitename)//'_'//climateyear_char//'.txt'
    ! filnam_dprec = 'sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dprec_'//trim(sitename)//'_'//climateyear_char//'.txt'
    ! filnam_dfsun = 'sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dfsun_'//trim(sitename)//'_'//climateyear_char//'.txt'
    ! filnam_dvapr = 'sitedata/climate/'//trim(sitename)//'/'//climateyear_char//'/'//'dvapr_'//trim(sitename)//'_'//climateyear_char//'.txt'
    
    write(0,*) 'GETCLIMATE_SITE: use climate data of year ', climateyear_char

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


  function getlanduse( runname, sitename, forcingyear, do_grharvest_forcing_file, const_lu, firstyeartrend ) result( out_landuse )
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
    logical, intent(in) :: const_lu
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

    if (const_lu) then
      readyear = firstyeartrend
    else
      readyear = forcingyear
    end if    

    ! get harvest data for forcing year
    if (present(do_grharvest_forcing_file)) then
      if (readyear>2002) then
        write(0,*) 'GETLANDUSE: held harvest dates fixed after 2002'
        write(landuseyear_char,999) 2002
      else
        ! create 4-digit string for year  
        write(landuseyear_char,999) readyear
      end if
      filnam = 'sitedata/landuse/'//trim(sitename)//'/'//landuseyear_char//'/'//trim(do_grharvest_forcing_file)//'_'//trim(sitename)//'_'//landuseyear_char//'.txt'
      inquire( file='./input/'//trim(filnam), exist=file_exists )
      
      if ( file_exists ) then
        ! found data file
        write(0,*) 'GETLANDUSE: use harvest data for year ', readyear
        tmp(:) = read1year_daily( trim(filnam) )
      else
        ! find first year with data available
        findyear = readyear
        do while ( .not. file_exists )
          findyear = findyear + 1
          write(landuseyear_char,999) findyear
          filnam = 'sitedata/landuse/'//trim(sitename)//'/'//landuseyear_char//'/'//trim(do_grharvest_forcing_file)//'_'//trim(sitename)//'_'//landuseyear_char//'.txt'
          inquire( file='./input/'//trim(filnam), exist=file_exists )
          if ( file_exists ) tmp(:) = read1year_daily( trim(filnam) )
        end do   
        write(0,*) 'GETLANDUSE: found harvest data for first year  ', findyear
      end if

      ! write(0,*) 'GETLANDUSE: forced no harvest  ', findyear

      ! translate zeros and ones to boolean
      do doy=1,ndayyear
        if (tmp(doy)==1.0) then
          out_landuse%do_grharvest(doy) = .true.
          ! out_landuse%do_grharvest(doy) = .false.
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

    !! calculate VPD in units of kPa
    vpd = ( 0.611 * exp( (17.27 * my_tc)/(my_tc + 237.3) ) - 0.10 * vap )    

    ! if (vpd<0.0) then
    !   print*,'temp: ', my_tc
    !   print*,'vapr: ', vap
    !   print*,'vpd : ', vpd
    !   print*,'SETTING VPD TO ZERO'
    !   vpd = 0.0
    !   stop
    ! end if

    !! convert to Pa
    vpd = vpd * 1.0e3

  end function calc_vpd


  !===========================LOW-LEVEL============================

  ! function read1year_daily( filename )
  !   !////////////////////////////////////////////////////////////////
  !   ! Function reads a file that contains 365 lines, each line for
  !   ! a daily value. 
  !   !----------------------------------------------------------------
  !   use md_params_core, only: ndayyear
  !   implicit none

  !   ! arguments
  !   character(len=*), intent(in) :: filename

  !   ! local variables
  !   real, dimension(ndayyear) :: dval

  !   ! function return value
  !   real, dimension(ndayyear) :: read1year_daily

  !   open(20, file='./input/'//filename, status='old',  form='formatted', action='read', err=888)
  !   read(20,*) dval
  !   close(20)

  !   read1year_daily = dval

  !   return
  !   600 format (F9.7)
  !   888 write(0,*) 'READ1YEAR: error opening file '//trim(filename)//'. Abort. '
  !   stop

  ! end function read1year_daily


  ! function read1year_monthly( filename )
  !   !////////////////////////////////////////////////////////////////
  !   ! Function reads a file that contains 365 lines, each line for
  !   ! a daily value. 
  !   !----------------------------------------------------------------
  !   use md_params_core, only: nmonth
  !   implicit none

  !   ! arguments
  !   character(len=*), intent(in) :: filename

  !   ! local variables
  !   real, dimension(nmonth) :: mval

  !   ! function return value
  !   real, dimension(nmonth) :: read1year_monthly

  !   open(20, file='./input/'//trim(filename), status='old',  form='formatted', action='read', err=888)
  !   read(20,*) mval
  !   close(20)

  !   read1year_monthly = mval

  !   return
  !   600 format (F9.7)
  !   888 write(0,*) 'READ1YEAR: error opening file '//trim(filename)//'. Abort. '
  !   stop

  ! end function read1year_monthly


  ! function getvalreal( filename, realyear, day, dm, mo )
  !   !////////////////////////////////////////////////////////////////
  !   ! Function reads one (annual) value corresponding to the given 
  !   !  year from a time series ascii file. 
  !   !----------------------------------------------------------------

  !   implicit none
  !   ! arguments
  !   character(len=*), intent(in) :: filename
  !   integer, intent(in) :: realyear
  !   integer, intent(in), optional :: day ! day in year (1:365)
  !   integer, intent(in), optional :: dm  ! day in month (1:31)
  !   integer, intent(in), optional :: mo  ! month in year (1:12)

  !   ! function return value
  !   real :: getvalreal

  !   ! local variables
  !   integer :: l
  !   real :: tmp(3) ! 3 so that an additional value for this year could be read
  !   real :: realyear_decimal 

  !   if (present(day)) then
  !    ! convert day number into decimal number
  !    realyear_decimal = real(realyear) + real(day)/real(ndayyear)
  !   endif

  !   open(20, file=filename, status='old',  form='formatted', err=888)

  !   if (present(day)) then
  !    ! find corresponding day in first column and read 3 values on this line
  !    read(20, 100, err=999) (tmp(l), l=1,3)  
  !    do while (abs(realyear_decimal-tmp(1)).gt.1.0d-8)
  !      read(20, 100, err=999) (tmp(l), l=1,3)
  !    enddo

  !   else
  !    ! find corresponding year in first column and read 3 values on this line
  !    read(20, 100, err=999) (tmp(l), l=1,3)  
  !    do while (abs(realyear-tmp(1)).gt.1.0d-8)
  !      read(20, 100, err=999) (tmp(l), l=1,3)
  !    enddo

  !   endif

  !   getvalreal = tmp(2) 

  !   100     format (30d16.8)
  !   close(20)

  !   return

  !   888     write(0,*) 'GETVALREAL: error opening file '//trim(filename)//'. Abort. '
  !   stop
  !   999     write(0,*) 'GETVALREAL: error reading file '//trim(filename)//'. Abort. '
  !   stop 

  ! end function getvalreal


  ! function getvalreal_STANDARD( filename, realyear, mo, dm, day, realyear_decimal )
  !   !////////////////////////////////////////////////////////////////
  !   !  SR reads one (annual) value corresponding to the given year 
  !   ! from a time series ascii file. File has to be located in 
  !   !  ./input/ and has to contain only rows formatted like
  !   !  '2002  1  1 0.496632 0.054053', which represents 
  !   !  'YYYY MM DM      GPP GPP err.'. DM is the day within the month.
  !   !  If 'realyear' is dummy (-9999), then it's interpreted as to 
  !   !  represent a mean climatology for the course of a year.
  !   !----------------------------------------------------------------

  !   implicit none
  !   ! arguments
  !   character(len=*), intent(in) :: filename
  !   integer, intent(in), optional :: realyear ! year AD as integer
  !   integer, intent(in), optional :: mo  ! month in year (1:12)
  !   integer, intent(in), optional :: dm  ! day in month (1:31 or 1:31 or 1:28)
  !   integer, intent(in), optional :: day ! day in year (1:365)
  !   real,    intent(in), optional :: realyear_decimal ! year AD as decimal number corresponding to day in the year

  !   ! function return value
  !   real :: getvalreal_STANDARD

  !   ! local variables
  !   integer :: inyear
  !   integer :: inmo
  !   integer :: indm
  !   integer :: inday
  !   real    :: inyear_decimal
  !   real    :: inval1
  !   real    :: inval2

  !   !print*,'looking for realyear, mo, dm',realyear,mo,dm

  !   ! open file
  !   open(20, file='./input/'//filename, status='old', form='formatted', err=888)

  !   if (present(realyear)) then
  !      ! DATA FOR EACH YEAR
  !      if (present(mo)) then
  !          ! DATA FOR EACH MONTH
  !          if (present(dm)) then
  !              ! DATA FOR EACH DAY IN THE MONTH
  !              ! read the 2 values for this day in this year
  !              read(20, 100, err=999) inyear, inmo, indm, inval1, inval2
  !              do while ( (realyear-inyear).ne.0 .or. (mo-inmo).ne.0 .or. (dm-indm).ne.0 )
  !                read(20, 100, err=999) inyear, inmo, indm, inval1, inval2
  !              enddo
  !          else           
  !              ! read the 2 values for this month in this year
  !              read(20, 200, err=999) inyear, inmo, inval1, inval2
  !              do while ( (realyear-inyear).ne.0 .or. (mo-inmo).ne.0 )
  !                read(20, 200, err=999) inyear, inmo, inval1, inval2
  !              enddo
  !          end if
  !      else if (present(day)) then
  !          ! DATA FOR EACH DAY IN THE YEAR
  !          ! read the 2 values for this day in this year
  !          read(20, 700, err=999) inyear, inday, inval1, inval2
  !          do while ( (realyear-inyear).ne.0 .or. (day-inday).ne.0 )
  !            read(20, 700, err=999) inyear, inday, inval1, inval2
  !          enddo
  !      else
  !          ! read the 2 values for this year
  !          read(20, 300, err=999) inyear, inval1, inval2
  !          do while ( (realyear-inyear).ne.0 )
  !            read(20, 300, err=999) inyear, inval1, inval2
  !          enddo
  !      end if
  !   else if (present(realyear_decimal)) then
  !     ! DATA PROVIDED FOR EACH DAY AS A DECIMAL OF REALYEAR
  !     ! find corresponding day in first column and read 3 values on this line
  !     read(20, 900, err=999) inyear_decimal, inval1, inval2  
  !     do while (abs(realyear_decimal-inyear_decimal).gt.1.0d-8)
  !       read(20, 900, err=999) inyear_decimal, inval1, inval2  
  !     enddo
  !   else
  !      ! DATA AS AVERAGE OVER MULTIPLE YEARS (recycle climatology)
  !      ! FOR EACH MONTH, AND DAY-IN-THE-MONTH
  !      if (present(mo)) then
  !          if (present(dm)) then
  !              ! read the 2 values for this day
  !              read(20, 400, err=999) inmo, indm, inval1, inval2
  !              !print*,'inmo, indm, inval1, inval2', inmo, indm, inval1, inval2
  !              do while ( (mo-inmo).ne.0 .or. (dm-indm).ne.0 )
  !                read(20, 400, err=999) inmo, indm, inval1, inval2
  !                !print*,'inmo, indm, inval1, inval2', inmo, indm, inval1, inval2
  !              enddo
  !          else           
  !              ! read the 2 values for this month
  !              read(20, 500, err=999) inmo, inval1, inval2
  !              do while ( (mo-inmo).ne.0 )
  !                read(20, 500, err=999) inmo, inval1, inval2
  !              enddo

  !          end if
  !      else if (present(day)) then
  !          ! DATA FOR EACH DAY IN THE YEAR
  !          ! read the 2 values for this day
  !          read(20, 800, err=999) inday, inval1, inval2
  !          do while ( (day-inday).ne.0 )
  !            read(20, 800, err=999) inday, inval1, inval2
  !          enddo
  !      else
  !          ! read the 2 values in this input file
  !          read(20, 600, err=999) inval1, inval2
  !      end if
  !   endif

  !   !print*,'found realyear, mo, dm      ',inyear,inmo,indm,inval1

  !   getvalreal_STANDARD = inval1

  !   100     format (I4,I3,I3,F9.7,F9.7)
  !   200     format (I4,I3,F9.7,F9.7)
  !   300     format (I4,F9.7,F9.7)
  !   400     format (I3,I3,F9.7,F9.7)
  !   500     format (I3,F9.7,F9.7)
  !   600     format (F9.7,F9.7)
  !   700     format (I4,I4,F9.7,F9.7)
  !   800     format (I4,F9.7,F9.7)
  !   900     format (30d16.8,F9.7,F9.7)

  !   close(20)

  !   return

  !   888     write(0,*) 'GETVALREAL_STANDARD: error opening file '//trim(filename)//'. Abort. '
  !   stop
  !   999     write(0,*) 'GETVALREAL_STANDARD: error reading file '//trim(filename)//'. Abort. '
  !   stop 

  ! end function getvalreal_STANDARD

end module md_forcing_siterun

