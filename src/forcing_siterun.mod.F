module _forcing_siterun
  !////////////////////////////////////////////////////////////////
  ! Module contains forcing variables (climate, co2, ...), and
  ! subroutines used to read forcing input files.
  ! This module is only used on the level of 'sofun', but not
  ! within 'biosphere', as these variables are passed on to 'biosphere'
  ! as arguments.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core
  use _sofunutils, only: daily2monthly, read1year_daily, read1year_monthly, &
    getvalreal
  implicit none

  real, dimension(ndayyear,maxgrid) :: dtemp_field, dprec_field, &
    dfsun_field, dvpd_field
  real, dimension(ndayyear,maxgrid) :: dndep_field
  real, dimension(nlu,maxgrid)      :: lu_area
  real                              :: pco2

contains

  subroutine getco2( runname, realyear )
    !////////////////////////////////////////////////////////////////
    !  SR reads this year's atmospheric CO2 from input
    !----------------------------------------------------------------
    use _params_core
    use _params_siml

    implicit none
    character(len=*), intent(in) :: runname
    integer, intent(in) :: realyear

    pco2 = getvalreal( trim(input_dir)//trim(co2_forcing_file), realyear )

    return

  end subroutine getco2


  subroutine getclimate_fields( runname )
    !////////////////////////////////////////////////////////////////
    !  SR reads this year's daily temperature and precipitation.
    !----------------------------------------------------------------
    use _params_core, only: nmonth, ndaymonth
    use _params_siml, only: sitename
    implicit none

    ! arguments
    character(len=*), intent(in) :: runname     

    ! local variables
    integer :: day, mo, dm, yr
    real, dimension(nmonth) :: mtemp, mfsun, mvapr, mvpd


    ! PRESCRIBED DAILY CLIMATE (TEMP, PREC, FSUN) FOR ONE YEAR
    ! xxx deal with jpngr dimension only when using NetCDF
    write(0,*) 'prescribe daily climate (temp, prec, fsun) ...'

    dtemp_field(:,:) = 0.0
    dprec_field(:,:) = 0.0
    dfsun_field(:,:) = 0.0

    dtemp_field(:,1) = read1year_daily('dtemp_'//sitename//'_2002.txt')
    dprec_field(:,1) = read1year_daily('dprec_'//sitename//'_2002.txt')

    mfsun(:) = read1year_monthly('mfsun_'//sitename//'_2002.txt')
    mvapr(:) = read1year_monthly('mvapr_'//sitename//'_2002.txt')

    ! calculate monthly VPD based on monthly vapour pressure and temperature data
    mtemp(:) = daily2monthly( dtemp_field(:,1), "mean")
    do mo=1,nmonth
      mvpd(mo) = calc_vpd( mtemp(mo), mvapr(mo) )
    end do
    write(0,*) 'mvpd', mvpd(:)

    ! use monthly value for each day in month for dfsun, and dvpd
    day=0
    do mo=1,nmonth
      do dm=1,ndaymonth(mo)
        day=day+1
        dfsun_field(day,1) = mfsun(mo)
        dvpd_field(day,1)  = mvpd(mo)
      end do
    end do

    !insf  = (/0.21, 0.27, 0.30, 0.40, 0.39, 0.39, 0.40, 0.43, 0.36, 0.32, 0.23, 0.19/)
    !intc  = (/4.80, 4.85, 7.10, 9.10, 12.4, 15.3, 17.6, 17.3, 14.6, 11.2, 7.55, 5.05/)
    !inppt = (/61.0, 41.2, 44.5, 48.0, 46.4, 44.6, 46.0, 52.3, 50.3, 71.8, 66.3, 62.9/)
    !day=0
    !do mo=1,nmonth
    !  do dm=1,ndaymonth(mo)
    !    day=day+1
    !    dtemp_field(day,1) = intc(mo)
    !    dprec_field(day,1) = inppt(mo)/ndaymonth(mo)
    !    dfsun_field(day,1) = insf(mo)
    !  end do
    !enddo

    !day=0
    !do mo=1,nmonth
    !  do dm=1,ndaymonth(mo)
    !    day=day+1
    !    print*,'day, mo, dm ',day, mo, dm
    !    dtemp_field(day,1) = getvalreal_STANDARD( &
    !      sitename//'_dtemp_2000_STANDARD.txt', mo=mo, dm=dm &
    !      )
    !    dprec_field(day,1) = getvalreal_STANDARD( &
    !      sitename//'_dprec_2000_STANDARD.txt', mo=mo, dm=dm &
    !      )
    !    dfsun_field(day,1) = getvalreal_STANDARD( &
    !      sitename//'_dfsun_2000_STANDARD.txt', mo=mo, dm=dm &
    !      )
    !  end do
    !enddo
    
    write(0,*) '... done. Good job, beni.'

    return

  end subroutine getclimate_fields


  !subroutine getclimate( day, jpngr, dtemp, dprec, dfsun )
  !  !////////////////////////////////////////////////////////////////
  !  !  SR reads this year's daily temperature and precipitation.
  !  !----------------------------------------------------------------
  !  use _params_core
  !  implicit none

!  ! ! arguments
  !  integer, intent(in) :: day     
  !  integer, intent(in) :: jpngr     
  !  real, intent(out) :: dtemp, dprec, dfsun

!  !  do yr=1,nyeartrend
  !    day=0
  !    do mo=1,nmonth
  !      do dm=1,ndaymonth(mo)
  !        day=day+1
  !        do pft=1,npft
  !          dtemp = dtemp_field(day,jpngr)
  !          dprec = dprec_field(day,jpngr)
  !          dfsun = dfsun_field(day,jpngr)
  !        end do
  !      end do
  !    enddo
  !  enddo
  !  write(0,*) '... done'

!  !  return

!  !end subroutine getclimate


  subroutine getndep( runname, realyear )
    !////////////////////////////////////////////////////////////////
    !  SR reads this year's annual ndeposition and distributes it
    !  over days according to daily precipitation.
    !----------------------------------------------------------------
    use _params_core
    implicit none

    ! arguments
    character(len=*), intent(in) :: runname
    integer, intent(in) :: realyear

    ! local variables
    real :: in_andep
    real, dimension(ndayyear) :: dprec_rel
    integer :: jpngr

    ! Get annual N-deposition           
    in_andep = 1.0  

    ! Distribute annual Ndep to days by daily precipitation
    do jpngr=1,maxgrid
      dprec_rel(:) = dprec_field(:,jpngr)/sum(dprec_field(:,jpngr))
      dndep_field(:,jpngr) = in_andep * dprec_rel(:)
    end do

    return

  end subroutine getndep


  subroutine getlanduse( runname, realyear )
    !////////////////////////////////////////////////////////////////
    !  SR reads this year's annual landuse state
    !----------------------------------------------------------------
    use _params_core, only: lunat
    implicit none

    ! arguments
    character(len=*), intent(in)      :: runname
    integer, intent(in)               :: realyear

    lu_area(lunat,:) = 1.0

    return

  end subroutine getlanduse


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
    real, intent(out) :: vpd       !  mean monthly vapor pressure deficit, Pa

    if ( present(tmin) .and. present(tmax) ) then
      my_tc = 0.5 * (tmin + tmax)
    else
      my_tc = tc
    end if

    !! calculate VPD in units of kPa
    vpd = ( 0.611 * exp( (17.27 * my_tc)/(my_tc + 237.3) ) - 0.10 * vap )    

    !! convert to Pa
    vpd = vpd * 1.0e3

  end function calc_vpd


  !===========================LOW-LEVEL============================

  ! function read1year_daily( filename )
  !   !////////////////////////////////////////////////////////////////
  !   ! Function reads a file that contains 365 lines, each line for
  !   ! a daily value. 
  !   !----------------------------------------------------------------
  !   use _params_core, only: ndayyear
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
  !   use _params_core, only: nmonth
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
  !   !  Function reads one (annual) value corresponding to the given 
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
  !   !  from a time series ascii file. File has to be located in 
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

end module _forcing_siterun

