module _sofunutils
  !/////////////////////////////////////////////////////////////////////////
  ! Contains utility functions to deal with arrays containing 365 daily 
  ! values or 12 monthly values representing one year.
  ! - running: calculates running mean
  ! - daily2monthly: calculates monthly mean/sum of daily values for resp. mo.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !-------------------------------------------------------------------------
  implicit none

contains

  function running( presval, inow, lenval, lenper, method, prevval ) result( runningval )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns running sum or average over. 'prevval' is optional, if not pro-
    ! vided, sum/average is taken only over preceeding days/months of current
    ! year.
    !-------------------------------------------------------------------------
    ! arguments
    ! xxx instead of dimension declaration with 'lenval', use size(presval)
    real, dimension(lenval), intent(in) :: presval            ! vector containing 'lenvals' values for each timestep in this year
    integer, intent(in) :: inow                               ! index corresponding to "now" (day of year or month of year)  
    integer, intent(in) :: lenval                             ! number of timesteps per year
    integer, intent(in) :: lenper                             ! number of timesteps over which to average/sum
    character(len=*), intent(in) :: method                    ! either "sum" or "mean" for running sum or running mean
    real, dimension(lenval), intent(in), optional :: prevval  ! vector containing 'lenvals' values for each timestep in previous year

    ! local variables
    real, dimension(lenval) :: valbuf

    ! function return variable
    real, intent(out) :: runningval

    !print*,'day, lenper ',inow, lenper

    if (present(prevval)) then
      !print*,'A filling up valbuf from ',(lenval-(inow-1)),'to',lenval
      !print*,'A with values from        1 to     ',inow
      valbuf((lenval-(inow-1)):lenval) = presval(1:inow)
      !print*,'B filling up valbuf from  1 to',(lenval-inow)
      !print*,'B with values from       ',(inow+1),' to ',lenval
      valbuf(1:(lenval-inow)) = prevval((inow+1):lenval)
    else
      !print*,'A filling up valbuf from  1 to',inow
      !print*,'A with values from        1 to ',inow
      valbuf(1:inow) = presval(1:inow)
      !print*,'B filling up valbuf from  ',(inow+1),'to',lenval
      !print*,'B with values zero'
      valbuf((inow+1):lenval) = 0.0
    end if

    select case (method)
      case ("sum")
        runningval = sum(valbuf((lenval-lenper+1):lenval))
      case ("mean")
        if (present(prevval)) then
          runningval = sum(valbuf((lenval-lenper+1):lenval))/lenper
        else
          runningval = sum(valbuf((lenval-lenper+1):lenval))/inow
        end if
      case default
        stop 'RUNNING: declare valid method.'
    end select

  end function running


  function daily2monthly( dval, method ) result( mval )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns monthly values as a mean over daily values in each month.
    ! Arguments:
    ! - dval   : vector containing 365 (366 in case lapyear is TRUE) daily values
    ! - method : true of monthly values represent total of daily values in resp. month
    !-------------------------------------------------------------------------
    use _params_core, only: ndaymonth, cumdaymonth, ndayyear, nmonth

    ! arguments
    real, intent(in), dimension(ndayyear) :: dval
    character(len=*), intent(in) :: method

    ! function return variable
    real, dimension(nmonth), intent(out) :: mval

    ! local variables
    integer :: moy
    integer, dimension(nmonth) :: istart, iend

    istart = cumdaymonth - ndaymonth + 1
    iend   = cumdaymonth

    ! loop over months and take sum/mean of days in that month
    do moy=1,nmonth
      select case (method)
        case ("sum")
          mval(moy) = sum( dval( istart(moy) : iend(moy) ))
        case ("mean")
          mval(moy) = sum( dval( istart(moy) : iend(moy) )) / ndaymonth(moy)
        case default
          stop 'DAILY2MONTHLY: select valid method (sum, mean)' 
      end select
    end do

    return

  end function daily2monthly


  function monthly2daily_weather( mval_prec, mval_wet, prdaily_random ) result( dval_prec )
    !--------------------------------------------------------------------
    ! Distributes monthly total precipitation to days, given number of 
    ! monthly wet days. Adopted from LPX.
    !--------------------------------------------------------------------
    use _params_core, only: nmonth, ndayyear, ndaymonth

    ! arguments
    real, dimension(nmonth), intent(in)     :: mval_prec  ! monthly precipitation totals
    real, dimension(nmonth), intent(in)     :: mval_wet   ! monthly number of wet days
    real, dimension(ndayyear,2), intent(in) :: prdaily_random ! random number seed

    ! local vars
    integer :: doy, moy, dm, iloop
    integer :: daysum            !accumulated days at beginning of month
    integer :: nwet              !# of generated wet days
    
    logical :: dry

    real :: prob, vv, v1, rand
    real, parameter :: c1 = 1.0
    real, parameter :: c2 = 1.2
    real, dimension(nmonth) :: prob_rain
    real, dimension(nmonth) :: mprecave      !average precipitation on wet days
    real, dimension(nmonth) :: mprecip       !acc. monthly generated precipitation

    ! function return variable
    real, dimension(ndayyear), intent(out) :: dval_prec

    doy = 0
    prob = 0.0
    do moy=1,nmonth
      prob_rain(moy) = 0.0
      mprecave(moy) = 0.0
      mprecip(moy) = 0.0
    enddo
    daysum = 0

    ! write(0,*) 'A'

    ! Initialize 2nd random number generator
    ! call srand( int( prdaily_random(1,1) * 100000 ) )

    ! xxx this solves the problem (removing multiplication with 1e5. But number of 
    ! actual wet days is not perfectly consistent with prescribed number of wet days.
    call srand( int( prdaily_random(1,1) ) )
   
    do moy=1,nmonth
      if ( mval_wet(moy)<=1.0 ) mval_wet(moy) = 1.0
      prob_rain(moy) = mval_wet(moy) / real( ndaymonth(moy) )
      mprecave(moy) = mval_prec(moy) / mval_wet(moy)
      dry = .TRUE.
      iloop = 0

      ! write(0,*) 'B'

      do while( dry )
        ! write(0,*) 'aa'
        iloop = iloop + 1
        nwet = 0
        do dm=1,ndaymonth(moy)
          ! write(0,*) 'a'
          doy = doy + 1
   
          ! Transitional probabilities (Geng et al. 1986)
          if (doy>1) then
            if (dval_prec(doy-1) < 0.1) then
              prob = 0.75 * prob_rain(moy)
            else 
              prob = 0.25 + (0.75 * prob_rain(moy))
            endif
          endif
          ! write(0,*) 'b'
        
          ! Determine we randomly and use Krysanova / Cramer estimates of 
          ! parameter values (c1,c2) for an exponential distribution
          if (iloop==1) then 
            ! write(0,*) 'getting prdaily_random'
            vv = real(prdaily_random(doy,1))
            ! write(0,*) 'vv', vv
          else
            ! write(0,*) 'calling rand'
            ! xxx problem: rand() generates a random number that leads to floating point exception
            vv = rand()
            ! write(0,*) 'vv'
            ! write(0,*) vv
          endif
          ! write(0,*) 'c'
          ! write(0,*) 'prob', prob
          if (vv>prob) then
            ! write(0,*) 'd1'
            dval_prec(doy) = 0.0
          else
            ! write(0,*) 'c1'
            nwet = nwet + 1
            ! write(0,*) 'c2'
            v1 = real( prdaily_random(doy,2) )
            ! write(0,*) 'c3'
            dval_prec(doy) = ((-log(v1)) ** c2) * mprecave(moy) * c1 
            ! write(0,*) 'c4'
            if (dval_prec(doy) < 0.1) dval_prec(doy) = 0.0
            ! write(0,*) 'c5'
          endif
          ! write(0,*) 'd'
          mprecip(moy) = mprecip(moy) + dval_prec(doy)
        enddo
    
        ! If it never rained this month and mprec(moy)>0 and mval_wet(moy)>0, do
        ! again
        dry = (nwet==0 .and. iloop<50 .and. mval_prec(moy)>0.1)
        if (iloop>50) then
          write(0,*) 'Daily.F, prdaily: Warning stopped after 50 tries in cell'
        endif

        ! Reset counter to start of month          
        if (dry) then
          doy = doy-ndaymonth(moy)
        endif

      enddo !while
      
      ! write(0,*) 'C'


      ! normalise generated precipitation by monthly CRU values
      if ( moy > 1 ) daysum = daysum + ndaymonth(moy-1)
      if ( mprecip(moy) < 1.0 ) mprecip(moy) = 1.0
      do dm=1,ndaymonth(moy)
        doy = daysum + dm
        dval_prec(doy) = dval_prec(doy) * (mval_prec(moy) / mprecip(moy))
        if ( dval_prec(doy) < 0.1 ) dval_prec(doy) = 0.0
        ! dval_prec(doy) = mval_prec(moy) / ndaymonth(moy)  !no generator
      enddo
         
      ! Alternative: equal distribution of rain for fixed number of wet days
      ! prob = prob_rain(moy) + prob
      ! if (prob.ge.1.0) then   
      !   dval_prec(doy) = mprec(moy)
      !   prob = prob-1.0
      ! else
      !   dval_prec(doy) = 0.0
      !   prob = prob
      ! endif
                      
    enddo                     !month 
    

  end function monthly2daily_weather


  function monthly2daily( mval, method, monthistotal, mval_pvy, mval_nxy ) result( dval )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns daily values based on monthly values, using a defined method.
    !-------------------------------------------------------------------------
    use _params_core, only: middaymonth, ndayyear, ndaymonth, nmonth
    
    ! arguments
    real, dimension(nmonth), intent(in) :: mval  ! vector containing 12 monthly values
    character(len=*), intent(in) :: method
    logical, intent(in), optional :: monthistotal ! true of monthly values represent total of daily values in resp. month
    real, dimension(nmonth), intent(in), optional :: mval_pvy  ! vector containing 12 monthly values of the previous year
    real, dimension(nmonth), intent(in), optional :: mval_nxy  ! vector containing 12 monthly values of the next year

    ! function return variable
    real, dimension(ndayyear) :: dval
    
    ! local variables
    integer :: moy, doy, today, dm, iloop
    real :: dd, todaysval

    real, dimension(0:(nmonth+1))    :: mval_ext
    !integer, dimension(0:(nmonth+1)) :: middaymonth_ext
    real :: startt, endt, starttemp, endtemp, dt, d2t, d3t, dtold, &
      dtnew, lastmonthtemp, nextmonthtemp, deltatemp, polya, polyb, polyc
        
    ! xxx implement select case also in '_rates' module

    select case (method)

      case ("interpol")
        !--------------------------------------------------------------------
        ! LINEAR INTERPOLATION
        ! of monthly to quasi-daily values.
        ! If optional argument 'mval_pvy' is provided, take December-value
        ! of previous year to interpolate to first 15 days of January,
        ! otherwise, use the same year's December value to get first 15 days.
        ! corresponds to subroutine 'daily' in LPX
        !--------------------------------------------------------------------

        ! define extended vector with monthly values for previous Dec and next Jan added
        mval_ext(1:nmonth)  = mval(1:nmonth)

        !middaymonth_ext(1:nmonth) = middaymonth(1:nmonth)
        !middaymonth_ext(0) = middaymonth(nmonth)
        !middaymonth_ext(nmonth+1) = 381

        if (present(mval_pvy)) then
          mval_ext(0) = mval_pvy(nmonth)   ! Dec value of previous year
        else
          mval_ext(0) = mval(nmonth)       ! take Dec value of this year ==> leads to jump!
        end if

        if (present(mval_nxy)) then
          mval_ext(nmonth+1) = mval_nxy(1) ! Jan value of next year
        else
          mval_ext(nmonth+1) = mval(1)     ! take Jan value of this year ==> leads to jump!
        end if

        do moy = 1,nmonth
          dd = (mval_ext(moy+1)-mval_ext(moy)) / real(middaymonth(moy+1) - middaymonth(moy))
          todaysval = mval_ext(moy)
          do doy = middaymonth(moy),middaymonth(moy+1)-1
            if (doy<=ndayyear) then
              today = doy
            else
              today = doy-ndayyear
            endif
            dval(today) = todaysval
            todaysval = todaysval + dd
          enddo
        enddo

        if (monthistotal) then
          doy = 0
          do moy=1,nmonth
            do dm=1,ndaymonth(moy)
              doy = doy+1
              dval(doy) = dval(doy) / real(ndaymonth(moy))
            enddo
          enddo
        endif


        !doy=1
        !do moy=1,nmonth
        !  do dm=1,ndaymonth(moy)
        !    doy=doy+1
        !    if (doy>middaymonth(moy)) then
        !      ! interpolate to next month
        !      dval(doy) = mval_ext(moy) + (doy-middaymonth_ext(moy))/ndaymonth_ext(moy) * (mval_ext(moy+1)-mval_ext(moy))
        !    else if (doy<middaymonth(moy)) then
        !      ! interpolate to previous month
        !      dval(doy) = mval_ext(moy-1) + (doy-middaymonth_ext(moy-1))/ndaymonth_ext(moy-1) * (mval_ext(moy)-mval_ext(moy-1))
        !    else
        !      ! take value directly
        !      dval(doy) = mval_ext(moy)
        !    end if
        !  end do
        !end do

      !  !if (iftotals) then
        !  doy=0
        !  do moy=1,nmonth
        !    do doyofmonth=1,ndaymonth(moy)
        !      doy=doy+1
        !      dval(doy)=dval(doy)/dble(ndaymonth(moy))
        !    enddo
        !  enddo
        !endif


      case ("polynom")
        !--------------------------------------------------------------------
        ! In addition to tempdaily daily values are calculated using a polynom of second
        ! order through the middpoints between months. Additionally, average of daily 
        ! values is identical to the monthly input data. That's crucial for modelling
        ! soil heat diffusion and freezing/thawing processes. 
        !--------------------------------------------------------------------!
        if (monthistotal) &
          stop 'MONTHLY2DAILY: no monthly totals allowed for polynom method'
        
        ! Starting conditons of december in previous year
        startt = -30.5               ! midpoint between Nov-Dec of previous year
        endt = 0.5                   ! midpoint between Dec-Jan of this year
        dt = real(ndaymonth(nmonth)) ! number of Dec days
        if (present(mval_pvy)) then
          lastmonthtemp = mval_pvy(nmonth) ! Dec mean temperature
        else
          lastmonthtemp = mval(nmonth)     ! Dec mean temperature
        end if

        doy = 0                      ! initialisation of this years days
        
        do moy=1,nmonth
          dtold = dt
          startt = endt
          endt = endt + dt
          if (moy<nmonth) then
            dtnew = real(ndaymonth(moy+1))
            nextmonthtemp = mval(moy+1)
          else
            dtnew = real(ndaymonth(1))
            if (present(mval_nxy)) then
              nextmonthtemp = mval_nxy(1)
            else
              nextmonthtemp = mval(1)
            end if
          endif

          starttemp = (mval(moy)*dt+lastmonthtemp*dtold)/(dt+dtold)
          endtemp = (nextmonthtemp*dtnew+mval(moy)*dt)/(dtnew+dt)
          deltatemp = endtemp-starttemp
          
          ! Calculate vars for a,b,c coefficients in polynom y = ax^2 +bx + c
          d2t = endt**2.0 - startt**2.0
          d3t = endt**3.0 - startt**3.0

          ! Take a sheet of paper and try solve the polynom, well here is the outcome
          polya = (mval(moy)*dt - deltatemp*d2t/dt/2.0 - starttemp*dt + deltatemp*startt) / (d3t/3.0 - d2t**2.0/dt/2.0 - dt*startt**2.0 + startt*d2t)
          polyb = deltatemp/dt - polya*(startt+endt)
          polyc = starttemp - polya*startt**2.0 - polyb*startt

          ! Calculate daily values with the polynom function
          do dm=1,ndaymonth(moy)
            doy = doy + 1
            dval(doy) = polya * real(doy)**2.0 + polyb * real(doy) + polyc
          enddo
          lastmonthtemp = mval(moy)
        enddo

      case( "uniform" )
        !--------------------------------------------------------------------
        ! Each day in month has the same (monthly) value
        !--------------------------------------------------------------------!      
        doy=0
        do moy=1,nmonth
          do dm=1,ndaymonth(moy)
            doy=doy+1
            dval(doy) = mval(moy)
          end do
        end do

      case default

        stop 'MONTHLY2DAILY: select viable case.'

    end select

  end function monthly2daily


  function ftemp( temp, method, ref_temp )
    !////////////////////////////////////////////////////////////////
    !  Temperature response function
    !----------------------------------------------------------------

    ! arguments
    real, intent(in)             :: temp ! temperature [in Degrees Celsius]
    character(len=*), intent(in) :: method
    real, intent(in), optional   :: ref_temp

    ! local variables
    real                         :: ref_temp_local  ! local copy of ref_temp

    ! for lloyd and taylor method
    real, parameter :: E0 = 308.56      ! Activation Energy
    real, parameter :: T0 = 227.13      ! calibration temperature [K]
    real, parameter :: Tzero = 273.15   ! 0°C = 273.15 K 

    ! function return variable
    real :: ftemp

    ! set default reference temperature to 10 deg C
    if (present(ref_temp)) then
     ref_temp_local = ref_temp
    else
     ref_temp_local = 10.0
    endif

    select case (method)

      case ("lloyd_and_taylor")
        !----------------------------------------------------------------
        ! LLOYD AND TAYLOR
        ! Temperature response function is a modified Q10 relationship
        ! (Lloyd & Taylor 1994)
        !----------------------------------------------------------------
        if (temp.ge.-40.0) then 
          ! avoid numerical errors
          ftemp = exp(E0*((1.0/(ref_temp_local+Tzero-T0))-(1.0/(temp+Tzero-T0))))
        else
          ! set temperature response to a constant at value of -40°C
          ftemp = exp(E0*((1.0/(ref_temp_local+Tzero-T0))-(1.0/(-40.0+Tzero-T0))))
        end if

      case default

        stop 'FTEMP: select valid method'

    end select

    return

  end function ftemp


  function fmoist( moist, method )
    !////////////////////////////////////////////////////////////////
    !  Temperature response function
    !----------------------------------------------------------------
    ! arguments
    real, intent(in)             :: moist ! temperature [in Degrees Celsius]
    character(len=*), intent(in) :: method

    ! function return variable
    real :: fmoist

    select case (method)

      case ("foley")
        !----------------------------------------------------------------
        ! FOLEY
        ! Calculates decomposition rate modifier for a given water fraction
        ! according to Foley 1995
        !----------------------------------------------------------------
        fmoist = (0.25+(0.75*moist))

      case default

        stop 'FMOIST: select valid method'

    end select

    return

  end function fmoist


  function read1year_daily( filename )
    !////////////////////////////////////////////////////////////////
    ! Function reads a file that contains 365 lines, each line for
    ! a daily value. 
    !----------------------------------------------------------------
    use _params_core, only: ndayyear
    implicit none

    ! arguments
    character(len=*), intent(in) :: filename

    ! local variables
    real, dimension(ndayyear) :: dval

    ! function return value
    real, dimension(ndayyear) :: read1year_daily

    open(20, file='./input/'//filename, status='old',  form='formatted', action='read', err=888)
    read(20,*) dval
    close(20)

    read1year_daily = dval

    return
    600 format (F9.7)
    888 write(0,*) 'READ1YEAR_DAILY: error opening file '//trim(filename)//'. Abort. '
    stop

  end function read1year_daily


  function read1year_monthly( filename )
    !////////////////////////////////////////////////////////////////
    ! Function reads a file that contains 365 lines, each line for
    ! a daily value. 
    !----------------------------------------------------------------
    use _params_core, only: nmonth
    implicit none

    ! arguments
    character(len=*), intent(in) :: filename

    ! local variables
    real, dimension(nmonth) :: mval

    ! function return value
    real, dimension(nmonth) :: read1year_monthly

    open(20, file='./input/'//trim(filename), status='old',  form='formatted', action='read', err=888)
    read(20,*) mval
    close(20)

    read1year_monthly = mval

    return
    600 format (F9.7)
    888 write(0,*) 'READ1YEAR_MONTHLY: error opening file ./input/'//trim(filename)//'. Abort. '
    stop

  end function read1year_monthly


  function getvalreal( filename, realyear, day, dm, mo ) result( valreal )
    !////////////////////////////////////////////////////////////////
    !  Function reads one (annual) value corresponding to the given 
    !  year from a time series ascii file. 
    !----------------------------------------------------------------
    use _params_core, only: ndayyear

    ! arguments
    character(len=*), intent(in) :: filename
    integer, intent(in) :: realyear
    integer, intent(in), optional :: day ! day in year (1:365)
    integer, intent(in), optional :: dm  ! day in month (1:31)
    integer, intent(in), optional :: mo  ! month in year (1:12)

    ! function return value
    real :: valreal

    ! local variables
    integer :: l
    real :: tmp(3) ! 3 so that an additional value for this year could be read
    real :: realyear_decimal 

    if (present(day)) then
     ! convert day number into decimal number
     realyear_decimal = real(realyear) + real(day)/real(ndayyear)
    endif

    open(20, file='./input/'//filename, status='old',  form='formatted', err=888)

    if (present(day)) then
     ! find corresponding day in first column and read 3 values on this line
     read(20, 100, err=999) (tmp(l), l=1,3)  
     do while (abs(realyear_decimal-tmp(1))>1.0d-8)
       read(20, 100, err=999) (tmp(l), l=1,3)
     enddo

    else
     ! find corresponding year in first column and read 3 values on this line
     read(20, 100, err=999) (tmp(l), l=1,3)  
     do while (abs(realyear-tmp(1))>1.0d-8)
       read(20, 100, err=999) (tmp(l), l=1,3)
     enddo

    endif

    valreal = tmp(2) 

    100     format (30d16.8)
    close(20)

    return

    888     write(0,*) 'GETVALREAL: error opening file '//trim(filename)//'. Abort. '
    stop
    999     write(0,*) 'GETVALREAL: error reading file '//trim(filename)//'. Abort. '
    stop 

  end function getvalreal
  

  function getvalreal_STANDARD( filename, realyear, mo, dm, day, realyear_decimal ) result( valreal )
    !////////////////////////////////////////////////////////////////
    !  SR reads one (annual) value corresponding to the given year 
    !  from a time series ascii file. File has to be located in 
    !  ./input/ and has to contain only rows formatted like
    !  '2002  1  1 0.496632 0.054053', which represents 
    !  'YYYY MM DM      GPP GPP err.'. DM is the day within the month.
    !  If 'realyear' is dummy (-9999), then it's interpreted as to 
    !  represent a mean climatology for the course of a year.
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: realyear ! year AD as integer
    integer, intent(in), optional :: mo  ! month in year (1:12)
    integer, intent(in), optional :: dm  ! day in month (1:31 or 1:31 or 1:28)
    integer, intent(in), optional :: day ! day in year (1:365)
    real,    intent(in), optional :: realyear_decimal ! year AD as decimal number corresponding to day in the year

    ! function return value
    real :: valreal

    ! local variables
    integer :: inyear
    integer :: inmo
    integer :: indm
    integer :: inday
    real    :: inyear_decimal
    real    :: inval1
    real    :: inval2

    !print*,'looking for realyear, mo, dm',realyear,mo,dm

    ! open file
    open(20, file='./input/'//filename, status='old', form='formatted', err=888)

    if (present(realyear)) then
       ! DATA FOR EACH YEAR
       if (present(mo)) then
           ! DATA FOR EACH MONTH
           if (present(dm)) then
               ! DATA FOR EACH DAY IN THE MONTH
               ! read the 2 values for this day in this year
               read(20, 100, err=999) inyear, inmo, indm, inval1, inval2
               do while ( (realyear-inyear).ne.0 .or. (mo-inmo).ne.0 .or. (dm-indm).ne.0 )
                 read(20, 100, err=999) inyear, inmo, indm, inval1, inval2
               enddo
           else           
               ! read the 2 values for this month in this year
               read(20, 200, err=999) inyear, inmo, inval1, inval2
               do while ( (realyear-inyear).ne.0 .or. (mo-inmo).ne.0 )
                 read(20, 200, err=999) inyear, inmo, inval1, inval2
               enddo
           end if
       else if (present(day)) then
           ! DATA FOR EACH DAY IN THE YEAR
           ! read the 2 values for this day in this year
           read(20, 700, err=999) inyear, inday, inval1, inval2
           do while ( (realyear-inyear).ne.0 .or. (day-inday).ne.0 )
             read(20, 700, err=999) inyear, inday, inval1, inval2
           enddo
       else
           ! read the 2 values for this year
           read(20, 300, err=999) inyear, inval1, inval2
           do while ( (realyear-inyear).ne.0 )
             read(20, 300, err=999) inyear, inval1, inval2
           enddo
       end if
    else if (present(realyear_decimal)) then
      ! DATA PROVIDED FOR EACH DAY AS A DECIMAL OF REALYEAR
      ! find corresponding day in first column and read 3 values on this line
      read(20, 900, err=999) inyear_decimal, inval1, inval2  
      do while (abs(realyear_decimal-inyear_decimal).gt.1.0d-8)
        read(20, 900, err=999) inyear_decimal, inval1, inval2  
      enddo
    else
       ! DATA AS AVERAGE OVER MULTIPLE YEARS (recycle climatology)
       ! FOR EACH MONTH, AND DAY-IN-THE-MONTH
       if (present(mo)) then
           if (present(dm)) then
               ! read the 2 values for this day
               read(20, 400, err=999) inmo, indm, inval1, inval2
               !print*,'inmo, indm, inval1, inval2', inmo, indm, inval1, inval2
               do while ( (mo-inmo).ne.0 .or. (dm-indm).ne.0 )
                 read(20, 400, err=999) inmo, indm, inval1, inval2
                 !print*,'inmo, indm, inval1, inval2', inmo, indm, inval1, inval2
               enddo
           else           
               ! read the 2 values for this month
               read(20, 500, err=999) inmo, inval1, inval2
               do while ( (mo-inmo).ne.0 )
                 read(20, 500, err=999) inmo, inval1, inval2
               enddo

           end if
       else if (present(day)) then
           ! DATA FOR EACH DAY IN THE YEAR
           ! read the 2 values for this day
           read(20, 800, err=999) inday, inval1, inval2
           do while ( (day-inday).ne.0 )
             read(20, 800, err=999) inday, inval1, inval2
           enddo
       else
           ! read the 2 values in this input file
           read(20, 600, err=999) inval1, inval2
       end if
    endif

    !print*,'found realyear, mo, dm      ',inyear,inmo,indm,inval1

    valreal = inval1

    100     format (I4,I3,I3,F9.7,F9.7)
    200     format (I4,I3,F9.7,F9.7)
    300     format (I4,F9.7,F9.7)
    400     format (I3,I3,F9.7,F9.7)
    500     format (I3,F9.7,F9.7)
    600     format (F9.7,F9.7)
    700     format (I4,I4,F9.7,F9.7)
    800     format (I4,F9.7,F9.7)
    900     format (30d16.8,F9.7,F9.7)

    close(20)

    return

    888     write(0,*) 'GETVALREAL_STANDARD: error opening file '//trim(filename)//'. Abort. '
    stop
    999     write(0,*) 'GETVALREAL_STANDARD: error reading file '//trim(filename)//'. Abort. '
    stop 

  end function getvalreal_STANDARD


  function getparreal( filename, paraname ) result( paravalue )
    !////////////////////////////////////////////////////////////////
    !  "Low-level" function for reading parameter values from text file
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filename, paraname

    ! function return value
    real, intent(out) :: paravalue

    ! local variables
    integer :: filehandle
    character(len=40) :: readname, readvalue

    filehandle = 111
    open(filehandle,status='old',err=19,file=filename)
    9    read(filehandle,12,end=10)readname,readvalue
    if (trim(readname)==paraname) then
      read(readvalue,*) paravalue
      goto 11
    else
      goto 9
    endif
    10   continue
    write(0,*) 'getparreal: '//paraname//' of type real not found'
    stop

    11   continue
    12   format(2a40)

    close(filehandle)

    return

    19   write(0,*) 'getparreal: '//filename//' not found!'
    stop

  end function getparreal


  function getparint( filename, paraname ) result( paravalue )
    !////////////////////////////////////////////////////////////////
    ! "Low-level" function for reading parameter values from text file
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filename, paraname 
     
    ! function return variable    
    integer, intent(out) :: paravalue

    ! local variables
    integer :: filehandle
    character(len=40) :: readname, readvalue

    filehandle = 111
    open(filehandle,status='old',err=19,file=filename)
    9    read(filehandle,12,end=10)readname,readvalue
    if (trim(readname)==paraname) then
      read(readvalue,*) paravalue
      goto 11
    else
      goto 9
    endif
    10   continue
    write(0,*) 'getparint: in file '//filename//':'
    write(0,*) 'getparint: '//paraname//' of type integer not found'
    stop

    11   continue
    12   format(2a40)

    close(filehandle)

    return

    19   write(0,*) 'getparint: file '//filename//' not found:'
    write(0,*) filename
    stop

  end function getparint


  function getparlogical( filename, paraname ) result( paravalue )
    !////////////////////////////////////////////////////////////////
    ! "Low-level" function for reading parameter values from text file
    !----------------------------------------------------------------
    ! arguments
    character(len=*), intent(in) :: filename, paraname
    ! function return variable
    logical, intent(out) :: paravalue

    ! local variables
    integer :: filehandle
    character(len=40) :: readname, readvalue

    filehandle = 111
    open(filehandle,status='old',err=19,file=filename)
    9    read(filehandle,12,end=10)readname,readvalue
    if (trim(readname)==paraname) then
      read(readvalue,*) paravalue
      goto 11
    else
      goto 9
    endif
    10   continue
    write(0,*) 'GETPARLOGICAL: '//paraname//' of type logical not found'
    stop

    11   continue
    12   format(2a40)

    close(filehandle)

    return

    19   write(0,*) 'GETPARLOGICAL: '//filename//' not found!'
    stop

  end function getparlogical


  subroutine getparstring( filename, paraname, paravalue )
    !////////////////////////////////////////////////////////////////
    ! "Low-level" function for reading parameter values from text file
    !----------------------------------------------------------------
    character(len=*), intent(in) :: filename, paraname
    character(len=*), intent(out) :: paravalue

    integer :: filehandle,i
    character(len=40) :: readname
    character(len=1024) :: readvalue

    filehandle = 111
    open(filehandle,status='old',err=19,file=filename)
    9    read(filehandle,12,end=10)readname,readvalue
    if (trim(readname)==paraname) then
      ! Strip leading and trailing whitespace from readvalue
      i=1
      do while (readvalue(i:i)==' ' .and. i.lt.len(readvalue))
        i=i+1
      enddo
      paravalue = trim(readvalue(i:))
      goto 11
    else
      goto 9
    endif
    10   continue
    write(0,*) 'GETSTRING: '//paraname//' of type string not found'
    stop

    11   continue
    12   format(a40,a)

    close(filehandle)

    return

    19   write(0,*) 'GETSTRING: '//filename//' not found!'
    stop

  end subroutine getparstring


  function getparchar( filename, paraname ) result( paravalue )
    !////////////////////////////////////////////////////////////////
    ! "Low-level" function for reading parameter values from text file
    !----------------------------------------------------------------
    character(len=*), intent(in) :: filename, paraname
    character(len=1), intent(out) :: paravalue

    integer :: filehandle
    character(len=40) :: readname, readvalue

    filehandle = 111
    open(filehandle,status='old',err=19,file=filename)
    9    read(filehandle,12,end=10)readname,readvalue
    if (trim(readname)==paraname) then
      read(readvalue,*) paravalue
      goto 11
    else
      goto 9
    endif
    10   continue
    write(0,*) 'getparchar: '//paraname//' of type char not found'
    stop

    11   continue
    12   format(2a40)

    close(filehandle)

    return

    19   write(0,*) 'getparchar: '//filename//' not found!'
    stop

  end function getparchar

 !  subroutine setreal(filename,paraname,default,paravalue)
 !          !////////////////////////////////////////////////////////////////
 !          ! "Low-level" SR for reading/defining parameter values from text file
 !          !----------------------------------------------------------------

 !          implicit none

 !          integer filehandle
 !          character(len=40) :: readname, readvalue
 !          character(len=*), intent(in) :: filename, paraname

 !          real,intent(out) :: paravalue
 !          real,intent(in)  :: default

 !          filehandle = 111
 !          open(filehandle,status='old',err=19,file=filename)
 ! 9        read(filehandle,12,end=10)readname,readvalue
 !          if (trim(readname)==paraname) then
 !             read(readvalue,*) paravalue
 !             goto 11
 !          else
 !             goto 9
 !          endif
 ! 10       continue
 !          paravalue = default

 ! 11       continue
 ! 12       format(2a40)

 !          close(filehandle)

 !          return

 ! 19       stop 'GETPAR: '//filename//' not found!'

 !          end

 !          subroutine setint(filename,paraname,default,paravalue)
 !          !////////////////////////////////////////////////////////////////
 !          ! "Low-level" SR for reading/defining parameter values from text file
 !          !----------------------------------------------------------------

 !          implicit none

 !          integer filehandle
 !          character(len=40) :: readname, readvalue
 !          character(len=*), intent(in) :: filename, paraname

 !          integer, intent(out) :: paravalue
 !          integer :: default

 !          filehandle = 111
 !          open(filehandle,status='old',err=19,file=filename)
 ! 9        read(filehandle,12,end=10)readname,readvalue
 !          if (trim(readname)==paraname) then
 !             read(readvalue,*) paravalue
 !             goto 11
 !          else
 !             goto 9
 !          endif
 ! 10       continue
 !          paravalue = default

 ! 11       continue
 ! 12       format(2a40)

 !          close(filehandle)

 !          return

 ! 19       stop 'GETPAR: '//filename//' not found!'

 !          end

 !          subroutine setlogical(filename,paraname,default,paravalue)
 !          !////////////////////////////////////////////////////////////////
 !          ! "Low-level" SR for reading/defining parameter values from text file
 !          !----------------------------------------------------------------

 !          implicit none

 !          character(len=*), intent(in) :: filename, paraname
 !          logical, intent(out) :: paravalue

 !          integer filehandle
 !          character(len=40) :: readname, readvalue
 !          logical :: default

 !          filehandle = 111
 !          open(filehandle,status='old',err=19,file=filename)
 ! 9        read(filehandle,12,end=10)readname,readvalue
 !          if (trim(readname)==paraname) then
 !             read(readvalue,*) paravalue
 !             goto 11
 !          else
 !             goto 9
 !          endif
 ! 10       continue
 !          paravalue = default

 ! 11       continue
 ! 12       format(2a40)

 !          close(filehandle)

 !          return

 ! 19       stop 'GETPAR: '//filename//' not found!'

 !          end


 !          subroutine setstring(filename,paraname,default,paravalue)
 !          !////////////////////////////////////////////////////////////////
 !          ! "Low-level" SR for reading/defining parameter values from text file
 !          !----------------------------------------------------------------

 !          implicit none

 !          ! Arguments
 !          character(len=*), intent(in) :: filename, paraname
 !          character(len=*), intent(out) :: paravalue, default

 !          ! Local variables
 !          integer filehandle,i
 !          character readname*40,readvalue*1024

 !          filehandle = 111
 !          open(filehandle,status='old',err=19,file=filename)
 ! 9        read(filehandle,12,end=10)readname,readvalue
 !          if (trim(readname)==paraname) then
 !             ! Strip leading and trailing whitespace from readvalue
 !             i=1
 !             do while (readvalue(i:i)==' ' .and. i.lt.len(readvalue))
 !                i=i+1
 !             enddo
 !             paravalue = trim(readvalue(i:))
 !             goto 11
 !          else
 !             goto 9
 !          endif
 ! 10       continue
 !          paravalue = default

 ! 11       continue

 ! 12       format(a40,a)

 !          close(filehandle)

 !          return

 ! 19       stop 'GETPAR: '//filename//' not found!'

 !          end



  ! function icumsum( arr, seed )
  !   !/////////////////////////////////////////////////////////////////////////
  !   ! Returns the running (cumulative) sum of a vector of integers. 
  !   ! Adopted from: 
  !   ! Numerical Recipes in Fortran 90: Volume 2, edited by William H. Press
  !   !-------------------------------------------------------------------------
  !   implicit none 

  !   ! arguments
  !   integer, dimension(:), intent(in)  :: arr
  !   integer, intent(in), optional      :: seed
  !   integer, dimension(:), intent(out) :: icumsum

  !   ! local variables
  !   integer :: n, j
  !   integer :: sd

  !   n = size(arr)
  !   if (n==0) return
  !   sd = 0
  !   if ( present(seed) ) sd = seed
  !   icumsum(1) = arr(1) + sd
  !   do j=2,n
  !     icumsum(j) = icumsum(j-1) + arr(j)
  !   end do

  ! end function icumsum


end module _sofunutils
