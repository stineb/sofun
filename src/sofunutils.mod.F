module sofunutils
  !/////////////////////////////////////////////////////////////////////////
  ! Contains utility functions to deal with arrays containing 365 daily 
  ! values or 12 monthly values representing one year.
  ! - running: calculates running mean
  ! - daily2monthly: calculates monthly mean/sum of daily values for resp. mo.
  !-------------------------------------------------------------------------
  implicit none

contains

  function running( presval, inow, lenval, lenper, method, prevval )
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
    real, intent(out) :: running

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
        running = sum(valbuf((lenval-lenper+1):lenval))
      case ("mean")
        if (present(prevval)) then
          running = sum(valbuf((lenval-lenper+1):lenval))/lenper
        else
          running = sum(valbuf((lenval-lenper+1):lenval))/inow
        end if
      case default
        stop 'RUNNING: declare valid method.'
    end select

  end function running


  function daily2monthly( dval, method )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns monthly values as a mean over daily values in each month.
    ! Arguments:
    ! - dval   : vector containing 365 (366 in case lapyear is TRUE) daily values
    ! - method : true of monthly values represent total of daily values in resp. month
    !-------------------------------------------------------------------------
    use params_core, only: ndaymonth, cumdaymonth, ndayyear, nmonth

    ! arguments
    real, intent(in), dimension(ndayyear) :: dval
    character(len=*), intent(in) :: method

    ! function return variable
    real, dimension(nmonth), intent(out) :: daily2monthly

    ! local variables
    integer :: moy
    integer, dimension(nmonth) :: istart, iend

    istart = cumdaymonth - ndaymonth + 1
    iend   = cumdaymonth

    ! loop over months and take sum/mean of days in that month
    do moy=1,nmonth
      select case (method)
        case ("sum")
          daily2monthly(moy) = sum( dval( istart(moy) : iend(moy) ))
        case ("mean")
          daily2monthly(moy) = sum( dval( istart(moy) : iend(moy) )) / ndaymonth(moy)
        case default
          stop 'DAILY2MONTHLY: select valid method (sum, mean)' 
      end select
    end do

    return

  end function daily2monthly


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


  function monthly2daily( mval, method, monthistotal, mval_pvy, mval_nxy )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns daily values based on monthly values, using a defined method.
    !-------------------------------------------------------------------------
    use params_core, only: middaymonth, ndayyear, ndaymonth, nmonth
    
    ! arguments
    real, dimension(nmonth), intent(in) :: mval  ! vector containing 12 monthly values
    character(len=*), intent(in) :: method
    logical, intent(in) :: monthistotal ! true of monthly values represent total of daily values in resp. month
    real, dimension(nmonth), intent(in), optional :: mval_pvy  ! vector containing 12 monthly values of the previous year
    real, dimension(nmonth), intent(in), optional :: mval_nxy  ! vector containing 12 monthly values of the next year
    
    ! local variables
    integer :: moy, doy, today, dm
    real :: dd, todaysval

    real, dimension(0:(nmonth+1))    :: mval_ext
    !integer, dimension(0:(nmonth+1)) :: middaymonth_ext
    real :: startt, endt, starttemp, endtemp, dt, d2t, d3t, dtold, &
      dtnew, lastmonthtemp, nextmonthtemp, deltatemp, polya, polyb, polyc
    
    ! function return variable
    real, dimension(ndayyear) :: monthly2daily
    
    ! xxx implement select case also in 'rates' module

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
            monthly2daily(today) = todaysval
            todaysval = todaysval + dd
          enddo
        enddo

        if (monthistotal) then
          doy = 0
          do moy=1,nmonth
            do dm=1,ndaymonth(moy)
              doy = doy+1
              monthly2daily(doy) = monthly2daily(doy) / real(ndaymonth(moy))
            enddo
          enddo
        endif


        !doy=1
        !do moy=1,nmonth
        !  do dm=1,ndaymonth(moy)
        !    doy=doy+1
        !    if (doy>middaymonth(moy)) then
        !      ! interpolate to next month
        !      monthly2daily(doy) = mval_ext(moy) + (doy-middaymonth_ext(moy))/ndaymonth_ext(moy) * (mval_ext(moy+1)-mval_ext(moy))
        !    else if (doy<middaymonth(moy)) then
        !      ! interpolate to previous month
        !      monthly2daily(doy) = mval_ext(moy-1) + (doy-middaymonth_ext(moy-1))/ndaymonth_ext(moy-1) * (mval_ext(moy)-mval_ext(moy-1))
        !    else
        !      ! take value directly
        !      monthly2daily(doy) = mval_ext(moy)
        !    end if
        !  end do
        !end do

      !  !if (iftotals) then
        !  doy=0
        !  do moy=1,nmonth
        !    do doyofmonth=1,ndaymonth(moy)
        !      doy=doy+1
        !      monthly2daily(doy)=monthly2daily(doy)/dble(ndaymonth(moy))
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
            monthly2daily(doy) = polya * real(doy)**2.0 + polyb * real(doy) + polyc
          enddo
          lastmonthtemp = mval(moy)
        enddo

      case default

        stop 'MONTHLY2DAILY: select viable case.'

    end select

  end function monthly2daily


!  function monthly2daily_precip( mval_prec, mval_wet )!

!    use params_core, only: nmonth, ndayyear

!    ! arguments
!    real, dimension(nmonth), intent(in) :: mval_prec  ! monthly precipitation totals
!    real, dimension(nmonth), intent(in) :: mval_wet   ! monthly number of wet days

!    ! local variables
!    integer :: dm,day,d,m
!    integer :: daysum            !accumulated days at beginning of month
!    
!    real :: prob, vv, v1, random
!    real, parameter :: c1 = 1.0
!    real, parameter :: c2 = 1.2
!    real, dimension(nmonth) :: prob_rain
!    real, dimension(nmonth) :: mprecave     !average precipitation on wet days
!    real, dimension(nmonth) :: mprecip      !acc. monthly generated precipitation
!    real, dimension(ndayyear) :: dval_prec_save

!    ! function return variable
!    real, dimension(ndayyear), intent(out) :: dval_prec


!    day=0
!    prob=0.0
!    do m=1,nmonth
!      prob_rain(m)=0.0
!      mprecave(m)=0.0
!      mprecip(m)=0.0
!    enddo
!    daysum=0
!   
!    do m=1,nmonth
!      if(mval_wet(m)<=1.0) mval_wet(m)=1.0
!      prob_rain(m)=mval_wet(m)/real(ndaymonth(m))
!      mprecave(m)=mval_prec(m)/mval_wet(m)
!       
!      do dm=1,ndaymonth(m)
!        day=day+1
!   
!        ! Transitional probabilities (Geng et al. 1986)
!        if (day>1) then
!          if (dval_prec(day-1)<0.1) then
!            prob=0.75*prob_rain(m)
!          else 
!            prob=0.25+(0.75*prob_rain(m))
!          endif
!        endif
!        
!        ! Determine wet days randomly and use Krysanova/Cramer estimates of 
!        ! parameter values (c1,c2) for an exponential distribution
!        vv=real(prdaily_random(jpngr,day,1))
!        if(vv>prob) then
!          dval_prec(day)=0.0
!        else
!          v1=real(prdaily_random(jpngr,day,2))
!          dval_prec(day)=((-log(v1))**c2)*mprecave(m)*c1 
!          if(dval_prec(day)<0.1) dval_prec(day)=0.0
!        endif
!        mprecip(m)=mprecip(m)+dval_prec(day)
!        dval_prec_save(day)=dval_prec(day) !
!      enddo
!      
!      ! normalise generated precipitation by monthly values
!      if(m>1) daysum=daysum+ndaymonth(m-1)
!      if(mprecip(m)<1.0) mprecip(m)=1.0
!      do dm=1,ndaymonth(m)
!        day=daysum+dm
!        dval_prec(day)=dval_prec(day)*(mval_prec(m)/mprecip(m))
!        if (dval_prec(day)<0.1d0) dval_prec(day)=0.0
!      ! dval_prec(day)=mval_prec(m)/ndaymonth(m)  !no generator
!      enddo
!      
!    enddo                     !month 
!    
!    ! calculate monthly means after interpolation - not absolutely identical to input
!    day=0
!    do m=1,nmonth
!      mprecint(m)=0.0
!      do d=1,ndaymonth(m)
!        day = day + 1
!        mprecint(m)=mprecint(m)+dval_prec(day)!/dble(ndaymonth(m))
!      enddo
!    enddo
!          
!    return!

!  end function monthly2daily_precip!
!      
!  function random(k)!
!  
!    ! arguments
!    integer, dimension(4), intent(in) :: k
!    ! local variables
!    integer :: i

!    k(4)=3*k(4)+k(2)
!    k(3)=3*k(3)+k(1)
!    k(2)=3*k(2)
!    k(1)=3*k(1)
!    i=k(1)/1000
!    k(1)=k(1)-i*1000
!    k(2)=k(2)+i
!    i=k(2)/100
!    k(2)=k(2)-100*i
!    k(3)=k(3)+i
!    i=k(3)/1000
!    k(3)=k(3)-i*1000
!    k(4)=k(4)+i
!    i=k(4)/100
!    k(4)=k(4)-100*i
!    random=(((dble(k(1))*0.001d0+dble(k(2)))*0.01d0+dble(k(3))) * 0.001d0+dble(k(4)))*0.01d0
!    
!  end function random          
      

end module sofunutils