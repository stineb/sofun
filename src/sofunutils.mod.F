module sofunutils
    !/////////////////////////////////////////////////////////////////////////
    ! Contains utility functions to deal with arrays containing 365 daily 
    ! values or 12 monthly values representing one year.
    ! - running: calculates running mean
    ! - daily2monthly: calculates monthly mean/sum of daily values for resp. mo.
    !-------------------------------------------------------------------------

contains

  function running( presval, inow, lenval, lenper, method, prevval )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns running sum or average over. 'prevval' is optional, if not pro-
    ! vided, sum/average is taken only over preceeding days/months of current
    ! year.
    !-------------------------------------------------------------------------

    implicit none

    ! arguments
    real, dimension(lenval), intent(in) :: presval            ! vector containing 'lenvals' values for each timestep in this year
    integer, intent(in) :: inow                               ! index corresponding to "now" (day of year or month of year)  
    integer, intent(in) :: lenval                             ! number of timesteps per year
    integer, intent(in) :: lenper                             ! number of timesteps over which to average/sum
    character(len=*), intent(in) :: method                    ! either "sum" or "mean" for running sum or running mean
    real, dimension(lenval), intent(in), optional :: prevval  ! vector containing 'lenvals' values for each timestep in previous year

    ! local variables
    real, dimension(lenval) :: valbuf

    ! function return variable
    real, dimension(lenval), intent(out) :: running

    if (present(prevval)) then
      valbuf((lenval-(inow-1)):lenval) = presval(1:inow)
      valbuf(1:(lenval-inow)) = prevval((inow+1):lenval)
    else
      valbuf(1:inow) = presval(1:inow)
      valbuf((inow+1):lenval) = 0.0
    end if

    select case (method)
      case ("sum")
        running = sum(valbuf)
      case ("mean")
        if (present(prevval)) then
          running = sum(valbuf((lenval-lenper+1):lenval))/lenval
        else
          running = sum(valbuf((lenval-lenper+1):lenval))/inow
        end if
      case default
        stop 'RUNNING: declare valid method.'
    end select

  end function runningmean


  function daily2monthly( dval, method )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns monthly values as a mean over daily values in each month
    !-------------------------------------------------------------------------
    use params_core, only: ndayyear, ndaymonth, nmonth
    implicit none

    ! arguments
    real, dimension(ndayyear), intent(in) :: dval  ! vector containing 365 daily values
    character(len=*), intent(in) :: method ! true of monthly values represent total of daily values in resp. month

    ! function return variable
    real, dimension(nmonth), intent(out) :: daily2monthly

    do moy=1,nmonth
      select case (method)
        case ("sum")
          daily2monthly(moy) = sum(dval(sum(ndaymonth(1:(moy-1)))+1:sum(ndaymonth(1:moy))))
        case ("mean")
          daily2monthly(moy) = sum(dval(sum(ndaymonth(1:(moy-1)))+1:sum(ndaymonth(1:moy))))/ndaymonth(moy)
        case default
          stop 'DAILY2MONTHLY: select valid method' 
      end select
    end do

  end function daily2monthly


!  function monthly2daily( mval, method, monthistotal, mval_pvy, mval_nxy )
!    !/////////////////////////////////////////////////////////////////////////
!    ! Returns daily values based on monthly values, using a defined method.
!    !-------------------------------------------------------------------------
!    use params_core, only: middaymonth, ndayyear, ndaymonth!

!    ! arguments
!    real, dimension(nmonth), intent(in) :: mval  ! vector containing 12 monthly values
!    character(len=*), intent(in) :: method
!    logical, intent(in) :: monthistotal ! true of monthly values represent total of daily values in resp. month
!    real, dimension(nmonth), intent(in), optional :: mval_pvy  ! vector containing 12 monthly values of the previous year
!    real, dimension(nmonth), intent(in), optional :: mval_nxy  ! vector containing 12 monthly values of the next year!

!    ! local variables
!    real :: startt, endt, starttemp, endtemp, dt, d2t, d3t, dtold, &
!      dtnew, lastmonthtemp, nextmonthtemp, deltatemp, polya, polyb, polyc!

!    ! function return variable
!    real, dimension(ndayyear) :: monthly2daily !

!    ! implement select case also in 'rates' module
!    select case (method)!

!      case ("interpol")
!        !--------------------------------------------------------------------
!        ! LINEAR INTERPOLATION
!        ! of monthly to quasi-daily values.
!        ! If optional argument 'mval_pvy' is provided, take December-value
!        ! of previous year to interpolate to first 15 days of January,
!        ! otherwise, use the same year's December value to get first 15 days.
!        ! corresponds to subroutine 'daily' in LPX
!        !--------------------------------------------------------------------
!        doy=1
!        do moy=1,nmonth
!          do dm=1,ndaymonth
!            doy=doy+1
!            if (doy>middaymonth(moy)) then
!              monthly2daily(doy) = mval(moy) + (doy-middaymonth(moy))/ndaymonth(moy) * (mval(moy+1)-mval(moy))
!            else if (doy<middaymonth(moy)) then
!              if (moy==1) then
!                if (present(mval_pvy)) then
!                  monthly2daily(doy) = mval_pvy(nmonth) + (middaymonth(nmonth)+doy)/ndaymonth(nmonth) * (mval(moy)-mval_pvy(nmonth))
!                else
!                  monthly2daily(doy) = mval(nmonth) + (middaymonth(nmonth)+doy)/ndaymonth(nmonth) * (mval(moy)-mval(nmonth))
!                end if
!              else
!                monthly2daily(doy) = mval(moy-1) + (doy-middaymonth(moy-1))/ndaymonth(moy-1) * (mval(moy)-mval(moy-1))
!              end if
!            else
!              monthly2daily(doy) = mval(moy)
!            end if
!          end do
!        end do!

!        if (iftotals) then
!          doy=0
!          do moy=1,nmonth
!            do doyofmonth=1,ndaymonth(moy)
!              doy=doy+1
!              monthly2daily(doy)=monthly2daily(doy)/dble(ndaymonth(moy))
!            enddo
!          enddo
!        endif!
!

!      case ("polynom")
!        !--------------------------------------------------------------------
!        ! In addition to tempdaily daily values are calculated using a polynom of second
!        ! order through the middpoints between months. Additionally, average of daily 
!        ! values is identical to the monthly input data. That's crucial for modelling
!        ! soil heat diffusion and freezing/thawing processes. 
!        !--------------------------------------------------------------------!

!        if (monthistotal) &
!          stop 'MONTHLY2DAILY: no monthly totals allowed for polynom method'
!        
!        ! Starting conditons of december in previous year
!        startt = -30.5               ! midpoint between Nov-Dec of previous year
!        endt = 0.5                   ! midpoint between Dec-Jan of this year
!        dt = real(ndaymonth(month))  ! number of Dec days
!        if (present(mval_pvy)) then
!          lastmonthtemp = mval_pvy[nmonth] ! Dec mean temperature
!        else
!          lastmonthtemp = mval[nmonth]     ! Dec mean temperature
!        end if!

!        doy = 0                      ! initialisation of this years days
!        
!        do month=1,nmonth
!          dtold = dt
!          startt = endt
!          endt = endt + dt
!          if (month<nmonth) then
!            dtnew = real(ndaymonth(month+1))
!            nextmonthtemp = mval(month+1)
!          else
!            dtnew = real(ndaymonth(1))
!            if (present(mval_nxy)) then
!              nextmonthtemp = mval_nxy(1)
!            else
!              nextmonthtemp = mval(1)
!            end if
!          endif!

!          starttemp = (mval(month)*dt+lastmonthtemp*dtold)/(dt+dtold)
!          endtemp = (nextmonthtemp*dtnew+mval(month)*dt)/(dtnew+dt)
!          deltatemp = endtemp-starttemp
!          
!          ! Calculate vars for a,b,c coefficients in polynom y = ax^2 +bx + c
!          d2t = endt**2.0 - startt**2.0
!          d3t = endt**3.0 - startt**3.0!

!          ! Take a sheet of paper and try solve the polynom, well here is the outcome
!          polya = (mval(month)*dt - deltatemp*d2t/dt/2.0 - starttemp*dt + deltatemp*startt) / (d3t/3.0 - d2t**2.0/dt/2.0 - dt*startt**2.0 + startt*d2t)
!          polyb = deltatemp/dt - polya*(startt+endt)
!          polyc = starttemp - polya*startt**2.0 - polyb*startt!

!          ! Calculate daily values with the polynom function
!          do dm=1,ndaymonth(month)
!            doy = doy + 1
!            monthly2daily(doy) = polya * real(doy)**2.0 + polyb * real(doy) + polyc
!          enddo
!          lastmonthtemp = mval(month)
!        enddo!

!      case default!

!        stop 'MONTHLY2DAILY: select viable case.'!

!    end select!

!  end function monthly2daily!
!

!  function monthly2daily_precip( mval_prec, mval_wet )!

!    use params_core, only: nmonth, ndayyear

!    implicit none

!    ! arguments
!    real, dimension(nmonth), intent(in) :: mval_prec  ! monthly precipitation totals
!    real, dimension(nmonth), intent(in) :: mval_wet   ! monthly number of wet days

!    ! local variables
!    integer :: dayofmonth,day,d,m
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
!      do dayofmonth=1,ndaymonth(m)
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
!      do dayofmonth=1,ndaymonth(m)
!        day=daysum+dayofmonth
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
!    random=(((dble(k(1))*0.001d0+dble(k(2)))*0.01d0+dble(k(3))) * 0.001d0+dble(k(4)))*0.01d0!

!  end function random          
!      

end module sofunutils