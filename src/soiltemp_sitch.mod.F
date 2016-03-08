module _soiltemp
  !////////////////////////////////////////////////////////////////
  ! SITCH SOILTEMP MODULE
  ! Contains the "main" subroutine 'soiltemp' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'soiltemp' must contain this list 
  ! of subroutines (names that way), in order to be exchangeable
  ! with this module:
  !   - getpar_modl_soiltemp
  !   - initio_soiltemp
  !   - initoutput_soiltemp
  !   - getout_daily_soiltemp
  !   - getout_monthly_soiltemp
  !   - writeout_ascii_soiltemp
  !   - soiltemp
  ! Required module-independent model state variables (necessarily 
  ! updated by 'waterbal') are:
  !   - soil moisture ('dwtot')
  !   - soil temperature ('psoilphys%temp')
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _vars_core, only: psoilphys

  implicit none

contains

  subroutine soiltemp( jpngr, moy, day, dtemp ) 
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates soil temperature based on.
    !-------------------------------------------------------------------------
    use _params_core, only: ndayyear, nlu, nmonth, maxgrid, ndaymonth, pi
    use _params_siml, only: init
    use _sofunutils, only: running, daily2monthly
    use _vars_core, only: psoilphys, soilphys

    ! arguments
    integer, intent(in)                   :: jpngr
    integer, intent(in)                   :: moy
    integer, intent(in)                   :: day                            ! current day of year
    real, dimension(ndayyear), intent(in) :: dtemp        ! daily temperature (deg C)

    ! local variables
    real, dimension(ndayyear,maxgrid), save     :: dtemp_pvy    ! daily temperature of previous year (deg C)
    real, dimension(nlu,ndayyear,maxgrid), save :: wscal_pvy ! daily Cramer-Prentice-Alpha of previous year (unitless) 
    real, dimension(nlu,ndayyear), save         :: wscal_alldays

    !real, dimension(nmonth) :: mtemp                     ! monthly temperature as a mean of daily values in resp. month
    !real, dimension(nmonth) :: mtemp_pvy                 ! monthly temperature as a mean of daily values in resp. month previous year

    !real, dimension(ndayyear), save :: dtemp_buf        ! daily temperature vector containing values of the present day and the preceeding 364 days. Updated daily. (deg C)
    !real, dimension(ndayyear), save :: dwtot_buf        ! daily soil moisture content, containing values of the present day and the preceeding 364 days. Updated daily

    integer :: pm ,ppm, lu

    real :: avetemp, meanw1
    real :: tempthismonth, templastmonth
    real :: diffus_wp
    real :: diffus_15whc
    real :: diffus_fc
    real :: diffus
    real :: alag, amp, lag, lagtemp

    ! xxx read soil parameters from parameter input file
    diffus_wp    = 0.2   ! soilpar(5,jpngr) = store(3,soilcode(jpngr))
    diffus_15whc = 0.65  ! soilpar(6,jpngr) = store(4,soilcode(jpngr))
    diffus_fc    = 0.4   ! soilpar(7,jpngr) = store(5,soilcode(jpngr))

    ! in first year, use this years air temperature (available for all days in this year)
    if ( init .and. day==1 ) then
      dtemp_pvy(:,jpngr) = dtemp(:)
    end if

    wscal_alldays(:,day) = soilphys(:)%wscal

    avetemp = running( dtemp, day, ndayyear, ndayyear, "mean", dtemp_pvy(:,jpngr) ) 

    ! get monthly mean temperature vector from daily vector
    !mtemp     = daily2monthly( dtemp,     "mean" )
    !mtemp_pvy = daily2monthly( dtemp_pvy, "mean" )

    ! get average temperature of the preceeding N days in month (30/31/28 days)
    if (moy==1) then
      pm = 12
      ppm = 11
    else if (moy==2) then
      pm = 1
      ppm = 12
    else
      pm = moy - 1
      ppm = moy - 2
    end if
    tempthismonth = running( dtemp, day, ndayyear, ndaymonth(pm), "mean", dtemp_pvy(:,jpngr) )
    templastmonth = running( dtemp, modulo(day-ndaymonth(pm),ndayyear), ndayyear, ndaymonth(ppm), "mean", dtemp_pvy(:,jpngr) )


    do lu=1,nlu
      !-------------------------------------------------------------------------
      ! recalculate running mean of previous 12 month's temperature and soil moisture
      ! avetemp stores running mean temperature of previous 12 months.
      ! meanw1 stores running mean soil moisture in layer 1 of previous 12 months 
      !-------------------------------------------------------------------------
      if (init) then
        meanw1  = running( wscal_alldays(lu,:), day, ndayyear, ndayyear, "mean"  )
      else
        meanw1  = running( wscal_alldays(lu,:), day, ndayyear, ndayyear, "mean", wscal_pvy(lu,:,jpngr)  )
      end if

      ! In case of zero soil water, return with soil temp = air temp
      if (meanw1==0.0) then
        psoilphys(lu,jpngr)%temp = dtemp(day)
        return
      endif
          
      ! Interpolate thermal diffusivity function against soil water content
      if (meanw1<0.15) then
        diffus = (diffus_15whc-diffus_wp)/0.15*meanw1+diffus_wp
      else
        diffus = (diffus_fc-diffus_15whc)/0.85*(meanw1-0.15)+diffus_15whc
      endif
          
      ! Convert diffusivity from mm2/s to m2/month
      ! multiplication by 1e-6 (-> m2/s) * 2.628e6 (s/month)  =  2.628
      diffus = diffus*2.628
          
      ! Calculate amplitude fraction and lag at soil depth 0.25 m
      alag = 0.25/sqrt(12.0*diffus/pi)
      amp  = exp(-alag)
      lag  = alag*(6.0/pi)                                 !convert lag from angular units to months
          
      ! Calculate monthly soil temperatures for this year.  For each month,
      ! calculate average air temp for preceding 12 months (including this one)
          
      ! Estimate air temperature "lag" months ago by linear interpolation
      ! between air temperatures for this and last month
      lagtemp = (tempthismonth-templastmonth) * (1.0-lag) + templastmonth
          
      ! Adjust amplitude of lagged air temp to give estimated soil temp
      psoilphys(lu,jpngr)%temp = avetemp + amp * (lagtemp-avetemp)

    end do

    ! save temperature for next year
    if (day==ndayyear) then
      dtemp_pvy(:,jpngr) = dtemp(:)
      wscal_pvy(:,:,jpngr) = wscal_alldays(:,:)
    end if

    return

  end subroutine soiltemp

end module _soiltemp
