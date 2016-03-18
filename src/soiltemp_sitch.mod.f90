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
  !   - soil temperature ('soiltemp')
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: nlu, maxgrid

  implicit none

  !----------------------------------------------------------------
  ! Module-specific state variables
  !----------------------------------------------------------------
  real, dimension(nlu,maxgrid) :: dtemp_soil          ! soil temperature [deg C]

  !----------------------------------------------------------------
  ! Module-specific daily output variables
  !----------------------------------------------------------------
  real, allocatable, dimension(:,:,:) :: outdtemp_soil

contains

  subroutine soiltemp( jpngr, moy, day, dtemp, params_soil ) 
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates soil temperature based on.
    !-------------------------------------------------------------------------
    use _params_core, only: ndayyear, nlu, maxgrid, ndaymonth, pi
    use _params_siml, only: init
    use _sofunutils, only: running, daily2monthly
    use _waterbal, only: soilphys
    use _params_soil, only: paramtype_soil


    ! arguments
    integer, intent(in)                   :: jpngr
    integer, intent(in)                   :: moy
    integer, intent(in)                   :: day                            ! current day of year
    real, dimension(ndayyear), intent(in) :: dtemp        ! daily temperature (deg C)
    type( paramtype_soil ), intent(in)    :: params_soil

    ! local variables
    real, dimension(ndayyear,maxgrid), save     :: dtemp_pvy    ! daily temperature of previous year (deg C)
    real, dimension(nlu,ndayyear,maxgrid), save :: wscal_pvy ! daily Cramer-Prentice-Alpha of previous year (unitless) 
    real, dimension(nlu,ndayyear), save         :: wscal_alldays

    !real, dimension(ndayyear), save :: dtemp_buf        ! daily temperature vector containing values of the present day and the preceeding 364 days. Updated daily. (deg C)
    !real, dimension(ndayyear), save :: dwtot_buf        ! daily soil moisture content, containing values of the present day and the preceeding 364 days. Updated daily

    integer :: pm ,ppm, lu

    real :: avetemp, meanw1
    real :: tempthismonth, templastmonth
    real :: diffus
    real :: alag, amp, lag, lagtemp


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
    templastmonth = running( dtemp, modulo( day - ndaymonth(pm), ndayyear ), ndayyear, ndaymonth(ppm), "mean", dtemp_pvy(:,jpngr) )


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
        dtemp_soil(lu,jpngr) = dtemp(day)
        return
      endif
          
      ! Interpolate thermal diffusivity function against soil water content
      if (meanw1<0.15) then
        diffus = ( params_soil%thdiff_whc15 - params_soil%thdiff_wp ) / 0.15 * meanw1 + params_soil%thdiff_wp
      else
        diffus = ( params_soil%thdiff_fc - params_soil%thdiff_whc15 ) / 0.85 * ( meanw1 - 0.15 ) + params_soil%thdiff_whc15
      endif
          
      ! Convert diffusivity from mm2/s to m2/month
      ! multiplication by 1e-6 (-> m2/s) * 2.628e6 (s/month)  =  2.628
      diffus = diffus * 2.628
          
      ! Calculate amplitude fraction and lag at soil depth 0.25 m
      alag = 0.25 / sqrt( 12.0 * diffus / pi )
      amp  = exp(-alag)
      lag  = alag * ( 6.0 / pi )                                 !convert lag from angular units to months
          
      ! Calculate monthly soil temperatures for this year.  For each month,
      ! calculate average air temp for preceding 12 months (including this one)
          
      ! Estimate air temperature "lag" months ago by linear interpolation
      ! between air temperatures for this and last month
      lagtemp = ( tempthismonth - templastmonth ) * ( 1.0 - lag ) + templastmonth
          
      ! Adjust amplitude of lagged air temp to give estimated soil temp
      dtemp_soil(lu,jpngr) = avetemp + amp * ( lagtemp - avetemp )

    end do

    ! save temperature for next year
    if (day==ndayyear) then
      dtemp_pvy(:,jpngr) = dtemp(:)
      wscal_pvy(:,:,jpngr) = wscal_alldays(:,:)
    end if

    return

  end subroutine soiltemp


  subroutine initio_soiltemp()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use _params_siml, only: runname, loutdtemp_soil

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(runname)

    ! soil temperature
    if (loutdtemp_soil) then
      filnam=trim(prefix)//'.d.soiltemp.out'
      open(109,file=filnam,err=999,status='unknown')
    end if

    return

    999  stop 'INITIO_soiltemp: error opening output files'

  end subroutine initio_soiltemp


  subroutine initoutput_soiltemp()
    !////////////////////////////////////////////////////////////////
    !  Initialises soiltemp-specific output variables
    !----------------------------------------------------------------
    use _params_core, only: ndayyear
    use _params_siml, only: runname, loutdtemp_soil

    if (loutdtemp_soil) allocate( outdtemp_soil(nlu,ndayyear,maxgrid)  )

  end subroutine initoutput_soiltemp


  subroutine getout_daily_soiltemp( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    use _params_siml, only: loutdtemp_soil

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy    
    integer, intent(in) :: doy    

    if (loutdtemp_soil) outdtemp_soil(:,doy,jpngr) = dtemp_soil(:,jpngr)

  end subroutine getout_daily_soiltemp


  subroutine writeout_ascii_soiltemp( year )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE soiltemp-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use _params_core, only: ndayyear
    use _params_siml, only: spinup, daily_out_startyr, daily_out_endyr, outyear, loutdtemp_soil

    ! arguments
    integer, intent(in) :: year       ! simulation year

    ! Local variables
    real :: itime
    integer :: doy, jpngr
    
    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii_soiltemp: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    !-------------------------------------------------------------------------
    if ( .not. spinup .and. outyear>=daily_out_startyr .and. outyear<=daily_out_endyr ) then

      ! Write daily output only during transient simulation
      do doy=1,ndayyear

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(outyear) + real(doy-1)/real(ndayyear)

        if (nlu>1) stop 'writeout_ascii_soiltemp: write out lu-area weighted sum'

        if (loutdtemp_soil) write(109,999) itime, sum(outdtemp_soil(:,doy,jpngr))

      end do
    end if

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_soiltemp


end module _soiltemp
