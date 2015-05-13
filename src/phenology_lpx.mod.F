module _phenology
  !////////////////////////////////////////////////////////////////
  ! TEMPERATURE-DRIVEN PHENOLOGY 
  ! Adopted from LPX-Bern
  ! Contains the "main" subroutine 'summerphenology' and all 
  ! necessary subroutines for handling input/output. 
  ! Every module that implements 'summerphenology' must contain 
  ! this list of subroutines (names that way).
  !   - summerphenology
  !   - getpar_modl_summerphenology
  !   - initio_summerphenology, (optional)
  !   - initoutput_summerphenology, (optional)
  !   - getout_daily_summerphenology, (optional)
  !   - getout_monthly_summerphenology, (optional)
  !   - writeout_ascii_summerphenology, (optional)
  ! Required module-independent model state variables (necessarily 
  ! updated by 'waterbal') are:
  !   - (none)
  !----------------------------------------------------------------
  use _params_core, only: npft, ndayyear
  
  implicit none

  ! PHENOLOGY PARAMETERS
  real, dimension(npft) :: GDDBASE ! GDD base, for PFT11-14, a T0 is chosen to be 0deg C (Prentice et al. 1992, J.o.Biogeography), pftpar(pft,33) in LPX
  real, dimension(npft) :: RAMP    ! summergreen phenology ramp, GDD requirement to grow full leaf canopy

  logical, dimension(npft) :: SUMMERGREEN, EVERGREEN, RAINGREEN  ! phenology type

  ! MODULE-SPECIFIC VARIABLES
  real, dimension(ndayyear,npft) :: dtphen   ! daily temperature-driven phenology (=dphen_t in LPX)


contains

  subroutine summerphenology( jpngr, dtemp )
    !//////////////////////////////////////////////////////////
    ! Defines dtphen, the temperature-driven phenology
    !----------------------------------------------------------
    use _params_core, only: ndayyear, maxgrid, nmonth, middaymonth
    use _params_modl, only: tree
    use _sofunutils, only: daily2monthly, monthly2daily

    ! arguments
    integer, intent(in) :: jpngr
    real, dimension(ndayyear), intent(in) :: dtemp

    ! local variables
    integer :: warmest, coldest, month, midsummer, firstday, d, pft, day
    real    :: leafon_n, aphen, gdd
    real, dimension(nmonth)         :: mtemp       ! monthly temperature as a mean of daily values in resp. month
    real, dimension(nmonth,maxgrid) :: mtemp_pvy   ! monthly temperature as a mean of daily values in resp. month, previous year
    real, dimension(ndayyear)       :: dtemp_int   ! daily temperature as linearly interpolated from monthly temperature
    logical, save :: firstcall = .true.


    ! initialise
    dtphen(:,:) = 0.0

    ! Phenology is driven by monthly temperatures and daily temperatures
    ! as interpolated from monthly temperatures to remove day-to-day
    ! variability
    mtemp = daily2monthly( dtemp, "mean" )
    if (firstcall) then
      mtemp_pvy(:,jpngr) = mtemp(:)
      firstcall = .false.
    end if
    dtemp_int = monthly2daily( mtemp, "interpol", .false., mtemp_pvy )

    ! First find warmest and coldest month and mid-summer day
    warmest=1
    do month=1,nmonth
      if (mtemp(month)>mtemp(warmest)) warmest=month
    enddo
    coldest=1
    do month=1,nmonth
      if (mtemp(month)<mtemp(coldest)) coldest=month
    enddo
    midsummer = middaymonth( warmest )

    do pft=1,npft
      !----------------------------------------------------------
      ! Find day of leaf abscission ('firstday') at end of summer
      ! i.e. when daily temperature falls below GDDBASE.
      !----------------------------------------------------------
      firstday=midsummer+1
      do while (dtemp_int(firstday)>=GDDBASE(pft) .and. firstday/=midsummer)
        firstday=firstday+1
        if (firstday>ndayyear) firstday=1
      enddo
      
      write(0,*) 'summergreen', summergreen

      if (summergreen(pft)) then
        !----------------------------------------------------------
        ! SUMMERGREEN TAXA
        !----------------------------------------------------------
        if (firstday==midsummer) then 
          dtphen(:,pft)=1.0     ! no leaf abscission
        else
          gdd=0.0               ! accumulated growing degree days
          day=firstday+1
          if (day>ndayyear) day=1
          do while (day/=firstday)
            if (dtemp_int(day)>GDDBASE(pft)) then ! growing day
              gdd = gdd + dtemp_int(day) - GDDBASE(pft)
              if (RAMP(pft)>0.0) then
                dtphen(day,pft) = min( gdd / RAMP(pft), 1.0 )
              else
                dtphen(day,pft) = 1.0
              endif
            endif
            day=day+1
            if (day>ndayyear) day=1
            write(0,*) 'SUMMERPHENOLOGY: day, dtphen', day, dtphen(day,pft)
            stop 'in summerphenology'
          enddo
        endif
        
        if (tree(pft)) then
          !----------------------------------------------------------
          ! TREES
          !----------------------------------------------------------
          aphen=sum(dtphen(:,pft))
          if (aphen>210) then 
            do d=middaymonth(coldest),middaymonth(coldest)+75
              if (d<=ndayyear) then
                day=d
              else
                day=d-ndayyear      
              endif
              dtphen(day,pft)=0.0
            enddo
            do d=middaymonth(coldest)-75,middaymonth(coldest)
              if (d>=1) then
                day=d
              else
                day=ndayyear+d
              endif
              dtphen(day,pft)=0.0
            enddo
          endif
        endif

      else
        !----------------------------------------------------------
        ! NON-SUMMERGREEN TAXA
        !----------------------------------------------------------
        dtphen(:,pft)=1.0
      endif
      
    enddo                     !pft

    ! save monthly temperature for next year
    mtemp_pvy(:,jpngr) = mtemp(:)
      
    return

  end subroutine summerphenology


  subroutine getpar_pft_summerphenology( )
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    ! local variables
    integer, parameter    :: npar = 2
    integer               :: pft
    real, dimension(npar) :: params_array
    character(len=50)     :: paramfilnam
    
    paramfilnam = 'params_pft_summerphenology.dat'

    open(unit=07,file=trim(paramfilnam),status='OLD')      
    read(07,*) params_array
    close(07)

    ! xxx implement array readin (for each pft, then loop over pfts)
    pft = 1

    ! growing degree days base (usually 5 deg C)
    GDDBASE(pft) = params_array(1)

    ! ramp slope for phenology (1 for grasses: immediate phenology turning on)
    RAMP(pft)    = params_array(2)
    
    ! phenology type
    EVERGREEN(:)   = .false.
    SUMMERGREEN(:) = .false.
    RAINGREEN(:)   = .false.
    if (params_array(3)==1.0) EVERGREEN(pft)   = .true.
    if (params_array(3)==2.0) SUMMERGREEN(pft) = .true.
    if (params_array(3)==3.0) RAINGREEN(pft)   = .true.

    return

  end subroutine getpar_pft_summerphenology

end module _phenology





