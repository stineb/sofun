module _phenology
  !////////////////////////////////////////////////////////////////
  ! TEMPERATURE-DRIVEN PHENOLOGY 
  ! Adopted from LPX-Bern
  ! Contains the "main" subroutine 'gettempphenology and phenology' and all 
  ! necessary subroutines for handling input/output. 
  ! Every module that implements 'tempphenology' must contain 
  ! this list of subroutines (names that way).
  !   - tempphenology
  !   - getpar_modl_tempphenology
  ! Required module-independent model state variables (necessarily 
  ! updated by 'waterbal') are:
  !   - (none)
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: npft, ndayyear
  
  implicit none

  ! PHENOLOGY PARAMETERS
  real                  :: GDDBASE ! GDD base, for PFT11-14, a T0 is chosen to be 0deg C (Prentice et al. 1992, J.o.Biogeography), pftpar(pft,33) in LPX
  real, dimension(npft) :: RAMP    ! summergreen phenology RAMP, GDD requirement to grow full leaf canopy

  logical, dimension(npft) :: summergreen, evergreen, raingreen  ! phenology type

  ! MODULE-SPECIFIC VARIABLES
  real, dimension(ndayyear,npft) :: dtphen   ! daily temperature-driven phenology (=dphen_t in LPX)
  logical, dimension(ndayyear,npft) :: sprout   ! boolean whether PFT is present
  logical, dimension(ndayyear,npft) :: shedleaves   ! boolean whether PFT is present


contains

  subroutine gettempphenology( jpngr, dtemp )
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
    dtphen(:,:)     = 0.0
    ! sprout(:,:)     = .false.
    ! shedleaves(:,:) = .false.

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
      do while (dtemp_int(firstday)>=GDDBASE .and. firstday/=midsummer)
        firstday=firstday+1
        if (firstday>ndayyear) firstday=1
      enddo
      
      ! write(0,*) 'dtemp_int'
      ! write(0,*) dtemp_int
      ! write(0,*) 'midsummer', midsummer
      ! write(0,*) 'firstday', firstday
      ! write(0,*) 'summergreen', summergreen
      ! write(0,*) 'GDDBASE', GDDBASE

      if (summergreen(pft)) then
        !----------------------------------------------------------
        ! summergreen TAXA
        !----------------------------------------------------------
        if (firstday==midsummer) then 
          dtphen(:,pft)=1.0     ! no leaf abscission
        else
          gdd=0.0               ! accumulated growing degree days
          day=firstday+1
          if (day>ndayyear) day=1
          do while (day/=firstday)
            if (dtemp_int(day)>GDDBASE) then ! growing day
              gdd = gdd + dtemp_int(day) - GDDBASE
              if (RAMP(pft)>0.0) then
                dtphen(day,pft) = min( gdd / RAMP(pft), 1.0 )
              else
                dtphen(day,pft) = 1.0
              endif
            endif
            ! write(0,*) 'day, dtphen', day, dtphen(day,pft)
            day=day+1
            if (day>ndayyear) day=1
          enddo
          ! write(0,*) 'gettempphenology: dtphen(day,pft) '
          ! write(0,*) dtphen(:,pft)
          ! stop
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
        ! NON-summergreen TAXA
        !----------------------------------------------------------
        dtphen(:,pft)=1.0
      endif
      
    enddo                     !pft

    ! save monthly temperature for next year
    mtemp_pvy(:,jpngr) = mtemp(:)

    ! do day=1,ndayyear
    !   if (sprout(day,1)) write(0,*) 'sprouting on day',day
    !   if (shedleaves(day,1)) write(0,*) 'shedleavesing on day',day
    ! end do
    ! write(0,*) shedleaves
    ! stop

    ! XXX MOVE THIS TO SEPARATE ROUTINE
    do day=1,ndayyear
      do pft=1,npft

        if (summergreen(pft)) then
          !----------------------------------------------------------
          ! temperature-driven phenology summergreen
          !----------------------------------------------------------

          if ( dtphen(day,pft) > 0.0 .and. dtphen(day-1,pft) == 0.0 ) then
            !----------------------------------------------------------
            ! beginning of season (spring)
            !----------------------------------------------------------
            sprout(day,pft) = .true.
            shedleaves(day,pft) = .false.
            ! write(0,*) 'sprouting on day ', day 

          else if ( dtphen(day,pft) > 0.0 ) then
            !----------------------------------------------------------
            ! during season (after spring and before autumn)
            !----------------------------------------------------------
            sprout(day,pft) = .false.
            shedleaves(day,pft) = .false.

          else if ( dtphen(day,pft) == 0.0 .and. dtphen(day-1,pft) > 0.0 ) then
            !----------------------------------------------------------
            ! end of season (autumn)
            !----------------------------------------------------------
            sprout(day,pft) = .false.
            shedleaves(day,pft) = .true.
            ! write(0,*) 'shedding leaves on day ', day 

          else if ( dtphen(day,pft) == 0.0 ) then
            !----------------------------------------------------------
            ! during dormant season (after autumn and before spring)
            !----------------------------------------------------------
            sprout(day,pft) = .false.
            shedleaves(day,pft) = .false.

          end if

        else

          stop 'estab_daily not implemented for trees'

        end if

      end do
    end do
    ! XXXX

    return

  end subroutine gettempphenology


  subroutine getpar_phenology( )
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use _sofunutils, only: getparreal

    ! local variables
    character*2           :: char_pftno
    real, dimension(npft) :: PHENTYPE
    integer               :: pft

    ! initialise

    ! growing degree days base (usually 5 deg C)
    GDDBASE = getparreal( 'params_phenology.dat', 'GDDBASE' )

    do pft=1,npft

      ! define PFT-extension used for parameter names in parameter file
      write(char_pftno, 999) pft

      ! RAMP slope for phenology (1 for grasses: immediate phenology turning on)
      RAMP(pft) = getparreal( 'params_phenology.dat', trim('RAMP_PFT')//char_pftno )

      ! phenology type
      PHENTYPE(pft) = getparreal( 'params_phenology.dat', trim('PHENTYPE_PFT')//char_pftno )

      if (PHENTYPE(pft)==1.0) evergreen(pft)   = .true.
      if (PHENTYPE(pft)==2.0) summergreen(pft) = .true.
      if (PHENTYPE(pft)==3.0) raingreen(pft)   = .true.

    end do
 
    return
 
    999  format (I2.2)

  end subroutine getpar_phenology


end module _phenology





