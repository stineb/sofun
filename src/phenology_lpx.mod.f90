module _phenology
  !////////////////////////////////////////////////////////////////
  ! TEMPERATURE-DRIVEN PHENOLOGY 
  ! Adopted from LPX-Bern
  ! Contains the "main" subroutine 'gettempphenology and phenology' and all 
  ! necessary subroutines for handling input/output. 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: npft, ndayyear
  
  implicit none

  private
  public gettempphenology, sprout, shedleaves, params_pft_pheno

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  real, dimension(ndayyear,npft)    :: dtphen       ! daily temperature-driven phenology (=dphen_t in LPX)
  logical, dimension(ndayyear,npft) :: sprout       ! boolean whether PFT is present
  logical, dimension(ndayyear,npft) :: shedleaves   ! boolean whether PFT is present


  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------


  !----------------------------------------------------------------
  ! Parameters
  !----------------------------------------------------------------
  type paramstype_pheno
    real :: gddbase ! GDD base, for PFT11-14, a T0 is chosen to be 0deg C (Prentice et al. 1992, J.o.Biogeography), pftpar(pft,33) in LPX
  end type paramstype_pheno

  type( paramstype_pheno ) :: params_pheno

  type pftparamstype_pheno
    real    :: ramp    ! summergreen phenology ramp, GDD requirement to grow full leaf canopy
    logical :: evergreen
    logical :: summergreen
    logical :: raingreen  
  end type pftparamstype_pheno

  type( pftparamstype_pheno ), dimension(npft) :: params_pft_pheno


contains

  subroutine gettempphenology( jpngr, dtemp )
    !//////////////////////////////////////////////////////////
    ! Defines dtphen, the temperature-driven phenology
    !----------------------------------------------------------
    use _params_core, only: ndayyear, maxgrid, nmonth, middaymonth
    use _plant, only: params_pft_plant
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
      ! i.e. when daily temperature falls below gddbase.
      !----------------------------------------------------------
      firstday=midsummer+1
      do while (dtemp_int(firstday)>=params_pheno%gddbase .and. firstday/=midsummer)
        firstday=firstday+1
        if (firstday>ndayyear) firstday=1
      enddo
      
      ! write(0,*) 'dtemp_int'
      ! write(0,*) dtemp_int
      ! write(0,*) 'midsummer', midsummer
      ! write(0,*) 'firstday', firstday
      ! write(0,*) 'summergreen', params_pft_pheno%summergreen
      ! write(0,*) 'params_pheno%gddbase', params_pheno%gddbase

      if (params_pft_pheno(pft)%summergreen) then
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
            if (dtemp_int(day)>params_pheno%gddbase) then ! growing day
              gdd = gdd + dtemp_int(day) - params_pheno%gddbase
              if (params_pft_pheno%ramp(pft)>0.0) then
                dtphen(day,pft) = min( gdd / params_pft_pheno%ramp(pft), 1.0 )
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
        
        if (params_pft_plant(pft)%tree) then
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

        if (params_pft_pheno(pft)%summergreen) then
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
    
    return

  end subroutine gettempphenology


  subroutine getpar_phenology()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use _sofunutils, only: getparreal
    use _plant, only: params_pft_plant

    ! local variables
    real        :: phentype
    integer     :: pft

    ! growing degree days base (usually 5 deg C)
    params_pheno%gddbase = getparreal( 'params/params_phenology.dat', 'gddbase' )

    do pft=1,npft

      ! ramp slope for phenology (1 for grasses: immediate phenology turning on)
      params_pft_pheno(pft)%ramp = getparreal( 'params/params_phenology.dat', 'ramp_pft_'//params_pft_plant(pft)%pftname )

      ! phenology type
      phentype = getparreal( 'params/params_phenology.dat', 'phentype_pft_'//params_pft_plant(pft)%pftname )

      if (phentype==1.0) params_pft_pheno(pft)%evergreen   = .true.
      if (phentype==2.0) params_pft_pheno(pft)%summergreen = .true.
      if (phentype==3.0) params_pft_pheno(pft)%raingreen   = .true.

    end do
 
    return
 
    999  format (I2.2)

  end subroutine getpar_phenology


end module _phenology





