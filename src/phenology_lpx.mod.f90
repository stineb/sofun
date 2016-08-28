module md_phenology
  !////////////////////////////////////////////////////////////////
  ! TEMPERATURE-DRIVEN PHENOLOGY 
  ! Adopted from LPX-Bern
  ! Contains the "main" subroutine 'gettempphenology and phenology' and all 
  ! necessary subroutines for handling input/output. 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: npft, ndayyear
  
  implicit none

  private
  public phenology_type, gettempphenology, getpar_modl_phenology

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  type phenology_type
    real,    dimension(ndayyear) :: dtphen       ! daily temperature-driven phenology (=dphen_t in LPX)
    logical, dimension(ndayyear) :: sprout       ! boolean whether PFT is present
    logical, dimension(ndayyear) :: shedleaves   ! boolean whether PFT is present
  end type phenology_type

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

  function gettempphenology( jpngr, dtemp ) result( pheno )
    !//////////////////////////////////////////////////////////
    ! Defines dtphen, the temperature-driven phenology
    !----------------------------------------------------------
    use md_params_core, only: ndayyear, maxgrid, nmonth, middaymonth
    use md_plant, only: params_pft_plant
    use md_sofunutils, only: daily2monthly, monthly2daily

    ! arguments
    integer, intent(in) :: jpngr
    real, dimension(ndayyear), intent(in) :: dtemp

    ! function return variable
    type( phenology_type ), dimension(npft) :: pheno

    ! local variables
    integer :: warmest, coldest, month, midsummer, firstday, d, pft, day
    real    :: leafon_n, aphen, gdd
    real, dimension(nmonth)         :: mtemp       ! monthly temperature as a mean of daily values in resp. month
    real, dimension(nmonth,maxgrid) :: mtemp_pvy   ! monthly temperature as a mean of daily values in resp. month, previous year
    real, dimension(ndayyear)       :: dtemp_int   ! daily temperature as linearly interpolated from monthly temperature
    logical, save :: firstcall = .true.

    do pft=1,npft

      pheno(pft)%dtphen(:)     = 0.0
      pheno(pft)%sprout(:)     = .false.
      pheno(pft)%shedleaves(:) = .false.

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

      if (npft>1) stop 'in phenology: think of something nice'
      pft = 1

      !----------------------------------------------------------
      ! Find day of leaf abscission ('firstday') at end of summer
      ! i.e. when daily temperature falls below gddbase.
      !----------------------------------------------------------
      firstday=midsummer+1
      do while (dtemp_int(firstday)>=params_pheno%gddbase .and. firstday/=midsummer)
        firstday=firstday+1
        if (firstday>ndayyear) firstday=1
      enddo
      
      if (params_pft_pheno(pft)%summergreen) then
        !----------------------------------------------------------
        ! summergreen TAXA
        !----------------------------------------------------------
        if (firstday==midsummer) then 
          pheno(pft)%dtphen(:,pft)=1.0     ! no leaf abscission
        else
          gdd=0.0               ! accumulated growing degree days
          day=firstday+1
          if (day>ndayyear) day=1
          do while (day/=firstday)
            if (dtemp_int(day)>params_pheno%gddbase) then ! growing day
              gdd = gdd + dtemp_int(day) - params_pheno%gddbase
              if (params_pft_pheno(pft)%ramp>0.0) then
                pheno(pft)%dtphen(day) = min( gdd / params_pft_pheno(pft)%ramp, 1.0 )
              else
                pheno(pft)%dtphen(day) = 1.0
              endif
            endif
            day=day+1
            if (day>ndayyear) day=1
          enddo
        endif
        
      else
        !----------------------------------------------------------
        ! NON-summergreen TAXA
        !----------------------------------------------------------
        pheno(pft)%dtphen(:,pft)=1.0

      endif
        
      ! save monthly temperature for next year
      mtemp_pvy(:,jpngr) = mtemp(:)

      ! xxx try: really weird: when appplying a loop over pft, pheno(pft)%dtphen, sprout, 
      ! and shedleaves are all set to false after finishing each iteration
      ! therefore set to pft=1 here.
      if (npft>1) stop 'in phenology: think of something nice'

      do day=2,ndayyear

        if (params_pft_pheno(pft)%summergreen) then
          !----------------------------------------------------------
          ! temperature-driven phenology summergreen
          !----------------------------------------------------------

          if ( pheno(pft)%dtphen(day) > 0.0 .and. pheno(pft)%dtphen(day-1) == 0.0 ) then
            !----------------------------------------------------------
            ! beginning of season (spring)
            !----------------------------------------------------------
            pheno(pft)%sprout(day) = .true.
            pheno(pft)%shedleaves(day) = .false.
            ! print*, 'sprouting on day ', day 
            ! print*, sprout(38,pft)

          else if ( pheno(pft)%dtphen(day) > 0.0 ) then
            !----------------------------------------------------------
            ! during season (after spring and before autumn)
            !----------------------------------------------------------
            pheno(pft)%sprout(day) = .false.
            pheno(pft)%shedleaves(day) = .false.
            ! print*, 'active on day ', day

          else if ( pheno(pft)%dtphen(day) == 0.0 .and. pheno(pft)%dtphen(day-1) > 0.0 ) then
            !----------------------------------------------------------
            ! end of season (autumn)
            !----------------------------------------------------------
            pheno(pft)%sprout(day) = .false.
            pheno(pft)%shedleaves(day) = .true.
            ! print*, 'shedding leaves on day ', day 
            ! print*, shedleaves(345,pft)

          else if ( pheno(pft)%dtphen(day) == 0.0 ) then
            !----------------------------------------------------------
            ! during dormant season (after autumn and before spring)
            !----------------------------------------------------------
            pheno(pft)%sprout(day) = .false.
            pheno(pft)%shedleaves(day) = .false.
            ! print*, 'dormant on day ', day

          end if

        else

          stop 'estab_daily not implemented for trees'

        end if

      end do

      ! xxx debug
      ! print*,'PHENOLOGY: overriding shedleaves'
      pheno(pft)%shedleaves(:) = .false.

    end do

  end function gettempphenology


  subroutine getpar_modl_phenology()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal
    use md_plant, only: params_pft_plant

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

  end subroutine getpar_modl_phenology


end module md_phenology





