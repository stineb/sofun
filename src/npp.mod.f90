module md_npp
  !////////////////////////////////////////////////////////////////
  ! NPP_LPJ MODULE
  ! Contains the "main" subroutine 'npp' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'npp' must contain this list 
  ! of subroutines (names that way).
  !   - npp
  !   - initio_npp
  !   - initoutput_npp
  !   - getout_daily_npp
  !   - getout_monthly_npp
  !   - writeout_ascii_npp
  ! Required module-independent model state variables (necessarily 
  ! updated by 'waterbal') are:
  !   - daily NPP ('dnpp')
  !   - soil temperature ('xxx')
  !   - inorganic N _pools ('no3', 'nh4')
  !   - xxx 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_classdefs
  use md_plant
  use md_params_core, only: npft, maxgrid

  implicit none

  private
  public npp, calc_cexu, calc_resp_maint, initoutput_npp, &
    initio_npp, getout_daily_npp, writeout_ascii_npp

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, allocatable, dimension(:,:,:) :: outdrleaf
  real, allocatable, dimension(:,:,:) :: outdrroot
  real, allocatable, dimension(:,:,:) :: outdrgrow

  ! annual
  real, dimension(npft,maxgrid) :: outarleaf
  real, dimension(npft,maxgrid) :: outarroot
  real, dimension(npft,maxgrid) :: outargrow


contains

  subroutine npp( plant, plant_fluxes, dtemp ) !jpngr, dtemp, doy )
    !/////////////////////////////////////////////////////////////////////////
    ! NET PRIMARY PRODUCTIVITY
    ! Calculate maintenance and growth respiration and substract this from GPP 
    ! to get NPP before additional root respiration for nutrient uptake (see 
    ! SR nuptake). NPP is defined so as to include C allocated to growth as 
    ! well as to exudation. Thus, exudation is not part of autotrophic respir-
    ! ation, but diverted to the quickly decaying exudates pool. Exudates decay
    ! is calculated in SR 'littersom' and is kept track of as soil respiration 
    ! ('rsoil'). This implies that growth respiration is "paid" also on exu-
    ! dates. 
    !-------------------------------------------------------------------------
    use md_params_core, only: npft, ndayyear
    use md_soiltemp, only: dtemp_soil
    use md_gpp, only: drd
    use md_turnover, only: turnover_leaf, turnover_root, turnover_labl
    use md_phenology, only: sprout
    use md_plant, only: plant_type, plant_fluxes_type
    use md_interface

    ! arguments
    type( plant_type ), dimension(npft), intent(in) :: plant ! npft counts over PFTs in all land units (tiles)
    type( plant_fluxes_type ), dimension(npft), intent(inout) :: plant_fluxes
    real, intent(in) :: dtemp      ! air temperature at this day

    ! integer, intent(in) :: jpngr
    ! integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu

    real, parameter :: dleaf_die = 0.012
    real, parameter :: droot_die = 0.012
    real, parameter :: dlabl_die = 0.0

    logical, save :: check_sprout = .false.

    ! print*, '---- in npp:'

    !-------------------------------------------------------------------------
    ! PFT LOOP
    !-------------------------------------------------------------------------
    do pft=1,npft

      if (plant(pft)%plabl%c%c12<0.0) stop 'before npp labile C is neg.'
      if (plant(pft)%plabl%n%n14<0.0) stop 'before npp labile N is neg.'

      lu = params_pft_plant(pft)%lu_category
      
      !/////////////////////////////////////////////////////////////////////////
      ! MAINTENANCE RESPIRATION
      ! use function 'resp_main'
      !-------------------------------------------------------------------------
      ! fine roots should have a higher repsiration coefficient than other tissues (Franklin et al., 2007).
      plant_fluxes(pft)%drleaf = drd(pft)  ! leaf respiration is given by dark respiration as calculated in P-model.       
      plant_fluxes(pft)%drroot = calc_resp_maint( plant(pft)%proot%c%c12 * nind(pft,jpngr), params_plant%r_root, dtemp )
      if (params_pft_plant(pft)%tree) then
        drsapw(pft) = calc_resp_maint( psapw(pft,jpngr)%c%c12 * nind(pft,jpngr), params_plant%r_sapw, dtemp )
      endif


      !/////////////////////////////////////////////////////////////////////////
      ! DAILY NPP AND C EXPORT
      ! NPP is the sum of C available for growth and for N uptake 
      ! This is where isotopic signatures are introduced because only 'dbminc'
      ! is diverted to a pool and re-emission to atmosphere gets delayed. Auto-
      ! trophic respiration is immediate, it makes thus no sense to calculate 
      ! full isotopic effects of gross exchange _fluxes.
      ! Growth respiration ('drgrow') is deduced from 'dnpp' in allocation SR.
      !-------------------------------------------------------------------------
      plant_fluxes(pft)%dnpp = carbon( plant_fuxes(pft)%dgpp - plant_fluxes(pft)%drleaf - plant_fluxes(pft)%drroot )
      plant_fluxes(pft)%dcex = calc_cexu( plant(pft)%proot%c%c12, dtemp )   


      !/////////////////////////////////////////////////////////////////////////
      ! SAFETY AND DEATH
      ! If negative C balance results from GPP - Rleaf - Rroot - Cex then ...
      ! ... first, change allocation to 100% leaves
      ! ... second, when this still leads to a complete depletion of the labile
      !     pool (negative values), shut down organism (zero GPP, NPP, etc., 
      !     but continuing turnover).
      !-------------------------------------------------------------------------
      ! This option (deactivate_root) leads to good results, the alternative leads to on-off growth. Unclear why.
      if ( (plant(pft)%plabl%c%c12 + plant_fluxes(pft)%dnpp%c12 - plant_fluxes(pft)%dcex) < 0.0 ) then
        call deactivate_root( plant_fuxes(pft)%dgpp, plant_fluxes(pft)%drleaf, plant(pft)%plabl%c%c12, plant(pft)%proot, plant_fluxes(pft)%drroot, plant_fluxes(pft)%dnpp%c12, plant_fluxes(pft)%dcex, dtemp, plitt_bg(pft,jpngr) )
      end if

      ! ! -------------------------------------------------------------------------
      ! ! the alternative formulation with shutting all fluxes down and decaying
      ! ! -------------------------------------------------------------------------
      ! if ( (plant(pft)%plabl%c%c12 + plant_fluxes(pft)%dnpp%c12 - plant_fluxes(pft)%dcex) < 0.0 ) then
      !   ! stop exuding
      !   plant_fluxes(pft)%dcex = 0.0

      !   if ( ( plant(pft)%plabl%c%c12 + plant_fluxes(pft)%dnpp%c12 ) < 0.0 ) then

      !     ! ! after C balance has become negative wait until it gets positive again to trigger sprouting
      !     ! ! print*,'setting check_sprout = T ', doy
      !     ! check_sprout = .true.

      !     ! slow death
      !     ! print*,'slow death', doy
      !     plant_fuxes(pft)%dgpp   = 0.0
      !     plant_fluxes(pft)%drleaf = 0.0
      !     plant_fluxes(pft)%drroot = 0.0
      !     drd(pft)    = 0.0
      !     plant_fluxes(pft)%dcex   = 0.0
      !     plant_fluxes(pft)%dnpp   = carbon(0.0)

      !     call turnover_leaf( dleaf_die, pft, jpngr )
      !     call turnover_root( droot_die, pft, jpngr )
      !     ! call turnover_labl( dlabl_die, pft, jpngr )

      !   end if

      ! else
      !   ! normal growth

      !   ! ! trigger sprouting now that C balance is positive again
      !   ! if (check_sprout) then
      !   !   ! print*,'sprouting next day'
      !   !   sprout(doy+1,pft) = .true.
      !   ! end if
      !   ! check_sprout = .false.

      !   ! ! print*,'normal growth', doy
      !   ! if ( .not. interface%steering%dofree_alloc ) frac_leaf(pft) = 0.5
      
      ! end if


      !/////////////////////////////////////////////////////////////////////////
      ! TO LABILE POOL
      ! NPP available for growth first enters the labile pool ('plabl ').
      ! XXX Allocation is called here without "paying"  growth respir.?
      !-------------------------------------------------------------------------
      call ccp( carbon( plant_fluxes(pft)%dcex ), plant(pft)%pexud )
      call ccp( cminus( plant_fluxes(pft)%dnpp, carbon(plant_fluxes(pft)%dcex) ), plant(pft)%plabl%c )

      if (plant(pft)%plabl%c%c12< -1.0e-13) stop 'after npp labile C is neg.'
      if (plant(pft)%plabl%n%n14< -1.0e-13) stop 'after npp labile N is neg.'

      ! print*,'gpp, dclabl', doy, plant_fuxes(pft)%dgpp, cminus( plant_fluxes(pft)%dnpp, carbon(plant_fluxes(pft)%dcex) )

    end do

    ! print*, '---- finished npp'

  end subroutine npp


  subroutine deactivate_root( mygpp, mydrleaf, myplabl, myproot, rroot, npp, cexu, dtemp, myplitt )
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates amount of root mass supportable by (GPP-Rd+Clabl='avl'), so that
    ! NPP is zero and doesn't get negative. Moves excess from pool 'myproot' to
    ! pool 'myplitt'.
    !-------------------------------------------------------------------------
    ! argument
    real, intent(in) :: mygpp
    real, intent(in) :: mydrleaf
    real, intent(in) :: myplabl
    type( orgpool ), intent(inout) :: myproot
    real, intent(out) :: rroot
    real, intent(out) :: npp
    real, intent(out) :: cexu
    real, intent(in) :: dtemp
    type( orgpool ), intent(inout), optional :: myplitt
    
    ! local variables
    real :: croot_trgt
    real :: droot
    type( orgpool ) :: rm_turn

    real, parameter :: safety = 0.9999

    ! calculate target root mass
    croot_trgt = safety * ( mygpp - mydrleaf + myplabl) / ( params_plant%r_root + params_plant%exurate )
    droot      = ( 1.0 - croot_trgt / myproot%c%c12 )
    rm_turn    = orgfrac( droot, myproot )
    if (present(myplitt)) then
      call orgmv( rm_turn, myproot, myplitt )
    else
      myproot = orgminus( myproot, rm_turn )
    end if

    ! update fluxes based on corrected root mass
    rroot = calc_resp_maint( myproot%c%c12, params_plant%r_root, dtemp )
    npp   = mygpp - mydrleaf - rroot
    cexu  = calc_cexu( myproot%c%c12 , dtemp )     

  end subroutine deactivate_root


  function calc_resp_maint( cmass, rresp, dtemp ) result( resp_maint )
    !////////////////////////////////////////////////////////////////
    ! Returns maintenance respiration
    !----------------------------------------------------------------
    use md_rates, only: ftemp
    use md_gpp, only: ramp_gpp_lotemp     ! same ramp as for GPP 

    ! arguments
    real, intent(in) :: cmass   ! N mass per unit area [gN/m2]
    real, intent(in) :: rresp   ! respiration coefficient [gC gC-1 d-1]
    real, intent(in) :: dtemp   ! temperature (soil or air, deg C)

    ! function return variable
    real :: resp_maint                    ! return value: maintenance respiration [gC/m2]

    resp_maint = cmass * rresp * ramp_gpp_lotemp( dtemp )

    ! LPX-like temperature dependeneo of respiration rates
    ! resp_maint = cmass * rresp * ftemp( dtemp, "lloyd_and_taylor" ) * ramp_gpp_lotemp( dtemp )

  end function calc_resp_maint


  function calc_cexu( croot, dtemp ) result( cexu )
    !/////////////////////////////////////////////////////////////////
    ! Constant exudation rate
    !-----------------------------------------------------------------
    use md_gpp, only: ramp_gpp_lotemp     ! same ramp as for GPP 

    ! arguments
    real, intent(in)           :: croot
    real, intent(in), optional :: dtemp   ! temperature (soil or air, deg C)

    ! function return variable
    real :: cexu

    ! low-temperature ramp is included here to prevent negative C balance after exudation
    cexu = params_plant%exurate * croot * ramp_gpp_lotemp( dtemp )

  end function calc_cexu


  subroutine initoutput_npp()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_params_core, only: npft, ndayyear, maxgrid
    use md_interface

    if (interface%params_siml%loutnpp) then

      if (interface%steering%init) then
        allocate( outdrleaf(npft,ndayyear,maxgrid) )
        allocate( outdrroot(npft,ndayyear,maxgrid) )
        allocate( outdrgrow(npft,ndayyear,maxgrid) )
      end if

      outdrleaf(:,:,:) = 0.0
      outdrroot(:,:,:) = 0.0
      outdrgrow(:,:,:) = 0.0

      ! annual output variables
      outarleaf(:,:) = 0.0
      outarroot(:,:) = 0.0
      outargrow(:,:) = 0.0

    end if

  end subroutine initoutput_npp


  subroutine initio_npp()
    !////////////////////////////////////////////////////////////////
    ! Opens input/output files.
    !----------------------------------------------------------------
    use md_interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !////////////////////////////////////////////////////////////////
    ! DAILY OUTPUT: OPEN ASCII OUTPUT FILES 
    !----------------------------------------------------------------
    if (interface%params_siml%loutnpp) then 

      ! LEAF RESPIRATION
      filnam=trim(prefix)//'.d.rleaf.out'
      open(450,file=filnam,err=999,status='unknown')

      ! ROOT RESPIRATION
      filnam=trim(prefix)//'.d.rroot.out'
      open(451,file=filnam,err=999,status='unknown')

      ! GROWTH RESPIRATION
      filnam=trim(prefix)//'.d.rgrow.out'
      open(454,file=filnam,err=999,status='unknown')


      !////////////////////////////////////////////////////////////////
      ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
      !----------------------------------------------------------------

      ! LEAF RESPIRATION
      filnam=trim(prefix)//'.a.rleaf.out'
      open(452,file=filnam,err=999,status='unknown')

      ! ROOT RESPIRATION
      filnam=trim(prefix)//'.a.rroot.out'
      open(453,file=filnam,err=999,status='unknown')

      ! GRWOTH RESPIRATION
      filnam=trim(prefix)//'.a.rgrow.out'
      open(455,file=filnam,err=999,status='unknown')

    end if

    return

    999  stop 'INITIO: error opening output files'

  end subroutine initio_npp


  subroutine getout_daily_npp( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    if (interface%params_siml%loutnpp) then
      !----------------------------------------------------------------
      ! DAILY
      ! Collect daily output variables
      ! so far not implemented for isotopes
      !----------------------------------------------------------------
      outdrleaf(:,doy,jpngr) = drleaf(:)
      outdrroot(:,doy,jpngr) = drroot(:)
      outdrgrow(:,doy,jpngr) = drgrow(:)

      !----------------------------------------------------------------
      ! ANNUAL SUM OVER DAILY VALUES
      ! Collect annual output variables as sum of daily values
      !----------------------------------------------------------------
      outarleaf(:,jpngr) = outarleaf(:,jpngr) + drleaf(:)
      outarroot(:,jpngr) = outarroot(:,jpngr) + drroot(:)
      outargrow(:,jpngr) = outargrow(:,jpngr) + drgrow(:)

    end if

  end subroutine getout_daily_npp


  subroutine writeout_ascii_npp( year )
    !/////////////////////////////////////////////////////////////////////////
    ! Write daily ASCII output
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !-------------------------------------------------------------------------
    use md_params_core, only: ndayyear, nlu
    use md_interface

    ! arguments
    integer, intent(in) :: year       ! simulation year

    ! local variables
    real :: itime
    integer :: day, moy, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! Collect variables to output variables
    !-------------------------------------------------------------------------
    if (nlu>1) stop 'Output only for one LU category implemented.'

    if (interface%params_siml%loutnpp) then
      !-------------------------------------------------------------------------
      ! DAILY OUTPUT
      ! Write daily value, summed over all PFTs / LUs
      ! xxx implement taking sum over PFTs (and gridcells) in this land use category
      !-------------------------------------------------------------------------
      if ( .not. interface%steering%spinup &
        .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
        .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

        ! Write daily output only during transient simulation
        do day=1,ndayyear

          ! Define 'itime' as a decimal number corresponding to day in the year + year
          itime = real(interface%steering%outyear) + real(day-1)/real(ndayyear)
          
          write(450,999) itime, sum(outdrleaf(:,day,jpngr))
          write(451,999) itime, sum(outdrroot(:,day,jpngr))
          write(454,999) itime, sum(outdrgrow(:,day,jpngr))

        end do

        !-------------------------------------------------------------------------
        ! ANNUAL OUTPUT
        ! Write annual value, summed over all PFTs / LUs
        ! xxx implement taking sum over PFTs (and gridcells) in this land use category
        !-------------------------------------------------------------------------
        itime = real(interface%steering%outyear)

        write(452,999) itime, outarleaf(:,jpngr)
        write(453,999) itime, outarroot(:,jpngr)
        write(455,999) itime, outargrow(:,jpngr)

      end if

    end if

    return

    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_npp


end module md_npp
