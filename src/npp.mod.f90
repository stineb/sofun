module md_npp
  !////////////////////////////////////////////////////////////////
  ! NPP_LPJ MODULE
  ! Contains the "main" subroutine 'npp' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'npp' must contain this list 
  ! of subroutines (names that way).
  !   - npp
  !   - ((interface%steering%init))io_npp
  !   - ((interface%steering%init))output_npp
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
  public npp, calc_cexu, calc_resp_maint, deactivate_root, initoutput_npp, &
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

  subroutine npp( jpngr, dtemp, doy )
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
    use md_gpp, only: dgpp, drd
    use md_turnover, only: turnover_leaf, turnover_root
    use md_phenology, only: sprout

    ! arguments
    integer, intent(in) :: jpngr
    real, intent(in)    :: dtemp      ! air temperature at this day
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu
    
    real :: cbal
    real :: avl

    real, parameter :: dleaf_die = 0.005
    real, parameter :: droot_die = 0.005


    ! print*, '---- in npp:'

    !-------------------------------------------------------------------------
    ! PFT LOOP
    !-------------------------------------------------------------------------
    do pft=1,npft

      ! print*,'BEFORE: clabl    ', plabl(pft,jpngr)%c

      if (plabl(pft,jpngr)%c%c12<0.0) stop 'before npp labile C is neg.'
      if (plabl(pft,jpngr)%n%n14<0.0) stop 'before npp labile N is neg.'

      lu = params_pft_plant(pft)%lu_category
      
      !/////////////////////////////////////////////////////////////////////////
      ! MAINTENANCE RESPIRATION
      ! use function 'resp_main'
      !-------------------------------------------------------------------------
      ! fine roots should have a higher repsiration coefficient than other tissues (Franklin et al., 2007).
      drleaf(pft) = drd(pft)  ! leaf respiration is given by dark respiration as calculated in P-model.       
      drroot(pft) = calc_resp_maint( proot(pft,jpngr)%c%c12 * nind(pft,jpngr), params_plant%r_root, dtemp )
      if (params_pft_plant(pft)%tree) then
        drsapw(pft) = calc_resp_maint( psapw(pft,jpngr)%c%c12 * nind(pft,jpngr), params_plant%r_sapw, dtemp )
      endif


      !!<<<<<<<new:

      ! dnpp(pft) = carbon( dgpp(pft) - drleaf(pft) - drroot(pft) )
      ! dcex(pft) = calc_cexu( proot(pft,jpngr)%c%c12 , dtemp )
      ! if (dgpp(pft)>0.0) print*,'DOY, resp. as fraction of GPP:', doy, (drleaf(pft) + drroot(pft) + dcex(pft)) / dgpp(pft)

      ! cbal = dnpp(pft)%c12 - dcex(pft)

      ! if ( cbal>0.0 ) then
      !   ! positive C balance after respiration and C export => PFT continues growing
      !   ! cleaf + croot = 0.0 after initialisation of PFT in vegdynamics
      !   isgrowing(pft,jpngr) = .true.
      !   ! isdying(pft,jpngr)   = .false.

      ! else
      !   ! no positive C balance after respiration and C export => PFT stops growing
      !   isgrowing(pft,jpngr) = .false.
      !   ! isdying(pft,jpngr)   = .false.
      !   dcex(pft) = 0.0

      !   if ( (dnpp(pft)%c12 + plabl(pft,jpngr)%c%c12) < 0.0 ) then
      !     ! labile pool is depleted
      !     ! print*,'cbal  ', cbal
      !     ! print*,'clabl ', plabl(pft,jpngr)%c%c12
      !     ! isdying(pft,jpngr) = .true.

      !     call turnover_leaf( dleaf_die, pft, jpngr )
      !     call turnover_root( droot_die, pft, jpngr )

      !     dgpp(pft)   = 0.0
      !     drleaf(pft) = 0.0
      !     drroot(pft) = 0.0
      !     drd(pft)    = 0.0
      !     dnpp(pft)   = carbon(0.0)

      !     ! print*,'dcex ', dcex(pft)
      !     ! print*,'dnpp ', dnpp(pft)
      !     ! print*,'clabl', plabl(pft,jpngr)
      !     ! stop 'in npp'

      !   end if

      ! end if

      ! !/////////////////////////////////////////////////////////////////////////
      ! ! C TO/FROM LABILE POOL AND TO EXUDATES POOL
      ! !-------------------------------------------------------------------------
      ! call ccp( carbon( dcex(pft) ), pexud(pft,jpngr) )
      ! call ccp( cminus( dnpp(pft), carbon(dcex(pft)) ), plabl(pft,jpngr)%c )

      ! ! ! If C used for root respiration and export is not available, then reduce 
      ! ! ! root mass to match 
      ! ! if ( avl < 0.0 ) then
      ! !   print*,'resize_plant ...'
      ! !   call resize_plant( dgpp(pft), drleaf(pft), plabl(pft,jpngr)%c%c12, proot(pft,jpngr), pleaf(pft,jpngr), drroot(pft), dnpp(pft)%c12, dcex(pft), dtemp, plitt_af(pft,jpngr), plitt_bg(pft,jpngr) )
      ! !   print*,'... done'
      ! ! end if


      ! if (plabl(pft,jpngr)%c%c12< -1.0e-13) stop 'after npp labile C is neg.'
      ! if (plabl(pft,jpngr)%n%n14< -1.0e-13) stop 'after npp labile N is neg.'

      !!===========   

      !/////////////////////////////////////////////////////////////////////////
      ! DAILY NPP AND C EXPORT
      ! NPP is the sum of C available for growth and for N uptake 
      ! This is where isotopic signatures are introduced because only 'dbminc'
      ! is diverted to a pool and re-emission to atmosphere gets delayed. Auto-
      ! trophic respiration is immediate, it makes thus no sense to calculate 
      ! full isotopic effects of gross exchange _fluxes.
      ! Growth respiration ('drgrow') is deduced from 'dnpp' in allocation SR.
      !-------------------------------------------------------------------------
      dnpp(pft) = carbon( dgpp(pft) - drleaf(pft) - drroot(pft) )
      dcex(pft) = calc_cexu( proot(pft,jpngr)%c%c12 , dtemp )   

      ! if (dgpp(pft)>0.0) then
      !   depletionfrac(pft) = max( 0.0, 1.0 - (drleaf(pft) + drroot(pft) + dcex(pft)) / dgpp(pft) )
      ! else
      !   if (sprout(doy,pft)) then
      !     depletionfrac(pft) = 1.0
      !   else
      !     depletionfrac(pft) = 0.0

      !     ! enhance turnover when plant has negative C balance
      !     call turnover_leaf( params_pft_plant(pft)%k_decay_leaf_base, pft, jpngr )
      !     call turnover_root( params_pft_plant(pft)%k_decay_root, pft, jpngr )
      !     dgpp(pft)   = 0.0
      !     dnpp(pft)   = carbon(0.0)
      !     dcex(pft)   = 0.0
      !     drleaf(pft) = 0.0
      !     drroot(pft) = 0.0
      !     drd(pft)    = 0.0

      !   end if
      ! end if
      ! print*,'doy, clabl, depl. ', doy, plabl(pft,jpngr)%c%c12, depletionfrac(pft)

      ! ! xxx try
      ! if ( (drleaf(pft) + drroot(pft) + dcex(pft)) > dgpp(pft) ) then
        
      !   ! if (lai_ind(pft,jpngr)>1.0)  then
      !   !   isgrowing(pft,jpngr) = .false.
      !   ! else
      !   !   isgrowing(pft,jpngr) = .true.
      !   ! end if

      !   dcex(pft) = 0.0

      !   ! if ( (plabl(pft,jpngr)%c%c12 + dnpp(pft)%c12)<0.0 ) then

      !   !   call turnover_leaf( dleaf_die, pft, jpngr )
      !   !   call turnover_root( droot_die, pft, jpngr )

      !   !   dcex(pft)   = 0.0
      !   !   dgpp(pft)   = 0.0
      !   !   drleaf(pft) = 0.0
      !   !   drroot(pft) = 0.0
      !   !   drd(pft)    = 0.0
      !   !   dnpp(pft)   = carbon(0.0)

      !   !   ! stop 'sharp decline'

      !   ! end if

      ! ! else

      ! !   isgrowing(pft,jpngr) = .true.

      ! end if


      !/////////////////////////////////////////////////////////////////////////
      ! SAFETY AND DEATH
      ! If negative C balance results from GPP - Rleaf - Rroot - Cex then ...
      ! ... first, change allocation to 100% leaves
      ! ... second, when this still leads to a complete depletion of the labile
      !     pool (negative values), shut down organism (zero GPP, NPP, etc., 
      !     but continuing turnover).
      !-------------------------------------------------------------------------
      if ( (plabl(pft,jpngr)%c%c12 + dnpp(pft)%c12 - dcex(pft)) < 0.0 ) then
        ! slow death
        print*,'slow death', doy
        frac_leaf(pft) = 1.0
        dgpp(pft)   = 0.0
        drleaf(pft) = 0.0
        drroot(pft) = 0.0
        drd(pft)    = 0.0
        dcex(pft)   = 0.0
        dnpp(pft)   = carbon(0.0)

        call turnover_leaf( dleaf_die, pft, jpngr )
        call turnover_root( droot_die, pft, jpngr )

        ! print*,'deactivating root'
        ! call deactivate_root( dgpp(pft), drleaf(pft), plabl(pft,jpngr)%c%c12, proot(pft,jpngr), drroot(pft), dnpp(pft)%c12, dcex(pft), dtemp, plitt_bg(pft,jpngr) )
      else if ( dnpp(pft)%c12 - dcex(pft) < 0.0 ) then
        ! negative C balance -> no more allocation to roots (no growth anyways)
        print*,'put all to leaves', doy
        frac_leaf(pft) = 1.0
      else
        ! normal growth
        print*,'normal growth', doy
        frac_leaf(pft) = 0.5
      end if

      !/////////////////////////////////////////////////////////////////////////
      ! TO LABILE POOL
      ! NPP available for growth first enters the labile pool ('plabl ').
      ! XXX Allocation is called here without "paying"  growth respir.?
      !-------------------------------------------------------------------------
      call ccp( carbon( dcex(pft) ), pexud(pft,jpngr) )
      call ccp( cminus( dnpp(pft), carbon(dcex(pft)) ), plabl(pft,jpngr)%c )

      ! print*,'doy ', doy
      ! print*,'npp - cex', dnpp(pft)%c12 - dcex(pft)
      ! print*,'clabl    ', plabl(pft,jpngr)%c

      if (plabl(pft,jpngr)%c%c12< -1.0e-13) stop 'after npp labile C is neg.'
      if (plabl(pft,jpngr)%n%n14< -1.0e-13) stop 'after npp labile N is neg.'

      ! ! xxx try:
      ! if (plabl(pft,jpngr)%c%c12< -1.0e-13) plabl(pft,jpngr)%c%c12 = 0.0
      ! if (plabl(pft,jpngr)%n%n14< -1.0e-13) plabl(pft,jpngr)%n%n14 = 0.0

      !!>>>>>>>>:old

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
    real, intent(in)           :: cmass   ! N mass per unit area [gN/m2]
    real, intent(in)           :: rresp   ! respiration coefficient [gC gC-1 d-1]
    real, intent(in), optional :: dtemp   ! temperature (soil or air, deg C)

    ! function return variable
    real :: resp_maint                    ! return value: maintenance respiration [gC/m2]

    resp_maint = cmass * rresp ! * ramp_gpp_lotemp( dtemp )

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
    use md_interface
    use md_params_core, only: npft, ndayyear, maxgrid

    if (interface%steering%init .and. interface%params_siml%loutnpp) then
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

    print*,'interface%params_siml%loutnpp ', interface%params_siml%loutnpp

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

    end if


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

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    if (interface%params_siml%loutnpp) outdrleaf(:,doy,jpngr) = drleaf(:)
    if (interface%params_siml%loutnpp) outdrroot(:,doy,jpngr) = drroot(:)
    if (interface%params_siml%loutnpp) outdrgrow(:,doy,jpngr) = drgrow(:)

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    outarleaf(:,jpngr) = outarleaf(:,jpngr) + drleaf(:)
    outarroot(:,jpngr) = outarroot(:,jpngr) + drroot(:)
    outargrow(:,jpngr) = outargrow(:,jpngr) + drgrow(:)


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

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutnpp) then
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
      end if
    end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    itime = real(interface%steering%outyear)

    write(452,999) itime, outarleaf(:,jpngr)
    write(453,999) itime, outarroot(:,jpngr)
    write(455,999) itime, outargrow(:,jpngr)

    return

    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_npp


end module md_npp
