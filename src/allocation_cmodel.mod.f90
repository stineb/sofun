module md_allocation
  !////////////////////////////////////////////////////////////////
  ! ALLOCATION MODULE
  ! Binary allocation formulation: either to leaves or to roots.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_classdefs
  use md_params_core, only: npft, nlu, maxgrid, ndaymonth, ndayyear, &
    c_molmass, n_molmass, nmonth

  implicit none

  private 
  public allocation_daily, initio_allocation, initoutput_allocation, &
    getout_daily_allocation, writeout_ascii_allocation

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  real, dimension(npft) :: dcleaf
  real, dimension(npft) :: dnleaf
  real, dimension(npft) :: dcroot
  real, dimension(npft) :: dnroot

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! output variables
  real, dimension(npft,maxgrid) :: outaCalclm
  real, dimension(npft,maxgrid) :: outaNalclm
  real, dimension(npft,maxgrid) :: outaCalcrm
  real, dimension(npft,maxgrid) :: outaNalcrm

contains

  subroutine allocation_daily( plant, plant_fluxes, solar, out_pmodel, dtemp )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: plant_type, plant_fluxes_type, params_plant, params_pft_plant, &
      update_leaftraits, update_leaftraits_init, get_fapar
    use md_params_core, only: eps
    use md_waterbal, only: solartype
    use md_gpp, only: outtype_pmodel

    ! ! xxx debug
    ! use md_nuptake, only: calc_dnup, outtype_calc_dnup
    ! use md_waterbal, only: solar, evap
    ! use md_gpp, only: calc_dgpp, calc_drd
    ! use md_npp, only: calc_resp_maint, calc_cexu
    ! use md_gpp, only: drd 
    ! use md_plant, only: dgpp, dnpp, drleaf, drroot, dcex, dnup
    ! use md_interface

    ! arguments
    type( plant_type ), dimension(npft), intent(inout)        :: plant ! npft counts over PFTs in all land units (tiles)
    type( plant_fluxes_type ), dimension(npft), intent(inout) :: plant_fluxes ! npft counts over PFTs in all land units (tiles)
    type( outtype_pmodel ), dimension(npft,nmonth), intent(in):: out_pmodel
    type( solartype ), intent(in)                             :: solar
    real, intent(in)                                          :: dtemp   ! air temperaure, deg C

    ! local variables
    integer :: lu
    integer :: pft
    real :: avl
    real, parameter :: f_labl = 0.1
    real, parameter :: k_labl = 0.05
    real :: c_reserve_pot, c_to_growth, c_to_storage, c_to_growth_direct, c_to_growth_from_storage


    ! xxx debug
    type( orgpool ) :: bal1, bal2, bald

    ! Variables N balance test
    logical, parameter :: baltest_trans = .false.  ! set to true to do mass conservation test during transient simulation
    logical :: verbose = .false.  ! set to true to activate verbose mode
    logical :: baltest
    type( orgpool ) :: orgtmp1, orgtmp2, orgbal1
    real :: ctmp

    !------------------------------------------------------------------
    baltest = .false.
    verbose = .false.
    !------------------------------------------------------------------

    ! initialise
    dcleaf(:) = 0.0
    dnleaf(:) = 0.0
    dcroot(:) = 0.0
    dnroot(:) = 0.0
    plant_fluxes(:)%drgrow = 0.0

    do pft=1,npft

      lu = params_pft_plant(pft)%lu_category

      !/////////////////////////////////////////////////////////////////////////
      ! Allocation to labile pool
      !-------------------------------------------------------------------------
      ! Get the maximum potential NSC storage pool as a function of leaf and root mass
      c_reserve_pot = f_labl * ( plant(pft)%pleaf%c%c12 + plant(pft)%proot%c%c12 )

      ! Get "demand-factor" of NSC storage pool, proportional to its "emptiness"
      f_demand_labl = max( 0.0, 1.0 - plant(pft)%plabl%c%c12 / c_reserve_pot )

      ! Get C flux to fill up NSC storage
      c_to_storage = f_demand_labl * ( plant_fluxes(pft)%dnpp%c12 - plant_fluxes(pft)%dcex )


      ! Get C that allocated to fill up storage and used directly for growth (because storage is full)
      c_to_growth = 0.0
      c_to_storage = min( c_reserve_pot - plant(pft)%plabl%c%c12, ( plant_fluxes(pft)%dnpp%c12 - plant_fluxes(pft)%dcex ) )
      c_to_growth_direct = ( plant_fluxes(pft)%dnpp%c12 - plant_fluxes(pft)%dcex ) - c_to_storage

      ! Put C in storage
      plant(pft)%plabl%c%c12 = plant(pft)%plabl%c%c12 + c_to_storage

      ! Get C for growth from storage
      c_to_growth_from_storage = k_labl * plant(pft)%plabl%c%c12
      plant(pft)%plabl%c%c12 = plant(pft)%plabl%c%c12 - c_to_growth_from_storage

      ! total C for growth
      avl = c_to_growth_from_storage + c_to_growth_direct

      ! tmp = plant_fluxes(pft)%dnpp%c12 - plant_fluxes(pft)%dcex
      print*,'GPP, Rl, Rr, Cex, dC, Cl, LAI, Cb: ', plant_fluxes(pft)%dgpp, plant_fluxes(pft)%drleaf, plant_fluxes(pft)%drroot, plant_fluxes(pft)%dcex, avl, plant(pft)%pleaf, plant(pft)%lai_ind, plant(pft)%plabl

      ! call ccp( cminus( plant_fluxes(pft)%dnpp, carbon(plant_fluxes(pft)%dcex) ), plant(pft)%plabl%c )

      ! if (plant(pft)%plabl%c%c12 < (-1)*eps) stop 'after npp labile C is neg.'
      ! if (plant(pft)%plabl%n%n14 < (-1)*eps) stop 'after npp labile N is neg.'

      if (params_pft_plant(pft)%grass) then

        if ( plant(pft)%plabl%c%c12>0.0 .and. dtemp>0.0 ) then

          !------------------------------------------------------------------
          ! Calculate maximum C allocatable based on current labile pool size.
          ! Maximum is the lower of all labile C and the C to be matched by all labile N,
          ! discounted by the yield factor.
          !------------------------------------------------------------------
          if (plant(pft)%pleaf%c%c12==0.0) then
            call update_leaftraits_init( plant(pft), pft, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )
            ! leaftraits(pft) = get_leaftraits_init( pft, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )
          end if

          ! Determine allocation to roots and leaves, fraction given by 'frac_leaf'
          ! avl = max( 0.0, plant(pft)%plabl%c%c12 - f_labl * plant(pft)%pleaf%c%c12 )
          dcleaf(pft) = params_plant%frac_leaf * params_plant%growtheff * avl
          dcroot(pft) = (1.0 - params_plant%frac_leaf) * params_plant%growtheff * avl
          dnroot(pft) = dcroot(pft) * params_pft_plant(pft)%r_ntoc_root          

          ! print*,'         doy, pleaf ', doy,  pleaf

          !-------------------------------------------------------------------
          ! LEAF ALLOCATION
          !-------------------------------------------------------------------
          if (baltest) orgtmp1 = orgminus( orgplus( plant(pft)%pleaf, plant(pft)%proot, plant(pft)%plabl, orgpool( carbon(plant_fluxes(pft)%drgrow), nitrogen(0.0) ) ), orgpool(carbon(0.0),plant_fluxes(pft)%dnup) )
          if (verbose) print*, 'calling allocate_leaf() ... '
          if (verbose) print*, '              with state variables:'
          if (verbose) print*, '              pleaf = ', plant(:)%pleaf
          if (verbose) print*, '              proot = ', plant(:)%proot
          if (verbose) print*, '              plabl = ', plant(:)%plabl
          if (verbose) print*, '              drgrow= ', plant_fluxes(:)%drgrow
          if (verbose) print*, '              dnup  = ', plant_fluxes(1)%dnup%n14
          call allocate_leaf( &
            pft, dcleaf(pft), &
            plant(pft)%pleaf%c%c12, plant(pft)%pleaf%n%n14, &
            plant(pft)%plabl%c%c12, plant(pft)%plabl%n%n14, &
            solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs, &
            plant(pft)%lai_ind, dnleaf(pft), &
            plant_fluxes(pft) &
            )
          if (verbose) print*, '              ==> returned: '
          if (verbose) print*, '              pleaf = ', plant(:)%pleaf
          if (verbose) print*, '              proot = ', plant(:)%proot
          if (verbose) print*, '              plabl = ', plant(:)%plabl
          if (baltest) ctmp = ( 1.0 - params_plant%growtheff ) * ( dcleaf(pft) ) / params_plant%growtheff
          if (verbose) print*, '              drgrow= ', ctmp
          if (verbose) print*, '              dnup  = ', plant_fluxes(1)%dnup%n14
          if (baltest) orgtmp2 = orgminus( orgplus( plant(pft)%pleaf, plant(pft)%proot, plant(pft)%plabl, orgpool( carbon(ctmp), nitrogen(0.0) ) ), orgpool(carbon(0.0),plant_fluxes(pft)%dnup) )
          if (baltest) orgbal1 = orgminus( orgtmp2, orgtmp1 )
          if (baltest) print*, '       balance A =', orgbal1
          if (baltest .and. abs(orgbal1%c%c12)>eps) stop 'balance A not satisfied for C'
          if (baltest .and. abs(orgbal1%n%n14)>eps) stop 'balance A not satisfied for N'

          !-------------------------------------------------------------------  
          ! Update leaf traits
          !-------------------------------------------------------------------  
          call update_leaftraits( plant(pft), pft, plant(pft)%lai_ind, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )

          !-------------------------------------------------------------------  
          ! Update fpc_grid and fapar_ind (not lai_ind)
          !-------------------------------------------------------------------  
          plant(pft)%fapar_ind = get_fapar( plant(pft)%lai_ind )
          ! print*,'plant(pft)%fapar_ind: ', plant(pft)%fapar_ind

          !-------------------------------------------------------------------
          ! ROOT ALLOCATION
          !-------------------------------------------------------------------
          if (baltest) orgtmp1 = orgminus( orgplus( plant(pft)%pleaf, plant(pft)%proot, plant(pft)%plabl, orgpool( carbon(plant_fluxes(pft)%drgrow), nitrogen(0.0) ) ), orgpool(carbon(0.0),plant_fluxes(pft)%dnup) )
          if (verbose) print*, 'calling allocate_root() ... '
          if (verbose) print*, '              with state variables:'
          if (verbose) print*, '              pleaf = ', plant(:)%pleaf
          if (verbose) print*, '              proot = ', plant(:)%proot
          if (verbose) print*, '              plabl = ', plant(:)%plabl
          if (verbose) print*, '              drgrow= ', plant_fluxes(:)%drgrow
          if (verbose) print*, '              dnup  = ', plant_fluxes(1)%dnup%n14
          call allocate_root( &
            pft, dcroot(pft), dnroot(pft), &
            plant(pft)%proot%c%c12, plant(pft)%proot%n%n14, &
            plant(pft)%plabl%c%c12, plant(pft)%plabl%n%n14,  &
            plant_fluxes(pft) &
            )
          if (verbose) print*, '              ==> returned: '
          if (verbose) print*, '              pleaf = ', plant(:)%pleaf
          if (verbose) print*, '              proot = ', plant(:)%proot
          if (verbose) print*, '              plabl = ', plant(:)%plabl
          if (baltest) ctmp = ( 1.0 - params_plant%growtheff ) * ( dcroot(pft) ) / params_plant%growtheff
          if (verbose) print*, '              drgrow= ', ctmp
          if (verbose) print*, '              dnup  = ', plant_fluxes(1)%dnup%n14
          if (baltest) orgtmp2 = orgminus( orgplus( plant(pft)%pleaf, plant(pft)%proot, plant(pft)%plabl, orgpool( carbon(ctmp), nitrogen(0.0) ) ), orgpool(carbon(0.0),plant_fluxes(pft)%dnup) )
          if (baltest) orgbal1 = orgminus( orgtmp2, orgtmp1 )
          if (baltest) print*, '       balance B =', orgbal1
          if (baltest .and. abs(orgbal1%c%c12)>eps) stop 'balance B not satisfied for C'
          if (baltest .and. abs(orgbal1%n%n14)>eps) stop 'balance B not satisfied for N'

          !-------------------------------------------------------------------
          ! GROWTH RESPIRATION, NPP
          !-------------------------------------------------------------------
          ! add growth respiration to autotrophic respiration and substract from NPP
          ! (note that NPP is added to plabl in and growth resp. is implicitly removed
          ! from plabl above)
          plant_fluxes(pft)%drgrow   = ( 1.0 - params_plant%growtheff ) * ( dcleaf(pft) + dcroot(pft) ) / params_plant%growtheff

          if ( plant(pft)%plabl%n%n14<0.0 ) plant(pft)%plabl%n%n14 = 0.0

        end if

      else

        stop 'allocation_daily not implemented for trees'

      end if

    end do

    ! print*, '--- END allocation_daily:'

  end subroutine allocation_daily


  subroutine allocate_leaf( pft, mydcleaf, cleaf, nleaf, clabl, nlabl, meanmppfd, nv, lai, mydnleaf, plant_fluxes )
    !///////////////////////////////////////////////////////////////////
    ! LEAF ALLOCATION
    ! Sequence of steps:
    ! - increment foliage C pool
    ! - update LAI
    ! - calculate canopy-level foliage N as a function of LAI 
    ! - reduce labile pool by C and N increments
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: plant_fluxes_type, params_plant, get_leaf_n_canopy, get_lai
    use md_params_core, only: eps

    ! arguments
    integer, intent(in)                    :: pft
    real, intent(in)                       :: mydcleaf
    real, intent(inout)                    :: cleaf, nleaf
    real, intent(inout)                    :: clabl, nlabl
    real, dimension(nmonth), intent(in)    :: meanmppfd
    real, dimension(nmonth), intent(in)    :: nv
    real, intent(out)                      :: lai
    real, intent(out)                      :: mydnleaf
    type(plant_fluxes_type), intent(inout) :: plant_fluxes

    ! local variables
    real :: nleaf0
    real :: dclabl, dnlabl

    ! Calculate LAI as a function of leaf C
    cleaf  = cleaf + mydcleaf
    lai = get_lai( pft, cleaf, meanmppfd(:), nv(:) )

    ! calculate canopy-level leaf N as a function of LAI and implied 
    ! depletion of labile N pool
    nleaf0   = nleaf      
    nleaf    = get_leaf_n_canopy( pft, lai, meanmppfd(:), nv(:) )
    mydnleaf = nleaf - nleaf0

    ! depletion of labile C pool is enhanced by growth respiration
    dclabl = 1.0 / params_plant%growtheff * mydcleaf

    ! substract from labile pools
    clabl  = clabl - dclabl
    nlabl  = nlabl - mydnleaf

    if ( clabl < (-1)*eps) then
      stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf C'
    else if ( clabl < 0.0 ) then
      ! numerical imprecision
      print*,'numerical imprecision?'
      clabl = 0.0
    end if

    ! If labile N gets negative, account gap as N fixation
    if ( nlabl < 0.0 ) then
      plant_fluxes%dnup%n14 = plant_fluxes%dnup%n14 - nlabl
      plant_fluxes%dnup_fix = plant_fluxes%dnup_fix - nlabl
      nlabl = 0.0
    end if

    ! XXX this is the N mass conserving way:
    ! if ( nlabl < -1e-8 ) then
    !   stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf N'
    ! else if ( nlabl < 0.0 ) then
    !   ! numerical imprecision
    !   nlabl = 0.0
    ! end if
    
  end subroutine allocate_leaf


  subroutine allocate_root( pft, mydcroot, mydnroot, croot, nroot, clabl, nlabl, plant_fluxes )
    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: plant_fluxes_type, params_plant, params_pft_plant
    use md_params_core, only: eps

    ! arguments
    integer, intent(in) :: pft
    real, intent(in)    :: mydcroot
    real, intent(in)    :: mydnroot
    real, intent(inout) :: croot, nroot
    real, intent(inout) :: clabl, nlabl
    type(plant_fluxes_type), intent(inout) :: plant_fluxes

    ! local variables
    real :: dclabl
    real :: dnlabl

    ! update root pools
    croot = croot + mydcroot
    nroot = nroot + mydnroot

    ! depletion of labile C pool is enhanced by growth respiration
    dclabl = 1.0 / params_plant%growtheff * mydcroot

    ! substract from labile pools
    clabl  = clabl - dclabl
    nlabl  = nlabl - mydnroot

    if ( clabl < (-1)*eps ) then
      stop 'ALLOCATE_ROOT: trying to remove too much from labile pool: leaf C'
    else if ( clabl < 0.0 ) then
      ! numerical imprecision
      clabl = 0.0
    end if

    ! If labile N gets negative, account gap as N fixation
    if ( nlabl < 0.0 ) then
      plant_fluxes%dnup%n14 = plant_fluxes%dnup%n14 - nlabl
      plant_fluxes%dnup_fix = plant_fluxes%dnup_fix - nlabl
      nlabl = 0.0
    end if

    ! XXX this is the N mass conserving way:
    ! if ( nlabl < -1e-8 ) then
    !   stop 'ALLOCATE_ROOT: trying to remove too much from labile pool: leaf N'
    ! else if ( nlabl < 0.0 ) then
    !   ! numerical imprecision
    !   nlabl = 0.0
    ! end if
  
  end subroutine allocate_root


  subroutine initio_allocation()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use md_interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------
    if (interface%params_siml%loutalloc) then

      ! C ALLOCATED TO LEAF GROWTH 
      filnam=trim(prefix)//'.a.calclm.out'
      open(350,file=filnam,err=999,status='unknown')

      ! N ALLOCATED TO LEAF GROWTH 
      filnam=trim(prefix)//'.a.nalclm.out'
      open(351,file=filnam,err=999,status='unknown')

      ! C ALLOCATED TO ROOT GROWTH 
      filnam=trim(prefix)//'.a.calcrm.out'
      open(352,file=filnam,err=999,status='unknown')

      ! N ALLOCATED TO ROOT GROWTH 
      filnam=trim(prefix)//'.a.nalcrm.out'
      open(353,file=filnam,err=999,status='unknown')

    end if

    return

    999  stop 'INITIO_ALLOCATION: error opening output files'

  end subroutine initio_allocation


  subroutine initoutput_allocation()
    !////////////////////////////////////////////////////////////////
    !  Initialises nuptake-specific output variables
    !----------------------------------------------------------------
    use md_interface

    ! xxx remove their day-dimension
    if (interface%params_siml%loutalloc) then
      outaCalclm(:,:) = 0.0
      outaNalclm(:,:) = 0.0
      outaCalcrm(:,:) = 0.0
      outaNalcrm(:,:) = 0.0
    end if
    
    ! print*, 'initialising outaCalloc',outaCalloc

  end subroutine initoutput_allocation


  subroutine getout_daily_allocation( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    if (interface%params_siml%loutalloc) then
      outaCalclm(:,jpngr) = outaCalclm(:,jpngr) + dcleaf(:) 
      outaNalclm(:,jpngr) = outaNalclm(:,jpngr) + dnleaf(:)
      outaCalcrm(:,jpngr) = outaCalcrm(:,jpngr) + dcroot(:) 
      outaNalcrm(:,jpngr) = outaNalcrm(:,jpngr) + dnroot(:)
    end if

    ! print*, 'collecting outaCalloc',outaCalloc

  end subroutine getout_daily_allocation


  subroutine writeout_ascii_allocation( year )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: year       ! simulation year

    ! local variables
    real    :: itime
    integer :: jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutalloc) then

      itime = real(interface%steering%outyear)

      write(350,999) itime, sum(outaCalclm(:,jpngr))
      write(351,999) itime, sum(outaNalclm(:,jpngr))
      write(352,999) itime, sum(outaCalcrm(:,jpngr))
      write(353,999) itime, sum(outaNalcrm(:,jpngr))
    end if

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_allocation

end module md_allocation
