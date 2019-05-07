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
    real, parameter :: f_labl_max = 0.5
    real, parameter :: f_labl_min = 0.05
    real, parameter :: k_labl = 0.1
    real :: c_reserve_max, c_to_storage, f_to_storage, f_demand_labl
    type( orgpool ) :: pgrow, org_reserve_min

    ! xxx debug
    type( orgpool ) :: bal1, bal2, bald

    ! Variables N balance test
    logical, parameter :: baltest_trans = .true.  ! set to true to do mass conservation test during transient simulation
    logical :: verbose
    logical :: baltest
    type( orgpool ) :: orgtmp1, orgtmp2, orgbal1, orgbal2
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
      ! C available for growth in labile, leaves or roots
      pgrow = orgpool( carbon( plant_fluxes(pft)%dnpp%c12 - plant_fluxes(pft)%dcex ), plant_fluxes(pft)%dnup )
      if (verbose) print*,'pgrow 1: ', pgrow
      if (verbose) print*,'plabl   1: ', plant(pft)%plabl

      if (pgrow%c%c12<0.0) then

        if (verbose) print*,'Negative growth:', pgrow%c%c12

        ! deplete labile pool and no growth this day
        f_to_storage = 1.0
        call orgmv( orgfrac( f_to_storage, pgrow ), pgrow, plant(pft)%plabl )

      else

        ! Get the maximum and minimum NSC storage pool as a function of leaf and root mass
        c_reserve_max   = f_labl_max * ( plant(pft)%pleaf%c%c12 + plant(pft)%proot%c%c12 )
        org_reserve_min = orgfrac( f_labl_min, orgplus( plant(pft)%pleaf, plant(pft)%proot ) )

        ! Get "demand-factor" of NSC storage pool, proportional to its "emptiness"
        if (c_reserve_max==0.0) then
          f_demand_labl = 0.0
        else
          f_demand_labl = max( 0.0, 1.0 - plant(pft)%plabl%c%c12 / c_reserve_max )
        end if

        ! Get C flux to fill up NSC storage
        c_to_storage = min( f_demand_labl * (c_reserve_max - plant(pft)%plabl%c%c12), pgrow%c%c12 )
        if (pgrow%c%c12==0.0) then
          f_to_storage = 0.0
        else
          f_to_storage = c_to_storage / pgrow%c%c12
        end if
        if (verbose) print*,'fraction to storage:', f_to_storage

        ! C and N balance 1: addition to storage
        call orgmv( orgfrac( f_to_storage, pgrow ), pgrow, plant(pft)%plabl )
        ! if (verbose) print*,'pgrow 2: ', pgrow
        ! if (verbose) print*,'plabl   2: ', plant(pft)%plabl

        ! C and N balance 2: add C and N available for growth from turnover of NSC storage pool
        call orgmv( orgfrac( k_labl, orgminus( plant(pft)%plabl, org_reserve_min ) ), plant(pft)%plabl, pgrow )
        if (verbose) print*,'pgrow 3: ', pgrow
        if (verbose) print*,'plabl   3: ', plant(pft)%plabl

        if (params_pft_plant(pft)%grass) then

          if (pgrow%c%c12>0.0) then
            !------------------------------------------------------------------
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by all labile N,
            ! discounted by the yield factor.
            !------------------------------------------------------------------
            if (plant(pft)%pleaf%c%c12==0.0) then
              call update_leaftraits_init( plant(pft), pft, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )
              ! leaftraits(pft) = get_leaftraits_init( pft, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )
            end if

            !-------------------------------------------------------------------
            ! GROWTH RESPIRATION, NPP
            !-------------------------------------------------------------------
            plant_fluxes(pft)%drgrow = (1.0 - params_plant%growtheff) * pgrow%c%c12
            pgrow%c%c12 = pgrow%c%c12 - plant_fluxes(pft)%drgrow       

            !-------------------------------------------------------------------
            ! Determine allocation to roots and leaves, fraction given by 'frac_leaf'
            !-------------------------------------------------------------------
            dcleaf(pft) = params_plant%frac_leaf * pgrow%c%c12
            dcroot(pft) = (1.0 - params_plant%frac_leaf) * pgrow%c%c12

            !-------------------------------------------------------------------
            ! LEAF ALLOCATION
            !-------------------------------------------------------------------
            ! if (verbose) print*, 'calling allocate_leaf() ... '
            ! if (verbose) print*, '              with state variables:'
            ! if (verbose) print*, '              pleaf = ', plant(:)%pleaf
            ! if (verbose) print*, '              proot = ', plant(:)%proot
            ! if (verbose) print*, '              pgrow = ', pgrow
            ! if (baltest) orgbal1 = orgplus( plant(pft)%pleaf, plant(pft)%proot, pgrow )
            call allocate_leaf( pft, &
                                dcleaf(pft), &
                                plant(pft)%pleaf, &
                                plant(pft)%lai_ind, &
                                pgrow, &
                                solar%meanmppfd(:), &
                                out_pmodel(pft,:)%actnv_unitiabs, &
                                dnleaf(pft), &
                                plant_fluxes(pft) &
                                )
            ! if (baltest) orgbal2 = orgplus( plant(pft)%pleaf, plant(pft)%proot, pgrow )
            ! if (verbose) print*, '              ==> returned from allocate_leaf(): '
            ! if (verbose) print*, '              pleaf = ', plant(:)%pleaf
            ! if (verbose) print*, '              proot = ', plant(:)%proot
            ! if (verbose) print*, '              pgrow = ', pgrow
            ! if (baltest) orgbal1 = orgminus( orgbal2, orgbal1 )
            ! if (baltest) print*, '            balance =', orgbal1
            ! if (baltest .and. abs(orgbal1%c%c12)>eps*10) stop 'balance A not satisfied for C'

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
            ! if (verbose) print*, 'calling allocate_root() ... '
            ! if (verbose) print*, '              with state variables:'
            ! if (verbose) print*, '              pleaf = ', plant(:)%pleaf
            ! if (verbose) print*, '              proot = ', plant(:)%proot
            ! if (verbose) print*, '              pgrow = ', pgrow
            ! if (baltest) orgbal1 = orgplus( plant(pft)%pleaf, plant(pft)%proot, pgrow )
            call allocate_root( &
                                pft, &
                                dcroot(pft), &
                                plant(pft)%proot, &
                                pgrow, &
                                dnroot(pft), &
                                plant_fluxes(pft) &
                                )
            ! if (baltest) orgbal2 = orgplus( plant(pft)%pleaf, plant(pft)%proot, pgrow )
            ! if (verbose) print*, '              ==> returned from allocate_root(): '
            ! if (verbose) print*, '              pleaf = ', plant(:)%pleaf
            ! if (verbose) print*, '              proot = ', plant(:)%proot
            ! if (verbose) print*, '              pgrow = ', pgrow
            ! if (baltest) orgbal1 = orgminus( orgbal2, orgbal1 )
            ! if (baltest) print*, '            balance =', orgbal1
            ! if (baltest .and. abs(orgbal1%c%c12)>eps*10) stop 'after allocate_root() balance not satisfied for C'
            ! if (baltest .and. pgrow%c%c12>eps*10) stop 'C for growth not used up' 

            if ( pgrow%n%n14<0.0 ) pgrow%n%n14 = 0.0

          end if

        else

          stop 'allocation_daily not implemented for trees'

        end if

      end if

    end do

    ! print*, '--- END allocation_daily:'

  end subroutine allocation_daily


  subroutine allocate_leaf( pft, dcleaf, pleaf, lai, pgrow, meanmppfd, nv, dnleaf, plant_fluxes )
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
    real, intent(in)                       :: dcleaf
    type(orgpool), intent(inout)           :: pleaf
    real, intent(out)                      :: lai
    type(orgpool), intent(inout)           :: pgrow
    real, dimension(nmonth), intent(in)    :: meanmppfd
    real, dimension(nmonth), intent(in)    :: nv
    real, intent(out)                      :: dnleaf
    type(plant_fluxes_type), intent(inout) :: plant_fluxes

    ! local variables
    real :: nleaf0

    ! Calculate LAI as a function of leaf C
    pleaf%c%c12  = pleaf%c%c12 + dcleaf
    lai = get_lai( pft, pleaf%c%c12, meanmppfd(:), nv(:) )

    ! calculate canopy-level leaf N as a function of LAI and implied 
    nleaf0      = pleaf%n%n14      
    pleaf%n%n14 = get_leaf_n_canopy( pft, lai, meanmppfd(:), nv(:) )
    dnleaf = pleaf%n%n14 - nleaf0

    ! substract from growth pool
    call orgsub( orgpool( carbon(dcleaf), nitrogen(dnleaf) ), pgrow )

    if ( pgrow%c%c12 < (-1)*eps) then
      stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf C'
    else if ( pgrow%c%c12 < 0.0 ) then
      ! numerical imprecision
      print*,'numerical imprecision?'
      pgrow%c%c12 = 0.0
    end if

    ! If labile N gets negative, account gap as N fixation
    if ( pgrow%n%n14 < 0.0 ) then
      plant_fluxes%dnup%n14 = plant_fluxes%dnup%n14 - pgrow%n%n14
      plant_fluxes%dnup_fix = plant_fluxes%dnup_fix - pgrow%n%n14
      pgrow%n%n14 = 0.0
    end if

    ! XXX this is the N mass conserving way:
    ! if ( nlabl < -1e-8 ) then
    !   stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf N'
    ! else if ( nlabl < 0.0 ) then
    !   ! numerical imprecision
    !   nlabl = 0.0
    ! end if
    
  end subroutine allocate_leaf


  subroutine allocate_root( pft, dcroot, proot, pgrow, dnroot, plant_fluxes )
    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: plant_fluxes_type, params_plant, params_pft_plant
    use md_params_core, only: eps

    ! arguments
    integer, intent(in) :: pft
    real, intent(in)    :: dcroot
    type(orgpool), intent(inout)           :: proot
    type(orgpool), intent(inout)           :: pgrow
    real, intent(out)                      :: dnroot
    type(plant_fluxes_type), intent(inout) :: plant_fluxes

    !------------------------------------------------------------------
    logical, parameter :: baltest = .true.
    type( orgpool ) :: orgbal1, orgbal2
    !------------------------------------------------------------------

    ! get change in root N from prescribed root C:N ratio
    dnroot = dcroot * params_pft_plant(pft)%r_ntoc_root

    ! if (baltest) orgbal1 = orgplus( pgrow, proot )

    ! ! update root pools
    ! proot%c%c12 = proot%c%c12 + dcroot
    ! proot%n%n14 = proot%n%n14 + dnroot

    ! ! update growth pool
    ! pgrow%c%c12 = pgrow%c%c12 - dcroot
    ! pgrow%n%n14 = pgrow%n%n14 - dnroot

    ! if (baltest) orgbal2 = orgplus( pgrow, proot )
    ! if (baltest) orgbal1 = orgminus( orgbal2, orgbal1 )
    ! if (baltest) print*, 'allocate_root balance =', orgbal1
    ! if (baltest .and. abs(orgbal1%c%c12)>eps*10) stop 'allocate_root: balance not satisfied for C'

    ! if (baltest) orgbal1 = orgplus( pgrow, proot )
    call orgmv( orgpool( carbon(dcroot), nitrogen(dnroot) ), pgrow, proot )
    ! if (baltest) orgbal2 = orgplus( pgrow, proot )
    ! if (baltest) orgbal1 = orgminus( orgbal2, orgbal1 )
    ! if (baltest) print*, 'allocate_root balance =', orgbal1
    ! if (baltest .and. abs(orgbal1%c%c12)>eps*10) stop 'allocate_root: balance not satisfied for C'

    if ( pgrow%c%c12 < (-1)*eps ) then
      stop 'ALLOCATE_ROOT: trying to remove too much from labile pool: leaf C'
    else if ( pgrow%c%c12 < 0.0 ) then
      ! numerical imprecision
      print*,'allocate_root: warning: numerical imprecision'
      pgrow%c%c12 = 0.0
    end if

    ! If labile N gets negative, account gap as N fixation
    if ( pgrow%n%n14 < 0.0 ) then
      plant_fluxes%dnup%n14 = plant_fluxes%dnup%n14 - pgrow%n%n14
      plant_fluxes%dnup_fix = plant_fluxes%dnup_fix - pgrow%n%n14
      pgrow%n%n14 = 0.0
      ! print*,'allocate_root: warning: filling up missing N'
    end if
  
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
