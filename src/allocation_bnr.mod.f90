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

  subroutine allocation_daily( jpngr, doy, moy, dtemp )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant, params_pft_plant, pleaf, proot, &
      plabl, drgrow, lai_ind, nind, canopy, leaftraits, &
      get_canopy, get_leaftraits, get_leaftraits_init, frac_leaf
    use md_waterbal, only: solar
    use md_gpp, only: mlue, mrd_unitiabs, mactnv_unitiabs
    use md_soiltemp, only: dtemp_soil
    use md_ntransform, only: pninorg
    use md_params_core, only: eps

    ! xxx debug
    use md_nuptake, only: calc_dnup, outtype_calc_dnup
    use md_waterbal, only: solar, evap
    use md_gpp, only: calc_dgpp, calc_drd
    use md_npp, only: calc_resp_maint, calc_cexu, deactivate_root
    use md_gpp, only: dgpp, drd 
    use md_plant, only: dnpp, drleaf, drroot, dcex, dnup
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy     ! day of year
    integer, intent(in) :: moy     ! month of year
    real,    intent(in) :: dtemp   ! air temperaure, deg C

    ! local variables
    integer :: lu
    integer :: pft
    integer :: usemoy, usedoy
    real    :: cavl, navl, avl
    real, parameter :: freserve = 0.0

    ! xxx debug
    type( orgpool ) :: bal1, bal2, bald
    logical, save :: toleaves = .true.       ! boolean determining whether C and N in this time step are allocated to leaves or roots
    logical, save :: nignore = .true.

    ! Variables N balance test
    logical, parameter :: baltest_trans = .false.  ! set to true to do mass conservation test during transient simulation
    logical :: verbose = .false.  ! set to true to activate verbose mode
    logical :: baltest
    type( orgpool ) :: orgtmp1, orgtmp2, orgbal1
    real :: ctmp

    !------------------------------------------------------------------
    ! Turn mass conservation tests on and off
    !------------------------------------------------------------------
    baltest = .false.
    verbose = .false.

    !------------------------------------------------------------------
    ! initialise
    !------------------------------------------------------------------
    dcleaf(:) = 0.0
    dnleaf(:) = 0.0
    dcroot(:) = 0.0
    dnroot(:) = 0.0
    drgrow(:) = 0.0


    do pft=1,npft

      lu = params_pft_plant(pft)%lu_category

      if (params_pft_plant(pft)%grass) then

        if ( interface%steering%dofree_alloc ) then
          !------------------------------------------------------------------
          ! Free allocation
          !------------------------------------------------------------------

          if ( plabl(pft,jpngr)%c%c12>0.0 .and. plabl(pft,jpngr)%n%n14>0.0 .and. dtemp>0.0 ) then
            !------------------------------------------------------------------
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by all labile N,
            ! discounted by the yield factor.
            !------------------------------------------------------------------
            if (pleaf(pft,jpngr)%c%c12==0.0) then
              leaftraits(pft) = get_leaftraits_init( pft, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
            end if

            ! first day of free allocation, put all to leaves
            if (frac_leaf(pft)==0.5) frac_leaf(pft) = 1.0

            ! abort when labile N pool gets negative upon allocation
            nignore = .false.

            ! get allocatable C to roots and leaves 
            if (frac_leaf(pft)==1.0) then

              dcleaf(pft) = min( &
                params_plant%growtheff * (plabl(pft,jpngr)%c%c12), &
                (plabl(pft,jpngr)%n%n14) * leaftraits(pft)%r_cton_leaf &
                )
              dcroot(pft) = 0.0
              dnroot(pft) = 0.0

            else if (frac_leaf(pft)==0.0) then

              dcroot(pft) = min( &
                params_plant%growtheff * (plabl(pft,jpngr)%c%c12), &
                (plabl(pft,jpngr)%n%n14) * params_pft_plant(pft)%r_cton_root &
                )
              dnroot(pft) = dcroot(pft) * params_pft_plant(pft)%r_ntoc_root
              dcleaf(pft) = 0.0

            else
              print*,'frac_leaf ', frac_leaf
              stop 'frac_leaf/=0 or 1'

            end if

            print*,'----', doy, '----'
            print*,'plabl before ', plabl(pft,jpngr)
            print*,'frac_leaf    ', frac_leaf(pft)

            !-------------------------------------------------------------------
            ! LEAF ALLOCATION
            !-------------------------------------------------------------------
            if (dcleaf(pft)>0.0) then

              call allocate_leaf( &
                pft, dcleaf(pft), &
                pleaf(pft,jpngr)%c%c12, pleaf(pft,jpngr)%n%n14, &
                plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
                solar%meanmppfd(:), mactnv_unitiabs(pft,:), &
                lai_ind(pft,jpngr), dnleaf(pft), nignore=nignore &
                )

              !-------------------------------------------------------------------  
              ! Update leaf traits
              !-------------------------------------------------------------------  
              leaftraits(pft) = get_leaftraits( pft, lai_ind(pft,jpngr), solar%meanmppfd(:), mactnv_unitiabs(pft,:) )

              !-------------------------------------------------------------------  
              ! Update fpc_grid and fapar_ind (not lai_ind)
              !-------------------------------------------------------------------  
              canopy(pft) = get_canopy( lai_ind(pft,jpngr) )

            end if

            !-------------------------------------------------------------------
            ! ROOT ALLOCATION
            !-------------------------------------------------------------------
            if (dcroot(pft)>0.0) then

              call allocate_root( &
                pft, dcroot(pft), dnroot(pft), &
                proot(pft,jpngr)%c%c12, proot(pft,jpngr)%n%n14, &
                plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
                nignore=nignore &
                )

            end if

            !-------------------------------------------------------------------
            ! GROWTH RESPIRATION, NPP
            !-------------------------------------------------------------------
            ! add growth respiration to autotrophic respiration and substract from NPP
            ! (note that NPP is added to plabl in and growth resp. is implicitly removed
            ! from plabl above)
            drgrow(pft)   = ( 1.0 - params_plant%growtheff ) * ( dcleaf(pft) + dcroot(pft) ) / params_plant%growtheff
            
            print*,'plabl after ', plabl(pft,jpngr)

            if ( plabl(pft,jpngr)%n%n14 == 0.0 ) then

              frac_leaf(pft) = 0.0 
            
            else if ( cton( plabl(pft,jpngr) ) > max( params_pft_plant(pft)%r_cton_root, leaftraits(pft)%r_cton_leaf ) ) then
            
              frac_leaf(pft) = 0.0 
            
            else
            
              frac_leaf(pft) = 1.0
            
            end if  

          !   if ( plabl(pft,jpngr)%c%c12 == 0.0 .or. plabl(pft,jpngr)%n%n14 == 0.0 &
          !     .or. cton( plabl(pft,jpngr)) > max( params_pft_plant(pft)%r_cton_root, leaftraits(pft)%r_cton_leaf ) ) then

          !     !-------------------------------------------------------------------  
          !     ! If C is left in labile pool, then N wasn't sufficient => get more
          !     ! N by allocating to roots next day.
          !     !-------------------------------------------------------------------  
          !     if ( plabl(pft,jpngr)%c%c12 > 0.0 ) frac_leaf = 0.0

          !     !-------------------------------------------------------------------  
          !     ! If N is left in labile pool, then C wasn't sufficient => get more
          !     ! C by allocating to leaves next day.
          !     !-------------------------------------------------------------------  
          !     if ( plabl(pft,jpngr)%n%n14 > 0.0 ) frac_leaf = 1.0

          !   ! else 

          !   !   stop 'middle way'
          !   !   frac_leaf = 0.5

          !   end if

          end if

        else
          !------------------------------------------------------------------
          ! Fixed allocation 
          !------------------------------------------------------------------
          if ( plabl(pft,jpngr)%c%c12>0.0 .and. dtemp>0.0 ) then

            !------------------------------------------------------------------
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by all labile N,
            ! discounted by the yield factor.
            !------------------------------------------------------------------
            if (pleaf(pft,jpngr)%c%c12==0.0) then
              leaftraits(pft) = get_leaftraits_init( pft, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
            end if

            ! if ( interface%steering%dofree_alloc ) then
            !   if (frac_leaf(pft)==0.5) frac_leaf(pft) = 1.0
            !   nignore = .false.
            !   if (frac_leaf(pft)==1.0) then
            !     dcleaf(pft) = min( &
            !       params_plant%growtheff * (plabl(pft,jpngr)%c%c12), &
            !       (plabl(pft,jpngr)%n%n14) * leaftraits(pft)%r_cton_leaf &
            !       )
            !     dcroot(pft) = 0.0
            !     dnroot(pft) = 0.0
            !   else if (frac_leaf(pft)==0.0) then
            !     dcroot(pft) = min( &
            !       params_plant%growtheff * (plabl(pft,jpngr)%c%c12), &
            !       (plabl(pft,jpngr)%n%n14) * params_pft_plant(pft)%r_cton_root &
            !       )
            !     dnroot(pft) = dcroot(pft) * params_pft_plant(pft)%r_ntoc_root
            !     dcleaf(pft) = 0.0
            !   else
            !     print*,'frac_leaf ', frac_leaf
            !     stop 'frac_leaf/=0 or 1'
            !   end if
            ! else
              ! Determine allocation to roots and leaves, fraction given by 'frac_leaf'
              nignore = .true.
              avl = max( 0.0, plabl(pft,jpngr)%c%c12 - freserve * pleaf(pft,jpngr)%c%c12 )
              dcleaf(pft) = frac_leaf(pft) * params_plant%growtheff * avl
              dcroot(pft) = (1.0 - frac_leaf(pft)) * params_plant%growtheff * avl
              dnroot(pft) = dcroot(pft) * params_pft_plant(pft)%r_ntoc_root          
            ! end if

            !-------------------------------------------------------------------
            ! LEAF ALLOCATION
            !-------------------------------------------------------------------
            if (dcleaf(pft)>0.0) then

              call allocate_leaf( &
                pft, dcleaf(pft), &
                pleaf(pft,jpngr)%c%c12, pleaf(pft,jpngr)%n%n14, &
                plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
                solar%meanmppfd(:), mactnv_unitiabs(pft,:), &
                lai_ind(pft,jpngr), dnleaf(pft), nignore=nignore &
                )

              !-------------------------------------------------------------------  
              ! Update leaf traits
              !-------------------------------------------------------------------  
              leaftraits(pft) = get_leaftraits( pft, lai_ind(pft,jpngr), solar%meanmppfd(:), mactnv_unitiabs(pft,:) )

              !-------------------------------------------------------------------  
              ! Update fpc_grid and fapar_ind (not lai_ind)
              !-------------------------------------------------------------------  
              canopy(pft) = get_canopy( lai_ind(pft,jpngr) )

            end if

            !-------------------------------------------------------------------
            ! ROOT ALLOCATION
            !-------------------------------------------------------------------
            if (dcroot(pft)>0.0) then

              call allocate_root( &
                pft, dcroot(pft), dnroot(pft), &
                proot(pft,jpngr)%c%c12, proot(pft,jpngr)%n%n14, &
                plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
                nignore=nignore &
                )

            end if

            !-------------------------------------------------------------------
            ! GROWTH RESPIRATION, NPP
            !-------------------------------------------------------------------
            ! add growth respiration to autotrophic respiration and substract from NPP
            ! (note that NPP is added to plabl in and growth resp. is implicitly removed
            ! from plabl above)
            drgrow(pft)   = ( 1.0 - params_plant%growtheff ) * ( dcleaf(pft) + dcroot(pft) ) / params_plant%growtheff

          end if

        end if

      else

        stop 'allocation_daily not implemented for trees'

      end if

    end do

    ! print*, '--- END allocation_daily:'

  end subroutine allocation_daily


  subroutine allocate_leaf( pft, mydcleaf, cleaf, nleaf, clabl, nlabl, meanmppfd, nv, lai, mydnleaf, nignore )
    !///////////////////////////////////////////////////////////////////
    ! LEAF ALLOCATION
    ! Sequence of steps:
    ! - increment foliage C pool
    ! - update LAI
    ! - calculate canopy-level foliage N as a function of LAI 
    ! - reduce labile pool by C and N increments
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant, get_leaf_n_canopy, get_lai, dnup
    use md_nuptake, only: dnup_fix
    use md_params_core, only: eps

    ! arguments
    integer, intent(in)                 :: pft
    real, intent(in)                    :: mydcleaf
    real, intent(inout)                 :: cleaf, nleaf
    real, intent(inout)                 :: clabl, nlabl
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv
    real, intent(out)                   :: lai
    real, intent(out)                   :: mydnleaf
    logical, intent(in)                 :: nignore

    ! local variables
    real :: nleaf0
    real :: dclabl, dnlabl

    ! Calculate LAI as a function of leaf C
    cleaf  = cleaf + mydcleaf
    lai = get_lai( pft, cleaf, meanmppfd(:), nv(:) )

    ! calculate canopy-level leaf N as a function of LAI
    nleaf0   = nleaf      
    nleaf    = get_leaf_n_canopy( pft, lai, meanmppfd(:), nv(:) )
    mydnleaf = nleaf - nleaf0

    ! depletion of labile C pool is enhanced by growth respiration
    dclabl = 1.0 / params_plant%growtheff * mydcleaf

    ! substract from labile pools
    clabl  = clabl - dclabl
    nlabl  = nlabl - mydnleaf

    if ( clabl < -1.0*eps ) then
      stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf C'
    else if ( clabl < 0.0 ) then
      ! numerical imprecision
      ! print*,'numerical imprecision?'
      ! print*,'clabl ', clabl
      ! stop 'allocate leaf'
      clabl = 0.0
    end if

    if (nignore) then
      ! If labile N gets negative, account gap as N fixation
      if ( nlabl < 0.0 ) then
        dnup(pft)%n14 = dnup(pft)%n14 - nlabl
        dnup_fix(pft) = dnup_fix(pft) - nlabl
        nlabl = 0.0
      end if
    else
      if ( nlabl < -1.0*eps ) then
        stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf N'
      else if ( nlabl < 0.0 ) then
        ! numerical imprecision
        ! print*,'numerical imprecision?'
        ! print*,'nlabl ', nlabl
        ! stop 'allocate leaf'
        nlabl = 0.0
      end if
    end if  

    ! dnlabl = min( nlabl, mydnleaf )
    ! if ( (dclabl - clabl) > 1e-8 ) stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf C'
    ! if ( (dnlabl - nlabl) > 1e-8 ) stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf N'
    ! clabl  = clabl - dclabl
    ! nlabl  = nlabl - dnlabl

  end subroutine allocate_leaf


  subroutine allocate_root( pft, mydcroot, mydnroot, croot, nroot, clabl, nlabl, nignore )
    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant, params_pft_plant, dnup
    use md_nuptake, only: dnup_fix
    use md_params_core, only: eps

    ! arguments
    integer, intent(in) :: pft
    real, intent(in)    :: mydcroot
    real, intent(in)    :: mydnroot
    real, intent(inout) :: croot, nroot
    real, intent(inout) :: clabl, nlabl
    logical, intent(in) :: nignore

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

    if ( clabl < -1.0*eps ) then
      stop 'ALLOCATE_ROOT: trying to remove too much from labile pool: root C'
    else if ( clabl < 0.0 ) then
      ! numerical imprecision
      ! print*,'numerical imprecision?'
      ! stop 'allocate root'
      clabl = 0.0
    end if

    if (nignore) then
      ! If labile N gets negative, account gap as N fixation
      if ( nlabl < 0.0 ) then
        dnup(pft)%n14 = dnup(pft)%n14 - nlabl
        dnup_fix(pft) = dnup_fix(pft) - nlabl
        nlabl = 0.0
      end if
    else
      if ( nlabl < -1.0*eps ) then
        stop 'ALLOCATE_ROOT: trying to remove too much from labile pool: root N'
      else if ( nlabl < 0.0 ) then
        ! numerical imprecision
        ! print*,'numerical imprecision?'
        ! stop 'allocate leaf'
        nlabl = 0.0
      end if
    end if

  end subroutine allocate_root


  function get_rcton_init( pft, meanmppfd, nv ) result( rcton )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    ! Cleaf = c_molmass * params_pft_plant(pft)%r_ctostructn_leaf * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * params_pft_plant(pft)%r_n_cw_v + LAI * params_pft_plant(pft)%ncw_min ]
    ! Nleaf = n_molmass * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * (params_pft_plant(pft)%r_n_cw_v + 1) + LAI * params_pft_plant(pft)%ncw_min ]
    ! linearization around LAI = 0 ==> (1-exp(-k*L)) ~= k*L
    ! ==> Cleaf ~= LAI * c_molmass * params_pft_plant(pft)%r_ctostructn_leaf * ( meanmppfd * kbeer * nv * params_pft_plant(pft)%r_n_cw_v + params_pft_plant(pft)%ncw_min )
    ! ==> Nleaf ~= LAI * n_molmass * ( meanmppfd * kbeer * nv * (params_pft_plant(pft)%r_n_cw_v + 1) + params_pft_plant(pft)%ncw_min )
    ! r_cton = Cleaf / Nleaf
    !----------------------------------------------------------------
    ! use md_params_core, only: nmonth
    use md_plant, only: params_plant, params_pft_plant

    ! arguments
    integer, intent(in)                 :: pft
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real :: rcton

    ! local variables
    real :: maxnv
    real :: tmp1, tmp2, tmp3

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )

    ! tmp1 = c_molmass * params_pft_plant(pft)%r_ctostructn_leaf
    ! tmp2 = maxnv * params_plant%kbeer * params_pft_plant(pft)%r_n_cw_v + params_pft_plant(pft)%ncw_min
    ! tmp3 = n_molmass * ( maxnv * params_plant%kbeer * ( params_pft_plant(pft)%r_n_cw_v + 1.0 ) + params_pft_plant(pft)%ncw_min )
    ! rcton = tmp1 * tmp2 / tmp3

    rcton = ( c_molmass * params_pft_plant(pft)%r_ctostructn_leaf * &
      ( maxnv * params_plant%kbeer * params_pft_plant(pft)%r_n_cw_v + params_pft_plant(pft)%ncw_min ) &
      ) / ( n_molmass * ( maxnv * params_plant%kbeer * ( params_pft_plant(pft)%r_n_cw_v + 1.0 ) + params_pft_plant(pft)%ncw_min ) )

  end function get_rcton_init


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

    return

    999  stop 'INITIO_ALLOCATION: error opening output files'

  end subroutine initio_allocation


  subroutine initoutput_allocation()
    !////////////////////////////////////////////////////////////////
    !  Initialises nuptake-specific output variables
    !----------------------------------------------------------------
    ! xxx remove their day-dimension
    outaCalclm(:,:) = 0.0
    outaNalclm(:,:) = 0.0
    outaCalcrm(:,:) = 0.0
    outaNalcrm(:,:) = 0.0

    ! print*, 'initialising outaCalloc',outaCalloc

  end subroutine initoutput_allocation


  subroutine getout_daily_allocation( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    outaCalclm(:,jpngr) = outaCalclm(:,jpngr) + dcleaf(:) 
    outaNalclm(:,jpngr) = outaNalclm(:,jpngr) + dnleaf(:)
    outaCalcrm(:,jpngr) = outaCalcrm(:,jpngr) + dcroot(:) 
    outaNalcrm(:,jpngr) = outaNalcrm(:,jpngr) + dnroot(:)

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
    itime = real(year) + real(interface%params_siml%firstyeartrend) - real(interface%params_siml%spinupyears)

    ! print*, 'writing time, outaCalloc',itime, sum(outaCalloc(:,jpngr))

    write(350,999) itime, sum(outaCalclm(:,jpngr))
    write(351,999) itime, sum(outaNalclm(:,jpngr))
    write(352,999) itime, sum(outaCalcrm(:,jpngr))
    write(353,999) itime, sum(outaNalcrm(:,jpngr))

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_allocation

end module md_allocation
