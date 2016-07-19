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

  subroutine allocation_daily( jpngr, usedoy, usemoy, dtemp )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant, params_pft_plant, pleaf, proot, &
      plabl, drgrow, lai_ind, nind, canopy, leaftraits, &
      get_canopy, depletionfrac, isgrowing, get_leaftraits, get_leaftraits_init
    use md_waterbal, only: solar
    use md_gpp, only: mlue, mrd_unitiabs, mactnv_unitiabs
    use md_soiltemp, only: dtemp_soil
    use md_ntransform, only: pninorg

    ! xxx debug
    use md_nuptake, only: calc_dnup, outtype_calc_dnup
    use md_waterbal, only: solar, evap
    use md_gpp, only: calc_dgpp, calc_drd
    use md_npp, only: calc_resp_maint, calc_cexu, deactivate_root
    use md_gpp, only: dgpp, drd 
    use md_plant, only: dnpp, drleaf, drroot, dcex, dnup

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: usedoy     ! day of year
    integer, intent(in) :: usemoy     ! month of year
    real,    intent(in) :: dtemp   ! air temperaure, deg C

    ! local variables
    integer :: lu
    integer :: pft
    type( orgpool ) :: reserve
    real, parameter :: freserve = 0.0000

    ! xxx debug
    type( orgpool ) :: bal1, bal2, bald
    real :: eps = 9.999e-8                   ! numerical imprecision allowed in mass conservation tests

    logical, save :: toleaves = .true.       ! boolean determining whether C and N in this time step are allocated to leaves or roots

    ! initialise
    dcleaf(:) = 0.0
    dnleaf(:) = 0.0
    dcroot(:) = 0.0
    dnroot(:) = 0.0
    drgrow(:) = 0.0

    do pft=1,npft

      reserve = orgfrac( freserve, pleaf(pft,jpngr) )
      ! print*,'reserve ', reserve
      ! print*,'plabl   ', plabl
      ! print*,'avl     ', orgminus( plabl(pft,jpngr), reserve )

      lu = params_pft_plant(pft)%lu_category

      if (params_pft_plant(pft)%grass) then

        ! if ( plabl(pft,jpngr)%c%c12>reserve%c%c12 .and. plabl(pft,jpngr)%n%n14>reserve%n%n14 .and. dtemp>0.0 ) then
        if ( plabl(pft,jpngr)%c%c12>0.0 .and. plabl(pft,jpngr)%n%n14>0.0 .and. dtemp>0.0 ) then

          ! ! mass balance test 
          ! print*,'pleaf ', pleaf
          ! print*,'proot ', proot
          ! print*,'plabl ', plabl
          ! print*,'drgrow', drgrow
          ! bal1 = orgplus( pleaf(pft,jpngr), proot(pft,jpngr), plabl(pft,jpngr), orgpool( carbon(drgrow(pft)), nitrogen(0.0) ) )

          !------------------------------------------------------------------
          ! Calculate maximum C allocatable based on current labile pool size.
          ! Maximum is the lower of all labile C and the C to be matched by all labile N,
          ! discounted by the yield factor.
          !------------------------------------------------------------------
          if (pleaf(pft,jpngr)%c%c12==0.0) then
            ! print*, 'Calculating initial C:N ratio'
            ! initial guess based on Taylor approximation of Cleaf and Nleaf function around cleaf=0
            leaftraits(pft)%r_cton_leaf = get_rcton_init( pft, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
            ! print*, 'solar%meanmppfd(:)', solar%meanmppfd(:)  
            ! print*, 'mactnv_unitiabs(pft,:)', mactnv_unitiabs(pft,:)  
            ! print*, 'initial guess: r_cton_leaf(pft,jpngr) ', leaftraits(pft)%r_cton_leaf  
            ! ! stop

            ! leaftraits(pft) = get_leaftraits_init( pft, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
            ! print*, 'initial guess: r_cton_leaf(pft,jpngr) ', leaftraits(pft)%r_cton_leaf  
            ! stop

          end if

          ! mass balance test 
          bal1 = orgplus( pleaf(pft,jpngr), proot(pft,jpngr), plabl(pft,jpngr) )

          if (toleaves) then
            !-------------------------------------------------------------------
            ! LEAF ALLOCATION
            !-------------------------------------------------------------------
            print*,'doy, toleaves', usedoy, toleaves
            ! dcleaf(pft) = min( &
            !   params_plant%growtheff * (plabl(pft,jpngr)%c%c12 - reserve%c%c12), &
            !   (plabl(pft,jpngr)%n%n14 - reserve%n%n14) * leaftraits(pft)%r_cton_leaf &
            !   )

            ! This criterium never completely depletes N pool because C:N decreases when
            ! leaf C and LAI increase. Therefore, the actual N needed is less than what
            ! is determined based on current C:N (leaftraits(pft)%r_cton_leaf).
            dcleaf(pft) = min( &
              params_plant%growtheff * (plabl(pft,jpngr)%c%c12), &
              (plabl(pft,jpngr)%n%n14) * leaftraits(pft)%r_cton_leaf &
              )


            call allocate_leaf( &
              pft, dcleaf(pft), &
              pleaf(pft,jpngr)%c%c12, pleaf(pft,jpngr)%n%n14, &
              plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
              solar%meanmppfd(:), mactnv_unitiabs(pft,:), &
              lai_ind(pft,jpngr), dnleaf(pft) &
              )

            !-------------------------------------------------------------------  
            ! Update leaf traits
            !-------------------------------------------------------------------  
            leaftraits(pft) = get_leaftraits( pft, lai_ind(pft,jpngr), solar%meanmppfd(:), mactnv_unitiabs(pft,:) )

            !-------------------------------------------------------------------  
            ! Update fpc_grid and fapar_ind (not lai_ind)
            !-------------------------------------------------------------------  
            canopy(pft) = get_canopy( lai_ind(pft,jpngr) )

            !-------------------------------------------------------------------  
            ! If C is left in labile pool, then N wasn't sufficient => get more
            ! N by allocating to roots next day.
            !-------------------------------------------------------------------  
            if ( plabl(pft,jpngr)%c%c12 > 0.0 ) toleaves = .false.

          else
            print*,'doy, toleaves', usedoy, toleaves
            !-------------------------------------------------------------------
            ! ROOT ALLOCATION
            !-------------------------------------------------------------------
            ! dcroot(pft) = min( &
            !   params_plant%growtheff * (plabl(pft,jpngr)%c%c12 - reserve%c%c12), &
            !   (plabl(pft,jpngr)%n%n14 - reserve%n%n14) * params_pft_plant(pft)%r_cton_root &
            !   )
            dcroot(pft) = min( &
              params_plant%growtheff * (plabl(pft,jpngr)%c%c12), &
              (plabl(pft,jpngr)%n%n14) * params_pft_plant(pft)%r_cton_root &
              )
            dnroot(pft) = dcroot(pft) * params_pft_plant(pft)%r_ntoc_root

            call allocate_root( &
              pft, dcroot(pft), dnroot(pft), &
              proot(pft,jpngr)%c%c12, proot(pft,jpngr)%n%n14, &
              plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14  &
              )

            !-------------------------------------------------------------------  
            ! If N is left in labile pool, then C wasn't sufficient => get more
            ! C by allocating to leaves next day.
            !-------------------------------------------------------------------  
            if ( plabl(pft,jpngr)%n%n14 > 0.0 ) toleaves = .true.

          end if

          !-------------------------------------------------------------------
          ! GROWTH RESPIRATION, NPP
          !-------------------------------------------------------------------
          ! add growth respiration to autotrophic respiration and substract from NPP
          ! (note that NPP is added to plabl in and growth resp. is implicitly removed
          ! from plabl above)
          drgrow(pft)   = ( 1.0 - params_plant%growtheff ) * ( dcleaf(pft) + dcroot(pft) ) / params_plant%growtheff

          ! mass balance test 
          ! print*,'checking balance ... '
          bal2 = orgplus( pleaf(pft,jpngr), proot(pft,jpngr), plabl(pft,jpngr), orgpool( carbon(drgrow(pft)), nitrogen(0.0) ) ) 
          bald = orgminus( bal2, bal1 )
          ! print*,'pleaf ', pleaf
          ! print*,'proot ', proot
          ! print*,'plabl ', plabl
          ! print*,'drgrow', drgrow
          ! print*,'                 ... ', bald
          if (abs(bald%c%c12)>eps) stop 'balance not satisfied for C'
          if (abs(bald%n%n14)>eps) stop 'balance not satisfied for N'
          ! print*,'... done'

          ! if (usedoy==40) stop 

        end if

      else

        stop 'allocation_daily not implemented for trees'

      end if

    end do

    ! print*, '--- END allocation_daily:'

  end subroutine allocation_daily


  subroutine allocate_leaf( pft, mydcleaf, cleaf, nleaf, clabl, nlabl, meanmppfd, nv, lai, mydnleaf )
    !///////////////////////////////////////////////////////////////////
    ! LEAF ALLOCATION
    ! Sequence of steps:
    ! - increment foliage C pool
    ! - update LAI
    ! - calculate canopy-level foliage N as a function of LAI 
    ! - reduce labile pool by C and N increments
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant, get_leaf_n_canopy, get_lai

    ! arguments
    integer, intent(in)                 :: pft
    real, intent(in)                    :: mydcleaf
    real, intent(inout)                 :: cleaf, nleaf
    real, intent(inout)                 :: clabl, nlabl
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv
    real, intent(out)                   :: lai
    real, intent(out)                   :: mydnleaf

    ! local variables
    real :: nleaf0
    real :: dclabl, dnlabl

    ! find LAI, given new leaf mass. This is necessary to get leaf-N as 
    ! a function of LAI.
    cleaf  = cleaf + mydcleaf

    ! Calculate LAI as a function of leaf C
    lai = get_lai( pft, cleaf, meanmppfd(:), nv(:) )

    ! calculate canopy-level leaf N as a function of LAI
    nleaf0   = nleaf      
    nleaf    = get_leaf_n_canopy( pft, lai, meanmppfd(:), nv(:) )
    mydnleaf = nleaf - nleaf0

    ! xxx DO THIS SANITY CHECK!!! 
    ! ! subtract from labile pool, making sure pool does not get negative
    ! if ( (mydnleaf - nlabl) > 1e-8 ) then
    !   print*,'nlabl    ', nlabl 
    !   print*,'mydnleaf ', mydnleaf
    !   stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf N'
    ! end if

    dclabl = min( clabl, 1.0 / params_plant%growtheff * mydcleaf )
    dnlabl = min( nlabl, mydnleaf )
    if ( (dclabl - clabl) > 1e-8 ) stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf C'
    if ( (dnlabl - nlabl) > 1e-8 ) stop 'ALLOCATE_LEAF: trying to remove too much from labile pool: leaf N'
    clabl  = clabl - dclabl
    nlabl  = nlabl - dnlabl

  end subroutine allocate_leaf


  subroutine allocate_root( pft, mydcroot, mydnroot, croot, nroot, clabl, nlabl )
    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant, params_pft_plant

    ! arguments
    integer, intent(in)         :: pft
    real, intent(in)            :: mydcroot
    real, intent(in)            :: mydnroot
    real, intent(inout)         :: croot, nroot
    real, intent(inout)         :: clabl, nlabl

    ! local variables
    real :: dclabl
    real :: dnlabl

    dclabl = min( clabl, 1.0 / params_plant%growtheff * mydcroot )
    dnlabl = min( nlabl, mydnroot )
    if ( (dnlabl - nlabl) > 1e-8 ) stop 'ALLOCATE ROOT: trying to remove too much from labile pool: root N'
    if ( (dclabl - clabl) > 1e-8 ) stop 'ALLOCATE ROOT: trying to remove too much from labile pool: root C'
    clabl  = clabl - dclabl
    nlabl  = nlabl - dnlabl

    if (mydcroot<0.0) stop 'ALLOCATE ROOT: root allocation neg.: C'
    if (mydnroot<0.0) stop 'ALLOCATE ROOT: root allocation neg.: N'

    croot = croot + mydcroot
    nroot = nroot + mydnroot
  
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
