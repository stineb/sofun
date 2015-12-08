module _allocation
  !////////////////////////////////////////////////////////////////
  ! ALLOCATION MODULE
  ! Contains the "main" subroutine 'allocation_daily' and all 
  ! necessary subroutines for handling input/output, and auxiliary
  ! subroutines.
  ! Every module that implements 'allocation_daily' must contain 
  ! this list of subroutines (names that way).
  !   - allocation_daily
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: npft, nlu, maxgrid, ndaymonth, ndayyear, c_molmass, n_molmass
  use _classdefs

  implicit none
  !------------------------------------------------------------------
  ! Define derived type to store current state variables.
  ! Is defined within all subroutines in this module.
  ! This avoids having to pass them as arguments.
  !------------------------------------------------------------------
  type statetype_eval_imbalance
    type(orgpool)  :: pleaf
    type(orgpool)  :: proot
    type(orgpool)  :: plabl
    real           :: crownarea
    real           :: mlue
    real           :: dppfd
    real           :: mrd_unitiabs
    real           :: nv
    real           :: meanmppfd
    type(nitrogen) :: pninorg
    real           :: pft 
  end type statetype_eval_imbalance

  type statetype_mustbe_zero_for_lai
    real :: meanmppfd
    real :: cleaf
    real :: nv
  end type statetype_mustbe_zero_for_lai

  type(statetype_eval_imbalance)      :: state_eval_imbalance
  type(statetype_mustbe_zero_for_lai) :: state_mustbe_zero_for_lai

  ! MODULE-SPECIFIC VARIABLES
  real, dimension(npft) :: dcleaf
  real, dimension(npft) :: dnleaf
  real, dimension(npft) :: dcroot
  real, dimension(npft) :: dnroot

  ! OUTPUT VARIABLES
  real, dimension(npft,maxgrid) :: outaCalclm
  real, dimension(npft,maxgrid) :: outaNalclm
  real, dimension(npft,maxgrid) :: outaCalcrm
  real, dimension(npft,maxgrid) :: outaNalcrm


contains

  subroutine allocation_daily( jpngr, doy, moy, dom )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use _params_modl, only: r_cton_root, r_ntoc_root, growtheff, grass, &
      lu_category
    use _vars_core, only: pleaf, proot, plabl, pninorg, r_cton_leaf, &
      r_ntoc_leaf, crownarea, drauto, dnpp, lai_ind
    use _gpp, only: mlue, mrd_unitiabs, mactnv_unitiabs
    use _waterbal, only: dppfd, meanmppfd
    use _findroot_fzeroin

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy     ! day of year
    integer, intent(in) :: moy     ! month of year
    integer, intent(in) :: dom     ! day of month

    ! local variables
    integer :: lu
    integer :: pft
    integer :: usemoy
    integer :: usedoy
    real    :: dclabl
    real    :: dnlabl
    logical :: cont          ! true if allocation to leaves (roots) is not 100% and not 0%
    real    :: max_dcleaf_n_constraint
    real    :: max_dcroot_n_constraint
    real    :: max_dc_n_constraint
    real    :: max_dc
    real    :: eval_allleaves
    real    :: eval_allroots
    real    :: abserr
    real    :: relerr
    real    :: nleaf0
    real    :: lai0, lai1
    integer, parameter :: nmax = 100

    type(outtype_zeroin) :: out_zeroin

    integer, save      :: invocation = 0                 ! internally counted simulation year
    integer, parameter :: spinupyr_phaseinit_2 = 1   ! this is unnecessary: might as well do flexible allocation right from the start.
    logical            :: flexalloc
    real, parameter    :: frac_shoot = 0.5

    ! xxx debug
    real    :: test

    !-------------------------------------------------------------------------
    ! Count number of calls (one for each simulation year) and allow flexible
    ! allocation only after year 'spinupyr_phaseinit_2'.
    !-------------------------------------------------------------------------
    if (doy==1) then
      invocation = invocation + 1
      ! write(0,*) 'WARNING: FIXED ALLOCATION'
    end if
    ! if ( invocation > spinupyr_phaseinit_2 ) then
    !   flexalloc = .true.
    ! else
    !   flexalloc = .false.
    ! end if
    flexalloc = .true.

    abserr=100.0*XMACHEPS*10e5
    relerr=1000.0*XMACHEPS*10e5

    do pft=1,npft

      lu = lu_category(pft)

      if (grass(pft)) then

        write(0,*) '--- allocation_daily, doy:',doy
        write(0,*)' pninorg(lu,jpngr)%n14',pninorg
        write(0,*)' plabl(lu,jpngr)',plabl
        ! stop

        if ( plabl(pft,jpngr)%c%c12>0.0 .and. plabl(pft,jpngr)%n%n14>0.0 ) then

          write(0,*) 'growing ...'

          if (flexalloc) then
            !------------------------------------------------------------------
            ! Store state variables for optimisation
            !------------------------------------------------------------------
            ! P-model uses monthly input values, implying jumps in LUC etc. 
            ! Anticipate next day's GPP etc. by using next day's (monthly) LUE.
            if (dom==ndaymonth(moy)) then
              usemoy = moy + 1
            else
              usemoy = moy
            end if
            if (doy==ndayyear) then
              usedoy = 1
            else
              usedoy = doy + 1
            end if

            ! state variables used in function eval_imbalance
            state_eval_imbalance%pleaf        = pleaf(pft,jpngr)
            state_eval_imbalance%proot        = proot(pft,jpngr)
            state_eval_imbalance%plabl        = plabl(pft,jpngr)
            state_eval_imbalance%crownarea    = crownarea(pft,jpngr)
            state_eval_imbalance%mlue         = mlue(usemoy)
            state_eval_imbalance%dppfd        = dppfd(usedoy)
            state_eval_imbalance%mrd_unitiabs = mrd_unitiabs(usemoy)
            state_eval_imbalance%meanmppfd    = meanmppfd(usemoy)
            state_eval_imbalance%pninorg      = pninorg(lu,jpngr)  ! the only data that is not yet available - use today's value 
            state_eval_imbalance%pft          = pft
            state_eval_imbalance%nv           = mactnv_unitiabs(usemoy)

            !------------------------------------------------------------------
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by all labile N,
            ! discounted by the yield factor.
            !------------------------------------------------------------------
            if (pleaf(pft,jpngr)%c%c12.eq.0.0) then
              ! initial guess based on Taylor approximation of Cleaf and Nleaf function around cleaf=0
              r_cton_leaf(pft,jpngr) = get_rcton_init( meanmppfd(usemoy), mactnv_unitiabs(usemoy) )
              ! write(0,*) 'meanmppfd(usemoy)      ', meanmppfd(usemoy)  
              ! write(0,*) 'mactnv_unitiabs(usemoy)', mactnv_unitiabs(usemoy)  
              ! write(0,*) 'r_cton_leaf(pft,jpngr) ', r_cton_leaf(pft,jpngr)  
              ! stop
            end if
            max_dcleaf_n_constraint = plabl(pft,jpngr)%n%n14 * r_cton_leaf(pft,jpngr)
            max_dcroot_n_constraint = plabl(pft,jpngr)%n%n14 * r_cton_root(pft) ! should be obsolete as generally r_ntoc_leaf > r_ntoc_root
            max_dc = min( growtheff * plabl(pft,jpngr)%c%c12, max_dcleaf_n_constraint, max_dcroot_n_constraint )

            ! write(0,*) 'plabl(pft,jpngr)', plabl(pft,jpngr)  
            ! write(0,*) 'r_cton_leaf(pft,jpngr)',r_cton_leaf(pft,jpngr)
            ! write(0,*) 'r_cton_root(pft)',r_cton_root(pft)
            ! write(0,*) 'growtheff', growtheff  
            ! write(0,*) 'max_dcleaf_n_constraint', max_dcleaf_n_constraint  
            ! write(0,*) 'max_dcroot_n_constraint', max_dcroot_n_constraint  
            ! write(0,*) 'max_dc', max_dc  

            ! write(0,*) 'moy',moy
            ! write(0,*) 'lai',lai
            ! write(0,*) 'mlue',mlue(moy) ! ok
            ! write(0,*) 'dppfd',dppfd(doy)  ! ok
            ! write(0,*) 'mrd_unitiabs',mrd_unitiabs(moy) ! ok
            ! write(0,*) 'meanmppfd',meanmppfd(moy) ! ok

            !------------------------------------------------------------------
            ! Optimisation by balanced growth
            ! Test I: Evaluate balance if all is put to roots.
            ! If C:N ratio of return is still greater than whole-plant C:N 
            ! ratio, then put all to roots.
            !------------------------------------------------------------------
            cont = .true.
            write(0,*) 'check alloation: all to roots'
            eval_allroots  = eval_imbalance( 0.0 )
            write(0,*) 'eval_allroots', eval_allroots  
            if (eval_allroots > 0.0) then
              dcleaf(pft) = 0.0
              cont = .false.
              ! write(0,*) '* putting all to roots *'
            end if

            !------------------------------------------------------------------
            ! Test II: Evaluate balance if all is put to leaves.
            ! If C:N ratio of return is still lower than whole-plant C:N ratio, 
            ! then put all to leaves.
            !------------------------------------------------------------------
            if (cont) then
              write(0,*) 'check alloation: all to leaves with dcleaf =', max_dc
              eval_allleaves = eval_imbalance( max_dc )
              write(0,*) 'eval_allleaves', eval_allleaves  
              if (eval_allleaves < 0.0) then
                dcleaf(pft) = max_dc
                cont = .false.
                write(0,*) '* putting all to leaves *'
              end if
            end if

            !------------------------------------------------------------------
            ! Optimum is between 0.0 and max_dc. Find root of function 
            ! 'eval_imbalance()' in the interval [0.0, max_dc].
            !------------------------------------------------------------------
            if (cont) then
              write(0,*) '*** finding root of eval_imbalance ***'
              out_zeroin = zeroin( eval_imbalance, abserr, relerr, nmax, '', 0.0, max_dc )
              if ( out_zeroin%error /= 0 ) then
                dcleaf(pft) = 0.0
              else
                dcleaf(pft) = out_zeroin%root
              end if
              ! write(0,*) 'no. of iterations', out_zeroin%niter
            end if
            ! stop 'do beni'

            !------------------------------------------------------------------
            ! xxx debug: project next-day's fluxes with optimal dcleaf, derived now
            !------------------------------------------------------------------
            ! write(0,*) 'dcleaf(pft)', dcleaf(pft)
            ! test = eval_imbalance( dcleaf(pft), .true. )
            ! stop

            !-------------------------------------------------------------------
            ! LEAF ALLOCATION
            !-------------------------------------------------------------------
            ! find LAI, given new leaf mass. This is necessary to get leaf-N as 
            ! a function of LAI.
            pleaf(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12 + dcleaf(pft)

            ! Calculate LAI as a function of canopy leaf mass
            !-------------------------------------------------------------------
            ! state variables used in function mustbe_zero_for_lai
            lai_ind(pft,jpngr) = get_lai( pleaf(pft,jpngr)%c%c12, meanmppfd(usemoy), mactnv_unitiabs(usemoy) )

            ! write(0,*) 'dcleaf ', dcleaf 
            ! write(0,*) 'cleaf  ', cleaf 
            ! write(0,*) 'lai    ', lai_ind(pft,jpngr)
            ! ! stop

            ! calculate canopy-level leaf N as a function of LAI
            nleaf0 = pleaf(pft,jpngr)%n%n14
            ! write(0,*) 'nleaf0 ', nleaf0 
            pleaf(pft,jpngr)%n%n14 = get_canopy_leaf_n( lai_ind(pft,jpngr), meanmppfd(usemoy), mactnv_unitiabs(usemoy) )
            ! write(0,*) 'nleaf ', nleaf 
            dnleaf(pft) = pleaf(pft,jpngr)%n%n14 - nleaf0
            ! write(0,*) 'dnleaf ', mydnleaf 
            ! if (mydnleaf.gt.0.0) then
            !   write(0,*) 'dcleaf/dnleaf ', dcleaf/mydnleaf 
            ! end if

            ! Update canopy-level C:N ratio in leaves 
            r_cton_leaf(pft,jpngr) = cton( pleaf(pft,jpngr) )
            r_ntoc_leaf(pft,jpngr) = ntoc( pleaf(pft,jpngr) )

            ! subtract from labile pool
            ! to ensure plabl not getting negative (within numerical imprecision)
            dclabl = min( plabl(pft,jpngr)%c%c12, 1.0 / growtheff * dcleaf(pft) )
            dnlabl = min( plabl(pft,jpngr)%n%n14, dnleaf(pft) )

            plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - dclabl
            plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnlabl
            

            !-------------------------------------------------------------------
            ! ROOT ALLOCATION
            !-------------------------------------------------------------------
            ! use remainder for allocation to roots
            dcroot(pft) = min( growtheff * plabl(pft,jpngr)%c%c12, r_cton_root(pft) * plabl(pft,jpngr)%n%n14 )
            dnroot(pft) = dcroot(pft) * r_ntoc_root(pft)

           ! to ensure plabl not getting negative (within numerical imprecision)
            dclabl = min( plabl(pft,jpngr)%c%c12, 1.0 / growtheff * dcroot(pft) )
            dnlabl = min( plabl(pft,jpngr)%n%n14, dnroot(pft) )

            plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - dclabl
            plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnlabl

            proot(pft,jpngr)%c%c12 = proot(pft,jpngr)%c%c12 + dcroot(pft)
            proot(pft,jpngr)%n%n14 = proot(pft,jpngr)%n%n14 + dnroot(pft)

            ! if (present(verbose)) then
            !   write(0,*) 'after presumed allocation:'
            !   write(0,*) 'clabl', clabl ! OK
            !   write(0,*) 'nlabl', nlabl ! OK
            !   ! write(0,*) 'dcleaf', dcleaf
            !   ! write(0,*) 'mydcroot', mydcroot
            !   ! write(0,*) 'mydnleaf', mydnleaf
            !   ! write(0,*) 'mydnroot', mydnroot

            !   ! write(0,*) 'after presumed allocation:'
            !   ! write(0,*) 'B cleaf', cleaf
            !   ! write(0,*) 'B root', root
            ! end if
        
          else
            !------------------------------------------------------------------
            ! Prescribe constant root:shoot ratio (in terms of C mass).
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by 
            ! all labile N, discounted by the yield factor.
            !------------------------------------------------------------------
            if (pleaf(pft,jpngr)%c%c12.eq.0.0) then
              ! initial guess based on Taylor approximation of Cleaf and Nleaf function around cleaf=0
              r_cton_leaf(pft,jpngr) = get_rcton_init( meanmppfd(usemoy), mactnv_unitiabs(usemoy) )
              r_ntoc_leaf(pft,jpngr) = 1.0 / r_cton_leaf(pft,jpngr)
            end if
            max_dc              = growtheff * plabl(pft,jpngr)%c%c12
            max_dc_n_constraint = plabl(pft,jpngr)%n%n14 / ( frac_shoot * r_ntoc_leaf(pft,jpngr) + ( 1.0 - frac_shoot) * r_ntoc_root(pft) )
            ! write(0,*) "clabl, nlabl", plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14 
            ! write(0,*) "max_dc, max_dc_n_constraint", max_dc, max_dc_n_constraint
            max_dc = min( max_dc, max_dc_n_constraint )

            !-------------------------------------------------------------------
            ! LEAF ALLOCATION
            !-------------------------------------------------------------------
            ! allocation to leaves is prescribed
            dcleaf(pft) = max_dc * frac_shoot
            pleaf(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12 + dcleaf(pft)

            ! Calculate LAI as a function of canopy leaf mass
            !-------------------------------------------------------------------
            lai_ind(pft,jpngr) = get_lai( pleaf(pft,jpngr)%c%c12, meanmppfd(usemoy), mactnv_unitiabs(usemoy) )
            stop 'shouldnt be here'
            ! write(0,*) 'dcleaf ', dcleaf 
            ! write(0,*) 'cleaf  ', cleaf 
            ! write(0,*) 'lai    ', lai_ind(pft,jpngr)
            ! ! stop

            ! calculate canopy-level leaf N as a function of LAI
            nleaf0 = pleaf(pft,jpngr)%n%n14
            ! write(0,*) 'nleaf0 ', nleaf0 
            pleaf(pft,jpngr)%n%n14 = get_canopy_leaf_n( lai_ind(pft,jpngr), meanmppfd(usemoy), mactnv_unitiabs(usemoy) )
            ! write(0,*) 'nleaf ', nleaf 
            dnleaf(pft) = pleaf(pft,jpngr)%n%n14 - nleaf0
            ! write(0,*) 'dnleaf ', mydnleaf 
            ! if (mydnleaf.gt.0.0) then
            !   write(0,*) 'dcleaf/dnleaf ', dcleaf/mydnleaf 
            ! end if

            ! Update canopy-level C:N ratio in leaves 
            r_cton_leaf(pft,jpngr) = cton( pleaf(pft,jpngr) )
            r_ntoc_leaf(pft,jpngr) = ntoc( pleaf(pft,jpngr) )

            ! subtract from labile pool
            ! to ensure plabl not getting negative (within numerical imprecision)
            dclabl = min( plabl(pft,jpngr)%c%c12, 1.0 / growtheff * dcleaf(pft) )
            dnlabl = min( plabl(pft,jpngr)%n%n14, dnleaf(pft) )

            plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - dclabl
            plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnlabl
          
            !-------------------------------------------------------------------
            ! ROOT ALLOCATION
            !-------------------------------------------------------------------
            dcroot(pft) = max_dc * ( 1.0 - frac_shoot )
            dnroot(pft) = dcroot(pft) * r_ntoc_root(pft)

            ! to ensure plabl not getting negative (within numerical imprecision)
            proot(pft,jpngr)%c%c12 = proot(pft,jpngr)%c%c12 + dcroot(pft)
            proot(pft,jpngr)%n%n14 = proot(pft,jpngr)%n%n14 + dnroot(pft)

            dclabl = min( plabl(pft,jpngr)%c%c12, 1.0 / growtheff * dcroot(pft) )
            dnlabl = min( plabl(pft,jpngr)%n%n14, dnroot(pft) )

            plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - dclabl
            plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnlabl

            ! add growth respiration to autotrophic respiration and substract from NPP
            ! (note that NPP is added to plabl in and growth resp. is implicitly removed
            ! from plabl above)
            drauto(pft)   = drauto(pft)     + ( 1.0 - growtheff ) * ( dcleaf(pft) + dcroot(pft) )
            dnpp(pft)%c12 = dnpp(pft)%c12   - ( 1.0 - growtheff ) * ( dcleaf(pft) + dcroot(pft) ) 

          end if

          !-------------------------------------------------------------------
          ! GROWTH RESPIRATION
          !-------------------------------------------------------------------
          ! add growth respiration to autotrophic respiration and substract from NPP
          ! (note that NPP is added to plabl in and growth resp. is implicitly removed
          ! from plabl above)
          drauto(pft)   = drauto(pft)     + ( 1.0 - growtheff ) * ( dcleaf(pft) + dcroot(pft) )
          dnpp(pft)%c12 = dnpp(pft)%c12   - ( 1.0 - growtheff ) * ( dcleaf(pft) + dcroot(pft) ) 

        else

          dcleaf(pft) = 0.0
          dcroot(pft) = 0.0
          dnleaf(pft) = 0.0
          dnroot(pft) = 0.0

          ! write(0,*) 'not growing ...'

        end if

      else

        stop 'allocation_daily not implemented for trees'

      end if

    end do

    ! test_calloc = test_calloc + dcleaf + dcroot
    ! write(0,*) 'test_calloc', test_calloc
    
    ! write(0,*) '--- END allocation_daily:'

  end subroutine allocation_daily


  function eval_imbalance( mydcleaf, verbose ) result ( eval )
    !/////////////////////////////////////////////////////////
    ! Evaluates C:N ratio of new assimilation after allocation 
    ! versus whole-plant C:N ratio after allocation. Optimal 
    ! allocation is where the two are equal. 
    !---------------------------------------------------------
    use _classdefs, only: orgpool, nitrogen
    use _params_modl, only: growtheff, k_decay_leaf, k_decay_root,&
     r_cton_root, r_ntoc_root
    use _gpp, only: calc_dgpp, calc_drd, r_n_cw_v, ncw_min, r_ctostructn_leaf
    use _nuptake, only: calc_dnup, exurate
    use _npp, only: calc_resp_maint, r_root
    use _vegdynamics, only: get_fapar
    use _findroot_fzeroin

    ! arguments
    real, intent(in)              :: mydcleaf
    logical, intent(in), optional :: verbose

    ! function return variable
    real, intent(out) :: eval

    ! local variables
    real    :: fapar
    real    :: r_ntoc_leaf
    real    :: crownarea
    real    :: mlue
    real    :: dppfd
    real    :: mrd_unitiabs
    real    :: meanmppfd
    real    :: ninorg
    real    :: nv
    integer :: pft

    real :: mydcroot
    real :: mydnleaf
    real :: mydnroot
    real :: cleaf
    real :: nleaf
    real :: croot
    real :: nroot
    real :: clabl
    real :: nlabl
    real :: lai
    real :: gpp
    real :: rd
    real :: mresp_root
    real :: exu
    real :: dc
    real :: dn
    real :: kcleaf
    real :: knleaf
    real :: kcroot
    real :: knroot

    real :: nleaf0
    real :: lai0, lai1

    type(outtype_zeroin) :: out_zeroin

    ! write(0,*) '--- in eval_imbalance with mydcleaf=',mydcleaf

    ! Copy to local variables for shorter writing
    cleaf        = state_eval_imbalance%pleaf%c%c12
    nleaf        = state_eval_imbalance%pleaf%n%n14
    croot        = state_eval_imbalance%proot%c%c12
    nroot        = state_eval_imbalance%proot%n%n14
    clabl        = state_eval_imbalance%plabl%c%c12
    nlabl        = state_eval_imbalance%plabl%n%n14
    crownarea    = state_eval_imbalance%crownarea
    mlue         = state_eval_imbalance%mlue         
    dppfd        = state_eval_imbalance%dppfd        
    mrd_unitiabs = state_eval_imbalance%mrd_unitiabs 
    meanmppfd    = state_eval_imbalance%meanmppfd    
    ninorg       = state_eval_imbalance%pninorg%n14 
    pft          = state_eval_imbalance%pft
    nv           = state_eval_imbalance%nv

    ! write(0,*) 'before presumed allocation:'
    ! write(0,*) 'cleaf', cleaf
    ! write(0,*) 'nleaf', nleaf
    ! write(0,*) 'croot', croot
    ! write(0,*) 'nroot', nroot
    ! write(0,*) 'clabl', clabl
    ! write(0,*) 'nlabl', nlabl

    !-------------------------------------------------------------------
    ! LEAF ALLOCATION
    !-------------------------------------------------------------------
    ! find LAI, given new leaf mass. This is necessary to get leaf-N as 
    ! a function of LAI.
    cleaf  = cleaf + mydcleaf

    ! Calculate LAI as a function of leaf C
    lai = get_lai( cleaf, meanmppfd, nv )

    ! write(0,*) 'dcleaf ', mydcleaf 
    ! write(0,*) 'cleaf  ', cleaf 
    ! write(0,*) 'lai    ', lai
    ! ! stop

    ! calculate canopy-level leaf N as a function of LAI
    nleaf0   = nleaf
    write(0,*) 'nleaf0 ', nleaf0
    nleaf    = get_canopy_leaf_n( lai, meanmppfd, nv )
    write(0,*) 'nleaf ', nleaf 
    mydnleaf = nleaf - nleaf0
    write(0,*) 'dnleaf ', mydnleaf 
    ! if (mydnleaf.gt.0.0) then
    !   write(0,*) 'dcleaf/dnleaf ', mydcleaf/mydnleaf 
    ! end if

    ! subtract from labile pool
    clabl  = clabl - 1.0 / growtheff * mydcleaf
    nlabl  = nlabl - mydnleaf
    ! write(0,*) 'clabl  ', clabl 
    ! write(0,*) 'nlabl  ', nlabl 
    ! write(0,*) 'r_ntoc_root(pft)  ', r_ntoc_root(pft)
    
    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    ! use remainder for allocation to roots
    mydcroot = min( growtheff * clabl, r_cton_root(pft) * nlabl )
    mydnroot = mydcroot * r_ntoc_root(pft)

    ! write(0,*) 'dcroot ', mydcroot 
    ! write(0,*) 'dnroot ', mydnroot 
    if (mydcroot.lt.0.0) stop 'root allocation neg.: C'
    if (mydnroot.lt.0.0) stop 'root allocation neg.: N'

    clabl    = clabl - 1.0 / growtheff * mydcroot
    nlabl    = nlabl - mydnroot
    croot    = croot + mydcroot
    nroot    = nroot + mydnroot

    ! if (present(verbose)) then
    !   write(0,*) 'after presumed allocation:'
    !   write(0,*) 'clabl', clabl ! OK
    !   write(0,*) 'nlabl', nlabl ! OK
    !   ! write(0,*) 'mydcleaf', mydcleaf
    !   ! write(0,*) 'mydcroot', mydcroot
    !   ! write(0,*) 'mydnleaf', mydnleaf
    !   ! write(0,*) 'mydnroot', mydnroot

    !   ! write(0,*) 'after presumed allocation:'
    !   ! write(0,*) 'B cleaf', cleaf
    !   ! write(0,*) 'B root', root
    ! end if

    !-------------------------------------------------------------------
    ! PROJECT NEXT DAY'S C AND N BALANCE:
    ! decay, GPP, respiration, N uptake
    !-------------------------------------------------------------------
    ! continuous decay
    kcleaf = cleaf * k_decay_leaf(pft)
    knleaf = nleaf * k_decay_leaf(pft)
    kcroot = croot * k_decay_root(pft)
    knroot = nroot * k_decay_root(pft)
    
    cleaf  = cleaf - kcleaf
    nleaf  = nleaf - knleaf
    croot  = croot - kcroot
    nroot  = nroot - knroot

    ! Calculate next day's C and N return after assumed allocation (tissue turnover happens before!)
    fapar = get_fapar( lai )

    ! if (present(verbose)) then
    !   write(0,*) 'cleaf',cleaf
    !   write(0,*) 'croot',croot
    !   write(0,*) 'fapar',fapar ! nicht ganz ok!
    !   write(0,*) 'mlue',mlue ! ok
    !   write(0,*) 'dppfd',dppfd  ! ok
    !   write(0,*) 'ninorg',ninorg
    !   ! write(0,*) 'mrd_unitiabs',mrd_unitiabs ! ok
    !   ! write(0,*) 'meanmppfd',meanmppfd ! ok
    !   ! write(0,*) 'exurate',exurate
    ! end if

    gpp        = calc_dgpp( fapar, mlue, dppfd ) ! XXX this gives incorrect results (next day's GPP is much smaller!)
    rd         = calc_drd(  fapar, mrd_unitiabs, meanmppfd )
    exu        = croot * exurate
    mresp_root = calc_resp_maint( croot, r_root )
    
    dc         = gpp - rd - mresp_root - exu
    dn         = calc_dnup( exu, ninorg )

    ! if (present(verbose)) then
      ! write(0,*) 'gpp                ', gpp
      ! write(0,*) 'rd                 ', rd
      ! write(0,*) 'exu                ', exu
      ! write(0,*) 'mresp_root         ', mresp_root
      ! write(0,*) 'dc                 ', dc
      ! write(0,*) 'dn                 ', dn
    !   write(0,*) 'current   plant C:N', ( cleaf + croot ) / ( nleaf + nroot )
    !   write(0,*) 'tomorrows plant C:N', growtheff * (dc + clabl) / (dn + nlabl)
    ! end if

    !-------------------------------------------------------------------
    ! EVALUATION QUANTITY - IS MINIMISED BY OPTIMISATION
    !-------------------------------------------------------------------
    ! Evaluation quantity is the difference between the 
    ! C:N ratio of new assimilates and the C:N ratio 
    ! of the whole plant after allocation.
    eval = growtheff * (dc + clabl) / (dn + nlabl) - ( cleaf + croot ) / ( nleaf + nroot )

    ! write(0,*) 'new assimilates stoichiometry: ', growtheff * (dc + clabl) / (dn + nlabl)
    ! write(0,*) 'current stoichiometry:         ', ( cleaf + croot ) / ( nleaf + nroot )
    ! write(0,*) 'eval                           ', eval
    ! stop

  end function eval_imbalance


  function get_lai( cleaf, meanmppfd, nv ) result( lai )
    !////////////////////////////////////////////////////////////////
    ! Calculates LAI as a function of canopy-level leaf-C:
    ! Cleaf = Mc * c * ( I0 * ( 1 - exp( -kL ) * nv * b + L * a ) )
    ! Cannot be solved analytically for L = f(Cleaf). Therefore, 
    ! numerical root-searching algorithm is applied so that
    ! Cleaf / ( Mc * c ) - ( I0 * ( 1 - exp( -kL ) * nv * b + L * a ) ) = 0
    ! This is implemented in function 'mustbe_zero_for_lai()'.
    !----------------------------------------------------------------
    use _findroot_fzeroin

    ! arguments 
    real, intent(in) :: cleaf
    real, intent(in) :: meanmppfd
    real, intent(in) :: nv 

    ! local variables
    real                 :: abserr
    real                 :: relerr
    real, parameter      :: lai0 = 0.0
    real, parameter      :: lai1 = 20.0
    integer, parameter :: nmax = 100
    type(outtype_zeroin) :: out_zeroin

    ! xxx debug
    real :: test

    ! function return value
    real, intent(out) :: lai

    ! xxx debug
    cleaf = 200.0

    write(0,*) 'what is LAI for Cleaf=', cleaf, '?'

    write(0,*) 'meanmppfd', meanmppfd
    write(0,*) 'nv', nv

    ! Update state. This derived-type variable is "global" within this module
    state_mustbe_zero_for_lai%cleaf     = cleaf
    state_mustbe_zero_for_lai%meanmppfd = meanmppfd
    state_mustbe_zero_for_lai%nv        = nv

    test = mustbe_zero_for_lai( 0.0 ) 
    write(0,*) 'LAI =  0 =>', test

    test = mustbe_zero_for_lai( 20.0 ) 
    write(0,*) 'LAI = 20 =>', test 

    ! call function zeroin to find root (value of LAI for which evaluation expression is zero)
    abserr=100.0*XMACHEPS*10e5
    relerr=1000.0*XMACHEPS*10e5

    write(0,*) 'abserr', abserr
    write(0,*) 'relerr', relerr
    ! stop 'here'

    write(0,*) '*** finding root of mustbe_zero_for_lai ***'
    out_zeroin = zeroin( mustbe_zero_for_lai, abserr, relerr, nmax, "fzeroin.log", 0.0, 20.0 )
    if ( out_zeroin%error /= 0 ) then
      lai = 0.0
    else
      lai = out_zeroin%root
    end if

    write(0,*) 'out_zeroin', out_zeroin
    write(0,*) 'cleaf', cleaf
    write(0,*) 'lai', lai
    stop

  end function get_lai


  function mustbe_zero_for_lai( mylai ) result( mustbe_zero )
    !/////////////////////////////////////////////////////////
    ! Returns LAI for a given leaf-C. Cannot solve this 
    ! analytically because of inhomogenous equation.
    !---------------------------------------------------------
    use _params_modl, only: kbeer
    use _gpp, only: r_n_cw_v, ncw_min, r_ctostructn_leaf
    ! use _vegdynamics, only: get_fapar

    ! arguments
    real, intent(in) :: mylai

    ! function return value
    real, intent(out) :: mustbe_zero

    ! local variables
    real :: meanmppfd
    real :: cleaf
    real :: nv

    ! Read from updated state. This derived-type variable is "global" within this module
    meanmppfd = state_mustbe_zero_for_lai%meanmppfd
    cleaf     = state_mustbe_zero_for_lai%cleaf
    nv        = state_mustbe_zero_for_lai%nv

    write(0,*) '----------'
    write(0,*) 'inside mustbe_zero_for_lai: '
    write(0,*) 'mylai', mylai
    write(0,*) 'cleaf', cleaf
    ! write(0,*) 'meanmppfd', meanmppfd
    ! write(0,*) 'nv', nv

    !---------------------------------------------------------
    ! Cleaf = c_molmass * r_ctostructn_leaf * N_canop_cellwall
    ! N_canop_cellwall = LAI * ncw_min + nv * Iabs * r_n_cw_v
    ! Iabs = meanmppfd * (1-exp( -kbeer * LAI))
    !---------------------------------------------------------
    mustbe_zero = cleaf / ( c_molmass * r_ctostructn_leaf ) - meanmppfd * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * nv * r_n_cw_v - mylai * ncw_min

    write(0,*) 'mustbe_zero', mustbe_zero


  end function mustbe_zero_for_lai


  function get_rcton_init( meanmppfd, nv ) result( rcton )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    !----------------------------------------------------------------
    use _params_modl, only: kbeer
    use _gpp, only: r_n_cw_v, ncw_min, r_ctostructn_leaf

    ! arguments
    real, intent(in) :: meanmppfd
    real, intent(in) :: nv

    ! function return variable
    real, intent(out) :: rcton

    rcton = ( c_molmass * r_ctostructn_leaf * ( meanmppfd * kbeer * nv * r_n_cw_v + ncw_min )) / ( n_molmass * ( meanmppfd * kbeer * nv * (r_n_cw_v + 1) + ncw_min ) )

  end function get_rcton_init


  function get_canopy_leaf_n( mylai, meanmppfd, nv ) result( mynleaf )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    !----------------------------------------------------------------
    use _params_modl, only: kbeer
    use _gpp, only: r_n_cw_v, ncw_min
    use _vegdynamics, only: get_fapar

    ! arguments
    real, intent(in) :: mylai
    real, intent(in) :: meanmppfd
    real, intent(in) :: nv

    ! function return variable
    real, intent(out) :: mynleaf

    mynleaf = n_molmass * ( meanmppfd * get_fapar( mylai ) * nv * ( 1.0 + r_n_cw_v ) + mylai * ncw_min )

  end function get_canopy_leaf_n


  subroutine initio_allocation( )
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use _params_siml, only: runname

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(runname)

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


  ! subroutine getpar_modl_allocation()
  !   !////////////////////////////////////////////////////////////////
  !   ! Subroutine reads nuptake module-specific parameters 
  !   ! from input file
  !   !----------------------------------------------------------------
  !   use _sofunutils, only: getreal

  ! end subroutine getpar_modl_allocation


  subroutine initoutput_allocation()
    !////////////////////////////////////////////////////////////////
    !  Initialises nuptake-specific output variables
    !----------------------------------------------------------------
    ! xxx remove their day-dimension
    outaCalclm(:,:) = 0.0
    outaNalclm(:,:) = 0.0
    outaCalcrm(:,:) = 0.0
    outaNalcrm(:,:) = 0.0

    ! write(0,*) 'initialising outaCalloc',outaCalloc

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

    ! write(0,*) 'collecting outaCalloc',outaCalloc

  end subroutine getout_daily_allocation


  subroutine writeout_ascii_allocation( year, spinup )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use _params_siml, only: firstyeartrend, spinupyears

    ! arguments
    integer, intent(in) :: year       ! simulation year
    logical, intent(in) :: spinup     ! true during spinup years

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
    itime = real(year) + real(firstyeartrend) - real(spinupyears)

    ! write(0,*) 'writing time, outaCalloc',itime, sum(outaCalloc(:,jpngr))

    write(350,999) itime, sum(outaCalclm(:,jpngr))
    write(351,999) itime, sum(outaNalclm(:,jpngr))
    write(352,999) itime, sum(outaCalcrm(:,jpngr))
    write(353,999) itime, sum(outaNalcrm(:,jpngr))

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_allocation

end module _allocation
