module md_allocation
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

  !-----------------------------------------------------------------------
  ! Uncertain (unknown) parameters. Runtime read-in
  !-----------------------------------------------------------------------
  type paramstype_alloc
    ! Parameters determining the relationship of structural N and C to metabolic N
    ! From regressing Narea to metabolic Narea in Hikosaka data
    ! real, parameter    :: r_n_cw_v = 1.23223            ! slope in the relationship of non-metabolic versus metabolic N per leaf area
    ! real, parameter    :: ncw_min = 0.056               ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area
    ! real :: r_n_cw_v = 1.23223            ! slope in the relationship of non-metabolic versus metabolic N per leaf area
    real :: ncw_min  = 0.056              ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area

    real :: r_n_cw_v = 0.4            ! slope in the relationship of non-metabolic versus metabolic N per leaf area
    ! real :: ncw_min           = 0.1       ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area
    
    real :: r_ctostructn_leaf = 50        ! constant ratio of C to structural N
  
  end type paramstype_alloc

  type( paramstype_alloc ) :: params_alloc

  real, parameter :: frac_shoot = 0.5

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC, PRIVATE VARIABLES
  !----------------------------------------------------------------
  type statetype_eval_imbalance
    type(orgpool)           :: pleaf
    type(orgpool)           :: proot
    type(orgpool)           :: plabl
    real                    :: usepft 
    integer                 :: usemoy
    integer                 :: usedoy
    integer                 :: usejpngr
    real                    :: airtemp
  end type statetype_eval_imbalance

  type(statetype_eval_imbalance)      :: state_eval_imbalance

  type statetype_mustbe_zero_for_lai
    real :: cleaf
    real :: maxnv
  end type statetype_mustbe_zero_for_lai

  ! states area global within module (instead of being passed on as arguments)
  type(statetype_mustbe_zero_for_lai) :: state_mustbe_zero_for_lai

  logical, parameter :: write_logfile_eval_imbalance = .false.
  real :: test

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
      plabl, drgrow, dnup, lai_ind, nind, canopy, leaftraits, &
      get_canopy, break_after_alloc
    use md_waterbal, only: solar
    use md_gpp, only: mlue, mrd_unitiabs, mactnv_unitiabs
    use md_phenology, only: dtphen
    use md_findroot_fzeroin
    use md_soiltemp, only: dtemp_soil
    use md_ntransform, only: pninorg

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: usedoy     ! day of year
    integer, intent(in) :: usemoy     ! month of year
    real,    intent(in) :: dtemp   ! air temperaure, deg C

    ! local variables
    integer :: lu
    integer :: pft
    logical :: cont          ! true if allocation to leaves (roots) is not 100% and not 0%
    real    :: max_dcleaf_n_constraint
    real    :: max_dcroot_n_constraint
    real    :: max_dc_buffr_constraint
    real    :: max_dc_n_constraint
    real    :: max_dc
    real    :: min_dc
    real    :: eval_allleaves
    real    :: eval_allroots
    real    :: abserr
    real    :: relerr
    real    :: nleaf0
    real    :: lai0, lai1
    integer, parameter :: nmax = 100

    type(outtype_zeroin)  :: out_zeroin

    integer, save      :: invocation = 0             ! internally counted simulation year
    integer, parameter :: spinupyr_phaseinit_2 = 1   ! this is unnecessary: might as well do flexible allocation right from the start.

    ! xxx verbose
    logical, parameter :: verbose = .false.

    abserr=100.0*XMACHEPS !*10e5
    relerr=1000.0*XMACHEPS !*10e5

    !-------------------------------------------------------------------------
    ! Count number of calls (one for each simulation year) and allow flexible
    ! allocation only after year 'spinupyr_phaseinit_2'.
    !-------------------------------------------------------------------------
    ! if (doy==1) then
    !   ! invocation = invocation + 1
    !   ! print*, 'WARNING: FIXED ALLOCATION'
    ! end if

    do pft=1,npft

      lu = params_pft_plant(pft)%lu_category

      if (params_pft_plant(pft)%grass) then

        ! print*, '--- allocation_daily, doy:',doy
        ! print*, 'doy, plabl  ', doy, plabl(pft,jpngr)
        ! print*, 'doy, dtphen ', dtphen(:,pft)
        ! stop 'do beni'

        ! print*, 'cleaf          ', pleaf(pft,jpngr)%c%c12
        ! print*, 'croot          ', proot(pft,jpngr)%c%c12
        ! print*, 'C:N in leaves  ', cton( pleaf(pft,jpngr), default=0.0 )


        if ( plabl(pft,jpngr)%c%c12>0.0 .and. plabl(pft,jpngr)%n%n14>0.0 .and. dtemp>0.0 ) then

          !------------------------------------------------------------------
          ! Store state variables for optimisation
          !------------------------------------------------------------------
          ! ! P-model uses monthly input values, implying jumps in LUC etc. 
          ! ! Anticipate next day's GPP etc. by using next day's (monthly) LUE.
          ! if (dom==ndaymonth(moy)) then
          !   usemoy = moy + 1
          !   if (usemoy==13) usemoy = 1
          ! else
          !   usemoy = moy
          ! end if
          ! if (doy==ndayyear) then
          !   usedoy = 1
          ! else
          !   usedoy = doy + 1
          ! end if

          ! ! state variables used in function eval_imbalance
          ! state_eval_imbalance%pleaf        = pleaf(pft,jpngr)
          ! state_eval_imbalance%proot        = proot(pft,jpngr)
          ! state_eval_imbalance%plabl        = plabl(pft,jpngr)
          ! state_eval_imbalance%nind         = nind(pft,jpngr)
          ! state_eval_imbalance%mlue         = mlue(pft,usemoy)
          ! state_eval_imbalance%dppfd        = solar%dppfd(usedoy)
          ! state_eval_imbalance%mrd_unitiabs = mrd_unitiabs(pft,usemoy)
          ! state_eval_imbalance%meanmppfd(:) = solar%meanmppfd(:)
          ! state_eval_imbalance%pninorg      = pninorg(lu,jpngr)  ! the only data that is not yet available - use today's value 
          ! state_eval_imbalance%pft          = pft
          ! state_eval_imbalance%nv(:)        = mactnv_unitiabs(pft,:)
          ! state_eval_imbalance%usemoy       = usemoy
          ! state_eval_imbalance%soiltemp     = dtemp_soil(lu,jpngr)
          ! state_eval_imbalance%airtemp      = dtemp

          ! state variables used in function eval_imbalance
          state_eval_imbalance%pleaf    = pleaf(pft,jpngr)
          state_eval_imbalance%proot    = proot(pft,jpngr)
          state_eval_imbalance%plabl    = plabl(pft,jpngr)
          state_eval_imbalance%usepft   = pft
          state_eval_imbalance%usemoy   = usemoy
          state_eval_imbalance%usedoy   = usedoy
          state_eval_imbalance%usejpngr = jpngr
          state_eval_imbalance%airtemp  = dtemp

          !------------------------------------------------------------------
          ! Calculate maximum C allocatable based on current labile pool size.
          ! Maximum is the lower of all labile C and the C to be matched by all labile N,
          ! discounted by the yield factor.
          !------------------------------------------------------------------
          if (pleaf(pft,jpngr)%c%c12==0.0) then
            ! print*, 'Calculating initial C:N ratio'
            ! initial guess based on Taylor approximation of Cleaf and Nleaf function around cleaf=0
            leaftraits(pft)%r_cton_leaf = get_rcton_init( solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
            ! print*, 'solar%meanmppfd(:)', solar%meanmppfd(:)  
            ! print*, 'mactnv_unitiabs(pft,:)', mactnv_unitiabs(pft,:)  
            ! print*, 'initial guess: r_cton_leaf(pft,jpngr) ', leaftraits(pft)%r_cton_leaf  
            ! stop
          end if

          max_dcleaf_n_constraint = plabl(pft,jpngr)%n%n14 * leaftraits(pft)%r_cton_leaf
          max_dcroot_n_constraint = plabl(pft,jpngr)%n%n14 * params_pft_plant(pft)%r_cton_root ! should be obsolete as generally r_ntoc_leaf > r_ntoc_root
          max_dc_buffr_constraint = max( 0.0, plabl(pft,jpngr)%c%c12 - ( params_plant%r_root + params_plant%exurate ) * proot(pft,jpngr)%c%c12 )
          ! if (plabl(pft,jpngr)%c%c12<( params_plant%r_root + params_plant%exurate ) * proot(pft,jpngr)%c%c12 ) stop 'not enough clabl'
          ! max_dc = min( params_plant%growtheff * plabl(pft,jpngr)%c%c12, max_dcleaf_n_constraint, max_dcroot_n_constraint )
          max_dc = min( params_plant%growtheff * max_dc_buffr_constraint, max_dcleaf_n_constraint, max_dcroot_n_constraint )
          min_dc = 0.0

          ! print*, 'plabl(pft,jpngr)', plabl(pft,jpngr)  
          ! print*, 'available after buffer c. ', max_dc_buffr_constraint
          ! print*, 'available for allocation: ', max_dc

          ! print*, 'r_cton_leaf(pft,jpngr)',r_cton_leaf(pft,jpngr)
          ! print*, 'r_cton_root(pft)',r_cton_root(pft)
          ! print*, 'growtheff', growtheff  
          ! print*, 'max_dcleaf_n_constraint', max_dcleaf_n_constraint  
          ! print*, 'max_dcroot_n_constraint', max_dcroot_n_constraint  
          ! print*, 'max_dc', max_dc  
          ! stop

          ! print*, 'moy',moy
          ! print*, 'lai',lai
          ! print*, 'mlue',mlue(moy) ! ok
          ! print*, 'dppfd',dppfd(doy)  ! ok
          ! print*, 'mrd_unitiabs',mrd_unitiabs(moy) ! ok
          ! print*, 'meanmppfd',meanmppfd(moy) ! ok

          ! print*, 'calculating LAI for Cleaf = 2.857124'
          ! test = get_lai( 2.857124, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
          ! print*, 'lai                       =', test
          ! print*, '----------------------------------'

          ! print*, 'calculating LAI for Cleaf = 100'
          ! test = get_lai( 100.0, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
          ! print*, 'lai                       =', test
          ! print*, '----------------------------------'

          ! print*, 'calculating LAI for Cleaf = 2.857124'
          ! test = get_lai( 2.857124, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
          ! print*, 'lai                       =', test
          ! print*, '----------------------------------'

          ! print*, 'calculating LAI for Cleaf = 0.5'
          ! test = get_lai( 0.5  , solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
          ! print*, 'lai                       =', test
          ! print*, '----------------------------------'

          ! print*, 'calculating LAI for Cleaf = 100'
          ! test = get_lai( 100.0, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )
          ! print*, 'lai                       =', test
          ! ! print*, '----------------------------------'

          ! print*, 'calculating eval_imbalance for dc = 0.0'
          ! test = eval_imbalance( 0.0 )
          ! print*, 'eval                              =', test
          ! print*, '----------------------------------'

          ! print*, 'calculating eval_imbalance for dc = ', 2.7
          ! test = eval_imbalance( 2.7 )
          ! print*, 'eval                              =', test
          ! print*, '----------------------------------'

          ! print*, 'calculating eval_imbalance for dc = ', 2.75
          ! test = eval_imbalance( 2.75 )
          ! print*, 'eval                              =', test
          ! print*, '----------------------------------'

          ! print*, 'calculating eval_imbalance for dc = ', 2.8
          ! test = eval_imbalance( 2.8 )
          ! print*, 'eval                              =', test
          ! print*, '----------------------------------'
          
          ! print*, 'calculating eval_imbalance for dc = ', 2.794742 
          ! test = eval_imbalance( 2.794742 )
          ! print*, 'eval                              =', test
          ! print*, '----------------------------------'

          ! print*, 'calculating eval_imbalance for dc = 0.0'
          ! test = eval_imbalance( 0.0 )
          ! print*, 'eval                              =', test
          ! print*, '----------------------------------'

          !------------------------------------------------------------------
          ! Optimisation by balanced growth
          ! Test I: Evaluate balance if all is put to roots.
          ! If C:N ratio of return is still greater than whole-plant C:N 
          ! ratio, then put all to roots.
          !------------------------------------------------------------------
          cont = .true.
          if (verbose) print*, 'check alloation: all to roots'
          eval_allroots  = eval_imbalance( min_dc )
          if (verbose) print*, 'eval_allroots', eval_allroots  
          if (eval_allroots > 0.0) then
            dcleaf(pft) = 0.0
            cont = .false.
            if (verbose) print*, '* putting all to roots *'
          end if

          !------------------------------------------------------------------
          ! Test II: Evaluate balance if all is put to leaves.
          ! If C:N ratio of return is still lower than whole-plant C:N ratio, 
          ! then put all to leaves.
          !------------------------------------------------------------------
          if (cont) then
            if (verbose) print*, 'check alloation: all to leaves with dcleaf =', max_dc
            eval_allleaves = eval_imbalance( max_dc )
            if (verbose) print*, 'eval_allleaves', eval_allleaves  
            if (eval_allleaves < 0.0) then
              dcleaf(pft) = max_dc
              cont = .false.
              if (verbose) print*, '* putting all to leaves *'
            end if
          end if

          !------------------------------------------------------------------
          ! Optimum is between 0.0 (=min_dc) and max_dc. Find root of function 
          ! 'eval_imbalance()' in the interval [0.0, max_dc].
          !------------------------------------------------------------------
          if (cont) then
            if (verbose) print*, '*** finding root of eval_imbalance ***'
            if (write_logfile_eval_imbalance) open(unit=666,file='eval_imbalance.log',status='unknown')
            out_zeroin = zeroin( eval_imbalance, abserr, relerr, nmax, min_dc, max_dc )
            if ( out_zeroin%error /= 0 ) then
              print*, 'error code ', out_zeroin%error
              stop 'zeroin for eval_imbalance() failed'
              dcleaf(pft) = 0.0
            else
              dcleaf(pft) = out_zeroin%root
            end if
            if (write_logfile_eval_imbalance) close(unit=666)
            if (verbose) print*, 'no. of iterations   ', out_zeroin%niter
            if (verbose) print*, 'dcleaf(pft) is root ', dcleaf(pft)
            test = eval_imbalance( dcleaf(pft), .true. )
            if (verbose) print*, 'eval                              =', test
            if (verbose) print*, '----------------------------------'
            break_after_alloc = .true.
            ! if (doy==200) stop 'do beni'
          else
            break_after_alloc = .false.
          end if

          ! !------------------------------------------------------------------
          ! ! xxx debug: project next-day's fluxes with optimal dcleaf, derived now
          ! !------------------------------------------------------------------
          ! test = eval_imbalance( dcleaf(pft), verbose )
          ! stop

          !-------------------------------------------------------------------
          ! LEAF ALLOCATION
          !-------------------------------------------------------------------
          ! print*, 'pleaf before ', pleaf(pft,jpngr)
          ! print*, 'plabl before ', plabl(pft,jpngr)
          call allocate_leaf( &
            dcleaf(pft), &
            pleaf(pft,jpngr)%c%c12, pleaf(pft,jpngr)%n%n14, &
            plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
            solar%meanmppfd(:), mactnv_unitiabs(pft,:), &
            lai_ind(pft,jpngr), dnleaf(pft) &
            )
          ! print*, 'pleaf after  ', pleaf(pft,jpngr)
          ! print*, 'plabl after  ', plabl(pft,jpngr)
          ! stop

          !-------------------------------------------------------------------  
          ! Update leaf traits
          !-------------------------------------------------------------------  
          leaftraits(pft) = get_leaftraits( lai_ind(pft,jpngr), solar%meanmppfd(:), mactnv_unitiabs(pft,:) )

          !-------------------------------------------------------------------  
          ! Update fpc_grid and fapar_ind (not lai_ind)
          !-------------------------------------------------------------------  
          canopy(pft) = get_canopy( lai_ind(pft,jpngr) )

          ! call update_fpc_grid( pft, jpngr )

          !-------------------------------------------------------------------
          ! ROOT ALLOCATION
          !-------------------------------------------------------------------
          ! print*, 'proot before ', proot(pft,jpngr)
          ! print*, 'plabl before ', plabl(pft,jpngr)
          call allocate_root( &
            proot(pft,jpngr)%c%c12, proot(pft,jpngr)%n%n14, &
            plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
            pft, dcroot(pft), dnroot(pft) &
            )
          ! print*, 'proot after  ', proot(pft,jpngr)
          ! print*, 'plabl after  ', plabl(pft,jpngr)
          ! stop 

          !-------------------------------------------------------------------
          ! GROWTH RESPIRATION, NPP
          !-------------------------------------------------------------------
          ! add growth respiration to autotrophic respiration and substract from NPP
          ! (note that NPP is added to plabl in and growth resp. is implicitly removed
          ! from plabl above)
          drgrow(pft)   = ( 1.0 - params_plant%growtheff ) * ( dcleaf(pft) + dcroot(pft) ) / params_plant%growtheff

        else

          dcleaf(pft) = 0.0
          dcroot(pft) = 0.0
          dnleaf(pft) = 0.0
          dnroot(pft) = 0.0
          drgrow(pft) = 0.0

          ! print*, 'not growing ...'

        end if

      else

        stop 'allocation_daily not implemented for trees'

      end if

    end do

    ! if (doy==39) stop 'on day 39'

    ! test_calloc = test_calloc + dcleaf + dcroot
    ! print*, 'test_calloc', test_calloc
    
    ! print*, '--- END allocation_daily:'

  end subroutine allocation_daily


  function eval_imbalance( mydcleaf, verbose ) result ( eval )
    !/////////////////////////////////////////////////////////
    ! Evaluates C:N ratio of new assimilation after allocation 
    ! versus whole-plant C:N ratio after allocation. Optimal 
    ! allocation is where the two are equal. 
    !---------------------------------------------------------
    use md_classdefs, only: orgpool, nitrogen
    use md_plant, only: params_pft_plant, params_plant, get_fapar, &
      canopy_type, get_canopy
    use md_gpp, only: calc_dgpp, calc_drd, mactnv_unitiabs, mlue, mrd_unitiabs
    use md_nuptake, only: calc_dnup, outtype_calc_dnup
    use md_npp, only: calc_resp_maint, calc_cexu, deactivate_root
    use md_findroot_fzeroin
    use md_waterbal, only: solar, evap
    use md_soiltemp, only: dtemp_soil
    use md_ntransform, only: pninorg

    ! arguments
    real, intent(in)              :: mydcleaf
    logical, intent(in), optional :: verbose

    ! function return variable
    real :: eval

    ! local variables
    real    :: cleaf
    real    :: nleaf
    real    :: croot
    real    :: nroot
    real    :: clabl
    real    :: nlabl
    integer :: usepft
    integer :: usemoy
    integer :: usedoy
    integer :: usejpngr
    real    :: airtemp

    integer :: lu

    real :: mydcroot
    real :: mydnleaf
    real :: mydnroot
    real :: mylai
    real :: gpp
    real :: npp
    real :: rd
    real :: mresp_root
    real :: cexu
    real :: avl
    real :: dc
    real :: dn
    real :: kcleaf
    real :: knleaf
    real :: kcroot
    real :: knroot

    real :: nleaf0
    real :: lai0, lai1

    type( outtype_zeroin )    :: out_zeroin
    type( outtype_calc_dnup ) :: out_calc_dnup
    type( canopy_type )       :: mycanopy

    ! write(0,*) '--- in eval_imbalance with mydcleaf=', mydcleaf

    ! Copy to local variables for shorter writing
    cleaf    = state_eval_imbalance%pleaf%c%c12
    nleaf    = state_eval_imbalance%pleaf%n%n14
    croot    = state_eval_imbalance%proot%c%c12
    nroot    = state_eval_imbalance%proot%n%n14
    clabl    = state_eval_imbalance%plabl%c%c12
    nlabl    = state_eval_imbalance%plabl%n%n14
    usepft   = state_eval_imbalance%usepft
    usemoy   = state_eval_imbalance%usemoy
    usedoy   = state_eval_imbalance%usedoy
    usejpngr = state_eval_imbalance%usejpngr
    airtemp  = state_eval_imbalance%airtemp

    !-------------------------------------------------------------------
    ! LEAF ALLOCATION
    !-------------------------------------------------------------------
    call allocate_leaf( mydcleaf, cleaf, nleaf, clabl, nlabl, solar%meanmppfd(:), mactnv_unitiabs(usepft,:), mylai, mydnleaf )

    !-------------------------------------------------------------------  
    ! Update fpc_grid and fapar_ind (not lai_ind)
    !-------------------------------------------------------------------  
    mycanopy = get_canopy( mylai )

    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    call allocate_root( croot, nroot, clabl, nlabl, usepft, mydcroot, mydnroot )

    !-------------------------------------------------------------------
    ! PROJECT NEXT DAY'S C AND N BALANCE:
    ! decay, GPP, respiration, N uptake
    !-------------------------------------------------------------------
    ! ! continuous decay
    ! kcleaf = cleaf * params_pft_plant(usepft)%k_decay_leaf
    ! knleaf = nleaf * params_pft_plant(usepft)%k_decay_leaf
    ! kcroot = croot * params_pft_plant(usepft)%k_decay_root
    ! knroot = nroot * params_pft_plant(usepft)%k_decay_root
    
    ! cleaf  = cleaf - kcleaf
    ! nleaf  = nleaf - knleaf
    ! croot  = croot - kcroot
    ! nroot  = nroot - knroot

    ! Calculate next day's C and N return after assumed allocation (tissue turnover happens before!)

    lu = params_pft_plant(usepft)%lu_category

    gpp           = calc_dgpp( mycanopy%fapar_ind, solar%dppfd(usedoy), mlue(usepft,usemoy), airtemp, evap(lu)%cpa )
    rd            = calc_drd(  mycanopy%fapar_ind, solar%meanmppfd(usemoy), mrd_unitiabs(usepft,usemoy), airtemp, evap(lu)%cpa  )
    mresp_root    = calc_resp_maint( croot, params_plant%r_root, airtemp )
    npp           = gpp - rd - mresp_root
    cexu          = calc_cexu( croot, airtemp ) 
    avl           = clabl + npp - cexu
    if (avl<0.0) call deactivate_root( gpp, rd, clabl, orgpool(carbon(croot),nitrogen(nroot)), mresp_root, npp, cexu, airtemp )
    dc            = npp - cexu
    out_calc_dnup = calc_dnup( cexu, pninorg(lu,usejpngr)%n14 )
    dn            = out_calc_dnup%fix + out_calc_dnup%act

    if (present(verbose)) then
      print*, '-------------'
      print*, 'in eval_imbalance, projected quantities:'
      print*, 'fapar              ', mycanopy%fapar_ind
      print*, 'dppfd              ', solar%dppfd(usedoy)
      print*, 'lue                ', mlue(usepft,usemoy)
      print*, 'temp               ', airtemp
      print*, 'CPA                ', evap(lu)%cpa
      print*, 'gpp                ', gpp
      print*, 'npp                ', npp
      print*, 'rd                 ', rd
      print*, 'cexu               ', cexu
      print*, 'mresp_root         ', mresp_root
      print*, 'dc                 ', dc
      print*, 'dn                 ', dn
      print*, 'current   plant C:N', ( cleaf + croot ) / ( nleaf + nroot )
      print*, '(dn + nlabl)       ', (dn + nlabl)
      print*, 'tomorrows plant C:N', params_plant%growtheff * (dc + clabl) / (dn + nlabl)
      print*, '-------------'
      ! stop
    end if

    !-------------------------------------------------------------------
    ! EVALUATION QUANTITY - IS MINIMISED BY OPTIMISATION
    ! Evaluation quantity is the difference between the 
    ! C:N ratio of new assimilates and the C:N ratio 
    ! of the whole plant after allocation.
    !-------------------------------------------------------------------
    if ((dn + nlabl)==0.0) then
      eval = -999.0
    else
      !     |---------------------------------------------------|  |------------------------------------|
      eval = params_plant%growtheff * (dc + clabl) / (dn + nlabl) - ( cleaf + croot ) / ( nleaf + nroot )
      !     |---------------------------------------------------|  |------------------------------------|
      !     |lab. pool C:N ratio after acq. nxt. day            |  | current whole-plant C:N ratio      |
      !     |---------------------------------------------------|  |------------------------------------|
    end if

    ! write(0,*) 'new assimilates stoichiometry: ', growtheff * (dc + clabl) / (dn + nlabl)
    ! write(0,*) 'current stoichiometry:         ', ( cleaf + croot ) / ( nleaf + nroot )
    if (write_logfile_eval_imbalance) write(666,*) mydcleaf, ",", eval

    ! write(0,*) mydcleaf, ",", eval
    ! write(0,*) 'eval                                       ', eval
    ! write(0,*) '-------------'
    ! stop

  end function eval_imbalance


  subroutine allocate_leaf( mydcleaf, cleaf, nleaf, clabl, nlabl, meanmppfd,    nv,    lai,   mydnleaf )
    !///////////////////////////////////////////////////////////////////
    ! LEAF ALLOCATION
    ! Sequence of steps:
    ! - increment foliage C pool
    ! - update LAI
    ! - calculate canopy-level foliage N as a function of LAI 
    ! - reduce labile pool by C and N increments
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant

    ! arguments
    real, intent(in)                    :: mydcleaf
    real, intent(inout)                 :: cleaf, nleaf
    real, intent(inout)                 :: clabl, nlabl
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv
    real, intent(out)                   :: lai
    real, optional, intent(out)         :: mydnleaf

    ! local variables
    real :: nleaf0
    real :: dclabl, dnlabl

    ! find LAI, given new leaf mass. This is necessary to get leaf-N as 
    ! a function of LAI.
    if (mydcleaf>0.0) then

      ! print*, 'cleaf before ', cleaf 
      cleaf  = cleaf + mydcleaf
      ! print*, 'mydcleaf     ', mydcleaf 

      ! print*, 'cleaf = 0.5', get_lai( 0.5  , meanmppfd(:), nv(:) )
      ! print*, 'cleaf = 100', get_lai( 100.0, meanmppfd(:), nv(:))
      ! stop

      ! Calculate LAI as a function of leaf C
      lai = get_lai( cleaf, meanmppfd(:), nv(:) )
      ! print*, 'in allocate_leaf: cleaf, lai :', cleaf, lai
      ! print*, 'mydcleaf ', mydcleaf 
      ! ! stop

      ! calculate canopy-level leaf N as a function of LAI
      nleaf0   = nleaf      
      nleaf    = get_leaf_n_canopy( lai, meanmppfd(:), nv(:) )
      mydnleaf = nleaf - nleaf0
      ! if (mydnleaf>0.0) then
      !   print*, 'mydcleaf/dnleaf ', mydcleaf/mydnleaf 
      ! end if

      ! subtract from labile pool, making sure pool does not get negative
      dclabl = min( clabl, 1.0 / params_plant%growtheff * mydcleaf )
      dnlabl = min( nlabl, mydnleaf )
      if ( (dclabl - clabl) > 1e-8 ) stop 'trying to remove too much from labile pool: leaf C'
      if ( (dnlabl - nlabl) > 1e-8 ) stop 'trying to remove too much from labile pool: leaf N'
      clabl  = clabl - dclabl
      nlabl  = nlabl - dnlabl

      ! write(0,*) 'r_ntoc_root(pft)  ', r_ntoc_root(pft)

    else

      lai      =  get_lai( cleaf, meanmppfd(:), nv(:) )
      mydnleaf = 0.0

    end if

  end subroutine allocate_leaf


  subroutine allocate_root( croot, nroot, clabl, nlabl, pft, mydcroot, mydnroot )
    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant, params_pft_plant

    ! arguments
    real, intent(inout)         :: croot, nroot
    real, intent(inout)         :: clabl, nlabl
    integer, intent(in)         :: pft
    real, optional, intent(out) :: mydcroot
    real, optional, intent(out) :: mydnroot

    ! local variables
    real :: dclabl
    real :: dnlabl

    if (clabl>0.0 .and. nlabl>0.0) then
      ! use remainder for allocation to roots
      mydcroot = min( params_plant%growtheff * clabl, params_pft_plant(pft)%r_cton_root * nlabl )
      mydnroot = min( mydcroot * params_pft_plant(pft)%r_ntoc_root, nlabl )

      dclabl = min( clabl, 1.0 / params_plant%growtheff * mydcroot )
      dnlabl = min( nlabl, mydnroot )
      if ( (dnlabl - nlabl) > 1e-8 ) stop 'trying to remove too much from labile pool: root N'
      if ( (dclabl - clabl) > 1e-8 ) stop 'trying to remove too much from labile pool: root C'
      clabl  = clabl - dclabl
      nlabl  = nlabl - dnlabl

      ! dclabl = min( plabl(pft,jpngr)%c%c12, 1.0 / growtheff * dcleaf(pft) )
      ! dnlabl = min( plabl(pft,jpngr)%n%n14, dnleaf(pft) )

      ! plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - dclabl
      ! plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnlabl

      ! write(0,*) 'mydcroot ', mydcroot 
      ! write(0,*) 'mydnroot ', mydnroot 
      if (mydcroot<0.0) stop 'root allocation neg.: C'
      if (mydnroot<0.0) stop 'root allocation neg.: N'

      croot = croot + mydcroot
      nroot = nroot + mydnroot

      ! if (present(verbose)) then
      !   print*, 'after presumed allocation:'
      !   print*, 'clabl', clabl ! OK
      !   print*, 'nlabl', nlabl ! OK
      !   ! print*, 'mydcleaf', mydcleaf
      !   ! print*, 'mydcroot', mydcroot
      !   ! print*, 'mydnleaf', mydnleaf
      !   ! print*, 'mydnroot', mydnroot

      !   ! print*, 'after presumed allocation:'
      !   ! print*, 'B cleaf', cleaf
      !   ! print*, 'B root', root
      ! end if
    end if

  end subroutine allocate_root


  function get_lai( cleaf, meanmppfd, nv ) result( lai )
    !////////////////////////////////////////////////////////////////
    ! Calculates LAI as a function of canopy-level leaf-C:
    ! Cleaf = Mc * c * ( I0 * ( 1 - exp( -kL ) * nv * b + L * a ) )
    ! Cannot be solved analytically for L = f(Cleaf). Therefore, 
    ! numerical root-searching algorithm is applied so that
    ! Cleaf / ( Mc * c ) - ( I0 * ( 1 - exp( -kL ) * nv * b + L * a ) ) = 0
    ! This is implemented in function 'mustbe_zero_for_lai()'.
    !----------------------------------------------------------------
    ! use md_params_core, only: nmonth
    use md_findroot_fzeroin

    ! arguments 
    real, intent(in)                    :: cleaf
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv 

    ! local variables
    real                 :: abserr
    real                 :: relerr
    real                 :: lower
    real                 :: upper
    integer, parameter   :: nmax = 100
    type(outtype_zeroin) :: out_zeroin

    ! xxx debug
    real :: test

    ! function return value
    real :: lai

    ! local variables
    real :: maxnv

    if (cleaf>0.0) then
      ! Metabolic N is predicted and is optimised at a monthly time scale. 
      ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
      ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost.
      ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
      ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
      maxnv = maxval( meanmppfd(:) * nv(:) )
      ! print*, 'maxnv ', maxnv
      ! stop

      ! print*, 'what is LAI for Cleaf=', cleaf
      ! print*, 'maxnv ', maxnv

      ! Update state. This derived-type variable is "global" within this module
      state_mustbe_zero_for_lai%cleaf = cleaf
      state_mustbe_zero_for_lai%maxnv = maxnv

      ! Calculate initial guess for LAI (always larger than actual LAI)
      lower = 0.0 ! uninformed lower bound
      upper = 1.0 / ( c_molmass * params_alloc%ncw_min * params_alloc%r_ctostructn_leaf ) * cleaf
      ! upper = 20.0
      ! print*, 'upper ', upper

      ! print*, 'lower =', lower
      ! test = mustbe_zero_for_lai( lower ) 
      ! print*, '=> test =', test

      ! print*, 'upper =', upper
      ! test = mustbe_zero_for_lai( upper ) 
      ! print*, '=> test =', test

      ! call function zeroin to find root (value of LAI for which evaluation expression is zero)
      abserr=100.0*XMACHEPS !*10e5
      relerr=1000.0*XMACHEPS !*10e5

      ! print*, 'abserr', abserr
      ! print*, 'relerr', relerr
      ! stop 'here'

      ! print*, '*** finding root of mustbe_zero_for_lai ***'
      out_zeroin = zeroin( mustbe_zero_for_lai, abserr, relerr, nmax, lower, upper )
      if ( out_zeroin%error /= 0 ) then
        lai = 0.0
        print*, 'error code', out_zeroin%error
        stop 'zeroin for mustbe_zero_for_lai() failed'
      else
        lai = out_zeroin%root
      end if

      ! print*, 'out_zeroin', out_zeroin
      ! print*, 'cleaf', cleaf
      ! print*, 'lai', lai
      ! stop
    
    else

      lai = 0.0

    end if

  end function get_lai


  function mustbe_zero_for_lai( mylai ) result( mustbe_zero )
    !/////////////////////////////////////////////////////////
    ! This function returns value of the expression 'mustbe_zero'. 
    ! If expression is zero, then 'mylai' is a root and is the 
    ! LAI for a given Cleaf (meanmppfd, cleaf, and nv are 
    ! passed on to this function as a derived type state.)
    !---------------------------------------------------------
    ! Cleaf = c_molmass * params_alloc%r_ctostructn_leaf * N_canop_cellwall
    ! N_canop_cellwall = LAI * params_alloc%ncw_min + nv * Iabs * params_alloc%r_n_cw_v
    ! Iabs = meanmppfd * (1-exp( -kbeer * LAI))
    ! ==> Cleaf = f(LAI) = c_molmass * params_alloc%r_ctostructn_leaf * [ meanmppfd * (1-exp( -kbeer * LAI)) * nv * params_alloc%r_n_cw_v + LAI * params_alloc%ncw_min ]
    ! ==> LAI = f(Cleaf) leads to inhomogenous equation. Therefore apply root finding algorithm so that:
    ! 0 = cleaf / ( c_molmass * params_alloc%r_ctostructn_leaf ) - meanmppfd * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * nv * params_alloc%r_n_cw_v - mylai * params_alloc%ncw_min
    !---------------------------------------------------------
    use md_plant, only: params_plant, get_fapar

    ! arguments
    real, intent(in) :: mylai

    ! function return value
    real :: mustbe_zero

    ! local variables
    real :: mycleaf
    real :: mymaxnv

    ! print*, '--- in mustbe_zero_for_lai with mydcleaf=', mylai

    ! Read from updated state. This derived-type variable is "global" within this module
    mycleaf = state_mustbe_zero_for_lai%cleaf
    mymaxnv = state_mustbe_zero_for_lai%maxnv

    ! print*, '----------'
    ! print*, 'inside mustbe_zero_for_lai: '
    ! print*, 'mylai', mylai
    ! print*, 'mycleaf', mycleaf
    ! print*, 'mymaxnv', mymaxnv

    ! mustbe_zero = cleaf / ( c_molmass * params_alloc%r_ctostructn_leaf ) - meanmppfd * nv * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * params_alloc%r_n_cw_v - mylai * params_alloc%ncw_min
    ! mustbe_zero = cleaf / ( c_molmass * params_alloc%r_ctostructn_leaf ) - maxnv * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * params_alloc%r_n_cw_v - mylai * params_alloc%ncw_min
    mustbe_zero = mycleaf / ( c_molmass * params_alloc%r_ctostructn_leaf ) - mymaxnv * get_fapar( mylai ) * params_alloc%r_n_cw_v - mylai * params_alloc%ncw_min

    ! print*, 'mustbe_zero                           ', mustbe_zero
    ! print*, '-------------'


  end function mustbe_zero_for_lai


  function get_rcton_init( meanmppfd, nv ) result( rcton )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    ! Cleaf = c_molmass * params_alloc%r_ctostructn_leaf * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * params_alloc%r_n_cw_v + LAI * params_alloc%ncw_min ]
    ! Nleaf = n_molmass * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * (params_alloc%r_n_cw_v + 1) + LAI * params_alloc%ncw_min ]
    ! linearization around LAI = 0 ==> (1-exp(-k*L)) ~= k*L
    ! ==> Cleaf ~= LAI * c_molmass * params_alloc%r_ctostructn_leaf * ( meanmppfd * kbeer * nv * params_alloc%r_n_cw_v + params_alloc%ncw_min )
    ! ==> Nleaf ~= LAI * n_molmass * ( meanmppfd * kbeer * nv * (params_alloc%r_n_cw_v + 1) + params_alloc%ncw_min )
    ! r_cton = Cleaf / Nleaf
    !----------------------------------------------------------------
    ! use md_params_core, only: nmonth
    use md_plant, only: params_plant

    ! arguments
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real :: rcton

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )
    rcton = ( c_molmass * params_alloc%r_ctostructn_leaf * ( maxnv * params_plant%kbeer * params_alloc%r_n_cw_v + params_alloc%ncw_min )) / ( n_molmass * ( maxnv * params_plant%kbeer * (params_alloc%r_n_cw_v + 1.0) + params_alloc%ncw_min ) )

    ! rcton = ( c_molmass * params_alloc%r_ctostructn_leaf * ( meanmppfd * kbeer * nv * params_alloc%r_n_cw_v + params_alloc%ncw_min )) / ( n_molmass * ( meanmppfd * kbeer * nv * (params_alloc%r_n_cw_v + 1) + params_alloc%ncw_min ) )

  end function get_rcton_init


  function get_leaf_n_metabolic_canopy( mylai, meanmppfd, nv ) result( mynleaf_metabolic )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! LAI * n_metabolic = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    !----------------------------------------------------------------
    ! use md_params_core, only: nmonth
    use md_plant, only: get_fapar

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real :: mynleaf_metabolic  ! mol N 

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )

    mynleaf_metabolic = maxnv * get_fapar( mylai )

  end function get_leaf_n_metabolic_canopy


  function get_leaf_n_structural_canopy( mylai, mynleaf_metabolic ) result( mynleaf_structural )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! LAI * n_structural = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    !----------------------------------------------------------------
    ! use md_params_core, only: nmonth

    ! arguments
    real, intent(in) :: mylai
    real, intent(in) :: mynleaf_metabolic

    ! function return variable
    real :: mynleaf_structural  ! mol N 

    mynleaf_structural = mynleaf_metabolic * params_alloc%r_n_cw_v + mylai * params_alloc%ncw_min

    ! print*, '--- in get_leaf_n_structural_canopy'
    ! print*, 'mylai ', mylai
    ! print*, 'mynleaf_metabolic ', mynleaf_metabolic
    ! print*, 'params_alloc%r_n_cw_v ', params_alloc%r_n_cw_v
    ! print*, 'params_alloc%ncw_min ', params_alloc%ncw_min
    ! print*, 'mynleaf_structural ', mynleaf_structural
    ! print*, '-------------------------------'

  end function get_leaf_n_structural_canopy


  function get_leaf_n_canopy( mylai, meanmppfd, nv ) result( mynleaf )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    ! Nleaf = LAI * (n_metabolic + n_cellwall) * n_molmass
    ! LAI * n_metabolic = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    ! LAI * n_cellwall = LAI * (params_alloc%ncw_min + params_alloc%r_n_cw_v * n_metabolic)
    ! ==> Nleaf = n_molmass * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * (params_alloc%r_n_cw_v + 1) + LAI * params_alloc%ncw_min ]
    !----------------------------------------------------------------
    ! use md_params_core, only: nmonth

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real :: mynleaf ! g N

    ! local variables
    real :: nleaf_metabolic   ! mol N m-2
    real :: nleaf_structural  ! mol N m-2

    nleaf_metabolic  = get_leaf_n_metabolic_canopy(  mylai, meanmppfd, nv )
    nleaf_structural = get_leaf_n_structural_canopy( mylai, nleaf_metabolic )
    mynleaf          = n_molmass * ( nleaf_metabolic + nleaf_structural )

    ! print*, '--- in get_leaf_n_canopy'
    ! print*, 'nleaf_metabolic ', nleaf_metabolic
    ! print*, 'nleaf_structural ', nleaf_structural
    ! print*, 'mynleaf ', mynleaf
    ! print*, '-------------------------------'

    ! mynleaf = n_molmass * ( maxnv * get_fapar( mylai ) * ( 1.0 + params_alloc%r_n_cw_v ) + mylai * params_alloc%ncw_min )

  end function get_leaf_n_canopy


  function get_leaftraits( mylai, meanmppfd, nv ) result( out_traits )
    !////////////////////////////////////////////////////////////////
    ! Calculates leaf traits based on (predicted) metabolic Narea and
    ! (prescribed) parameters that relate structural to metabolic
    ! Narea and Carea to structural Narea:
    ! Narea_metabolic  = predicted
    ! Narea_structural = a + b * Narea_metabolic
    ! Carea            = c * Narea_structural
    !----------------------------------------------------------------
    use md_params_core, only: c_content_of_biomass
    use md_plant, only: leaftraits_type

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    type( leaftraits_type ) :: out_traits

    ! local variables
    real :: mynarea_metabolic_canop   ! mol N m-2-ground
    real :: mynarea_structural_canop  ! mol N m-2-ground

    mynarea_metabolic_canop  = get_leaf_n_metabolic_canopy(  mylai, meanmppfd(:), nv(:) )     ! mol N m-2-ground    
    mynarea_structural_canop = get_leaf_n_structural_canopy( mylai, mynarea_metabolic_canop ) ! mol N m-2-ground
    
    out_traits%narea_metabolic  = n_molmass * mynarea_metabolic_canop / mylai   ! g N m-2-leaf
    out_traits%narea_structural = n_molmass * mynarea_structural_canop / mylai  ! g N m-2-leaf
    out_traits%narea            = n_molmass * ( mynarea_metabolic_canop + mynarea_structural_canop ) / mylai ! g N m-2-leaf
    out_traits%lma              = c_molmass * params_alloc%r_ctostructn_leaf * mynarea_structural_canop / mylai 
    out_traits%nmass            = out_traits%narea / ( out_traits%lma / c_content_of_biomass )
    out_traits%r_cton_leaf      = out_traits%lma / out_traits%narea
    out_traits%r_ntoc_leaf      = 1.0 / out_traits%r_cton_leaf

    ! print*, '--- in get_leaftraits'
    ! print*, 'mylai                  ', mylai
    ! print*, 'traits%narea_metabolic ', traits%narea_metabolic 
    ! print*, 'traits%narea_structural', traits%narea_structural
    ! print*, 'traits%narea           ', traits%narea           
    ! print*, 'traits%lma             ', traits%lma             
    ! print*, 'traits%nmass           ', traits%nmass           
    ! print*, 'traits%r_cton_leaf     ', traits%r_cton_leaf     
    ! print*, 'traits%r_ntoc_leaf     ', traits%r_ntoc_leaf     
    ! print*, '-------------------------------'
    ! stop 

  end function get_leaftraits


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
