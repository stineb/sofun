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
  use _params_core, only: npft, nlu, maxgrid, ndaymonth, ndayyear, &
    c_molmass, n_molmass, nmonth
  use _classdefs

  implicit none

  ! Parameters determining the relationship of structural N and C to metabolic N
  ! From regressing Narea to metabolic Narea in Hikosaka data
  ! real, parameter    :: r_n_cw_v = 1.23223            ! slope in the relationship of non-metabolic versus metabolic N per leaf area
  ! real, parameter    :: ncw_min = 0.056               ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area
  real, parameter    :: r_n_cw_v = 0.1            ! slope in the relationship of non-metabolic versus metabolic N per leaf area
  real, parameter    :: ncw_min = 0.1               ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area
  real, parameter    :: r_ctostructn_leaf = 40        ! constant ratio of C to structural N
  
  logical, parameter :: write_logfile_eval_imbalance = .false.

  !------------------------------------------------------------------
  ! Define derived type to store current state variables.
  ! Is defined within all subroutines in this module.
  ! This avoids having to pass them as arguments.
  !------------------------------------------------------------------

  type statetype_eval_imbalance
    type(orgpool)           :: pleaf
    type(orgpool)           :: proot
    type(orgpool)           :: plabl
    real                    :: crownarea
    real                    :: nind
    real                    :: mlue
    real                    :: dppfd
    real                    :: mrd_unitiabs
    real, dimension(nmonth) :: nv
    real, dimension(nmonth) :: meanmppfd
    type(nitrogen)          :: pninorg
    real                    :: pft 
    integer                 :: usemoy
    real                    :: soiltemp
  end type statetype_eval_imbalance

  type statetype_mustbe_zero_for_lai
    real :: cleaf
    real :: maxnv
  end type statetype_mustbe_zero_for_lai

  type leaftraits_type
    real :: narea
    real :: narea_metabolic
    real :: narea_structural
    real :: lma
    real :: nmass
    real :: r_cton_leaf
    real :: r_ntoc_leaf
  end type leaftraits_type

  ! states area global within module (instead of being passed on as arguments)
  type(statetype_eval_imbalance)      :: state_eval_imbalance
  type(statetype_mustbe_zero_for_lai) :: state_mustbe_zero_for_lai

  ! module-specific variables
  ! type(orgpool), dimension(npft) :: dpleaf
  ! type(orgpool), dimension(npft) :: dproot

  real, dimension(npft) :: dcleaf
  real, dimension(npft) :: dnleaf
  real, dimension(npft) :: dcroot
  real, dimension(npft) :: dnroot

  ! output variables
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
      r_ntoc_leaf, crownarea, drauto, dnpp, lai_ind, nind, fapar_ind, fpc_grid, &
      narea, narea_metabolic, narea_structural, lma, nmass, psoilphys, solar
    use _gpp, only: mlue, mrd_unitiabs, mactnv_unitiabs
    use _vegdynamics, only: update_fpc_grid
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
    real    :: min_dc
    real    :: eval_allleaves
    real    :: eval_allroots
    real    :: abserr
    real    :: relerr
    real    :: nleaf0
    real    :: lai0, lai1
    integer, parameter :: nmax = 100

    type(outtype_zeroin)  :: out_zeroin
    type(leaftraits_type) :: traits

    integer, save      :: invocation = 0             ! internally counted simulation year
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

    !-------------------------------------------------------------------------
    ! xxx try
    !-------------------------------------------------------------------------
    flexalloc = .true.
    !-------------------------------------------------------------------------

    abserr=100.0*XMACHEPS !*10e5
    relerr=1000.0*XMACHEPS !*10e5

    do pft=1,npft

      lu = lu_category(pft)

      if (grass(pft)) then

        ! write(0,*) '--- allocation_daily, doy:',doy
        ! write(0,*)' pninorg(lu,jpngr)%n14',pninorg
        ! write(0,*) 'plabl(lu,jpngr)', plabl
        ! write(0,*) 'cleaf          ', pleaf(pft,jpngr)%c%c12
        ! write(0,*) 'croot          ', proot(pft,jpngr)%c%c12
        ! write(0,*) 'C:N in leaves  ', cton( pleaf(pft,jpngr), default=0.0 )

        if ( plabl(pft,jpngr)%c%c12>0.0 .and. plabl(pft,jpngr)%n%n14>0.0 ) then

          write(0,*) 'starting to grow on day', doy
          write(0,*) 'with plabl =', plabl
          ! stop

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
            state_eval_imbalance%nind         = nind(pft,jpngr)
            state_eval_imbalance%mlue         = mlue(usemoy)
            state_eval_imbalance%dppfd        = solar%dppfd(usedoy)
            state_eval_imbalance%mrd_unitiabs = mrd_unitiabs(usemoy)
            state_eval_imbalance%meanmppfd(:) = solar%meanmppfd(:)
            state_eval_imbalance%pninorg      = pninorg(lu,jpngr)  ! the only data that is not yet available - use today's value 
            state_eval_imbalance%pft          = pft
            state_eval_imbalance%nv(:)        = mactnv_unitiabs(:)
            state_eval_imbalance%usemoy       = usemoy
            state_eval_imbalance%soiltemp     = psoilphys(lu,jpngr)%temp

            !------------------------------------------------------------------
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by all labile N,
            ! discounted by the yield factor.
            !------------------------------------------------------------------
            if (pleaf(pft,jpngr)%c%c12==0.0) then
              write(0,*) 'Calculating initial C:N ratio'
              ! initial guess based on Taylor approximation of Cleaf and Nleaf function around cleaf=0
              r_cton_leaf(pft,jpngr) = get_rcton_init( solar%meanmppfd(:), mactnv_unitiabs(:) )
              write(0,*) 'solar%meanmppfd(:)      ', solar%meanmppfd(:)  
              write(0,*) 'mactnv_unitiabs(:)', mactnv_unitiabs(:)  
              write(0,*) 'initial guess: r_cton_leaf(pft,jpngr) ', r_cton_leaf(pft,jpngr)  
              stop
            end if
            max_dcleaf_n_constraint = plabl(pft,jpngr)%n%n14 * r_cton_leaf(pft,jpngr)
            max_dcroot_n_constraint = plabl(pft,jpngr)%n%n14 * r_cton_root(pft) ! should be obsolete as generally r_ntoc_leaf > r_ntoc_root
            max_dc = min( growtheff * plabl(pft,jpngr)%c%c12, max_dcleaf_n_constraint, max_dcroot_n_constraint )
            min_dc = 0.0

            ! write(0,*) 'plabl(pft,jpngr)', plabl(pft,jpngr)  
            ! write(0,*) 'r_cton_leaf(pft,jpngr)',r_cton_leaf(pft,jpngr)
            ! write(0,*) 'r_cton_root(pft)',r_cton_root(pft)
            ! write(0,*) 'growtheff', growtheff  
            ! write(0,*) 'max_dcleaf_n_constraint', max_dcleaf_n_constraint  
            ! write(0,*) 'max_dcroot_n_constraint', max_dcroot_n_constraint  
            ! write(0,*) 'max_dc', max_dc  
            ! stop

            ! write(0,*) 'moy',moy
            ! write(0,*) 'lai',lai
            ! write(0,*) 'mlue',mlue(moy) ! ok
            ! write(0,*) 'dppfd',dppfd(doy)  ! ok
            ! write(0,*) 'mrd_unitiabs',mrd_unitiabs(moy) ! ok
            ! write(0,*) 'meanmppfd',meanmppfd(moy) ! ok

            ! write(0,*) 'calculating LAI for Cleaf = 2.857124'
            ! test = get_lai( 2.857124, solar%meanmppfd(:), mactnv_unitiabs(:) )
            ! write(0,*) 'lai                       =', test
            ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating LAI for Cleaf = 100'
            ! test = get_lai( 100.0, solar%meanmppfd(:), mactnv_unitiabs(:) )
            ! write(0,*) 'lai                       =', test
            ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating LAI for Cleaf = 2.857124'
            ! test = get_lai( 2.857124, solar%meanmppfd(:), mactnv_unitiabs(:) )
            ! write(0,*) 'lai                       =', test
            ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating LAI for Cleaf = 0.5'
            ! test = get_lai( 0.5  , solar%meanmppfd(:), mactnv_unitiabs(:) )
            ! write(0,*) 'lai                       =', test
            ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating LAI for Cleaf = 100'
            ! test = get_lai( 100.0, solar%meanmppfd(:), mactnv_unitiabs(:) )
            ! write(0,*) 'lai                       =', test
            ! ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating eval_imbalance for dc = 0.0'
            ! test = eval_imbalance( 0.0 )
            ! write(0,*) 'eval                              =', test
            ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating eval_imbalance for dc = ', 2.7
            ! test = eval_imbalance( 2.7 )
            ! write(0,*) 'eval                              =', test
            ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating eval_imbalance for dc = ', 2.75
            ! test = eval_imbalance( 2.75 )
            ! write(0,*) 'eval                              =', test
            ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating eval_imbalance for dc = ', 2.8
            ! test = eval_imbalance( 2.8 )
            ! write(0,*) 'eval                              =', test
            ! write(0,*) '----------------------------------'
            
            ! write(0,*) 'calculating eval_imbalance for dc = ', 2.794742 
            ! test = eval_imbalance( 2.794742 )
            ! write(0,*) 'eval                              =', test
            ! write(0,*) '----------------------------------'

            ! write(0,*) 'calculating eval_imbalance for dc = 0.0'
            ! test = eval_imbalance( 0.0 )
            ! write(0,*) 'eval                              =', test
            ! write(0,*) '----------------------------------'

            !------------------------------------------------------------------
            ! Optimisation by balanced growth
            ! Test I: Evaluate balance if all is put to roots.
            ! If C:N ratio of return is still greater than whole-plant C:N 
            ! ratio, then put all to roots.
            !------------------------------------------------------------------
            cont = .true.
            write(0,*) 'check alloation: all to roots'
            eval_allroots  = eval_imbalance( min_dc )
            write(0,*) 'eval_allroots', eval_allroots  
            if (eval_allroots > 0.0) then
              dcleaf(pft) = 0.0
              cont = .false.
              write(0,*) '* putting all to roots *'
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

            stop

            !------------------------------------------------------------------
            ! Optimum is between 0.0 (=min_dc) and max_dc. Find root of function 
            ! 'eval_imbalance()' in the interval [0.0, max_dc].
            !------------------------------------------------------------------
            if (cont) then
              write(0,*) '*** finding root of eval_imbalance ***'
              if (write_logfile_eval_imbalance) open(unit=666,file='eval_imbalance.log',status='unknown')
              ! write(0,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
              out_zeroin = zeroin( eval_imbalance, abserr, relerr, nmax, min_dc, max_dc )
              if ( out_zeroin%error /= 0 ) then
                write(0,*) 'error code ', out_zeroin%error
                stop 'zeroin for eval_imbalance() failed'
                dcleaf(pft) = 0.0
              else
                dcleaf(pft) = out_zeroin%root
              end if
              if (write_logfile_eval_imbalance) close(unit=666)
              write(0,*) 'no. of iterations   ', out_zeroin%niter
              write(0,*) 'dcleaf(pft) is root ', dcleaf(pft)
              write(0,*) 'found root of eval_imbalance for dc = ', dcleaf(pft)
              test = eval_imbalance( dcleaf(pft), .true. )
              write(0,*) 'eval                              =', test
              write(0,*) '----------------------------------'
              ! if (doy==200) stop 'do beni'
            end if

            !------------------------------------------------------------------
            ! xxx debug: project next-day's fluxes with optimal dcleaf, derived now
            !------------------------------------------------------------------
            ! test = eval_imbalance( dcleaf(pft), .true. )
            ! stop

            !-------------------------------------------------------------------
            ! LEAF ALLOCATION
            !-------------------------------------------------------------------  
            write(0,*) 'pleaf before ', pleaf(pft,jpngr)
            write(0,*) 'plabl before ', plabl(pft,jpngr)
            call allocate_leaf( &
              dcleaf(pft), &
              pleaf(pft,jpngr)%c%c12, pleaf(pft,jpngr)%n%n14, &
              plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
              solar%meanmppfd(:), mactnv_unitiabs(:), &
              lai_ind(pft,jpngr), dnleaf(pft) &
              )
            write(0,*) 'pleaf after  ', pleaf(pft,jpngr)
            write(0,*) 'plabl after  ', plabl(pft,jpngr)
            stop

            !-------------------------------------------------------------------  
            ! Update leaf traits
            !-------------------------------------------------------------------  
            traits                      = get_leaftraits( lai_ind(pft,jpngr), solar%meanmppfd(:), mactnv_unitiabs(:) )
            narea(pft)                  = traits%narea
            narea_metabolic(pft)        = traits%narea_metabolic
            narea_structural(pft)       = traits%narea_structural
            lma(pft,jpngr)              = traits%lma
            nmass(pft)                  = traits%nmass

            r_cton_leaf(pft,jpngr)      = traits%r_cton_leaf
            r_ntoc_leaf(pft,jpngr)      = traits%r_ntoc_leaf

            !-------------------------------------------------------------------  
            ! Update fpc_grid and fapar_ind (not lai_ind)
            !-------------------------------------------------------------------  
            call update_fpc_grid( pft, jpngr )

            ! write(0,*) '--- in main allocation'
            ! write(0,*) 'dcleaf           ', dcleaf(pft)
            ! write(0,*) 'lai_ind          ', lai_ind(pft,jpngr)
            ! write(0,*) 'cleaf            ', pleaf(pft,jpngr)%c%c12

            ! ! write(0,*) 'actual nleaf   ', pleaf(pft,jpngr)%n%n14
            ! ! write(0,*) 'traits nleaf   ', traits%narea * lai_ind(pft,jpngr)

            ! ! write(0,*) 'traits, narea', traits%narea
            ! ! write(0,*) 'actual, narea', pleaf(pft,jpngr)%n%n14 / lai_ind(pft,jpngr)

            ! ! write(0,*) 'lma                 ', traits%lma
            ! ! write(0,*) 'actual, Carea (=LMA)', pleaf(pft,jpngr)%c%c12 / lai_ind(pft,jpngr)
            ! write(0,*) 'actual nleaf   ', pleaf(pft,jpngr)%n%n14
            ! write(0,*) 'traits nleaf   ', traits%narea * lai_ind(pft,jpngr)

            ! write(0,*) 'traits, narea', traits%narea 
            ! write(0,*) 'actual, narea', pleaf(pft,jpngr)%n%n14 / lai_ind(pft,jpngr)

            ! write(0,*) 'traits lma          ', traits%lma
            ! write(0,*) 'actual, Carea (=LMA)', pleaf(pft,jpngr)%c%c12 / lai_ind(pft,jpngr)


            ! write(0,*) 'narea_metabolic  ', narea_metabolic(pft,jpngr)
            ! write(0,*) 'narea_structural ', narea_structural(pft,jpngr)
            ! write(0,*) 'lma              ', lma(pft,jpngr) 
            ! write(0,*) 'nmass(%)         ', nmass(pft,jpngr) * 100

            ! write(0,*) 'traits, cton_leaf', traits%r_cton_leaf
            ! write(0,*) 'actual, cton_leaf', cton( pleaf(pft,jpngr), default=0.0 )

            ! write(0,*) 'traits, lma      ', traits%lma
            ! write(0,*) 'actual, lma      ', pleaf(pft,jpngr)%c%c12 / lai_ind(pft,jpngr)

            ! stop

            !-------------------------------------------------------------------
            ! ROOT ALLOCATION
            !-------------------------------------------------------------------
            ! write(0,*) 'proot before ', proot(pft,jpngr)
            ! write(0,*) 'plabl before ', plabl(pft,jpngr)
            call allocate_root( &
              proot(pft,jpngr)%c%c12, proot(pft,jpngr)%n%n14, &
              plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
              pft, dcroot(pft), dnroot(pft) &
              )
            ! write(0,*) 'proot after  ', proot(pft,jpngr)
            ! write(0,*) 'plabl after  ', plabl(pft,jpngr)
            ! stop 
            ! write(0,*) 'dcleaf, dnleaf', dcleaf(pft), dnleaf(pft)
            ! write(0,*) 'dcroot, dnleaf', dcroot(pft), dnroot(pft)

            ! SUMMARY OF THIS DAY'S ALLOCATION
            ! write(0,*) 'dcleaf, dnleaf ', dcleaf(pft), dnleaf(pft)
            ! write(0,*) 'cleaf, nleaf   ', pleaf(pft,jpngr)
            ! write(0,*) 'croot, nroot   ', proot(pft,jpngr)
            ! write(0,*) 'lai_ind        ', lai_ind(pft,jpngr)
            ! write(0,*) 'C:N in leaves  ', cton( pleaf(pft,jpngr), default=0.0 )
            ! write(0,*) 'clabl, nlabl   ', plabl(pft,jpngr)
            ! if (doy>250) stop 'do beni, end of year'

          else
            !------------------------------------------------------------------
            ! Prescribe constant root:shoot ratio (in terms of C mass).
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by 
            ! all labile N, discounted by the yield factor.
            !------------------------------------------------------------------
            if (pleaf(pft,jpngr)%c%c12.eq.0.0) then
              ! initial guess based on Taylor approximation of Cleaf and Nleaf function around cleaf=0
              r_cton_leaf(pft,jpngr) = get_rcton_init( solar%meanmppfd(:), mactnv_unitiabs(:) )
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

            call allocate_leaf( &
              dcleaf(pft), &
              pleaf(pft,jpngr)%c%c12, pleaf(pft,jpngr)%n%n14, &
              plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
              solar%meanmppfd(:), mactnv_unitiabs(:), &
              lai_ind(pft,jpngr), dnleaf(pft) &
              )

            !-------------------------------------------------------------------  
            ! Update leaf traits
            !-------------------------------------------------------------------  
            traits                      = get_leaftraits( lai_ind(pft,jpngr), solar%meanmppfd(:), mactnv_unitiabs(:) )
            narea(pft)                  = traits%narea
            narea_metabolic(pft)        = traits%narea_metabolic
            narea_structural(pft)       = traits%narea_structural
            lma(pft,jpngr)              = traits%lma
            nmass(pft)                  = traits%nmass

            r_cton_leaf(pft,jpngr)      = traits%r_cton_leaf
            r_ntoc_leaf(pft,jpngr)      = traits%r_ntoc_leaf

            !-------------------------------------------------------------------  
            ! Update fpc_grid and fapar_ind (not lai_ind)
            !-------------------------------------------------------------------  
            call update_fpc_grid( pft, jpngr )

            !-------------------------------------------------------------------
            ! ROOT ALLOCATION
            !-------------------------------------------------------------------
            call allocate_root( &
              proot(pft,jpngr)%c%c12, proot(pft,jpngr)%n%n14, &
              plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14, &
              pft, dcroot(pft), dnroot(pft) &
              )

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

    ! if (doy==39) stop 'on day 39'

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
    use _gpp, only: calc_dgpp, calc_drd
    use _nuptake, only: calc_dnup, outtype_calc_dnup, calc_cexu
    use _npp, only: calc_resp_maint, r_root
    use _vegdynamics, only: get_fapar, get_fpc_grid
    use _findroot_fzeroin

    ! arguments
    real, intent(in)              :: mydcleaf
    logical, intent(in), optional :: verbose

    ! function return variable
    real, intent(out) :: eval

    ! local variables
    real                    :: dclabl, dnlabl 
    real                    :: mycrownarea
    real                    :: mynind
    real                    :: mlue
    real                    :: dppfd
    real                    :: mrd_unitiabs
    real, dimension(nmonth) :: meanmppfd
    real                    :: ninorg
    real, dimension(nmonth) :: nv
    integer                 :: pft
    integer                 :: usemoy
    real                    :: soiltemp

    real :: mydcroot
    real :: mydnleaf
    real :: mydnroot
    real :: cleaf
    real :: nleaf
    real :: croot
    real :: nroot
    real :: clabl
    real :: nlabl
    real :: mylai
    real :: myfapar_ind
    real :: myfpc_grid
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

    type(outtype_zeroin)    :: out_zeroin
    ! type(outtype_calc_dnup) :: out_calc_dnup
    real, dimension(2) :: out_calc_dnup

    write(0,*) '--- in eval_imbalance with mydcleaf=', mydcleaf

    ! Copy to local variables for shorter writing
    cleaf        = state_eval_imbalance%pleaf%c%c12
    nleaf        = state_eval_imbalance%pleaf%n%n14
    croot        = state_eval_imbalance%proot%c%c12
    nroot        = state_eval_imbalance%proot%n%n14
    clabl        = state_eval_imbalance%plabl%c%c12
    nlabl        = state_eval_imbalance%plabl%n%n14
    mycrownarea  = state_eval_imbalance%crownarea
    mynind       = state_eval_imbalance%nind
    mlue         = state_eval_imbalance%mlue         
    dppfd        = state_eval_imbalance%dppfd        
    mrd_unitiabs = state_eval_imbalance%mrd_unitiabs 
    meanmppfd(:) = state_eval_imbalance%meanmppfd(:)    
    ninorg       = state_eval_imbalance%pninorg%n14 
    pft          = state_eval_imbalance%pft
    nv(:)        = state_eval_imbalance%nv(:)
    usemoy       = state_eval_imbalance%usemoy
    soiltemp     = state_eval_imbalance%soiltemp

    write(0,*) '-------------'
    write(0,*) 'BEFORE presumed allocation:'
    write(0,*) 'cleaf', cleaf
    write(0,*) 'nleaf', nleaf
    write(0,*) 'croot', croot
    write(0,*) 'nroot', nroot
    write(0,*) 'clabl', clabl
    write(0,*) 'nlabl', nlabl
    write(0,*) '-------------'

    !-------------------------------------------------------------------
    ! LEAF ALLOCATION
    !-------------------------------------------------------------------
    call allocate_leaf( mydcleaf, cleaf, nleaf, clabl, nlabl, meanmppfd(:), nv(:), mylai, mydnleaf )

    ! Update foliar projective cover (=fAPAR) and fractional plant cover
    ! to calculate GPP below
    myfapar_ind = get_fapar( mylai )
    myfpc_grid = get_fpc_grid( mycrownarea, mynind, myfapar_ind ) 

    write(0,*) '-------------'
    write(0,*) 'AFTER presumed leaf allocation:'
    write(0,*) 'cleaf', cleaf
    write(0,*) 'nleaf', nleaf
    write(0,*) 'clabl', clabl
    write(0,*) 'nlabl', nlabl
    write(0,*) 'mylai', mylai
    write(0,*) 'mydnleaf', mydnleaf
    write(0,*) 'myfapar_ind', myfapar_ind
    write(0,*) 'myfpc_grid', myfpc_grid
    write(0,*) '-------------'

    !-------------------------------------------------------------------
    ! ROOT ALLOCATION
    !-------------------------------------------------------------------
    call allocate_root( croot, nroot, clabl, nlabl, pft, mydcroot, mydnroot )

    write(0,*) '-------------'
    write(0,*) 'AFTER presumed root allocation:'
    write(0,*) 'croot', croot
    write(0,*) 'nroot', nroot
    write(0,*) '-------------'

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

    ! if (present(verbose)) then
    !   write(0,*) '--- in eval_imbalance:'
    !   write(0,*) 'cleaf      ', cleaf
    !   write(0,*) 'mylai      ', mylai
    !   write(0,*) 'mycrownarea', mycrownarea
    !   write(0,*) 'mynind     ', mynind
    !   write(0,*) 'myfapar_ind  ', myfapar_ind
    !   write(0,*) 'myfpc_grid ', myfpc_grid

    !   write(0,*) 'cleaf',cleaf
    !   write(0,*) 'croot',croot
    !   write(0,*) 'fapar',fapar ! nicht ganz ok!
    !   write(0,*) 'mlue',mlue ! ok
    !   write(0,*) 'dppfd',dppfd  ! ok
    !   write(0,*) 'ninorg',ninorg
    !   ! write(0,*) 'mrd_unitiabs',mrd_unitiabs ! ok
    !   ! write(0,*) 'meanmppfd',meanmppfd ! ok
    ! end if

    gpp           = calc_dgpp( myfpc_grid, mlue, dppfd )
    rd            = calc_drd(  myfpc_grid, mrd_unitiabs, meanmppfd(usemoy) )
    exu           = calc_cexu( croot, cleaf )
    mresp_root    = calc_resp_maint( croot, r_root )
    
    dc            = gpp - rd - mresp_root - exu
    
    out_calc_dnup = calc_dnup( exu, ninorg, soiltemp )
    dn            = sum( out_calc_dnup )


    ! if (mydcroot>0.0) write(0,*) 'dcleaf/dcroot, dn ', mydcleaf/mydcroot, dn


    ! dn            = out_calc_dnup%dnup_act +  out_calc_dnup%dnup_fix

    ! write(0,*) 'exu ', exu
    ! write(0,*) 'ninorg', ninorg
    ! write(0,*) 'soiltemp', soiltemp
    ! write(0,*) 'out_calc_dnup%dnup_act', out_calc_dnup%dnup_act
    ! write(0,*) 'out_calc_dnup%dnup_fix', out_calc_dnup%dnup_fix

    ! if (out_calc_dnup%dnup_fix>0.0) stop 'finally fixing'

    ! write(0,*) 'AFTER presumed allocation:'
    ! write(0,*) 'cleaf      ', cleaf
    ! write(0,*) 'croot      ', croot
    ! write(0,*) 'myfpc_grid ', myfpc_grid
    ! write(0,*) 'dc         ', dc
    ! stop
    ! write(0,*) 'nleaf', nleaf
    ! write(0,*) 'nroot', nroot
    ! write(0,*) 'clabl', clabl
    ! write(0,*) 'nlabl', nlabl


    ! if (present(verbose)) then
      write(0,*) '-------------'
      write(0,*) 'gpp                ', gpp
      write(0,*) 'rd                 ', rd
      write(0,*) 'exu                ', exu
      write(0,*) 'mresp_root         ', mresp_root
      write(0,*) 'dc                 ', dc
      write(0,*) 'dn                 ', dn
      write(0,*) 'current   plant C:N', ( cleaf + croot ) / ( nleaf + nroot )
      write(0,*) 'tomorrows plant C:N', growtheff * (dc + clabl) / (dn + nlabl)
      write(0,*) '-------------'
      ! stop
    ! end if

    !-------------------------------------------------------------------
    ! EVALUATION QUANTITY - IS MINIMISED BY OPTIMISATION
    !-------------------------------------------------------------------
    ! Evaluation quantity is the difference between the 
    ! C:N ratio of new assimilates and the C:N ratio 
    ! of the whole plant after allocation.

    if ((dn + nlabl)==0.0) then
      eval = -999.0
    else
      !     |---------------------------------------|  |------------------------------------|
      eval = growtheff * (dc + clabl) / (dn + nlabl) - ( cleaf + croot ) / ( nleaf + nroot )
      !     |---------------------------------------|  |------------------------------------|
      !     |lab. pool C:N ratio after acq. nxt. day|  | current whole-plant C:N ratio      |
      !     |---------------------------------------|  |------------------------------------|
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
    use _classdefs
    ! use _params_core, only: nmonth
    use _params_modl, only: growtheff

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

      ! write(0,*) 'cleaf before ', cleaf 
      cleaf  = cleaf + mydcleaf
      ! write(0,*) 'mydcleaf     ', mydcleaf 

      ! write(0,*) 'cleaf = 0.5', get_lai( 0.5  , meanmppfd(:), nv(:) )
      ! write(0,*) 'cleaf = 100', get_lai( 100.0, meanmppfd(:), nv(:))
      ! stop

      ! Calculate LAI as a function of leaf C
      lai = get_lai( cleaf, meanmppfd(:), nv(:) )
      ! write(0,*) 'in allocate_leaf: cleaf, lai :', cleaf, lai
      ! write(0,*) 'mydcleaf ', mydcleaf 
      ! ! stop

      ! calculate canopy-level leaf N as a function of LAI
      nleaf0   = nleaf      
      nleaf    = get_canopy_leaf_n( lai, meanmppfd(:), nv(:) )
      mydnleaf = nleaf - nleaf0
      ! if (mydnleaf>0.0) then
      !   write(0,*) 'mydcleaf/dnleaf ', mydcleaf/mydnleaf 
      ! end if

      ! subtract from labile pool, making sure pool does not get negative
      dclabl = min( clabl, 1.0 / growtheff * mydcleaf )
      dnlabl = min( nlabl, mydnleaf )
      ! write(0,*) 'C balance after reducing labile:', (clabl - dclabl )
      ! write(0,*) 'N balance after reducing labile:', (nlabl - dnlabl )
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
    use _classdefs
    ! use _params_core, only: nmonth
    use _params_modl, only: growtheff, r_cton_root, r_ntoc_root

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
      mydcroot = min( growtheff * clabl, r_cton_root(pft) * nlabl )
      mydnroot = min( mydcroot * r_ntoc_root(pft), nlabl )

      dclabl = min( clabl, 1.0 / growtheff * mydcroot )
      dnlabl = min( nlabl, mydnroot)
      ! write (0,*) 'C balance after reducing labile:', (clabl - dclabl)
      ! write (0,*) 'N balance after reducing labile:', (nlabl - dnlabl)
      if ( (dclabl - clabl) > 1e-8 ) stop 'trying to remove too much from labile pool: root C'
      if ( (dnlabl - nlabl) > 1e-8 ) stop 'trying to remove too much from labile pool: root N'
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
    ! use _params_core, only: nmonth
    use _findroot_fzeroin

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
    real, intent(out) :: lai

    ! local variables
    real :: maxnv

    if (cleaf>0.0) then
      ! Metabolic N is predicted and is optimised at a monthly time scale. 
      ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
      ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost.
      ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
      ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
      maxnv = maxval( meanmppfd(:) * nv(:) )
      ! write(0,*) 'maxnv ', maxnv
      ! stop

      ! write(0,*) 'what is LAI for Cleaf=', cleaf
      ! write(0,*) 'maxnv ', maxnv

      ! Update state. This derived-type variable is "global" within this module
      state_mustbe_zero_for_lai%cleaf = cleaf
      state_mustbe_zero_for_lai%maxnv = maxnv

      ! Calculate initial guess for LAI (always larger than actual LAI)
      lower = 0.0 ! uninformed lower bound
      upper = 1.0 / ( c_molmass * ncw_min * r_ctostructn_leaf ) * cleaf
      ! upper = 20.0
      ! write(0,*) 'upper ', upper

      ! write(0,*) 'lower =', lower
      ! test = mustbe_zero_for_lai( lower ) 
      ! write(0,*) '=> test =', test

      ! write(0,*) 'upper =', upper
      ! test = mustbe_zero_for_lai( upper ) 
      ! write(0,*) '=> test =', test

      ! call function zeroin to find root (value of LAI for which evaluation expression is zero)
      abserr=100.0*XMACHEPS !*10e5
      relerr=1000.0*XMACHEPS !*10e5

      ! write(0,*) 'abserr', abserr
      ! write(0,*) 'relerr', relerr
      ! stop 'here'

      ! write(0,*) '*** finding root of mustbe_zero_for_lai ***'
      out_zeroin = zeroin( mustbe_zero_for_lai, abserr, relerr, nmax, lower, upper )
      if ( out_zeroin%error /= 0 ) then
        lai = 0.0
        write(0,*) 'error code', out_zeroin%error
        stop 'zeroin for mustbe_zero_for_lai() failed'
      else
        lai = out_zeroin%root
      end if

      ! write(0,*) 'out_zeroin', out_zeroin
      ! write(0,*) 'cleaf', cleaf
      ! write(0,*) 'lai', lai
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
    ! Cleaf = c_molmass * r_ctostructn_leaf * N_canop_cellwall
    ! N_canop_cellwall = LAI * ncw_min + nv * Iabs * r_n_cw_v
    ! Iabs = meanmppfd * (1-exp( -kbeer * LAI))
    ! ==> Cleaf = f(LAI) = c_molmass * r_ctostructn_leaf * [ meanmppfd * (1-exp( -kbeer * LAI)) * nv * r_n_cw_v + LAI * ncw_min ]
    ! ==> LAI = f(Cleaf) leads to inhomogenous equation. Therefore apply root finding algorithm so that:
    ! 0 = cleaf / ( c_molmass * r_ctostructn_leaf ) - meanmppfd * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * nv * r_n_cw_v - mylai * ncw_min
    !---------------------------------------------------------
    use _params_modl, only: kbeer
    use _vegdynamics, only: get_fapar

    ! arguments
    real, intent(in) :: mylai

    ! function return value
    real, intent(out) :: mustbe_zero

    ! local variables
    real :: mycleaf
    real :: mymaxnv

    ! write(0,*) '--- in mustbe_zero_for_lai with mydcleaf=', mylai

    ! Read from updated state. This derived-type variable is "global" within this module
    mycleaf = state_mustbe_zero_for_lai%cleaf
    mymaxnv = state_mustbe_zero_for_lai%maxnv

    ! write(0,*) '----------'
    ! write(0,*) 'inside mustbe_zero_for_lai: '
    ! write(0,*) 'mylai', mylai
    ! write(0,*) 'mycleaf', mycleaf
    ! write(0,*) 'mymaxnv', mymaxnv

    ! mustbe_zero = cleaf / ( c_molmass * r_ctostructn_leaf ) - meanmppfd * nv * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * r_n_cw_v - mylai * ncw_min
    ! mustbe_zero = cleaf / ( c_molmass * r_ctostructn_leaf ) - maxnv * ( 1.0 - exp( -1.0 * kbeer * mylai ) ) * r_n_cw_v - mylai * ncw_min
    mustbe_zero = mycleaf / ( c_molmass * r_ctostructn_leaf ) - mymaxnv * get_fapar( mylai ) * r_n_cw_v - mylai * ncw_min

    ! write(0,*) 'mustbe_zero                           ', mustbe_zero
    ! write(0,*) '-------------'


  end function mustbe_zero_for_lai


  function get_rcton_init( meanmppfd, nv ) result( rcton )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    ! Cleaf = c_molmass * r_ctostructn_leaf * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * r_n_cw_v + LAI * ncw_min ]
    ! Nleaf = n_molmass * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * (r_n_cw_v + 1) + LAI * ncw_min ]
    ! linearization around LAI = 0 ==> (1-exp(-k*L)) ~= k*L
    ! ==> Cleaf ~= LAI * c_molmass * r_ctostructn_leaf * ( meanmppfd * kbeer * nv * r_n_cw_v + ncw_min )
    ! ==> Nleaf ~= LAI * n_molmass * ( meanmppfd * kbeer * nv * (r_n_cw_v + 1) + ncw_min )
    ! r_cton = Cleaf / Nleaf
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth
    use _params_modl, only: kbeer

    ! arguments
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real, intent(out) :: rcton

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )
    rcton = ( c_molmass * r_ctostructn_leaf * ( maxnv * kbeer * r_n_cw_v + ncw_min )) / ( n_molmass * ( maxnv * kbeer * (r_n_cw_v + 1.0) + ncw_min ) )

    ! rcton = ( c_molmass * r_ctostructn_leaf * ( meanmppfd * kbeer * nv * r_n_cw_v + ncw_min )) / ( n_molmass * ( meanmppfd * kbeer * nv * (r_n_cw_v + 1) + ncw_min ) )

  end function get_rcton_init


  function get_canopy_leaf_n_metabolic( mylai, meanmppfd, nv ) result( mynleaf_metabolic )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! LAI * n_metabolic = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth
    use _vegdynamics, only: get_fapar

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real, intent(out) :: mynleaf_metabolic  ! mol N 

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )

    mynleaf_metabolic = maxnv * get_fapar( mylai )

  end function get_canopy_leaf_n_metabolic


  function get_canopy_leaf_n_structural( mylai, mynleaf_metabolic ) result( mynleaf_structural )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! LAI * n_structural = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth

    ! arguments
    real, intent(in) :: mylai
    real, intent(in) :: mynleaf_metabolic

    ! function return variable
    real, intent(out) :: mynleaf_structural  ! mol N 

    mynleaf_structural = mynleaf_metabolic * r_n_cw_v + mylai * ncw_min

    ! write(0,*) '--- in get_canopy_leaf_n_structural'
    ! write(0,*) 'mylai ', mylai
    ! write(0,*) 'mynleaf_metabolic ', mynleaf_metabolic
    ! write(0,*) 'r_n_cw_v ', r_n_cw_v
    ! write(0,*) 'ncw_min ', ncw_min
    ! write(0,*) 'mynleaf_structural ', mynleaf_structural
    ! write(0,*) '-------------------------------'

  end function get_canopy_leaf_n_structural


  function get_canopy_leaf_n( mylai, meanmppfd, nv ) result( mynleaf )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial guess based on Taylor approximation of 
    ! Cleaf and Nleaf function around cleaf=0.
    ! Nleaf = LAI * (n_metabolic + n_cellwall) * n_molmass
    ! LAI * n_metabolic = nv * Iabs
    ! Iabs = meanmppfd * (1-exp(-kbeer*LAI))
    ! LAI * n_cellwall = LAI * (ncw_min + r_n_cw_v * n_metabolic)
    ! ==> Nleaf = n_molmass * [ meanmppfd * (1-exp(-kbeer*LAI)) * nv * (r_n_cw_v + 1) + LAI * ncw_min ]
    !----------------------------------------------------------------
    ! use _params_core, only: nmonth

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real, intent(out) :: mynleaf ! g N

    ! local variables
    real :: nleaf_metabolic   ! mol N m-2
    real :: nleaf_structural  ! mol N m-2

    nleaf_metabolic  = get_canopy_leaf_n_metabolic(  mylai, meanmppfd, nv )
    nleaf_structural = get_canopy_leaf_n_structural( mylai, nleaf_metabolic )
    mynleaf          = n_molmass * ( nleaf_metabolic + nleaf_structural )

    ! write(0,*) '--- in get_canopy_leaf_n'
    ! write(0,*) 'nleaf_metabolic ', nleaf_metabolic
    ! write(0,*) 'nleaf_structural ', nleaf_structural
    ! write(0,*) 'mynleaf ', mynleaf
    ! write(0,*) '-------------------------------'

    ! mynleaf = n_molmass * ( maxnv * get_fapar( mylai ) * ( 1.0 + r_n_cw_v ) + mylai * ncw_min )

  end function get_canopy_leaf_n


  function get_leaftraits( mylai, meanmppfd, nv ) result( traits )
    !////////////////////////////////////////////////////////////////
    ! Calculates leaf traits based on (predicted) metabolic Narea and
    ! (prescribed) parameters that relate structural to metabolic
    ! Narea and Carea to structural Narea:
    ! Narea_metabolic  = predicted
    ! Narea_structural = a + b * Narea_metabolic
    ! Carea            = c * Narea_structural
    !----------------------------------------------------------------
    use _params_core, only: c_content_of_biomass
    use _vegdynamics, only: get_fapar

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    type(leaftraits_type) :: traits

    ! local variables
    real :: mynarea_metabolic_canop   ! mol N m-2-ground
    real :: mynarea_structural_canop  ! mol N m-2-ground

    mynarea_metabolic_canop  = get_canopy_leaf_n_metabolic(  mylai, meanmppfd(:), nv(:) )     ! mol N m-2-ground    
    mynarea_structural_canop = get_canopy_leaf_n_structural( mylai, mynarea_metabolic_canop ) ! mol N m-2-ground
    
    traits%narea_metabolic  = n_molmass * mynarea_metabolic_canop / mylai   ! g N m-2-leaf
    traits%narea_structural = n_molmass * mynarea_structural_canop / mylai  ! g N m-2-leaf
    traits%narea            = n_molmass * ( mynarea_metabolic_canop + mynarea_structural_canop ) / mylai ! g N m-2-leaf
    traits%lma              = c_molmass * r_ctostructn_leaf * mynarea_structural_canop / mylai 
    traits%nmass            = traits%narea / ( traits%lma / c_content_of_biomass )
    traits%r_cton_leaf      = traits%lma / traits%narea
    traits%r_ntoc_leaf      = 1.0 / traits%r_cton_leaf

    ! write(0,*) '--- in get_leaftraits'
    ! write(0,*) 'mylai                  ', mylai
    ! write(0,*) 'traits%narea_metabolic ', traits%narea_metabolic 
    ! write(0,*) 'traits%narea_structural', traits%narea_structural
    ! write(0,*) 'traits%narea           ', traits%narea           
    ! write(0,*) 'traits%lma             ', traits%lma             
    ! write(0,*) 'traits%nmass           ', traits%nmass           
    ! write(0,*) 'traits%r_cton_leaf     ', traits%r_cton_leaf     
    ! write(0,*) 'traits%r_ntoc_leaf     ', traits%r_ntoc_leaf     
    ! write(0,*) '-------------------------------'
    ! stop 

  end function get_leaftraits


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
