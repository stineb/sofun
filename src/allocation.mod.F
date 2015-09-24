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
  use _params_core, only: npft, maxgrid
  use _classdefs

  !------------------------------------------------------------------
  ! Define derived type to store current state variables.
  ! Is defined within all subroutines in this module.
  ! This avoids having to pass them as arguments.
  !------------------------------------------------------------------
  type statetype_eval_imbalance
    type(orgpool)  :: pleaf
    type(orgpool)  :: proot
    type(orgpool)  :: plabl
    real           :: r_ntoc_leaf
    real           :: sla
    real           :: crownarea
    real           :: mlue
    real           :: dppfd
    real           :: mrd_unitiabs
    real           :: meanmppfd
    type(nitrogen) :: pninorg
    real           :: pft 
  end type statetype_eval_imbalance

  type(statetype_eval_imbalance) :: state_eval_imbalance

  ! MODULE-SPECIFIC VARIABLES
  real, dimension(npft) :: dcleaf
  real, dimension(npft) :: dnleaf
  real, dimension(npft) :: dcroot
  real, dimension(npft) :: dnroot

  ! OUTPUT VARIABLES
  real, dimension(npft,maxgrid) :: outaCalloc
  real, dimension(npft,maxgrid) :: outaNalloc


contains

  subroutine allocation_daily( jpngr, doy, moy, dom )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use _params_core, only: ndaymonth, ndayyear
    use _params_modl, only: r_cton_root, r_ntoc_root, growtheff, grass, &
      lu_category
    use _vars_core, only: pleaf, proot, plabl, pninorg, r_cton_leaf, &
      r_ntoc_leaf, sla, crownarea, drauto, dnpp
    use _gpp, only: mlue, mrd_unitiabs
    use _waterbal, only: dppfd, meanmppfd
    use _findroot_fzeroin

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy     ! day of year
    integer, intent(in) :: moy     ! month of year
    integer, intent(in) :: dom     ! day of month

    ! local variables
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
    real    :: eval_allroots
    real    :: abserr
    real    :: relerr
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
    if ( invocation > spinupyr_phaseinit_2 ) then
      flexalloc = .true.
    else
      flexalloc = .false.
    end if

    abserr=100.0*XMACHEPS
    relerr=1000.0*XMACHEPS

    do pft=1,npft

      lu = lu_category(pft)

      if (grass(pft)) then

        ! write(0,*) '--- allocation_daily, doy:',doy
        ! write(0,*)' pninorg(lu,jpngr)%n14',pninorg
        ! write(0,*)' plabl(lu,jpngr)',plabl

        if ( plabl(pft,jpngr)%c%c12>0.0 .and. plabl(pft,jpngr)%n%n14>0.0 ) then

          ! write(0,*) 'growing ...'
          if (flexalloc) then
            !------------------------------------------------------------------
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by all labile N,
            ! discounted by the yield factor.
            !------------------------------------------------------------------
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
            ! write(0,*) 'sla',sla
            ! write(0,*) 'mlue',mlue(moy) ! ok
            ! write(0,*) 'dppfd',dppfd(doy)  ! ok
            ! write(0,*) 'mrd_unitiabs',mrd_unitiabs(moy) ! ok
            ! write(0,*) 'meanmppfd',meanmppfd(moy) ! ok

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
            state_eval_imbalance%pleaf        = pleaf(pft,jpngr)
            state_eval_imbalance%proot        = proot(pft,jpngr)
            state_eval_imbalance%plabl        = plabl(pft,jpngr)
            state_eval_imbalance%r_ntoc_leaf  = r_ntoc_leaf(pft,jpngr)
            state_eval_imbalance%sla          = sla(pft,jpngr)
            state_eval_imbalance%crownarea    = crownarea(pft,jpngr)
            state_eval_imbalance%mlue         = mlue(usemoy)
            state_eval_imbalance%dppfd        = dppfd(usedoy)
            state_eval_imbalance%mrd_unitiabs = mrd_unitiabs(usemoy)
            state_eval_imbalance%meanmppfd    = meanmppfd(usemoy)
            state_eval_imbalance%pninorg      = pninorg(lu,jpngr)  ! the only data that is not yet available - use today's value 
            state_eval_imbalance%pft          = pft


            !------------------------------------------------------------------
            ! Optimisation by balanced growth
            ! Test I: Evaluate balance if all is put to roots.
            ! If C:N ratio of return is still greater than whole-plant C:N 
            ! ratio, then put all to roots.
            !------------------------------------------------------------------
            cont = .true.
            eval_allroots  = eval_imbalance( 0.0 )
            ! write(0,*) 'eval_allroots', eval_allroots  
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
              eval_allleaves = eval_imbalance( max_dc )
              ! write(0,*) 'eval_allleaves', eval_allleaves  
              if (eval_allleaves < 0.0) then
                dcleaf(pft) = max_dc
                cont = .false.
                ! write(0,*) '* putting all to leaves *'
              end if
            end if

            !------------------------------------------------------------------
            ! Optimum is between 0.0 and max_dc. Find root of function 
            ! 'eval_imbalance()' in the interval [0.0, max_dc].
            !------------------------------------------------------------------
            if (cont) then
              ! write(0,*) '* finding root *'
              out_zeroin = zeroin( eval_imbalance, abserr, relerr, 100, 0.0 ,max_dc )
              if ( out_zeroin%error /= 0 ) then
                dcleaf(pft) = 0.0
              else
                dcleaf(pft) = out_zeroin%root
              end if
              ! write(0,*) 'no. of iterations', out_zeroin%niter
            end if

            !------------------------------------------------------------------
            ! xxx debug: project next-day's fluxes with optimal dcleaf, derived now
            !------------------------------------------------------------------
            ! write(0,*) 'dcleaf(pft)', dcleaf(pft)
            test = eval_imbalance( dcleaf(pft), .true. )
            ! stop

            !------------------------------------------------------------------
            ! Allocate to roots and leaves and extract from labile pool based 
            ! on 'dcleaf(pft)' determined above.
            !------------------------------------------------------------------
            dnleaf(pft) = dcleaf(pft) * r_ntoc_leaf(pft,jpngr)

            ! write(0,*) 'dcleaf(pft)', dcleaf(pft)  
            ! write(0,*) 'dnleaf(pft)', dnleaf(pft)  

            ! to ensure plabl not getting negative (within numerical imprecision)
            dclabl = min( plabl(pft,jpngr)%c%c12, 1.0 / growtheff * dcleaf(pft) )
            dnlabl = min( plabl(pft,jpngr)%n%n14, dnleaf(pft) )

            plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - dclabl
            plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnlabl

            pleaf(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12 + dcleaf(pft)
            pleaf(pft,jpngr)%n%n14 = pleaf(pft,jpngr)%n%n14 + dnleaf(pft)

            dcroot(pft) = min( growtheff * plabl(pft,jpngr)%c%c12, r_cton_root(pft) * plabl(pft,jpngr)%n%n14 )
            dnroot(pft) = dcroot(pft) * r_ntoc_root(pft)

            ! write(0,*) 'dcroot(pft):dcleaf(pft)', dcroot(pft) / dcleaf(pft)

            ! to ensure plabl not getting negative (within numerical imprecision)
            dclabl = min( plabl(pft,jpngr)%c%c12, 1.0 / growtheff * dcroot(pft) )
            dnlabl = min( plabl(pft,jpngr)%n%n14, dnroot(pft) )

            plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - dclabl
            plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnlabl

            proot(pft,jpngr)%c%c12 = proot(pft,jpngr)%c%c12 + dcroot(pft)
            proot(pft,jpngr)%n%n14 = proot(pft,jpngr)%n%n14 + dnroot(pft)

            ! add growth respiration to autotrophic respiration and substract from NPP
            ! (note that NPP is added to plabl in and growth resp. is implicitly removed
            ! from plabl above)
            drauto(pft)   = drauto(pft)     + ( 1.0 - growtheff ) * ( dcleaf(pft) + dcroot(pft) )
            dnpp(pft)%c12 = dnpp(pft)%c12   - ( 1.0 - growtheff ) * ( dcleaf(pft) + dcroot(pft) ) 

            ! write(0,*) 'dcroot(pft)'
            ! write(0,*) 'plabl(pft,jpngr)', plabl(pft,jpngr)
            ! if (doy==48) stop 'here'
        
          else
            !------------------------------------------------------------------
            ! Prescribe constant root:shoot ratio (in terms of C mass).
            ! Calculate maximum C allocatable based on current labile pool size.
            ! Maximum is the lower of all labile C and the C to be matched by 
            ! all labile N, discounted by the yield factor.
            !------------------------------------------------------------------
            max_dc = growtheff * plabl(pft,jpngr)%c%c12
            max_dc_n_constraint = plabl(pft,jpngr)%n%n14 / ( frac_shoot * r_ntoc_leaf(pft,jpngr) + ( 1.0 - frac_shoot) * r_ntoc_root(pft) )
            ! write(0,*) "clabl, nlabl", plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14 
            ! write(0,*) "max_dc, max_dc_n_constraint", max_dc, max_dc_n_constraint
            max_dc = min( max_dc, max_dc_n_constraint )

            !------------------------------------------------------------------
            ! Allocate to roots and leaves and extract from labile pool based 
            ! on 'dcleaf(pft)' determined above.
            !------------------------------------------------------------------
            ! allocation to leaves is prescribed
            dcleaf(pft) = max_dc * frac_shoot
            dnleaf(pft) = dcleaf(pft) * r_ntoc_leaf(pft,jpngr)

            dcroot(pft) = max_dc * ( 1.0 - frac_shoot )
            dnroot(pft) = dcroot(pft) * r_ntoc_root(pft)

            ! to ensure plabl not getting negative (within numerical imprecision)
            dclabl = min( plabl(pft,jpngr)%c%c12, 1.0 / growtheff * ( dcleaf(pft) + dcroot(pft) ) )
            dnlabl = min( plabl(pft,jpngr)%n%n14, ( dnleaf(pft) + dnroot(pft) ) )

            ! write(0,*) 'plabl  ', plabl(pft,jpngr) 
            ! write(0,*) 'dcleaf ', dcleaf(pft)
            ! write(0,*) 'dcroot ', dcroot(pft)
            ! write(0,*) 'dclabl ', dclabl
            ! write(0,*) 'dnlabl ', dnlabl

            plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - dclabl
            plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnlabl

            pleaf(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12 + dcleaf(pft)
            pleaf(pft,jpngr)%n%n14 = pleaf(pft,jpngr)%n%n14 + dnleaf(pft)

            proot(pft,jpngr)%c%c12 = proot(pft,jpngr)%c%c12 + dcroot(pft)
            proot(pft,jpngr)%n%n14 = proot(pft,jpngr)%n%n14 + dnroot(pft)

            ! add growth respiration to autotrophic respiration and substract from NPP
            ! (note that NPP is added to plabl in and growth resp. is implicitly removed
            ! from plabl above)
            drauto(pft)   = drauto(pft)     + ( 1.0 - growtheff ) * ( dcleaf(pft) + dcroot(pft) )
            dnpp(pft)%c12 = dnpp(pft)%c12   - ( 1.0 - growtheff ) * ( dcleaf(pft) + dcroot(pft) ) 

          end if

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
    use _gpp, only: calc_dgpp, calc_drd
    use _nuptake, only: calc_dnup, exurate
    use _npp, only: calc_resp_maint, r_root
    use _vegdynamics, only: get_fapar

    ! arguments
    real, intent(in)  :: mydcleaf
    logical, intent(in), optional :: verbose

    ! function return variable
    real, intent(out) :: eval

    ! local variables
    real    :: r_ntoc_leaf
    real    :: sla
    real    :: crownarea
    real    :: mlue
    real    :: dppfd
    real    :: mrd_unitiabs
    real    :: meanmppfd
    real    :: ninorg
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

    ! write(0,*) '--- in eval_imbalance with mydcleaf=',mydcleaf

    ! shorter ...
    cleaf        = state_eval_imbalance%pleaf%c%c12
    nleaf        = state_eval_imbalance%pleaf%n%n14
    croot        = state_eval_imbalance%proot%c%c12
    nroot        = state_eval_imbalance%proot%n%n14
    clabl        = state_eval_imbalance%plabl%c%c12
    nlabl        = state_eval_imbalance%plabl%n%n14
    r_ntoc_leaf  = state_eval_imbalance%r_ntoc_leaf
    sla          = state_eval_imbalance%sla          
    crownarea    = state_eval_imbalance%crownarea
    mlue         = state_eval_imbalance%mlue         
    dppfd        = state_eval_imbalance%dppfd        
    mrd_unitiabs = state_eval_imbalance%mrd_unitiabs 
    meanmppfd    = state_eval_imbalance%meanmppfd    
    ninorg       = state_eval_imbalance%pninorg%n14 
    pft          = state_eval_imbalance%pft          

    ! write(0,*) 'before presumed allocation:'
    ! write(0,*) 'cleaf', cleaf
    ! write(0,*) 'nleaf', nleaf
    ! write(0,*) 'croot', croot
    ! write(0,*) 'nroot', nroot
    ! write(0,*) 'clabl', clabl
    ! write(0,*) 'nlabl', nlabl

    ! Allocate
    mydnleaf = mydcleaf * r_ntoc_leaf
    clabl  = clabl - 1.0 / growtheff * mydcleaf
    nlabl  = nlabl - mydnleaf
    cleaf  = cleaf + mydcleaf
    nleaf  = nleaf + mydnleaf
    
    mydcroot = min( growtheff * clabl, r_cton_root(pft) * nlabl )
    mydnroot = mydcroot * r_ntoc_root(pft)
    clabl  = clabl - 1.0 / growtheff * mydcroot
    nlabl  = nlabl - mydnroot
    croot  = croot + mydcroot
    nroot  = nroot + mydnroot

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

    ! continuous decay
    kcleaf = cleaf * k_decay_leaf(pft)
    knleaf = nleaf * k_decay_leaf(pft)
    kcroot = croot * k_decay_root(pft)
    knroot = nroot * k_decay_root(pft)

    cleaf = cleaf - kcleaf
    nleaf = nleaf - knleaf
    croot = croot - kcroot
    nroot = nroot - knroot

    ! write(0,*) 'B cleaf', cleaf

    ! Calculate next day's C and N return after assumed allocation (tissue turnover happens before!)
    lai   = cleaf * sla / crownarea
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
    rd         = calc_drd( fapar, mrd_unitiabs, meanmppfd )
    exu        = croot * exurate
    mresp_root = calc_resp_maint( croot, r_root)
    
    dc         = gpp - rd - mresp_root - exu
    dn         = calc_dnup( exu, ninorg )

    ! if (present(verbose)) then
    !   write(0,*) 'gpp                ', gpp
    !   write(0,*) 'dc                 ', dc
    !   write(0,*) 'dn                 ', dn
    !   write(0,*) 'current   plant C:N', ( cleaf + croot ) / ( nleaf + nroot )
    !   write(0,*) 'tomorrows plant C:N', growtheff * (dc + clabl) / (dn + nlabl)
    ! end if

    ! Evaluation quantity is the difference between the 
    ! C:N ratio of new assimilates and the C:N ratio 
    ! of the whole plant after allocation.
    eval = growtheff * (dc + clabl) / (dn + nlabl) - ( cleaf + croot ) / ( nleaf + nroot )

    ! write(0,*) 'new assimilates stoichiometry: ', growtheff * (dc + clabl) / (dn + nlabl)
    ! write(0,*) 'current stoichiometry:         ', ( cleaf + croot ) / ( nleaf + nroot )
    ! write(0,*) 'eval                           ', eval

  end function eval_imbalance


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

    ! C ALLOCATED TO GROWTH 
    filnam=trim(prefix)//'.a.calloc.out'
    open(303,file=filnam,err=999,status='unknown')

    ! N ALLOCATED TO GROWTH 
    filnam=trim(prefix)//'.a.nalloc.out'
    open(309,file=filnam,err=999,status='unknown')

    write(0,*) 'opening files for outaCalloc'

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
    outaCalloc(:,:) = 0.0
    outaNalloc(:,:) = 0.0

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

    outaCalloc(:,jpngr) = outaCalloc(:,jpngr) + dcleaf(:) + dcroot(:)
    outaNalloc(:,jpngr) = outaNalloc(:,jpngr) + dnleaf(:) + dnroot(:)

    ! write(0,*) 'collecting outaCalloc',outaCalloc

  end subroutine getout_daily_allocation


  subroutine writeout_ascii_allocation( year, spinup )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use _params_core, only: ndayyear, nlu
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

    write(303,999) itime, sum(outaCalloc(:,jpngr))
    write(309,999) itime, sum(outaNalloc(:,jpngr))

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_allocation

end module _allocation
