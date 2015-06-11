module _allocation

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


contains

  subroutine allocation_daily( jpngr, doy, moy, dom )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use _params_core, only: npft, ndaymonth
    use _params_modl, only: r_cton_root, r_ntoc_root, growtheff, grass, &
      lu_category
    use _vars_core, only: pleaf, proot, plabl, pninorg, r_cton_leaf, &
      r_ntoc_leaf, sla, crownarea
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
    real    :: dcleaf
    real    :: dnleaf
    logical :: cont          ! true if allocation to leaves (roots) is not 100% and not 0%
    real    :: max_dcleaf_n_constraint
    real    :: max_dcroot_n_constraint
    real    :: max_dc
    real    :: eval_allroots
    real    :: abserr
    real    :: relerr
    type(outtype_zeroin) :: out_zeroin

    ! xxx debug
    real    :: test

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

          !------------------------------------------------------------------
          ! Calculate maximum C allocatable based on current labile pool size.
          ! Maximum is the lower of all labile C and the C to be matched by all labile N,
          ! discounted by the yield factor.
          !------------------------------------------------------------------
          max_dcleaf_n_constraint = plabl(pft,jpngr)%n%n14 * r_cton_leaf(pft,jpngr)
          max_dcroot_n_constraint = plabl(pft,jpngr)%n%n14 * r_cton_root(pft)
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
            dcleaf = 0.0
            ! dcroot = max_dc
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
              dcleaf = max_dc
              ! dcroot = 0.0
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
              dcleaf = 0.0
            else
              dcleaf = out_zeroin%root
            end if
            ! write(0,*) 'no. of iterations', out_zeroin%niter
          end if

          !------------------------------------------------------------------
          ! xxx debug: project next-day's fluxes with optimal dcleaf, derived now
          !------------------------------------------------------------------
          ! write(0,*) 'dcleaf', dcleaf
          test = eval_imbalance( dcleaf, .true. )
          ! stop

          !------------------------------------------------------------------
          ! Allocate to roots and leaves and extract from labile pool based 
          ! on 'dcleaf' determined above.
          !------------------------------------------------------------------
          dnleaf = dcleaf * r_ntoc_leaf(pft,jpngr)

          ! write(0,*) 'dcleaf', dcleaf  
          ! write(0,*) 'dnleaf', dnleaf  

          plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - 1.0 / growtheff * dcleaf
          plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnleaf
          pleaf(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12 + dcleaf
          pleaf(pft,jpngr)%n%n14 = pleaf(pft,jpngr)%n%n14 + dnleaf

          ! write(0,*) 'a dcroot', dcroot  
          dcroot = min( growtheff * plabl(pft,jpngr)%c%c12, r_cton_root(pft) * plabl(pft,jpngr)%n%n14 )
          dnroot = dcroot * r_ntoc_root(pft)

          ! write(0,*) 'b dcroot', dcroot  
          ! write(0,*) 'dcleaf', dcleaf  
          ! write(0,*) 'dnroot', dnroot  
          ! write(0,*) 'dcroot:dcleaf', dcroot / dcleaf

          plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - 1.0 / growtheff * dcroot
          plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnroot
          proot(pft,jpngr)%c%c12 = proot(pft,jpngr)%c%c12 + dcroot
          proot(pft,jpngr)%n%n14 = proot(pft,jpngr)%n%n14 + dnroot

          ! write(0,*) 'dcroot'
          ! write(0,*) 'plabl(pft,jpngr)', plabl(pft,jpngr)
          ! if (doy==48) stop 'here'

        else

          ! write(0,*) 'not growing ...'

        end if

      else

        stop 'allocation_daily not implemented for trees'

      end if

    end do

    ! write(0,*) '--- END allocation_daily:'

  end subroutine allocation_daily


  function eval_imbalance( dcleaf, verbose ) result ( eval )
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
    real, intent(in)  :: dcleaf
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

    real :: dnleaf
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

    ! write(0,*) '--- in eval_imbalance with dcleaf=',dcleaf

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
    dnleaf = dcleaf * r_ntoc_leaf
    clabl  = clabl - 1.0 / growtheff * dcleaf
    nlabl  = nlabl - dnleaf
    cleaf  = cleaf + dcleaf
    nleaf  = nleaf + dnleaf
    
    dcroot = min( growtheff * clabl, r_cton_root(pft) * nlabl )
    dnroot = dcroot * r_ntoc_root(pft)
    clabl  = clabl - 1.0 / growtheff * dcroot
    nlabl  = nlabl - dnroot
    croot  = croot + dcroot
    nroot  = nroot + dnroot

    ! if (present(verbose)) then
    !   write(0,*) 'after presumed allocation:'
    !   write(0,*) 'clabl', clabl ! OK
    !   write(0,*) 'nlabl', nlabl ! OK
    !   ! write(0,*) 'dcleaf', dcleaf
    !   ! write(0,*) 'dcroot', dcroot
    !   ! write(0,*) 'dnleaf', dnleaf
    !   ! write(0,*) 'dnroot', dnroot

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

end module _allocation
