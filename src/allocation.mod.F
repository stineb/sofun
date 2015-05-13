module _allocation

contains

  subroutine allocation_daily( jpngr, doy )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use _params_core, only: npft
    use _params_modl, only: r_cton_root, growtheff, grass
    use _pools, only: pleaf, proot, plabl, pninorg, r_cton_leaf, r_ntoc_leaf, sla
    use _gpp, only: mlue, mrd_unitiabs
    use _waterbal, only: dppfd, meanmppfd
    use _findroot_fzeroin

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    real    :: dcleaf
    real    :: dnleaf
    logical :: findroot          ! true if allocation to leaves (roots) is not 100% and not 0%
    real    :: max_dcleaf_n_constraint
    real    :: max_dcroot_n_constraint
    real    :: max_dcleaf
    real    :: eval_allroots
    real    :: abserr
    real    :: relerr
    type(outtype_zeroin) :: out_zeroin


    abserr=100.0*XMACHEPS
    relerr=1000.0*XMACHEPS

    do pft=1,npft

      if (grass(pft)) then

        ! Maximum is the lower of all labile C and the C to be matched by all labile N,
        ! discounted by the yield factor.
        max_dcleaf_n_constraint = plabl(pft,jpngr)%n%n14 * r_cton_leaf(pft,jpngr)
        max_dcroot_n_constraint = plabl(pft,jpngr)%n%n14 * r_cton_root(pft)
        max_dcleaf = min( growtheff * plabl(pft,jpngr)%c%c12, max_dcleaf_n_constraint )

        ! Optimisation by balanced growth
        ! Test I: Evaluate balance if all is put to roots.
        ! If C:N ratio of return is still greater than whole-plant C:N ratio, then put all to roots.
        findroot = .true.
        eval_allroots  = eval_imbalance( 0.0, pleaf(pft,jpngr), proot(pft,jpngr), plabl(pft,jpngr), r_ntoc_leaf(pft,jpngr), sla(pft,jpngr), mlue(moy), dppfd(doy), mrd_unitiabs(moy), meanmppfd(moy), pninorg(lu,jpngr), pft )
        if (eval_allroots > 0.0) then
          dcleaf = 0.0
          findroot = .false.
        end if

        ! Test II: Evaluate balance if all is put to leaves.
        ! If C:N ratio of return is still lower than whole-plant C:N ratio, then put all to leaves.
        eval_allleaves = eval_imbalance( max_dcleaf, pleaf(pft,jpngr), proot(pft,jpngr), plabl(pft,jpngr), r_ntoc_leaf(pft,jpngr), sla(pft,jpngr), mlue(moy), dppfd(doy), mrd_unitiabs(moy), meanmppfd(moy), pninorg(lu,jpngr), pft )
        if (eval_allleaves < 0.0) then
          dcleaf = max_dcleaf
          findroot = .false.
        end if

        if (findroot) then
          ! Find root of function 'eval_imbalance()' in the interval [0.0, max_dcleaf]
          print("finding root")
          out_zeroin = zeroin( eval_imbalance, abserr, relerr, 100, 0.0 ,max_dcleaf )
          if ( out_zeroin%error /= 0 ) then
            dcleaf = 0.0
          else
            dcleaf = out_zeroin%root
          end if
        end if

        ! Allocate based on 'dcleaf' determined above
        dnleaf = dcleaf * r_ntoc_leaf(pft,jpngr)

        plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - 1.0 / growtheff * dcleaf
        plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnleaf
        pleaf(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12 + dcleaf
        pleaf(pft,jpngr)%n%n14 = pleaf(pft,jpngr)%n%n14 + dnleaf

        dcroot = min( growtheff * clabl, r_cton_root(pft) * nlabl )
        dnroot = dcroot * r_ntoc_root

        plabl(pft,jpngr)%c%c12 = plabl(pft,jpngr)%c%c12 - 1.0 / growtheff * dcroot
        plabl(pft,jpngr)%n%n14 = plabl(pft,jpngr)%n%n14 - dnroot
        pleaf(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12 + dcroot
        pleaf(pft,jpngr)%n%n14 = pleaf(pft,jpngr)%n%n14 + dnroot

      else

        stop 'allocation_daily not implemented for trees'

      end if

    end do

  end subroutine allocation_daily


  function eval_imbalance( &
      dcleaf &
    , pleaf &
    , proot &
    , plabl &
    , r_ntoc_leaf &
    , sla &
    , mlue &
    , dppfd &
    , mrd_unitiabs &
    , meanmppfd &
    , ninorg &
    , pft &
    ) result ( eval )
    !/////////////////////////////////////////////////////////
    ! Evaluates C:N ratio of new assimilation after allocation 
    ! versus whole-plant C:N ratio after allocation. Optimal 
    ! allocation is where the two are equal. 
    !---------------------------------------------------------
    use _classdefs, only: orgpool, nitrogen
    use _params_modl, only: growtheff, k_decay_leaf, k_decay_root, r_cton_root
    use _gpp, only: calc_dgpp, calc_drd
    use _nuptake, only: calc_dnup, exurate
    use _npp, only: calc_resp_maint, r_root

    ! arguments
    real, intent(in)          :: dcleaf
    type(orgpool), intent(in) :: pleaf
    type(orgpool), intent(in) :: proot
    type(orgpool), intent(in) :: plabl
    real, intent(in)          :: r_ntoc_leaf
    real, intent(in)          :: sla
    real, intent(in)          :: mlue
    real, intent(in)          :: dppfd
    real, intent(in)          :: mrd_unitiabs
    real, intent(in)          :: meanmppfd
    type(nitrogen), intent(in):: ninorg
    integer                   :: pft

    ! function return variable
    real, intent(out) :: eval

    ! local variables
    real :: dnleaf
    real :: cleaf, nleaf, croot, nroot
    real :: lai
    real :: gpp
    real :: rd
    real :: mresp_root
    real :: exu
    real :: dc
    real :: dn

    ! shorter ...
    cleaf = pleaf%c%c12
    nleaf = pleaf%n%n14
    croot = proot%c%c12
    nroot = proot%n%n14
    clabl = plabl%c%c12
    nlabl = plabl%n%n14

    ! Allocate
    dnleaf = dcleaf * r_ntoc_leaf
    clabl  = clabl - 1.0 / growtheff * dcleaf
    nlabl  = nlabl - dnleaf
    cleaf  = cleaf + dcleaf
    nleaf  = nleaf + dnleaf
    
    dcroot = min( growtheff * clabl, r_cton_root(pft) * nlabl )
    dnroot = dcroot * r_ntoc_root
    clabl  = clabl - 1.0 / growtheff * dcroot
    nlabl  = nlabl - dnroot
    croot  = croot + dcroot
    nroot  = nroot + dnroot

    ! Allocation and decay
    cleaf = cleaf * (1.0 - k_decay_leaf(pft))
    croot = croot * (1.0 - k_decay_root(pft))
    nleaf = nleaf * (1.0 - k_decay_leaf(pft))
    nroot = nroot * (1.0 - k_decay_root(pft))

    ! Calculate next day's C and N return after assumed allocation (tissue turnover happens before!)
    lai = cleaf * sla
    gpp = calc_dgpp( lai, mlue, dppfd )
    rd  = calc_drd( lai, mrd_unitiabs, meanmppfd )
    exu = croot * exurate
    mresp_root = calc_resp_maint( croot, r_root)

    dc  = gpp - rd - mresp_root - exu
    dn  = calc_dnup( croot, ninorg%n14 )

    ! Evaluation quantity is the difference between the 
    ! C:N ratio of new assimilates and the C:N ratio 
    ! of the whole plant after allocation.
    eval = growtheff * (dc + clabl) / (dn + nlabl) - ( cleaf + croot ) / ( nleaf + nroot )

 end function eval_imbalance

end module _allocation
