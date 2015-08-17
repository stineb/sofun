module _turnover
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk

  implicit none

contains

  subroutine turnover( jpngr, doy )
    !////////////////////////////////////////////////////////////////
    !  Annual vegetation biomass turnover, called at the end of the
    !  year.
    !----------------------------------------------------------------
    use _classdefs
    use _params_modl, only: lu_category, tree, grass, &
      k_decay_leaf, k_decay_sapw, k_decay_root, npft, F_NRETAIN
    use _phenology, only: summergreen, shedleaves
    use _vars_core, only: pleaf, proot, psapw, pwood, plitt_af, &
      plitt_bg, plabl, pexud, pninorg
    use _vars_core , only: nind, ispresent
    use _outvars, only: outaCveg2lit, outaNveg2lit

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu
    real :: dleaf
    real :: dsapw
    real :: droot
    real :: dlabl
    type(orgpool) :: lm_turn
    type(orgpool) :: sm_turn
    type(orgpool) :: rm_turn
    type(orgpool) :: lb_turn

    do pft=1,npft

      if (ispresent(pft,jpngr)) then

        lu = lu_category(pft)

        ! Turnover _rates are reciprocals of tissue longevity
        ! dleaf=1.0/long_leaf(pft)
        ! assuming no continuous leaf turnover
        if (grass(pft)) then
          if (summergreen(pft)) then
            if (shedleaves(doy,pft)) then
              write(0,*) 'end of season: labile C, labile N', plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14
              dleaf = 1.0
              droot = 1.0
              dlabl = 1.0
              ! write(0,*) 'complete die-off on day', doy
            else
              dleaf = k_decay_leaf(pft)  
              droot = k_decay_root(pft)
              dlabl = 0.0
            end if
          else
            dleaf = k_decay_leaf(pft)
            droot = k_decay_root(pft)
            dlabl = 0.0
          end if   
        else
          stop 'turnover only implemented for grasses'        
        end if

        ! Calculate the biomass turnover in this year          
        lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) )
        rm_turn = orgfrac( droot, proot(pft,jpngr) )
        if (tree(pft)) then
          sm_turn = orgfrac( dsapw, psapw(pft,jpngr) )
        endif
        lb_turn = orgfrac( dlabl, plabl(pft,jpngr) )

        ! xxx try
        F_NRETAIN = 0.0

        ! reduce leaf mass and root mass
        call orgsub( lm_turn, pleaf(pft,jpngr) )
        call orgsub( rm_turn, proot(pft,jpngr) )
        if (dlabl>0.0) then
          call orgsub( lb_turn, plabl(pft,jpngr) )
        end if

        ! add all organic (fixed) C to litter
        call cmvRec( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))
        call cmvRec( rm_turn%c, rm_turn%c, plitt_bg(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))

        ! add labile C and N to litter as well
        if (dlabl>0.0) then
          write(0,*) 'TURNOVER: WARNING LABILE C AND N ARE SET TO ZERO WITHOUT MASS CONSERVATION'
          lb_turn%c%c12 = 0.0
          lb_turn%n%n14 = 0.0
          ! call cmvRec( lb_turn%c, lb_turn%c, plitt_bg(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
          ! call nmvRec( lb_turn%n, lb_turn%n, plitt_bg(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
          ! write(0,*) 'end of growing season-plabl:', plabl
          ! write(0,*) 'moving C and N to exudates and ninorg:', lb_turn
        end if

        ! retain fraction of N
        call nmv( nfrac( F_NRETAIN, lm_turn%n ), lm_turn%n, plabl(pft,jpngr)%n )
        call nmv( nfrac( F_NRETAIN, rm_turn%n ), rm_turn%n, plabl(pft,jpngr)%n )

        ! rest goes to litter
        call nmvRec( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
        call nmvRec( rm_turn%n, rm_turn%n, plitt_bg(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )

        if (tree(pft)) then
          ! Sapwood turns over into heartwood (non-living wood)
          call orgmv( sm_turn, psapw(pft,jpngr), pwood(pft,jpngr) )
        endif
        
        ! done here: Set isotopes to MISSING_VALUE if _pools are empty

        !print*,'TURNOVER B: C:N of whole plant ',cton( orgplus(pleaf(pft,jpngr),proot(pft,jpngr)), default=0.0 )
        !print*,'TURNOVER B: C:N of leaf mass   ',cton( pleaf(pft,jpngr), default=0.0 )
        !print*,'TURNOVER B: C:N of root mass   ',cton( proot(pft,jpngr), default=0.0 )

      endif                   !present
    enddo                     !pft

  end subroutine turnover

end module _turnover
