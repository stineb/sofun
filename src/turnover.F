subroutine turnover( jpngr )
!////////////////////////////////////////////////////////////////
!  Annual vegetation biomass turnover, called at the end of the
!  year.
!----------------------------------------------------------------
use _classdefs
use _params_modl, only: lu_category, tree, long_leaf, long_root, long_sapw, npft, F_NRETAIN
use _pools, only: pleaf, proot, psapw, pwood, plitt_af, plitt_bg, plabl
use _vegdynamics , only: nind, ispresent
use _outvars, only: outaCveg2lit, outaNveg2lit

implicit none

! ARGUMENTS
integer, intent(in) :: jpngr

! LOCAL VARIABLES:
integer :: pft
integer :: lu
real :: dleaf
real :: dsapw
real :: droot
type(orgpool) :: lm_turn
type(orgpool) :: sm_turn
type(orgpool) :: rm_turn

do pft=1,npft

  if (ispresent(pft,jpngr)) then

    lu=lu_category(pft)

    ! Turnover _rates are reciprocals of tissue longevity
    dleaf=1.0/long_leaf(pft)
    dsapw=1.0/long_sapw(pft)
    droot=1.0/long_root(pft)

    ! xxx debug
    !print*,'TURNOVER: dleaf set to droot'  
    !dleaf = droot  
    !long_leaf(pft) = 0.5
    !long_root(pft) = 0.5
    !long_sapw(pft) = 0.5

    !print*,'TURNOVER A: C:N of whole plant ',cton( orgplus(pleaf(pft,jpngr),proot(pft,jpngr)), default=0.0 )
    !print*,'TURNOVER A: C:N of leaf mass   ',cton( pleaf(pft,jpngr), default=0.0 )
    !print*,'TURNOVER A: C:N of root mass   ',cton( proot(pft,jpngr), default=0.0 )
      
    ! Calculate the biomass turnover in this year          
    lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) )
    rm_turn = orgfrac( droot, proot(pft,jpngr) )
    if (tree(pft)) then
      sm_turn = orgfrac( dsapw, psapw(pft,jpngr) )
    endif

#if _allocation_makela
    ! C: plant -> litter
    call cmvRec( lm_turn%c, pleaf(pft,jpngr)%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    call cmvRec( rm_turn%c, proot(pft,jpngr)%c, plitt_bg(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    if (tree(pft)) then
      call cmvRec( sm_turn%c, psapw(pft,jpngr)%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    endif
 
    ! N: plant -> storage pool
    call nmv( lm_turn%n, pleaf(pft,jpngr)%n, plabl(pft,jpngr)%n, scale=F_NRETAIN )
    call nmv( rm_turn%n, proot(pft,jpngr)%n, plabl(pft,jpngr)%n, scale=F_NRETAIN )
    if (tree(pft)) then
      call nmv( sm_turn%n, psapw(pft,jpngr)%n, plabl(pft,jpngr)%n, scale=F_NRETAIN )
    endif

    ! N: plant -> litter
    call nmvRec( lm_turn%n, pleaf(pft,jpngr)%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr)*(1.0-F_NRETAIN) )
    call nmvRec( rm_turn%n, proot(pft,jpngr)%n, plitt_bg(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr)*(1.0-F_NRETAIN) )
    if (tree(pft)) then
      call nmvRec( sm_turn%n, psapw(pft,jpngr)%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr)*(1.0-F_NRETAIN) )
    endif

#else
    ! xxx try
    F_NRETAIN = 0.0

    ! reduce leaf mass and root mass
    call orgsub( lm_turn, pleaf(pft,jpngr) )
    call orgsub( rm_turn, proot(pft,jpngr) )

    ! add all C to litter
    call cmvRec( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))
    call cmvRec( rm_turn%c, rm_turn%c, plitt_bg(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))

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
#endif

    !print*,'TURNOVER B: C:N of whole plant ',cton( orgplus(pleaf(pft,jpngr),proot(pft,jpngr)), default=0.0 )
    !print*,'TURNOVER B: C:N of leaf mass   ',cton( pleaf(pft,jpngr), default=0.0 )
    !print*,'TURNOVER B: C:N of root mass   ',cton( proot(pft,jpngr), default=0.0 )

  endif                   !present
enddo                     !pft

end subroutine turnover
