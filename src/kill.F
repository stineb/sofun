subroutine killpft( pft, jpngr )
  !////////////////////////////////////////////////////////////////
  !  Move organic mass of killed PFT to litter. Note that units 
  !  of pleaf and proot is gC/m2/ind., whereas units of plitt is
  !  gC/m2!
  !  b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _vars_core, only: pleaf, plitt_af, plitt_bg, psapw, plitt_as, pwood, proot, plabl
  use _params_modl, only: tree
  use _vegdynamics, only: nind
  use _outvars, only: outaCveg2lit, outaNveg2lit

  implicit none

  ! arguments
  integer, intent(in) :: pft
  integer, intent(in) :: jpngr

  ! sanity check
  if (pleaf(pft,jpngr)%c%c12<0.0) stop 'negative leaf mass'
  if (proot(pft,jpngr)%c%c12<0.0) stop 'negative root mass'

  ! transfer killed mass to litter
  call orgmvRec( pleaf(pft,jpngr), pleaf(pft,jpngr), plitt_af(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr), nind(pft,jpngr) )
  call orgmvRec( proot(pft,jpngr), proot(pft,jpngr), plitt_bg(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr), nind(pft,jpngr) )
  call orgmvRec( plabl(pft,jpngr), plabl(pft,jpngr), plitt_af(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr), nind(pft,jpngr) )  ! xxx or add this to below-ground litter?

  if (tree(pft)) then
    call orgmvRec( psapw(pft,jpngr), psapw(pft,jpngr) , plitt_af(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr), nind(pft,jpngr) )
    call orgmvRec( pwood(pft,jpngr), pwood(pft,jpngr) , plitt_as(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr), nind(pft,jpngr) )
  end if

  call initpft(pft,jpngr)

end subroutine killpft
