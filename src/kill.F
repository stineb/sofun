      subroutine killpft( pft, jpngr )
        !////////////////////////////////////////////////////////////////
        !  Move organic mass of killed PFT to litter. Note that units 
        !  of pleaf and proot is gC/m2/ind., whereas units of plitt is
        !  gC/m2!
        !  b.stocker@imperial.ac.uk
        !----------------------------------------------------------------
        use pools
        use params_modl, only: tree
        use treegeometry, only: nind

        implicit none

        ! ARGUMENTS
        integer, intent(in) :: pft
        integer, intent(in) :: jpngr

        ! sanity check
        if (pleaf(pft,jpngr)%c%c12<0.0) stop 'negative leaf mass'
        if (proot(pft,jpngr)%c%c12<0.0) stop 'negative root mass'

        ! transfer killed mass to litter
        call orgmv( pleaf(pft,jpngr), pleaf(pft,jpngr) &
          , plitt_af(pft,jpngr), scale=nind(pft,jpngr) & 
          )
        call orgmv( proot(pft,jpngr), proot(pft,jpngr) &
              , plitt_bg(pft,jpngr), scale=nind(pft,jpngr) & 
              )
        call orgmv( plabl(pft,jpngr), plabl(pft,jpngr) &
          , plitt_af(pft,jpngr), scale=nind(pft,jpngr) & 
          )  ! xxx or add this to below-ground litter?

        if (tree(pft)) then
          call orgmv( psapw(pft,jpngr), psapw(pft,jpngr) &
            , plitt_af(pft,jpngr), scale=nind(pft,jpngr) & 
            )
          call orgmv( pwood(pft,jpngr), pwood(pft,jpngr) &
            , plitt_as(pft,jpngr), scale=nind(pft,jpngr) & 
            )
        end if

        call initpft(pft,jpngr)

        return
      
      end subroutine killpft
