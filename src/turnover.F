      subroutine turnover( jpngr )
      !////////////////////////////////////////////////////////////////
      !  Annual vegetation biomass turnover, called at the end of the
      !  year.
      !----------------------------------------------------------------
      use classdefs
      use params_modl, only: lu_category, long_leaf, long_root, long_sapw, npft, F_NRETAIN
      use pools, only: pleaf, proot, plitt_af, plitt_bg, plabl
      use statevars, only: nind, ispresent

      implicit none

      ! ARGUMENTS
      integer, intent(in) :: jpngr
      
      ! LOCAL VARIABLES:
      integer :: pft
      integer :: lu
      real :: dleaf
      !real :: dsapw
      real :: droot
      type(orgpool) :: lm_turn
      type(orgpool) :: sm_turn
      type(orgpool) :: rm_turn

      do pft=1,npft

        if (ispresent(pft,jpngr)) then

          lu=lu_category(pft)

          ! Turnover rates are reciprocals of tissue longevity
          dleaf=1.0/long_leaf(pft)
          !dsapw=1.0/long_sapw(pft)
          droot=1.0/long_root(pft)

          ! Calculate the biomass turnover in this year          
          lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) )
          !lm_turn = orgfrac( dsapw, psapw(pft,jpngr) )
          rm_turn = orgfrac( droot, proot(pft,jpngr) )
          
          ! xxx try
          F_NRETAIN = 0.0

          ! Move turnover to litter
          ! All C from leaf turnover goes to litter
          call cmv( lm_turn%c, pleaf(pft,jpngr)%c, plitt_af(pft,jpngr)%c, scale=nind(pft,jpngr) )
          ! A fraction (F_NRETAIN) of leaf N turnover is retained and moved to temporary storage (plabl) ...
          call nmv( lm_turn%n, pleaf(pft,jpngr)%n, plabl(pft,jpngr)%n,    scale=F_NRETAIN )
          ! ... rest goes to litter
          call nmv( lm_turn%n, pleaf(pft,jpngr)%n, plitt_af(pft,jpngr)%n, scale=nind(pft,jpngr)*(1.-F_NRETAIN) )

          ! No root N is retained (xxx is this viable? xxx)
          call orgmv( rm_turn, proot(pft,jpngr), plitt_bg(pft,jpngr), scale=nind(pft,jpngr) )

          ! done here: Set isotopes to MISSING_VALUE if pools are empty
          
        endif                   !present
      enddo                     !pft

      end subroutine turnover
