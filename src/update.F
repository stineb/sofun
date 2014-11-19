      subroutine update_fpc(pft,jpngr)

      use params_modl, only: kbeer, sla
      use statevars
      use pools, only: pleaf

      implicit none

      ! ARGUMENTS
      integer,intent(in) :: pft
      integer,intent(in) :: jpngr

      ! LOCAL VARIABLES
      real :: fpc_ind
      
      if (crownarea(pft,jpngr)>0.0) then
        lai_ind(pft,jpngr) = (pleaf(pft,jpngr)%c%c12 * sla(pft)) / crownarea(pft,jpngr)
      else
        lai_ind(pft,jpngr) = 0.0
      endif
      
      fpc_ind = ( 1.0 - exp( -1.0 * kbeer * lai_ind(pft,jpngr) ) )
      fpc_grid(pft,jpngr) = crownarea(pft,jpngr) * nind(pft,jpngr) * fpc_ind
      
      end subroutine update_fpc
