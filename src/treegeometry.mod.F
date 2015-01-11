      module treegeometry
        !////////////////////////////////////////////////////////////////
        ! Module contains tree geometry variables and functions.
        ! Variables have 'jpngr' dimension, thus have "memory" from one
        ! year to the next like pool variables.
        !----------------------------------------------------------------
        use params_core, only: npft, maxgrid
        implicit none

        logical, dimension(npft,maxgrid) :: ispresent      ! boolean whether PFT is present
        real, dimension(npft,maxgrid)    :: height      ! tree height (m)
        real, dimension(npft,maxgrid)    :: crownarea   ! individual's tree crown area
        real, dimension(npft,maxgrid)    :: nind        ! number of individuals [1/m2]
        real, dimension(npft,maxgrid)    :: lai_ind
        real, dimension(npft,maxgrid)    :: fpc_grid
        real, dimension(npft,maxgrid)    :: fpc_ind 

      contains

        subroutine update_fpc( pft, jpngr )
          !////////////////////////////////////////////////////////////////
          ! Updates fpc_grid, fpc_ind, and lai_ind
          !----------------------------------------------------------------
          use pools, only: pleaf
          
          implicit none

          ! ARGUMENTS
          integer, intent(in) :: pft, jpngr

          lai_ind(pft,jpngr)  = get_lai_ind( pleaf(pft,jpngr)%c%c12, crownarea(pft,jpngr), pft )   
          fpc_ind(pft,jpngr)  = get_fpc_ind( lai_ind(pft,jpngr) )
          fpc_grid(pft,jpngr) = get_fpc_grid( crownarea(pft,jpngr), nind(pft,jpngr), fpc_ind(pft,jpngr) )

        end subroutine update_fpc


        function get_nind( fpc_grid, crownarea, fpc_ind )
          !////////////////////////////////////////////////////////////////
          !----------------------------------------------------------------
          use pools, only: pleaf

          implicit none

          ! ARGUMENTS
          real,intent(in) :: fpc_grid
          real,intent(in) :: crownarea
          real,intent(in) :: fpc_ind

          ! FUNCTION RETURN VARIABLE
          real, intent(out) :: get_nind

          get_nind = fpc_grid / ( crownarea * fpc_ind ) 

        end function get_nind


        function get_fpc_grid( crownarea, nind, fpc_ind )
          !////////////////////////////////////////////////////////////////
          ! Function returns total fractional plant cover of a PFT
          ! Eq. 8 in Sitch et al., 2003
          !----------------------------------------------------------------

          ! arguments
          real, intent(in) :: crownarea
          real, intent(in) :: nind
          real, intent(in) :: fpc_ind

          ! function return variable
          real, intent(out) :: get_fpc_grid

          get_fpc_grid = crownarea * nind * fpc_ind

        end function get_fpc_grid


        function get_fpc_ind( lai_ind )
          !////////////////////////////////////////////////////////////////
          ! Function returns fractional plant cover an individual
          ! Eq. 7 in Sitch et al., 2003
          !----------------------------------------------------------------
          use params_modl, only: kbeer

          ! arguments
          real, intent(in) :: lai_ind

          ! function return variable
          real, intent(out) :: get_fpc_ind

          get_fpc_ind = ( 1.0 - exp( -1.0 * kbeer * lai_ind) )

        end function get_fpc_ind


        function get_lai_ind( pleaf, crownarea, pft )
          !////////////////////////////////////////////////////////////////
          ! Function returns leaf area index of an individual
          ! Eq. 5 in Sitch et al., 2003
          !----------------------------------------------------------------
          use params_modl, only: sla
          implicit none

          ! arguments
          real, intent(in)    :: pleaf
          real, intent(in)    :: crownarea
          integer, intent(in) :: pft

          ! function return variable
          real :: get_lai_ind

#if _check_sanity
          if (crownarea<0.0) stop 'GET_LAI_IND: crownarea < 0'
#endif

          if (crownarea==0.0) then
            get_lai_ind = 0.0
          else
            get_lai_ind = pleaf * sla(pft) / crownarea
          endif

        end function get_lai_ind


        subroutine update_tree_geometry( pft, jpngr, sm_ind_temp, afntoc_sm )
          !////////////////////////////////////////////////////////////////
          ! Subroutine updating tree geometry.
          ! Requires global variables:
          !   - psapw, pwood
          ! Requires global model parameters:
          !   - wooddens, allom, crownarea_max, reinickerp, latosa, sla
          ! Updates global variables (geometry-related):
          !   - height, crownarea, psapw
          ! Summary of calculation:
          ! ( Cwood+Csap', Cleaf, Croot ) --> ( LA, SA, H, D, CA, Csap)
          ! Allometry is independent of sapwood/heartwood ratio,
          ! only sum of the two is used (VI)
          !----------------------------------------------------------------
          use classdefs
          use pools, only: pleaf, psapw, pwood
          use params_modl, only: wooddens, allom1, allom2, allom3, &
            crownarea_max, reinickerp, latosa, sla, pi

          ! arguments
          integer, intent(in)          :: pft, jpngr
          type(orgpool), intent(in)    :: sm_ind_temp
          real, intent(in)             :: afntoc_sm

          ! local variables
          real    :: stemdiam

          ! (I)   LA = latosa * SA
          ! (II)  Cleaf = lmtorm * Croot
          ! (III) H = allom2 * D**allom3
          ! (IV)  CA = min (allom1 * D**reinickerp, crownarea_max)
          ! (V)   LA = Cleaf * sla
          ! (VI)  wooddens = ( Cwood + Csap' ) / ( (D/2)**2 * pi * H )

          ! 1. solve (VI) for D using (III)
          stemdiam = (4.0*(sm_ind_temp%c%c12+pwood(pft,jpngr)%c%c12)/ &
            wooddens(pft)/pi/allom2(pft))**(1.0/(2.0+allom3))

          ! 2. solve (III) for H
          height(pft,jpngr) = allom2(pft)*stemdiam**allom3

          ! 3. (IV) for CA
          crownarea(pft,jpngr) = min(crownarea_max, allom1(pft)*stemdiam**reinickerp)

          ! 4. Recalculate sapwood mass, transferring excess sapwood to heartwood
          ! compartment, if necessary to satisfy (I)
          psapw(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12*height(pft,jpngr)*wooddens(pft)*sla(pft)/latosa(pft)
          call orgcp( orgminus( sm_ind_temp, psapw(pft,jpngr)), pwood(pft,jpngr) )


        end subroutine update_tree_geometry

      end module treegeometry

