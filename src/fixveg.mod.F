module fixveg
  !//////////////////////////////////////////////////////////////////////////////
  ! MODULE FIXVEG
  ! Contains all subroutines and functions to calculate vegetation cover with 
  ! fixed prescribed PFT coverate (fractional plant cover, FPC). These are:
  ! - establishment
  ! - light
  ! - mortality
  ! Variables and functions contained in module 'treegeometry' are not inter-
  ! changeable. This module contains subroutines that update/calculate variables
  ! contained in 'treegeometry'.
  ! This module requires optional argument 'fpc_grid_data' to be passed on to 
  ! 'biosphere'. 
  !-----------------------------------------------------------------------------
  use treegeometry

contains

  subroutine establishment( jpngr )
    !//////////////////////////////////////////////////////////////////////////////
    ! GET FIXED PRESCRIBED VEGETATION COVER
    ! Subroutine defines 'fpc_grid' and updates 'nind' following geo-
    ! metry constraints.
    !-----------------------------------------------------------------------------
    use classdefs
    use params_modl, only: tree
    use params_core, only: npft
    use pools, only: psapw, pleaf

    ! xxx debug
    use pools, only: plabl

    implicit none

    ! ARGUMENTS
    integer, intent(in)      :: jpngr

    ! LOCAL VARIABLES
    integer                  :: pft
    logical, dimension(npft) :: estab
    type(orgpool)            :: psapw_temp


    !-----------------------------------------------------------------------------
    ! derive which PFTs are present from fpc_grid (which is prescribed)
    !-----------------------------------------------------------------------------
    do pft=1,npft
      if (fpc_grid(pft,jpngr)>0.0) then
        ispresent(pft,jpngr) = .true.
      else
        ispresent(pft,jpngr) = .false.
      end if
    end do

    ! xxx debug
    !print*,'A fpc_grid  ',fpc_grid(:,1)
    !print*,'A crownarea ',crownarea(:,1)
    !print*,'A leaf mass ',pleaf(:,1)
    !print*,'A labl      ',plabl(:,1)

    do pft=1,npft

      ! xxx does this really also work for trees?

      !-----------------------------------------------------------------------------
      ! Re-calculate crownarea
      !-----------------------------------------------------------------------------
      psapw_temp = psapw(pft,jpngr)
      call update_tree_geometry( &
        pft, jpngr, psapw_temp, ntoc( psapw(pft,jpngr), default=0.0 ) &
        )

      !-----------------------------------------------------------------------------
      ! Update fpc_grid, fpc_ind, and lai_ind
      ! xxx shouldn't fpc_grid be prescribed?
      !-----------------------------------------------------------------------------
      lai_ind(pft,jpngr)  = get_lai_ind( pleaf(pft,jpngr)%c%c12, crownarea(pft,jpngr), pft )

      !-----------------------------------------------------------------------------
      ! Re-calculate 'nind' based on new leaf mass and crownarea
      !-----------------------------------------------------------------------------
      nind(pft,jpngr) = get_nind( tree(pft), fpc_grid(pft,jpngr), crownarea(pft,jpngr), fpc_ind(pft,jpngr) )

    end do

    ! xxx debug
    !print*,'B fpc_grid ',fpc_grid(:,1)
    !print*,'B crownarea ',crownarea(:,1)

    return

  end subroutine establishment


  subroutine light( jpngr )
    !////////////////////////////////////////////////////////////////
    ! LIGHT COMPETITION
    ! reduce FPC of trees to 'fpc_tree_max', killing according number of
    ! individuals (=reducing 'nind')
    ! only effective with dynamic vegetation
    ! b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use classdefs
    use params_core, only: nlu
    use params_modl, only: sla, islu, tree, grass, fpc_tree_max, kbeer
    use treegeometry, only: fpc_grid, ispresent, nind, update_fpc
    use pools
    implicit none

    ! ARGUMENTS
    integer, intent(in) :: jpngr


    ! LOCAL VARIABLES
    integer :: pft,lu

    real :: fpc_tree_total  ! total grid FPC for tree PFTs
    real :: fpc_grass_total ! total grid FPC for grass PFTs
    real :: fpc_total       ! total grid FPC for all PFTs
    real :: excess          ! tree FPC or grass cover to be reduced
    real :: nind_kill       ! reduction in individual density to reduce tree FPC to permitted maximum (indiv/m2)
    real :: cleaf_target    ! required grass leaf mass so that FPC constraints are satisfied
    real :: frac_kill       ! fraction of grass PFT pools removed so that FPC constraints are satisfied


    !do lu=1,nlu

 !  !   do pft=1,npft
   !     if (islu(pft,lu)) then
   !       if (ispresent(pft,jpngr)) then
   !         if (tree(pft)) then
   !           !-----------------------------------------------------------------------------
   !           ! TREES
   !           !-----------------------------------------------------------------------------
   !           ! xxx find a solution for this
   !           !if (fpc_tree_total> xxx prescribed fpc_grid of this pft xxx) then
   !           if (fpc_tree_total>fpc_tree_max) then
   !             ! test kuno 2005-01-7: split excess according to total FPC rather than inc
   !             ! note that fpc_inc is already included in fpc_grid, fpc_tree_total, since
   !             ! it is used to calculate the excess fpc (check origin of these vars).
   !             ! kuno 2005-09-23:  mixed scheme is abolished, use only area-based scheme
   !             excess = (fpc_tree_total-fpc_tree_max)*fpc_grid(pft,jpngr)/fpc_tree_total
   !                             
   !             ! Reduce individual density (and thereby gridcell-level biomass)
   !             ! so that total tree FPC reduced to 'fpc_tree_max'
   !             nind_kill = nind(pft,jpngr) * (excess/fpc_grid(pft,jpngr))
   !             nind(pft,jpngr) = nind(pft,jpngr) - nind_kill

 !  !             ! Transfer lost biomass to litter
   !             call orgcp( orgfrac( nind_kill, pleaf(pft,jpngr) ), plitt_af(pft,jpngr) )
   !             call orgcp( orgfrac( nind_kill, psapw(pft,jpngr) ), plitt_af(pft,jpngr) )
   !             call orgcp( orgfrac( nind_kill, pwood(pft,jpngr) ), plitt_as(pft,jpngr) )
   !             call orgcp( orgfrac( nind_kill, proot(pft,jpngr) ), plitt_bg(pft,jpngr) )

 !  !           end if               ! fpc_tree_total > fpc_tree_max

 !  !         else
   !           !-----------------------------------------------------------------------------
   !           ! GRASSES
   !           !-----------------------------------------------------------------------------
   !           ! If grasses only, fpc_grass_total ~= 1; fpc_tree_total = 0
   !           if ( fpc_grass_total > ( 1.0 - min( fpc_tree_max, fpc_tree_total ) ) ) then           

 !  !             excess = ( fpc_grass_total - 1.0 + min( fpc_tree_max, fpc_tree_total ) ) * ( fpc_grid(pft,jpngr) / fpc_grass_total ) + 1.0d-15 

 !  !             cleaf_target = -1.0 / kbeer * dlog( 1.0d0 - ( fpc_grid(pft,jpngr) - excess ) ) / sla(pft)
   !             frac_kill = (pleaf(pft,jpngr)%c%c12 - cleaf_target) / pleaf(pft,jpngr)%c%c12

 !  !             ! Transfer lost biomass to litter
   !             call orgmv( orgfrac( frac_kill, pleaf(pft,jpngr)), pleaf(pft,jpngr), plitt_af(pft,jpngr) )
   !             call orgmv( orgfrac( frac_kill, proot(pft,jpngr)), proot(pft,jpngr), plitt_bg(pft,jpngr) )

 !  !           end if
   !         
   !         end if ! tree/grass
   !    
   !         !-----------------------------------------------------------------------------
   !         ! Update fractional plant coverage. Should now comply with fpc_tree_max and 
   !         ! not exceed 1.
   !         !-----------------------------------------------------------------------------
   !         call update_fpc(pft,jpngr)

 !  !       end if  ! ispresent

 !  !     end if  ! islu

 !  !   end do  !pft

 !  ! end do

    return

  end subroutine light

end module fixveg