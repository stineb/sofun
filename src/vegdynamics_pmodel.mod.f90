module md_vegdynamics
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk

  implicit none

  private
  public vegdynamics

contains

  subroutine vegdynamics( tile, plant, solar, out_pmodel, fapar_prescr )
    !//////////////////////////////////////////////////////////////////
    ! Updates canopy and tile variables and calls 'estab' to 
    ! simulate establishment of new individuals
    !------------------------------------------------------------------
    use md_params_core, only: npft, nlu, nmonth, dummy
    use md_plant, only: initpft, get_leaftraits, plant_type, params_pft_plant
    use md_tile, only: tile_type
    use md_waterbal, only: solartype
    use md_gpp, only: outtype_pmodel

    ! arguments
    type( tile_type ), dimension(nlu), intent(inout)           :: tile
    type( plant_type ), dimension(npft), intent(inout)         :: plant ! npft counts over PFTs in all land units (tiles)
    type( solartype ), intent(in)                              :: solar
    type( outtype_pmodel ), dimension(npft,nmonth), intent(in) :: out_pmodel

    ! arguments (may be dummy)
    real, intent(in) :: fapar_prescr

    ! local variables
    integer :: pft, lu

    do lu=1,nlu
      !------------------------------------------------------------------
      ! Add individuals
      !------------------------------------------------------------------
      do pft=1,npft
        if (params_pft_plant(pft)%lu_category==lu) then

          ! Override interactively simulated fAPAR with data
          if (fapar_prescr/=dummy) plant(pft)%fapar_ind = fapar_prescr

          ! initialise all pools of this PFT with zero
          call initpft( plant(pft) )

          ! Set stuff that is required but obsolete for P-model
          plant(pft)%acrown        = 1.0
          tile(lu)%canopy%fpc_grid = 1.0

          ! get annually updated leaf traits (vary because of variations in light and CO2)
          call get_leaftraits( plant(pft), solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )

        end if
      end do
    end do

  end subroutine vegdynamics


end module md_vegdynamics