module md_vegdynamics
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk

  implicit none

  private
  public vegdynamics

  real, parameter :: diam_inc_init = 0.001 ! m

contains

  subroutine vegdynamics( tile, plant, solar, out_pmodel )
    !//////////////////////////////////////////////////////////////////
    ! Updates canopy and tile variables and calls 'estab' to 
    ! simulate establishment of new individuals
    !------------------------------------------------------------------
    use md_classdefs
    use md_params_core, only: npft, nlu, nmonth, ndayyear
    use md_plant, only: initpft, get_leaftraits, plant_type, params_pft_plant
    use md_allocation, only: update_tree
    use md_tile, only: tile_type
    use md_waterbal, only: solartype
    use md_gpp, only: outtype_pmodel

    ! arguments
    type( tile_type ), dimension(nlu), intent(inout)           :: tile
    type( plant_type ), dimension(npft), intent(in)            :: plant ! npft counts over PFTs in all land units (tiles)
    type( solartype ), intent(in)                              :: solar
    type( outtype_pmodel ), dimension(npft,nmonth), intent(in) :: out_pmodel

    ! local variables
    integer :: pft, lu

    do lu=1,nlu
      if (tile(lu)%fpc_grid == 0.0) then
        !------------------------------------------------------------------
        ! Add individuals
        !------------------------------------------------------------------
        do pft=1,npft
          if (params_pft_plant(pft)%lu_category==lu) then

            ! initialise all pools of this PFT with zero
            call initpft( plant(pft) )

            ! get annually updated leaf traits (vary because of variations in light and CO2)
            call get_leaftraits( plant(pft), solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )

            ! add a "seed" by forcing initial diameter increment
            call update_tree( plant(pft), diam_inc_init )

            ! give it some labile C and N to pay for turnover this year (until allocation at the end of year)
            plant(pft)%plabl = orgplus( &
              orgfrac( params_pft_plant(pft)%k_decay_leaf_base * ndayyear, plant(pft)%pleaf ), &
              orgfrac( params_pft_plant(pft)%k_decay_root * ndayyear, plant(pft)%proot ) &
              ) 

            ! xxx needs to be done: add implicit C and N fluxes to NPP and N-uptake

            ! Set number of individuals (XXX change this to m-2)
            tile(lu)%nind(pft)  = 1.0 
            tile%fpc_grid = tile(lu)%nind(pft) * plant(pft)%acrown

          end if
        end do
      end if
    end do

  end subroutine vegdynamics

end module md_vegdynamics
