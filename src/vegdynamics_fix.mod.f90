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
    use md_params_core, only: npft, nlu
    use md_plant, only: initpft, get_leaftraits
    use md_allocation, only: update_tree

    ! arguments
    type( tile_type ), dimension(nlu), intent(inout)           :: tile
    type( plant_type ), dimension(npft), intent(in)            :: plant ! npft counts over PFTs in all land units (tiles)
    type( solartype ), intent(in)                              :: solar
    type( outtype_pmodel ), dimension(npft,nmonth), intent(in) :: out_pmodel

    ! local variables
    integer :: pft

    do lu=1,nlu
      if (tile(lu)%fpc_grid = 0.0) then
        !------------------------------------------------------------------
        ! Add individuals
        !------------------------------------------------------------------
        do pft=1,npft
          if (params_pft_plant(pft)%lu_category==lu) then

            ! initialise all pools of this PFT with zero
            call initpft( plant(pft) )

            ! get annually updated leaf traits (vary because of variations in light and CO2)
            call get_leaftraits( plant(pft), solar%meanmppfd(:), out_pmodel%actnv_unitiabs(:) )

            ! add a "seed" by forcing initial diameter increment
            call update_tree( plant(pft), diam_inc_init )

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
