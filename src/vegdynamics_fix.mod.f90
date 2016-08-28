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
    use md_params_core, only: npft

    ! arguments
    type( tile_type ), intent(inout) :: tile
    type( solartype )                :: solar
    type( outtype_pmodel )           :: out_pmodel

    ! local variables
    integer :: pft

    if (tile%fpc_grid = 0.0) then

      do lu=1,nlu
        do pft=1,npft
          if (params_pft_plant(pft)%lu_category==lu) then

            ! initialise all pools of this PFT with zero
            call initpft( pft, jpngr )

            ! get annually updated leaf traits (vary because of variations in light and CO2)
            call get_leaftraits( tree, pft, solar%meanmppfd(:), out_pmodel%actnv_unitiabs(:) )

            ! add a "seed" by forcing initial diameter increment
            call update_tree( tree, diam_inc_init )

            ! xxx needs to be done: add implicit C and N fluxes to NPP and N-uptake

          end if
        end do
      end do

    end if

  end subroutine vegdynamics


end module md_vegdynamics
