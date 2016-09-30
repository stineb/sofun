module md_tile
  !////////////////////////////////////////////////////////////////
  ! Holds all tile-specific variables and procedurs
  ! --------------------------------------------------------------
  use md_params_core, only: npft

  implicit none

  private
  public tile_type, initglobal_tile

  !----------------------------------------------------------------
  ! Tile type
  !----------------------------------------------------------------
  type tile_type

    ! Index that goes along with this instance of 'tile'
    integer :: luno

    ! stand state variables
    real :: fpc_grid                ! fractional projective cover (sum of crownarea by canopy plants)
    real, dimension(npft) :: nind   ! root biomass [gC/ind.] (=rm_ind)

  end type tile_type

contains

  subroutine initglobal_tile( tile )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of all _pools on all gridcells at the beginning
    !  of the simulation.
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_params_core, only: nlu, maxgrid

    ! argument
    type( tile_type ), dimension(nlu,maxgrid), intent(inout) :: tile

    ! local variables
    integer :: lu
    integer :: jpngr

    !-----------------------------------------------------------------------------
    ! derive which PFTs are present from fpc_grid (which is prescribed)
    !-----------------------------------------------------------------------------
    do jpngr=1,maxgrid
      do lu=1,nlu
        call initlu( tile(lu,jpngr) )
        tile(lu,jpngr)%luno = lu
      end do
    end do

  end subroutine initglobal_tile


  subroutine initlu( tile )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of specified PFT on specified gridcell
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    ! argument
    type( tile_type ), intent(out) :: tile

    tile%fpc_grid = 0.0
    tile%nind(:) = 0.0

  end subroutine initlu


end module md_tile