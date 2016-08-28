module md_tile
  !////////////////////////////////////////////////////////////////
  ! Holds all tile-specific variables and procedurs
  ! --------------------------------------------------------------
  use md_classdefs

  implicit none

  private
  public tile_type

  !----------------------------------------------------------------
  ! Tile type
  !----------------------------------------------------------------
  type tile_type
    real :: fpc_grid   ! fractional projective cover (sum of crownarea by canopy plants)
    real :: nind       ! root biomass [gC/ind.] (=rm_ind)
  end type tile_type

contains

end module md_tile
