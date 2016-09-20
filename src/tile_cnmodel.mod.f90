module md_tile
  !////////////////////////////////////////////////////////////////
  ! Holds all tile-specific variables and procedurs
  ! --------------------------------------------------------------
  use md_params_core, only: npft
  use md_littersom, only: orgsoil_type
  use md_classdefs, only: orgpool, nitrogen

  implicit none

  private
  public tile_type, soil_type, orgsoil_type, inosoil_type, physoil_type, initglobal_tile

  !----------------------------------------------------------------
  ! organic soil pools
  !----------------------------------------------------------------
  type orgsoil_type
    type( orgpool ) :: sl ! soil organic matter, fast turnover [gC/m2]
    type( orgpool ) :: fs ! soil organic matter, slow turnover [gC/m2]
  end type orgsoil_type

  !----------------------------------------------------------------
  ! inorganic soil pools
  !----------------------------------------------------------------
  type inosoil_type
    type( nitrogen ) :: no3
    type( nitrogen ) :: nh4
  end type inosoil_type

  !----------------------------------------------------------------
  ! physical soil state variables with memory from year to year (~pools)
  !----------------------------------------------------------------
  type psoilphystype
    type( nitrogen ) :: temp
    type( nitrogen ) :: wcont
  end type psoilphystype

  !----------------------------------------------------------------
  ! Soil type
  !----------------------------------------------------------------
  type soil_type
    type( orgsoil_type )  :: org
    type( inosoil_type )  :: ino
    type( psoilphystype ) :: phy
  end type soil_type

  !----------------------------------------------------------------
  ! Litter type
  !----------------------------------------------------------------
  type litter_type
    type( orgpool ) :: af
    type( orgpool ) :: as
    type( orgpool ) :: bg
  end type litter_type


  !----------------------------------------------------------------
  ! Tile type
  !----------------------------------------------------------------
  type tile_type

    ! Index that goes along with this instance of 'tile'
    integer :: luno

    ! stand state variables
    real :: fpc_grid                ! fractional projective cover (sum of crownarea by canopy plants)
    real, dimension(npft) :: nind

    ! all organic, inorganic, and physical soil variables
    type( soil_type ) :: soil

    ! all organic litter pools
    type( litter_type ) :: litter

    ! exudates pool (very short turnover) [gC/m2]
    type(carbon) :: pexud            

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