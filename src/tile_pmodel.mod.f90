module md_tile
  !////////////////////////////////////////////////////////////////
  ! Holds all tile-specific variables and procedurs
  ! --------------------------------------------------------------
  use md_params_core, only: npft, nlu
  use md_params_soil, only: paramtype_soil
  use md_plant, only: plant_type, plant_fluxes_type, initglobal_plant

  implicit none

  private
  public tile_type, tile_fluxes_type, initglobal_tile, psoilphystype, soil_type, initdaily_tile


  !----------------------------------------------------------------
  ! Tile type with year-to-year memory
  !----------------------------------------------------------------
  type tile_type
    integer                           :: luno       ! Index that goes along with this instance of 'tile'
    type(soil_type)                   :: soil       ! all organic, inorganic, and physical soil variables
    type(canopy_type)                 :: canopy     ! mean canopy
    type(plant_type), dimension(npft) :: plant
  end type tile_type

  !----------------------------------------------------------------
  ! Soil type
  !----------------------------------------------------------------
  type soil_type
    type(psoilphystype)  :: phy      ! soil state variables
    type(paramtype_soil) :: params   ! soil parameters
  end type soil_type

  !----------------------------------------------------------------
  ! physical soil state variables with memory from year to year (~pools)
  !----------------------------------------------------------------
  type psoilphystype
    real :: temp        ! soil temperature [deg C]
    real :: wcont       ! liquid soil water mass [mm = kg/m2]
    real :: wscal       ! relative soil water content, between 0 and 1
    real :: snow        ! snow depth in liquid-water-equivalents [mm = kg/m2]
    real :: rlmalpha    ! rolling mean of annual mean alpha (AET/PET)
  end type psoilphystype

  !----------------------------------------------------------------
  ! Canopy type
  ! Contains tile-level aggregated variables related to the canopy
  !----------------------------------------------------------------
  type canopy_type
    real :: lai            ! leaf area index 
    real :: fapar          ! fraction of absorbed photosynthetically active radiation (unitless)
    real :: height         ! canopy height (m)
    real :: conductance    ! canopy conductance, upscaled from leaf-level stomatal conductance (m s-1)
    real :: fpc_grid       ! fractional projective cover (sum of crownarea by canopy plants)
  end type canopy_type


  !----------------------------------------------------------------
  ! Variables without memory (not necessarily just fluxes)
  !----------------------------------------------------------------
  type tile_fluxes_type
    type(canopy_fluxes_type) :: canopy
    type(plant_fluxes_type), dimension(npft) :: plant
  end type tile_fluxes_type


  type canopy_fluxes_type
    ! water
    real :: dro             ! daily runoff (mm d-1)
    real :: dfleach         ! daily fraction of soil water going to runoff (used for calculating leaching)
    real :: dwbal           ! daily water balance as precipitation and snow melt minus runoff and evapotranspiration (mm d-1)
    real :: econ            ! water-to-energy conversion factor (econ), m^3/J
    real :: rn              ! daily net radiation (J/m2/d)
    real :: rnn             ! nighttime net radiation (J m-1 d-1)
    real :: rnl             ! net longwave radiation (W m-2)
    real :: cn              ! daily condensation (mm d-1)
    real :: daet            ! daily total evapotranspiration (mm d-1)
    real :: daet_e          ! daily total evapotranspiration (J m-2 d-1)
    real :: daet_soil       ! daily soil evaporation (mm d-1)
    real :: daet_e_soil     ! daily soil evaporation (J m-2 d-1)
    real :: daet_canop      ! daily canopy transpiration (mm d-1)
    real :: daet_e_canop    ! daily canopy transpiration (J m-2 d-1)
    ! real :: sw        ! evaporative supply rate (mm/h)
    ! real :: rho_air         ! density of air (g m-3)
    ! real :: sat_slope       ! slope of saturation vap press temp curve, Pa/K 
    ! real :: lv              ! enthalpy of vaporization, J/kg
    ! real :: rho_water       ! density of water (g m-3)

    ! carbon 
    real :: dgpp
    real :: drd

    ! radiation
    real :: ppfd_splash
    real :: dayl
    real :: dra
    real :: nu
    real :: lambda

    ! real, dimension(ndayyear) :: dayl               ! day length (hours)
    ! real, dimension(ndayyear) :: dra                ! daily TOA solar irradiation (J/m2)
    ! real, dimension(ndayyear) :: dppfd_splash       ! daily total PPFD (mol m-2 d-1)
    ! real, dimension(nmonth)   :: mppfd_splash       ! monthly total PPFD (mol m-2 month-1)
    ! real, dimension(nmonth)   :: meanmppfd_splash   ! monthly mean PPFD, averaged over daylight seconds (mol m-2 s-1)

  end type canopy_fluxes_type

contains

  subroutine initglobal_tile( tile, ngridcells )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of all _pools on all gridcells at the beginning
    !  of the simulation.
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    ! argument
    integer, intent(in) :: ngridcells
    type(tile_type), dimension(nlu,ngridcells), intent(inout) :: tile

    ! local variables
    integer :: lu
    integer :: jpngr

    !-----------------------------------------------------------------------------
    ! derive which PFTs are present from fpc_grid (which is prescribed)
    !-----------------------------------------------------------------------------
    ! allocate( tile(nlu,ngridcells) )

    do jpngr=1,ngridcells
      do lu=1,nlu
        
        tile(lu,jpngr)%luno = lu

        ! initialise soil variables
        call initglobal_soil( tile(lu,jpngr)%soil )

        ! initialise canopy variables
        call initglobal_canopy( tile(lu,jpngr)%canopy )

        ! initialise plant variables
        call initglobal_plant( tile(lu,jpngr)%plant(:) )

      end do
    end do

  end subroutine initglobal_tile


  subroutine initglobal_canopy( canopy )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of specified PFT on specified gridcell
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    ! argument
    type(canopy_type), intent(inout) :: canopy

    lai         = 0.0
    fapar       = 0.0
    height      = 0.0
    conductance = 0.0
    fpc_grid    = 0.0

  end subroutine initglobal_canopy


  subroutine initglobal_soil( soil )
    !////////////////////////////////////////////////////////////////
    ! initialise soil variables globally
    !----------------------------------------------------------------
    use md_interface, only: interface

    ! argument
    type(soil_type), intent(inout) :: soil

    call initglobal_soil_phy( soil%phy )

    ! Copy soil parameters
    soil%params = interface%soilparams(jpngr)

  end subroutine initglobal_soil


  subroutine initglobal_soil_phy( phy )
    !////////////////////////////////////////////////////////////////
    ! initialise physical soil variables globally
    !----------------------------------------------------------------
    ! argument
    type(psoilphystype), intent(inout) :: phy

    ! initialise physical soil variables
    phy%wcont    = 50.0
    phy%temp     = 10.0
    phy%snow     = 0.0
    phy%rlmalpha = 0.0

  end subroutine initglobal_soil_phy


  subroutine initdaily_tile( tile_fluxes )
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables within derived type 'soilphys'.
    !----------------------------------------------------------------
    ! arguments
    type(tile_fluxes_type), dimension(nlu), intent(inout) :: tile_fluxes

    tile_fluxes(:)%canopy%dro = 0.0             
    tile_fluxes(:)%canopy%dfleach = 0.0         
    tile_fluxes(:)%canopy%dwbal = 0.0           
    tile_fluxes(:)%canopy%econ = 0.0            
    tile_fluxes(:)%canopy%rn = 0.0              
    tile_fluxes(:)%canopy%rnn = 0.0             
    tile_fluxes(:)%canopy%rnl = 0.0             
    tile_fluxes(:)%canopy%cn = 0.0              
    tile_fluxes(:)%canopy%daet = 0.0            
    tile_fluxes(:)%canopy%daet_e = 0.0          
    tile_fluxes(:)%canopy%daet_soil = 0.0       
    tile_fluxes(:)%canopy%daet_e_soil = 0.0     
    tile_fluxes(:)%canopy%daet_canop = 0.0      
    tile_fluxes(:)%canopy%daet_e_canop = 0.0    
    tile_fluxes(:)%canopy%dgpp = 0.0
    tile_fluxes(:)%canopy%drd = 0.0
    tile_fluxes(:)%canopy%ppfd_splash = 0.0
    tile_fluxes(:)%canopy%dayl = 0.0
    tile_fluxes(:)%canopy%dra = 0.0
    tile_fluxes(:)%canopy%nu = 0.0
    tile_fluxes(:)%canopy%lambda = 0.0

    call initdaily_plant( tile_fluxes(:)%plant(:) )

  end subroutine initdaily_tile

end module md_tile
