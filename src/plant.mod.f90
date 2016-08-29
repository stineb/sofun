module md_plant
  !////////////////////////////////////////////////////////////////
  ! Holds all PFT-specific pools, fluxes, and IO-handling prodecures
  ! --------------------------------------------------------------
  use md_classdefs
  use md_params_core, only: nlu, npft, maxgrid, ndayyear, lunat, nmonth

  implicit none

  private
  public plant_type, pexud, plitt_af, plitt_as, plitt_bg,                 &
    dgpp, dnpp, drgrow, drleaf, drroot, drsapw, dcex,                     &
    dnup, dnup_fix,                                                       &
    params_pft_plant, params_plant, initglobal_plant, initpft,            &
    initdaily_plant, outdnpp, outdnup, outdcleaf, outdcroot, outdclabl,   &
    outdnlabl, outdclitt, outdNlitt, outdCsoil, outdNsoil, outdlai,       &
    outdfapar, ddoc,                                                      &
    dnarea_mb, dnarea_cw, dlma, dcton_lm, get_fapar,                      &
    initoutput_plant, initio_plant, getout_daily_plant,                   &
    getout_annual_plant, writeout_ascii_plant, getpar_modl_plant,         &
    get_leaftraits,                                                       & 
    outacveg2lit, outanveg2lit

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  type plant_type

    ! PFT index that goes along with this instance of 'plant'
    integer :: pftno

    ! pools
    type(orgpool) :: pleaf      ! leaf biomass [gC/ind.] (=lm_ind)
    type(orgpool) :: proot      ! root biomass [gC/ind.] (=rm_ind)
    type(orgpool) :: psapw      ! sapwood biomass [gC/ind.] (=sm_ind)
    type(orgpool) :: pwood      ! heartwood (non-living) biomass [gC/ind.] (=hm_ind)
    type(orgpool) :: plabl      ! labile pool, temporary storage of N and C [gC/ind.] (=bm_inc but contains also N) 
    
    ! geometry
    real :: height              ! tree height (m)
    real :: diam                ! tree basal diameter (m)
    real :: acrown              ! crown area (m2)
    real :: fcrown              ! crown fraction of tree height (unitless)

    ! canopy
    real :: lai_ind             ! fraction of absorbed photosynthetically active radiation
    real :: fapar_ind           ! fraction of absorbed photosynthetically active radiation

    ! leaf traits
    real :: narea               ! total leaf N per unit leaf area (gN m-2)
    real :: narea_metabolic     ! metabolic leaf N per unit leaf area (gN m-2)
    real :: narea_structural    ! structural leaf N per unit leaf area (gN m-2)
    real :: lma                 ! leaf mass per area (gC m-2)
    real :: sla                 ! specific leaf area (m2 gC-1)
    real :: nmass               ! leaf N per unit leaf mass, g N / g-dry mass
    real :: r_cton_leaf         ! leaf C:N ratio [gC/gN] 
    real :: r_ntoc_leaf         ! leaf N:C ratio [gN/gC]

  end type plant_type


  type(carbon),  dimension(nlu,maxgrid)  :: pexud            ! exudates pool (very short turnover) [gC/m2]
  
  type(orgpool), dimension(npft,maxgrid) :: plitt_af         ! above-ground litter, fast turnover [gC/m2]
  type(orgpool), dimension(npft,maxgrid) :: plitt_as         ! above-ground litter, slow turnover [gC/m2]
  type(orgpool), dimension(npft,maxgrid) :: plitt_bg         ! below-ground litter [gC/m2]

  ! fluxes
  real, dimension(npft)                  :: dgpp             ! daily gross primary production [gC/m2/d]
  type(carbon), dimension(npft)          :: dnpp             ! net primary production [gC/m2/d]
  real, dimension(npft)                  :: drgrow           ! growth respiration (growth+maintenance resp. of all compartments), no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
  real, dimension(npft)                  :: drleaf           ! leaf maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
  real, dimension(npft)                  :: drroot           ! root maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
  real, dimension(npft)                  :: drsapw           ! sapwood maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
  real, dimension(npft)                  :: dcex             ! labile C exudation for N uptake, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]

  real, dimension(nlu)                   :: ddoc             ! surrogate for dissolved organic carbon used for denitrification rate (see ntransform)
  type(nitrogen), dimension(npft)        :: dnup             ! daily N uptake [gN/m2/d]
  real, dimension(npft)                  :: dnup_fix         ! daily N uptake by plant symbiotic N fixation [gN/m2/d]

  !-----------------------------------------------------------------------
  ! Uncertain (unknown) parameters. Runtime read-in
  !-----------------------------------------------------------------------
  ! NON PFT-DEPENDENT PARAMETERS
  type params_plant_type
    real :: r_root            ! Fine root-specific respiration rate (gC gC-1 d-1)
    real :: r_sapw            ! Sapwood-specific respiration rate (gC gC-1 d-1)
    real :: exurate           ! Fine root-specific C export rate (gC gC-1 d-1)
    real :: kbeer             ! canopy light extinction coefficient
    real :: f_nretain         ! fraction of N retained at leaf abscission 
    real :: fpc_tree_max      ! maximum fractional plant coverage of trees
    real :: growtheff         ! growth respiration coefficient = yield factor [unitless]
    real :: cton_soil         ! C:N ratio of soil organic matter (consider this to be equal to that of microbial biomass)
  end type params_plant_type

  type( params_plant_type ) :: params_plant

  ! PFT-DEPENDENT PARAMETERS
  type params_pft_plant_type
    character(len=4) :: pftname    ! standard PFT name with 4 characters length
    integer :: lu_category         ! land use category associated with PFT
    logical, dimension(nlu) :: islu! islu(ipft,ilu) is true if ipft belongs to ilu
    logical :: grass               ! boolean for growth form 'grass'
    logical :: tree                ! boolean for growth form 'tree'
    logical :: nfixer              ! whether plant is capable of symbiotically fixing N
    logical :: c3                  ! whether plant follows C3 photosynthesis
    logical :: c4                  ! whether plant follows C4 photosynthesis
    real    :: k_decay_leaf_base   ! base leaf decay constant [year-1]
    real    :: k_decay_leaf_width  ! shape parameter for turnover function if LAI
    real    :: k_decay_root        ! root decay constant [year-1]
    real    :: k_decay_labl        ! labile pool decay constant [year-1]
    real    :: r_cton_root         ! C:N ratio in roots (gC/gN)
    real    :: r_ntoc_root         ! N:C ratio in roots (inverse of 'r_cton_root', gN/gC)
    real    :: r_cton_wood         ! C:N ratio in wood (gC/gN)
    real    :: r_ntoc_wood         ! N:C ratio in wood (inverse of 'r_cton_root', gN/gC)
    real    :: sla                 ! specific leaf area (m2 gC-1)
    real    :: lma                 ! leaf mass per area (gC m-2)
    real    :: r_ntolma            ! constant ratio of structural N to C (LMA) (gN/gC)
    real    :: lai_ind             ! constant leaf area index within crown of an individual
  end type params_pft_plant_type

  type( params_pft_plant_type ), dimension(npft) :: params_pft_plant

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, allocatable, dimension(:,:,:) :: outdgpp    ! daily gross primary production [gC/m2/d]
  real, allocatable, dimension(:,:,:) :: outdnpp
  real, allocatable, dimension(:,:,:) :: outdnup
  real, allocatable, dimension(:,:,:) :: outdcex
  real, allocatable, dimension(:,:,:) :: outdcleaf
  real, allocatable, dimension(:,:,:) :: outdcroot
  real, allocatable, dimension(:,:,:) :: outdclabl
  real, allocatable, dimension(:,:,:) :: outdnlabl
  real, allocatable, dimension(:,:,:) :: outdclitt
  real, allocatable, dimension(:,:,:) :: outdNlitt
  real, allocatable, dimension(:,:,:) :: outdCsoil
  real, allocatable, dimension(:,:,:) :: outdNsoil
  real, allocatable, dimension(:,:,:) :: outdlai
  real, allocatable, dimension(:,:,:) :: outdfapar

  ! These are stored as dayly variables for annual output
  ! at day of year when LAI is at its maximum.
  real, dimension(npft,ndayyear) :: dnarea_mb
  real, dimension(npft,ndayyear) :: dnarea_cw
  real, dimension(npft,ndayyear) :: dlma
  real, dimension(npft,ndayyear) :: dcton_lm

  ! annual
  real, dimension(npft,maxgrid) :: outagpp
  real, dimension(npft,maxgrid) :: outanpp
  real, dimension(npft,maxgrid) :: outanup
  real, dimension(npft,maxgrid) :: outanup_fix
  real, dimension(npft,maxgrid) :: outacex
  real, dimension(npft,maxgrid) :: outacveg2lit
  real, dimension(npft,maxgrid) :: outanveg2lit
  real, dimension(npft,maxgrid) :: outanarea_mb
  real, dimension(npft,maxgrid) :: outanarea_cw
  real, dimension(npft,maxgrid) :: outalai
  real, dimension(npft,maxgrid) :: outalma
  real, dimension(npft,maxgrid) :: outacton_lm
  real, dimension(npft,maxgrid) :: outacroot
  real, dimension(npft,maxgrid) :: outacleaf
  real, dimension(npft,maxgrid) :: outaclabl
  real, dimension(npft,maxgrid) :: outanlabl

contains

  function get_fapar( lai ) result( fapar )
    !////////////////////////////////////////////////////////////////
    ! FOLIAGE PROJECTIVE COVER 
    ! = Fraction of Absorbed Photosynthetically Active Radiation
    ! Function returns fractional plant cover an individual
    ! Eq. 7 in Sitch et al., 2003
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: lai

    ! function return variable
    real :: fapar

    fapar = ( 1.0 - exp( -1.0 * params_plant%kbeer * lai ) )

  end function get_fapar


  function get_leaf_n_metabolic_canopy( mylai, meanmppfd, nv, myfapar ) result( mynleaf_metabolic )
    !////////////////////////////////////////////////////////////////
    ! Calculates metabolic leaf N at canopy-level, determined by 
    ! light conditions (meanmppfd) and the Rubisco-N per unit absorbed
    ! light.
    !----------------------------------------------------------------
    use md_params_core, only: nmonth

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv
    real, intent(in), optional          :: myfapar

    ! function return variable
    real :: mynleaf_metabolic  ! mol N m-2-ground

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )

    if (present(myfapar)) then
      mynleaf_metabolic = maxnv * myfapar
    else
      mynleaf_metabolic = maxnv * get_fapar( mylai )
    end if

  end function get_leaf_n_metabolic_canopy


  subroutine get_leaftraits( plant, meanmppfd, nv )
    !////////////////////////////////////////////////////////////////
    ! Calculates leaf traits based on (predicted) metabolic Narea and
    ! (prescribed) parameters that relate structural to metabolic
    ! Narea and Carea to structural Narea:
    ! Narea_metabolic  = predicted
    ! Narea_structural = rN:C_struct * LMA
    !----------------------------------------------------------------
    use md_params_core, only: c_content_of_biomass, nmonth, n_molmass, c_molmass

    ! arguments
    type( plant_type ), intent(inout)   :: plant
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! local variables
    real :: narea_metabolic_canopy   ! g N m-2-ground

    ! canopy-level, in units of gN / m2-ground 
    narea_metabolic_canopy  = n_molmass * get_leaf_n_metabolic_canopy(  params_pft_plant(plant%pftno)%lai_ind, meanmppfd(:), nv(:) )

    ! leaf-level, in units of gN / m2-leaf 
    plant%narea_metabolic  = narea_metabolic_canopy / params_pft_plant(plant%pftno)%lai_ind
    plant%narea_structural = params_pft_plant(plant%pftno)%r_ntolma * params_pft_plant(plant%pftno)%lma
    plant%narea            = plant%narea_metabolic + plant%narea_structural

    ! additional traits
    plant%nmass            = plant%narea / ( plant%lma / c_content_of_biomass )
    plant%r_cton_leaf      = params_pft_plant(plant%pftno)%lma / plant%narea
    plant%r_ntoc_leaf      = 1.0 / plant%r_cton_leaf

  end subroutine get_leaftraits


  subroutine getpar_modl_plant()
    !////////////////////////////////////////////////////////////////
    !  Subroutine reads model parameters from input file.
    !  It was necessary to separate this SR from module md_plant
    !  because this SR uses module md_waterbal, which also uses
    !  _plant.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------    
    use md_sofunutils, only: getparreal
    use md_interface, only: interface

    ! local variables
    integer :: pft
    integer :: npft_site

    !----------------------------------------------------------------
    ! NON-PFT DEPENDENT PARAMETERS
    !----------------------------------------------------------------
    ! canopy light extinction coefficient for Beer's Law
    params_plant%kbeer = getparreal( 'params/params_plant.dat', 'kbeer' )

    ! growth efficiency = yield factor, central value: 0.6, range: 0.5-0.7; Zhang et al. (2009), see Li et al., 2014
    params_plant%growtheff = getparreal( 'params/params_plant.dat', 'growtheff' )

    ! Fine-root mass specific respiration rate (gC gC-1 year-1)
    ! Central value: 0.913 year-1 (Yan and Zhao (2007); see Li et al., 2014)
    params_plant%r_root = getparreal( 'params/params_plant.dat', 'r_root' ) / ndayyear  ! conversion to rate per day

    ! Fine-root specific respiration rate (gC gC-1 year-1)
    ! Central value: 0.044 year-1 (Yan and Zhao (2007); see Li et al., 2014)
    ! (= 0.044 nmol mol-1 s-1; range: 0.5–10, 20 nmol mol-1 s-1 (Landsberg and Sands (2010))
    params_plant%r_sapw = getparreal( 'params/params_plant.dat', 'r_sapw' ) / ndayyear  ! conversion to rate per day

    ! C export rate per unit root mass
    params_plant%exurate = getparreal( 'params/params_plant.dat', 'exurate' )

    ! C:N ratio of soil organic matter [1]
    params_plant%cton_soil = getparreal( 'params/params_plant.dat', 'cton_soil' )

    !----------------------------------------------------------------
    ! PFT DEPENDENT PARAMETERS
    ! read parameter input file and store values in single array
    ! important: Keep this order of reading PFT parameters fixed.
    !----------------------------------------------------------------
    pft = 0
    if ( interface%params_siml%lTeBE ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TeBE' )

    else if ( interface%params_siml%lGrC3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GrC3' )

    else if ( interface%params_siml%lGNC3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GNC3' )

    else if ( interface%params_siml%lGrC4 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GrC4' )

    else
      stop 'PLANT:GETPAR_MODL_PLANT: PFT name not valid. See run/<simulationname>.sofun.parameter'
    end if
    npft_site = pft

  end subroutine getpar_modl_plant


  function getpftparams( pftname ) result( out_getpftparams )
    !----------------------------------------------------------------
    ! Read PFT parameters from respective file, given the PFT name
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    character(len=*), intent(in) :: pftname

    ! local variables
    real :: lu_category_prov    ! land use category associated with PFT (provisional)

    ! function return variable
    type( params_pft_plant_type ) :: out_getpftparams

    ! standard PFT name
    out_getpftparams%pftname = pftname

    ! PFT names
    ! GrC3 : C3 grass                          
    ! GrC4 : C4 grass     
    if (trim(pftname)=='TeBE') then
      out_getpftparams%grass   = .false.
      out_getpftparams%tree    = .true.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='GrC3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='GNC3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .true.
    else if (trim(pftname)=='GrC4') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .false.
      out_getpftparams%c4      = .true.
      out_getpftparams%nfixer  = .false.
    end if      
    
    ! land use category associated with PFT (provisional)
    lu_category_prov = getparreal( trim('params/params_plant_'//trim(pftname)//'.dat'), 'lu_category_prov' )
    if (lu_category_prov==1.0) then
      out_getpftparams%lu_category = lunat
      out_getpftparams%islu(lunat) = .true.
    else
      out_getpftparams%islu(lunat) = .false.
    end if

    ! leaf decay constant, read in as [years-1], central value: 0.0 yr-1 for deciduous plants
    out_getpftparams%k_decay_leaf_base = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_leaf_base' ) / ndayyear 

    ! shape parameter for turnover function if LAI
    out_getpftparams%k_decay_leaf_width = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_leaf_width' )

    ! root decay constant [days], read in as [years-1], central value: 1.04 (Shan et al., 1993; see Li et al., 2014)
    out_getpftparams%k_decay_root = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_root' ) / ndayyear 

    ! root decay constant [days], read in as [years-1]
    out_getpftparams%k_decay_labl = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_labl' ) / ndayyear 

    ! root C:N and N:C ratio (gC/gN and gN/gC)
    out_getpftparams%r_cton_root = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_cton_root' )
    out_getpftparams%r_ntoc_root = 1.0 / out_getpftparams%r_cton_root

    ! wood and sapwood C:N and N:C ratio (gC/gN and gN/gC)
    out_getpftparams%r_cton_wood = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_cton_wood' )
    out_getpftparams%r_ntoc_wood = 1.0 / out_getpftparams%r_cton_wood

    ! leaf mass per area (gC m-2)
    out_getpftparams%lma = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'lma' )
    out_getpftparams%sla = 1.0 / out_getpftparams%lma

    ! constant ratio of leaf structural N to LMA
    out_getpftparams%r_ntolma = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_ntolma' )

    ! constant leaf area index within the crown of an individual
    out_getpftparams%lai_ind = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'lai_ind' )

  end function getpftparams


  subroutine initglobal_plant( plant )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of all _pools on all gridcells at the beginning
    !  of the simulation.
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_params_core, only: npft, maxgrid

    ! argument
    type( plant_type ), dimension(npft,maxgrid), intent(inout) :: plant

    ! local variables
    integer :: pft
    integer :: jpngr

    !-----------------------------------------------------------------------------
    ! derive which PFTs are present from fpc_grid (which is prescribed)
    !-----------------------------------------------------------------------------
    do jpngr=1,maxgrid
      do pft=1,npft
        call initpft( plant(pft,jpngr) )
        plant(pft,jpngr)%pftno = pft
      end do
    end do
 
    ! initialise all PFT-specific non-plant pools with zero (not done in 'initpft')
    pexud(:,:)    = carbon(0.0)                           ! exudates in soil, carbon pool [gC/m2]

    plitt_af(:,:) = orgpool(carbon(0.0),nitrogen(0.0))    ! above-ground fine   litter, organic pool [gC/m2 and gN/m2]
    plitt_as(:,:) = orgpool(carbon(0.0),nitrogen(0.0))    ! above-ground coarse litter, organic pool [gC/m2 and gN/m2]
    plitt_bg(:,:) = orgpool(carbon(0.0),nitrogen(0.0))    ! below-ground fine   litter, organic pool [gC/m2 and gN/m2]


  end subroutine initglobal_plant


  subroutine initpft( plant )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of specified PFT on specified gridcell
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    ! argument
    type( plant_type ), intent(inout) :: plant

    ! initialise all pools with zero
    plant%pleaf = orgpool(carbon(0.0),nitrogen(0.0))
    plant%proot = orgpool(carbon(0.0),nitrogen(0.0))
    plant%psapw = orgpool(carbon(0.0),nitrogen(0.0))
    plant%pwood = orgpool(carbon(0.0),nitrogen(0.0))
    plant%plabl = orgpool(carbon(0.0),nitrogen(0.0))

    ! geometry
    plant%height = 0.0
    plant%diam   = 0.0
    plant%acrown = 0.0
    plant%fcrown = 0.0

    ! canpopy state variables
    plant%narea            = 0.0
    plant%narea_metabolic  = 0.0
    plant%narea_structural = 0.0
    plant%lma              = 0.0
    plant%sla              = 0.0
    plant%nmass            = 0.0
    plant%r_cton_leaf      = 0.0
    plant%r_ntoc_leaf      = 0.0

    ! canopy variables, fixed in T-model for an individual
    plant%lai_ind   = params_pft_plant(plant%pftno)%lai_ind
    plant%fapar_ind = get_fapar( plant%lai_ind )

  end subroutine initpft


  subroutine initdaily_plant()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    !----------------------------------------------------------------
    dgpp(:)     = 0.0
    dnpp(:)     = carbon(0.0)
    dcex(:)     = 0.0
    dnup(:)     = nitrogen(0.0)
    dnup_fix(:) = 0.0
    drgrow(:)   = 0.0
    drleaf(:)   = 0.0
    drroot(:)   = 0.0
    drsapw(:)   = 0.0

  end subroutine initdaily_plant


  subroutine initoutput_plant()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_interface, only: interface

    if (interface%steering%init .and. interface%params_siml%loutdgpp  ) allocate( outdgpp      (npft,ndayyear,maxgrid) )
    if (interface%steering%init .and. interface%params_siml%loutdnpp  ) allocate( outdnpp      (npft,ndayyear,maxgrid) )
    if (interface%steering%init .and. interface%params_siml%loutdnup  ) allocate( outdnup      (npft,ndayyear,maxgrid) )
    if (interface%steering%init .and. interface%params_siml%loutdcex  ) allocate( outdcex      (npft,ndayyear,maxgrid) )
    if (interface%steering%init)                                        allocate( outdcleaf    (npft,ndayyear,maxgrid) )
    if (interface%steering%init)                                        allocate( outdcroot    (npft,ndayyear,maxgrid) )
    if (interface%steering%init .and. interface%params_siml%loutdclabl) allocate( outdclabl    (npft,ndayyear,maxgrid) )
    if (interface%steering%init .and. interface%params_siml%loutdnlabl) allocate( outdnlabl    (npft,ndayyear,maxgrid) )
    if (interface%steering%init .and. interface%params_siml%loutdclitt) allocate( outdclitt    (npft,ndayyear,maxgrid) )
    if (interface%steering%init .and. interface%params_siml%loutdnlitt) allocate( outdNlitt    (npft,ndayyear,maxgrid) )
    if (interface%steering%init .and. interface%params_siml%loutdfapar) allocate( outdfapar    (npft,ndayyear,maxgrid) )

    ! this is needed also for other (annual) output variables
    allocate( outdlai(npft,ndayyear,maxgrid) )
    
    outdgpp  (:,:,:) = 0.0
    outdnpp  (:,:,:) = 0.0
    outdnup  (:,:,:) = 0.0
    outdcex  (:,:,:) = 0.0
    outdcleaf(:,:,:) = 0.0
    outdcroot(:,:,:) = 0.0
    outdclabl(:,:,:) = 0.0
    outdnlabl(:,:,:) = 0.0
    outdclitt(:,:,:) = 0.0
    outdNlitt(:,:,:) = 0.0
    outdlai  (:,:,:) = 0.0
    outdfapar(:,:,:) = 0.0

    ! annual output variables
    if (interface%params_siml%loutplant) then
      outagpp(:,:)      = 0.0
      outanpp(:,:)      = 0.0
      outanup(:,:)      = 0.0
      outanup_fix(:,:)  = 0.0
      outacex(:,:)      = 0.0
      outacveg2lit(:,:) = 0.0
      outanveg2lit(:,:) = 0.0
      outanarea_mb(:,:) = 0.0
      outanarea_cw(:,:) = 0.0
      outalai     (:,:) = 0.0
      outalma     (:,:) = 0.0
      outacton_lm (:,:) = 0.0
      outacleaf(:,:)    = 0.0
      outacroot(:,:)    = 0.0
      outaclabl(:,:)    = 0.0
      outanlabl(:,:)    = 0.0
    end if

  end subroutine initoutput_plant


  subroutine initio_plant()
    !////////////////////////////////////////////////////////////////
    ! Opens input/output files.
    !----------------------------------------------------------------
    use md_interface, only: interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !////////////////////////////////////////////////////////////////
    ! DAILY OUTPUT: OPEN ASCII OUTPUT FILES 
    !----------------------------------------------------------------
    ! GPP
    if (interface%params_siml%loutdgpp) then
      filnam=trim(prefix)//'.d.gpp.out'
      open(101,file=filnam,err=999,status='unknown')
    end if 

    ! NPP
    if (interface%params_siml%loutdnpp) then 
      filnam=trim(prefix)//'.d.npp.out'
      open(102,file=filnam,err=999,status='unknown')
    end if

    ! LEAF C
    if (interface%params_siml%loutdcleaf    ) then
      filnam=trim(prefix)//'.d.cleaf.out'
      open(103,file=filnam,err=999,status='unknown')
    end if 

    ! N UPTAKE
    if (interface%params_siml%loutdnup      ) then
      filnam=trim(prefix)//'.d.nup.out'
      open(104,file=filnam,err=999,status='unknown')
    end if

    ! C EXUDATION
    if (interface%params_siml%loutdcex      ) then
      filnam=trim(prefix)//'.d.cex.out'
      open(105,file=filnam,err=999,status='unknown')
    end if

    ! ! AIR TEMPERATURE
    ! filnam=trim(prefix)//'.d.temp.out'
    ! open(110,file=filnam,err=999,status='unknown')

    ! ROOT C
    if (interface%params_siml%loutdcroot    ) then
      filnam=trim(prefix)//'.d.croot.out'
      open(111,file=filnam,err=999,status='unknown')
    end if

    ! LABILE C
    if (interface%params_siml%loutdclabl    ) then
      filnam=trim(prefix)//'.d.clabl.out'
      open(112,file=filnam,err=999,status='unknown')
    end if

    ! LITTER C
    if (interface%params_siml%loutdclitt    ) then
      filnam=trim(prefix)//'.d.clitt.out'
      open(113,file=filnam,err=999,status='unknown')
    end if

    ! LABILE N
    if (interface%params_siml%loutdnlabl    ) then
      filnam=trim(prefix)//'.d.nlabl.out'
      open(115,file=filnam,err=999,status='unknown')
    end if

    ! LITTER N
    if (interface%params_siml%loutdnlitt    ) then
      filnam=trim(prefix)//'.d.nlitt.out'
      open(119,file=filnam,err=999,status='unknown')
    end if

    ! LAI
    if (interface%params_siml%loutdlai      ) then
      filnam=trim(prefix)//'.d.lai.out'
      open(121,file=filnam,err=999,status='unknown')
    end if

    ! FAPAR
    if (interface%params_siml%loutdfapar    ) then
      filnam=trim(prefix)//'.d.fapar.out'
      open(122,file=filnam,err=999,status='unknown')
    end if

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------
    if (interface%params_siml%loutplant) then

      ! GPP 
      filnam=trim(prefix)//'.a.gpp.out'
      open(310,file=filnam,err=999,status='unknown')

      ! NPP 
      filnam=trim(prefix)//'.a.npp.out'
      open(311,file=filnam,err=999,status='unknown')      

      ! C VEGETATION -> LITTER TRANSFER
      filnam=trim(prefix)//'.a.cveg2lit.out'
      open(307,file=filnam,err=999,status='unknown')

      ! N VEGETATION -> LITTER TRANSFER
      filnam=trim(prefix)//'.a.nveg2lit.out'
      open(308,file=filnam,err=999,status='unknown')

      ! LEAF C
      filnam=trim(prefix)//'.a.cleaf.out'
      open(312,file=filnam,err=999,status='unknown')

      ! ROOT C
      filnam=trim(prefix)//'.a.croot.out'
      open(324,file=filnam,err=999,status='unknown')

      ! N UPTAKE
      filnam=trim(prefix)//'.a.nup.out'
      open(317,file=filnam,err=999,status='unknown')

      ! N FIXATION
      filnam=trim(prefix)//'.a.nup_fix.out'
      open(327,file=filnam,err=999,status='unknown')

      ! C EXUDATION
      filnam=trim(prefix)//'.a.cex.out'
      open(405,file=filnam,err=999,status='unknown')

      ! LAI (ANNUAL MAXIMUM)
      filnam=trim(prefix)//'.a.lai.out'
      open(318,file=filnam,err=999,status='unknown')

      ! METABOLIC NAREA (AT ANNUAL LAI MAXIMUM)
      filnam=trim(prefix)//'.a.narea_mb.out'
      open(319,file=filnam,err=999,status='unknown')

      ! CELL WALL NAREA (AT ANNUAL LAI MAXIMUM)
      filnam=trim(prefix)//'.a.narea_cw.out'
      open(320,file=filnam,err=999,status='unknown')

      ! LEAF C:N RATIO (AT ANNUAL LAI MAXIMUM)
      filnam=trim(prefix)//'.a.cton_lm.out'
      open(321,file=filnam,err=999,status='unknown')

      ! LMA (AT ANNUAL LAI MAXIMUM)
      filnam=trim(prefix)//'.a.lma.out'
      open(322,file=filnam,err=999,status='unknown')

      ! LABILE C AT THE END OF THE YEAR
      if (interface%params_siml%loutdclabl) then
        filnam=trim(prefix)//'.a.clabl.out'
        open(325,file=filnam,err=999,status='unknown')
      end if

      ! LABILE N AT THE END OF THE YEAR
      if (interface%params_siml%loutdnlabl) then
        filnam=trim(prefix)//'.a.nlabl.out'
        open(326,file=filnam,err=999,status='unknown')
      end if

    end if

    return

    999  stop 'INITIO: error opening output files'

  end subroutine initio_plant


  subroutine getout_daily_plant( plant, jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft
    use md_interface, only: interface

    ! arguments
    type( plant_type ), dimension(npft), intent(in) :: plant
    integer, intent(in)                             :: jpngr
    integer, intent(in)                             :: moy
    integer, intent(in) :: doy

    ! LOCAL VARIABLES
    integer :: pft

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    if (interface%params_siml%loutdgpp   ) outdgpp(:,doy,jpngr)   = dgpp(:)
    if (interface%params_siml%loutdnpp   ) outdnpp(:,doy,jpngr)   = dnpp(:)%c12
    if (interface%params_siml%loutdnup   ) outdnup(:,doy,jpngr)   = dnup(:)%n14
    if (interface%params_siml%loutdcex   ) outdcex(:,doy,jpngr)   = dcex(:)
    if (interface%params_siml%loutdcleaf .or. interface%params_siml%loutplant ) outdcleaf(:,doy,jpngr) = plant(:)%pleaf%c%c12
    if (interface%params_siml%loutdcroot .or. interface%params_siml%loutplant ) outdcroot(:,doy,jpngr) = plant(:)%proot%c%c12
    if (interface%params_siml%loutdclabl ) outdclabl(:,doy,jpngr) = plant(:)%plabl%c%c12
    if (interface%params_siml%loutdnlabl ) outdnlabl(:,doy,jpngr) = plant(:)%plabl%n%n14
    if (interface%params_siml%loutdclitt ) outdclitt(:,doy,jpngr) = plitt_af(:,jpngr)%c%c12 + plitt_as(:,jpngr)%c%c12 + plitt_bg(:,jpngr)%c%c12
    if (interface%params_siml%loutdnlitt ) outdNlitt(:,doy,jpngr) = plitt_af(:,jpngr)%n%n14 + plitt_as(:,jpngr)%n%n14 + plitt_bg(:,jpngr)%n%n14
    if (interface%params_siml%loutdfapar ) outdfapar(:,doy,jpngr) = plant(:)%fapar_ind

    ! this is needed also for other (annual) output variables
    outdlai(:,doy,jpngr) = plant(:)%lai_ind
      
    if (interface%params_siml%loutplant) then
      dnarea_mb(:,doy) = plant(:)%narea_metabolic 
      dnarea_cw(:,doy) = plant(:)%narea_structural
      dcton_lm(:,doy)  = plant(:)%r_cton_leaf
      dlma(:,doy)      = plant(:)%lma
    end if

    !----------------------------------------------------------------
    ! ANNUAL SUM/MEAN OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    if (interface%params_siml%loutplant) then
      outacleaf(:,jpngr)   = outacleaf(:,jpngr) + plant(:)%pleaf%c%c12 / ndayyear
      outacroot(:,jpngr)   = outacroot(:,jpngr) + plant(:)%proot%c%c12 / ndayyear
      outagpp(:,jpngr)     = outagpp(:,jpngr) + dgpp(:)
      outanpp(:,jpngr)     = outanpp(:,jpngr) + dnpp(:)%c12
      outanup(:,jpngr)     = outanup(:,jpngr) + dnup(:)%n14
      outanup_fix(:,jpngr) = outanup_fix(:,jpngr) + dnup_fix(:)
      outacex(:,jpngr)     = outacex(:,jpngr) + dcex(:)
    end if

  end subroutine getout_daily_plant


  subroutine getout_annual_plant( plant, jpngr )
    !////////////////////////////////////////////////////////////////
    !  SR called once a year to gather annual output variables.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft
    use md_interface, only: interface

    ! arguments
    type( plant_type ), dimension(npft), intent(in) :: plant
    integer, intent(in)                             :: jpngr

    ! local variables
    integer :: pft
    integer :: doy

    ! Output annual value at day of peak LAI
    if (interface%params_siml%loutplant) then
      outanarea_mb(:,jpngr) = plant(:)%narea_metabolic 
      outanarea_cw(:,jpngr) = plant(:)%narea_structural
      outalai     (:,jpngr) = plant(:)%lai_ind
      outalma     (:,jpngr) = plant(:)%lma
      outacton_lm (:,jpngr) = plant(:)%r_cton_leaf
    end if

    if (interface%params_siml%loutdclabl) outaclabl(:,jpngr) = outdclabl(:,ndayyear,jpngr) ! taken at the end of the year
    if (interface%params_siml%loutdnlabl) outanlabl(:,jpngr) = outdnlabl(:,ndayyear,jpngr) ! taken at the end of the year

  end subroutine getout_annual_plant


  subroutine writeout_ascii_plant()
    !/////////////////////////////////////////////////////////////////////////
    ! Write daily ASCII output
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !-------------------------------------------------------------------------
    use md_params_core, only: ndayyear
    use md_interface, only: interface

    ! local variables
    real :: itime
    integer :: day, moy, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1


    !-------------------------------------------------------------------------
    ! Collect variables to output variables
    !-------------------------------------------------------------------------
    !do lu=1,nlu
    !  do pft=1,npft
    !    if (lu==lu_category(pft)) then
    !    end if
    !  end do
    !end do
    if (nlu>1) stop 'Output only for one LU category implemented.'


    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if ( .not. interface%steering%spinup &
      .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
      .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      ! Write daily output only during transient simulation
      do day=1,ndayyear

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(interface%steering%outyear) + real(day-1)/real(ndayyear)
        
        if (interface%params_siml%loutdgpp  ) write(101,999) itime, sum(outdgpp(:,day,jpngr))
        if (interface%params_siml%loutdnpp  ) write(102,999) itime, sum(outdnpp(:,day,jpngr))
        if (interface%params_siml%loutdcleaf) write(103,999) itime, sum(outdcleaf(:,day,jpngr))
        if (interface%params_siml%loutdnup  ) write(104,999) itime, sum(outdnup(:,day,jpngr))
        if (interface%params_siml%loutdcex  ) write(105,999) itime, sum(outdcex(:,day,jpngr))
        if (interface%params_siml%loutdcroot) write(111,999) itime, sum(outdcroot(:,day,jpngr))
        if (interface%params_siml%loutdclabl) write(112,999) itime, sum(outdclabl(:,day,jpngr))
        if (interface%params_siml%loutdclitt) write(113,999) itime, sum(outdclitt(:,day,jpngr))
        if (interface%params_siml%loutdnlabl) write(115,999) itime, sum(outdnlabl(:,day,jpngr))
        if (interface%params_siml%loutdnlitt) write(119,999) itime, sum(outdNlitt(:,day,jpngr))
        if (interface%params_siml%loutdlai  ) write(121,999) itime, sum(outdlai(:,day,jpngr))
        if (interface%params_siml%loutdfapar) write(122,999) itime, sum(outdfapar(:,day,jpngr))
          
      end do
    end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutplant) then

      itime = real(interface%steering%outyear)

      write(307,999) itime, sum(outacveg2lit(:,jpngr))
      write(308,999) itime, sum(outanveg2lit(:,jpngr))
      write(310,999) itime, sum(outagpp(:,jpngr))
      write(311,999) itime, sum(outanpp(:,jpngr))
      write(312,999) itime, sum(outacleaf(:,jpngr))
      write(324,999) itime, sum(outacroot(:,jpngr))
      write(317,999) itime, sum(outanup(:,jpngr))
      write(327,999) itime, sum(outanup_fix(:,jpngr))
      write(405,999) itime, sum(outacex(:,jpngr))
      write(318,999) itime, sum(outalai(:,jpngr))
      write(319,999) itime, sum(outanarea_mb(:,jpngr))
      write(320,999) itime, sum(outanarea_cw(:,jpngr))
      write(321,999) itime, sum(outacton_lm(:,jpngr))
      write(322,999) itime, sum(outalma(:,jpngr))
      if (interface%params_siml%loutdclabl) write(325,999) itime, sum(outaclabl(:,jpngr))
      if (interface%params_siml%loutdnlabl) write(326,999) itime, sum(outanlabl(:,jpngr))
    end if

    return

  999     format (F20.8,F20.8)

  end subroutine writeout_ascii_plant

end module md_plant
