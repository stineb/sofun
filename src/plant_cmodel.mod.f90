module md_plant
  !////////////////////////////////////////////////////////////////
  !  Module contains (constrainable) model parameters.
  !  Model parameters adopted here are from LPX C3 grass PFT
  !  Litter and soil turnover parameters are divided by 365 to 
  !  convert from [1/yr] to [1/d].
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core
  use md_classdefs, only: orgpool, carbon, nitrogen

  implicit none

  private
  public plant_type, plant_fluxes_type, getpar_modl_plant, params_pft_plant, &
    params_plant, &
    initdaily_plant, initoutput_plant, initio_plant, getout_daily_plant,     &
    writeout_ascii_plant, maxdoy, initglobal_plant, update_leaftraits, &
    update_leaftraits_init, initpft, getout_annual_plant, get_lai, seed

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  ! Pools and other variables with year-to-year memory
  !----------------------------------------------------------------
  type plant_type

    ! PFT index that goes along with this instance of 'plant'
    integer :: pftno

    ! canopy
    real :: nind                ! number of individuals
    real :: fpc_grid            ! fractional projective cover
    real :: lai_ind             ! fraction of absorbed photosynthetically active radiation
    real :: fapar_ind           ! fraction of absorbed photosynthetically active radiation
    real :: acrown              ! crown area

    ! leaf traits
    real :: leafc_canopy              ! g C m-2-ground, canopy-level
    real :: narea_canopy              ! g N m-2-ground, canopy-level
    real :: narea_metabolic_canopy    ! g N m-2-ground, canopy-level
    real :: narea_structural_canopy   ! g N m-2-ground, canopy-level
    real :: narea               ! total leaf N per unit leaf area (gN m-2)
    real :: narea_metabolic     ! metabolic leaf N per unit leaf area (gN m-2)
    real :: narea_structural    ! structural leaf N per unit leaf area (gN m-2)
    real :: lma                 ! leaf mass per area (gC m-2)
    real :: sla                 ! specific leaf area (m2 gC-1)
    real :: nmass               ! leaf N per unit leaf mass, g N / g-dry mass
    real :: r_cton_leaf         ! leaf C:N ratio [gC/gN] 
    real :: r_ntoc_leaf         ! leaf N:C ratio [gN/gC]

    ! pools
    type(orgpool) :: pleaf     ! leaf biomass [gC/ind.] (=lm_ind)
    type(orgpool) :: proot     ! root biomass [gC/ind.] (=rm_ind)
    type(orgpool) :: psapw     ! sapwood biomass [gC/ind.] (=sm_ind)
    type(orgpool) :: pwood     ! heartwood (non-living) biomass [gC/ind.] (=hm_ind)
    type(orgpool) :: plabl     ! labile pool, temporary storage of N and C [gC/ind.] (=bm_inc but contains also N) 
    
    type(carbon)  :: pexud     ! exudates pool (very short turnover) [gC/m2]
    
    type(orgpool) :: plitt_af  ! above-ground litter, fast turnover [gC/m2]
    type(orgpool) :: plitt_as  ! above-ground litter, slow turnover [gC/m2]
    type(orgpool) :: plitt_bg  ! below-ground litter [gC/m2]

  end type plant_type


  !----------------------------------------------------------------
  ! Fluxes and other variables with no memory
  !----------------------------------------------------------------
  type plant_fluxes_type

    real :: dgpp     ! daily gross primary production [gC/m2/d]           
    real :: drd      ! daily dark respiration [gC/m2/d]
    real :: drleaf   ! daily total leaf respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: drroot   ! root maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: drsapw   ! sapwood maintenance respiration, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: drgrow   ! growth respiration (growth+maintenance resp. of all compartments), no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: dcex     ! labile C exudation for N uptake, no explicit isotopic signature as it is identical to the signature of GPP [gC/m2/d]
    real :: dtransp  ! daily transpiration [mm]

    type(carbon) :: dnpp     ! daily net primary production (gpp-ra, npp=bp+cex) [gC/m2/d]

    type(nitrogen) :: dnup             ! daily N uptake [gN/m2/d]
    real           :: dnup_fix         ! daily N uptake by plant symbiotic N fixation [gN/m2/d]

  end type plant_fluxes_type

  !-----------------------------------------------------------------------
  ! Parameters. Runtime read-in
  !-----------------------------------------------------------------------
  ! NON PFT-DEPENDENT PARAMETERS
  type params_plant_type
    real :: kbeer             ! canopy light extinction coefficient
    real :: r_root            ! Fine root-specific respiration rate (gC gC-1 d-1)
    real :: r_sapw            ! Sapwood-specific respiration rate (gC gC-1 d-1)
    real :: exurate           ! Fine root-specific C export rate (gC gC-1 d-1)
    real :: f_nretain         ! fraction of N retained at leaf abscission 
    real :: fpc_tree_max      ! maximum fractional plant coverage of trees
    real :: growtheff         ! growth respiration coefficient = yield factor [unitless]
    real :: cton_soil         ! C:N ratio of soil organic matter (consider this to be equal to that of microbial biomass)
    real :: frac_leaf         ! fraction of allocatable C to leaf 
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
    ! real    :: sla                 ! specific leaf area (m2 gC-1)
    ! real    :: lma                 ! leaf mass per area (gC m-2)
    ! real    :: r_ntolma            ! constant ratio of structural N to C (LMA) (gN/gC)
    real    :: k_decay_leaf_base   ! base leaf decay constant [year-1]
    real    :: k_decay_leaf_width  ! shape parameter for turnover function if LAI
    real    :: k_decay_sapw        ! sapwood decay constant [year-1]
    real    :: k_decay_root        ! root decay constant [year-1]
    real    :: k_decay_labl        ! labile pool decay constant [year-1]
    real    :: r_cton_root         ! C:N ratio in roots (gC/gN)
    real    :: r_ntoc_root         ! N:C ratio in roots (inverse of 'r_cton_root', gN/gC)
    real    :: ncw_min             ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area    
    real    :: r_n_cw_v            ! slope in the relationship of non-metabolic versus metabolic N per leaf area              
    real    :: r_ctostructn_leaf   ! constant ratio of C to structural N (mol C / mol N)
  end type params_pft_plant_type

  type( params_pft_plant_type ), dimension(npft) :: params_pft_plant

  !-----------------------------------------------------------------------
  ! Fixed parameters
  !-----------------------------------------------------------------------
  ! type( orgpool ), parameter :: seed = orgpool( carbon(5.0), nitrogen(0.0) )
  type( orgpool ), parameter :: seed = orgpool( carbon(5.0), nitrogen(0.12) )
  ! type( orgpool ), parameter :: seed = orgpool( carbon(100.0), nitrogen(1 .0) )

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  ! annual
  real, dimension(:,:), allocatable :: outanarea_mb
  real, dimension(:,:), allocatable :: outanarea_cw
  real, dimension(:,:), allocatable :: outalai
  real, dimension(:,:), allocatable :: outalma
  real, dimension(:,:), allocatable :: outacton_lm

  ! required for outputting leaf trait variables in other modules
  integer, dimension(npft) :: maxdoy  ! DOY of maximum LAI

contains

  ! function get_canopy( lai ) result( out_canopy )
  !   !//////////////////////////////////////////////////////////////////
  !   ! Returs canopy variables as a function of LAI
  !   !------------------------------------------------------------------
  !   ! arguments
  !   real, intent(in) :: lai

  !   ! function return value
  !   type( canopy_type ) :: out_canopy

  !   out_canopy%fapar_ind = get_fapar( lai )

  ! end function get_canopy


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

    fapar = ( 1.0 - exp( -1.0 * params_plant%kbeer * lai) )

  end function get_fapar


  function get_lai( pft, cleaf, meanmppfd, nv ) result( lai )
    !////////////////////////////////////////////////////////////////
    ! Calculates LAI as a function of leaf-C. This is not so straight
    ! forward due to the dependence of canopy-metabolic leaf N on LAI,
    ! and the dependence of canopy-structural leaf N and C on canopy-
    ! metabolic leaf N.
    !----------------------------------------------------------------
    use md_params_core, only: nmonth, c_molmass
    use md_lambertw, only: calc_wapr

    ! arguments
    integer, intent(in)                 :: pft
    real, intent(in)                    :: cleaf
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv 

    ! function return variable
    real :: lai

    ! local variables
    real    :: alpha, beta, gamma ! variable substitutes
    real    :: maxnv
    real    :: arg_to_lambertw
    integer :: nerror


    if (cleaf>0.0) then

      ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
      ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
      maxnv = maxval( meanmppfd(:) * nv(:) )

      alpha = maxnv * params_pft_plant(pft)%r_n_cw_v
      beta  = params_pft_plant(pft)%ncw_min
      gamma = cleaf / ( c_molmass * params_pft_plant(pft)%r_ctostructn_leaf ) 

      arg_to_lambertw = alpha * params_plant%kbeer / beta * exp( (alpha - gamma) * params_plant%kbeer / beta )

      lai = 1.0 / (beta * params_plant%kbeer ) * ( -alpha * params_plant%kbeer + gamma * params_plant%kbeer + beta * calc_wapr( arg_to_lambertw, 0, nerror, 9999 ) )

    else

      lai = 0.0

    end if
    
  end function get_lai


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


  function get_leaf_n_structural_canopy( pft, mylai, mynleaf_metabolic ) result( mynleaf_structural )
    !////////////////////////////////////////////////////////////////
    ! Calculates structural leaf N at canopy-level, determined by 
    ! metabolic leaf N (linear relationship)
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: pft
    real, intent(in)    :: mylai
    real, intent(in)    :: mynleaf_metabolic

    ! function return variable
    real :: mynleaf_structural  ! mol N m-2-ground

    mynleaf_structural = mynleaf_metabolic * params_pft_plant(pft)%r_n_cw_v + mylai * params_pft_plant(pft)%ncw_min

  end function get_leaf_n_structural_canopy


  function get_leaf_n_canopy( pft, mylai, meanmppfd, nv ) result( mynleaf )
    !////////////////////////////////////////////////////////////////
    ! Calculates total leaf N at canopy-level, determined by 
    ! metabolic leaf N (linear relationship)
    ! Caution: this returns g N m-2-ground (not mol N m-2-ground)!
    !----------------------------------------------------------------
    use md_params_core, only: nmonth, n_molmass

    ! arguments
    integer, intent(in)                 :: pft
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! function return variable
    real :: mynleaf ! g N m-2-ground

    ! local variables
    real :: nleaf_metabolic   ! mol N m-2
    real :: nleaf_structural  ! mol N m-2

    nleaf_metabolic  = get_leaf_n_metabolic_canopy(  mylai, meanmppfd(:), nv(:) )
    nleaf_structural = get_leaf_n_structural_canopy( pft, mylai, nleaf_metabolic )
    mynleaf          = n_molmass * ( nleaf_metabolic + nleaf_structural )

  end function get_leaf_n_canopy


  subroutine update_leaftraits_init( plant, pft, meanmppfd, nv )
    !////////////////////////////////////////////////////////////////
    ! Calculates initial leaf traits (Taylor approximation for LAI -> 0)
    !----------------------------------------------------------------
    use md_params_core, only: c_content_of_biomass, nmonth, n_molmass, c_molmass

    ! arguments
    type( plant_type ), intent(inout)   :: plant
    integer, intent(in)                 :: pft
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! local variables
    real :: maxnv
    real :: mynarea_metabolic   ! mol N m-2-ground
    real :: mynarea_structural  ! mol N m-2-ground

    maxnv = maxval( meanmppfd(:) * nv(:) )

    mynarea_metabolic  = maxnv * params_plant%kbeer
    mynarea_structural = params_pft_plant(pft)%r_n_cw_v * maxnv * params_plant%kbeer + params_pft_plant(pft)%ncw_min

    ! leaf-level, in units of gN / m2-leaf 
    plant%narea_metabolic  = n_molmass * mynarea_metabolic  ! g N m-2-leaf
    plant%narea_structural = n_molmass * mynarea_structural ! g N m-2-leaf
    plant%narea            = n_molmass * ( mynarea_metabolic +  mynarea_structural ) ! g N m-2-leaf
    plant%lma              = c_molmass * params_pft_plant(pft)%r_ctostructn_leaf * mynarea_structural

    ! additional traits
    plant%nmass            = plant%narea / ( plant%lma / c_content_of_biomass )
    plant%r_cton_leaf      = plant%lma / plant%narea
    plant%r_ntoc_leaf      = 1.0 / plant%r_cton_leaf

    ! canopy-level, in units of gN / m2-ground 
    plant%narea_metabolic_canopy  = 0.0
    plant%narea_structural_canopy = 0.0
    plant%narea_canopy            = 0.0
    plant%leafc_canopy            = 0.0

  end subroutine update_leaftraits_init 


  subroutine update_leaftraits( plant, pft, mylai, meanmppfd, nv, myfapar )
    !////////////////////////////////////////////////////////////////
    ! Calculates leaf traits based on (predicted) metabolic Narea and
    ! (prescribed) parameters that relate structural to metabolic
    ! Narea and Carea to structural Narea:
    ! Narea_metabolic  = predicted
    ! Narea_structural = a + b * Narea_metabolic
    ! Carea            = c * Narea_structural
    !----------------------------------------------------------------
    use md_params_core, only: c_content_of_biomass, nmonth, n_molmass, c_molmass

    ! arguments
    type( plant_type ), intent(inout)   :: plant
    integer, intent(in)                 :: pft
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv
    real, intent(in), optional          :: myfapar

    ! local variables
    real :: mynarea_metabolic_canop   ! mol N m-2-ground
    real :: mynarea_structural_canop  ! mol N m-2-ground

    if (mylai==0.0) then
      ! canopy-level
      plant%narea_metabolic_canopy  = 0.0
      plant%narea_structural_canopy = 0.0
      plant%narea_canopy            = 0.0
      plant%leafc_canopy            = 0.0

      ! leaf-level
      plant%narea_metabolic  = 0.0
      plant%narea_structural = 0.0
      plant%narea            = 0.0
      plant%lma              = 0.0
      plant%nmass            = 0.0
      plant%r_cton_leaf      = 0.0
      plant%r_ntoc_leaf      = 0.0
    else
      ! calculate quantities in units of mol N
      mynarea_metabolic_canop  = get_leaf_n_metabolic_canopy(  mylai, meanmppfd(:), nv(:) )     ! mol N m-2-ground    
      mynarea_structural_canop = get_leaf_n_structural_canopy( pft, mylai, mynarea_metabolic_canop ) ! mol N m-2-ground

      ! canopy-level, in units of gN / m2-ground 
      plant%narea_metabolic_canopy  = n_molmass * mynarea_metabolic_canop ! g N m-2-ground 
      plant%narea_structural_canopy = n_molmass * mynarea_structural_canop ! g N m-2-ground
      plant%narea_canopy            = n_molmass * (mynarea_metabolic_canop + mynarea_structural_canop)  ! g N m-2-ground
      plant%leafc_canopy            = c_molmass * params_pft_plant(pft)%r_ctostructn_leaf * mynarea_structural_canop ! g C m-2-ground

      ! leaf-level, in units of gN / m2-leaf 
      plant%narea_metabolic  = plant%narea_metabolic_canopy / mylai   ! g N m-2-leaf
      plant%narea_structural = plant%narea_structural_canopy / mylai  ! g N m-2-leaf
      plant%narea            = plant%narea_canopy / mylai ! g N m-2-leaf
      plant%lma              = plant%leafc_canopy / mylai 

      ! additional traits
      plant%nmass            = plant%narea / ( plant%lma / c_content_of_biomass )
      plant%r_cton_leaf      = plant%lma / plant%narea
      plant%r_ntoc_leaf      = 1.0 / plant%r_cton_leaf
    end if

  end subroutine update_leaftraits


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
    use md_interface

    ! local variables
    integer :: pft
    integer :: npft_site

    !----------------------------------------------------------------
    ! NON-PFT DEPENDENT PARAMETERS
    !----------------------------------------------------------------
    ! canopy light extinction coefficient for Beer's Law
    params_plant%kbeer = getparreal( 'params/params_plant.dat', 'kbeer' )

    ! fraction of N retained at leaf abscission 
    params_plant%f_nretain = getparreal( 'params/params_plant.dat', 'f_nretain' )
    
    ! maximum fractional plant coverage of trees (sum of all tree PFTs)
    params_plant%fpc_tree_max = getparreal( 'params/params_plant.dat', 'fpc_tree_max' )

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
    if ( interface%params_siml%lTrE ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TrE' )
    end if

    if ( interface%params_siml%lTNE ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TNE' )
    end if

    if ( interface%params_siml%lTrD ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TrD' )
    end if

    if ( interface%params_siml%lTND ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TND' )
    end if

    if ( interface%params_siml%lGr3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'Gr3' )
    end if

    if ( interface%params_siml%lGN3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GN3' )
    end if

    if ( interface%params_siml%lGr4 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'Gr4' )
    end if

    npft_site = pft
    if (npft_site==0) stop 'PLANT:GETPAR_MODL_PLANT: PFT name not valid. See run/<simulationname>.sofun.parameter'

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
    real :: code_growthform
    real :: code_nfixer

    ! function return variable
    type( params_pft_plant_type ) :: out_getpftparams

    ! standard PFT name
    out_getpftparams%pftname = pftname

    ! PFT names
    ! Gr3 : C3 grass                          
    ! Gr4 : C4 grass     
    if (trim(pftname)=='Gr3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='GN3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .true.
    else if (trim(pftname)=='Gr4') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .false.
      out_getpftparams%c4      = .true.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='TrE') then
      out_getpftparams%grass   = .false.
      out_getpftparams%tree    = .true.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='TNE') then
      out_getpftparams%grass   = .false.
      out_getpftparams%tree    = .true.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .true.
    else if (trim(pftname)=='TND') then
      out_getpftparams%grass   = .false.
      out_getpftparams%tree    = .true.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .true.
    end if      

    ! land use category associated with PFT (provisional) 
    lu_category_prov = getparreal( trim('params/params_plant_'//trim(pftname)//'.dat'), 'lu_category_prov' )
    if (lu_category_prov==1.0) then
      out_getpftparams%lu_category = lunat
      out_getpftparams%islu(lunat) = .true.
    else
      out_getpftparams%islu(lunat) = .false.
    end if

    ! ! leaf mass per area (gC m-2)
    ! out_getpftparams%lma = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'lma' )
    ! out_getpftparams%sla = 1.0 / out_getpftparams%lma

    ! ! constant ratio of leaf structural N to LMA
    ! out_getpftparams%r_ntolma = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_ntolma' )

    ! leaf decay constant, read in as [years-1], central value: 0.0 yr-1 for deciduous plants
    out_getpftparams%k_decay_leaf_base = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_leaf_base' ) / ndayyear 

    ! shape parameter for turnover function if LAI
    out_getpftparams%k_decay_leaf_width = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_leaf_width' )

    ! sapwood decay constant [days], read in as [years-1], central value: xxx
    out_getpftparams%k_decay_sapw = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_sapw' ) / ndayyear 

    ! root decay constant [days], read in as [years-1], central value: 1.04 (Shan et al., 1993; see Li et al., 2014)
    out_getpftparams%k_decay_root = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_root' ) / ndayyear 

    ! root decay constant [days], read in as [years-1]
    out_getpftparams%k_decay_labl = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_labl' ) / ndayyear 

    ! root C:N and N:C ratio (gC/gN and gN/gC)
    out_getpftparams%r_cton_root = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_cton_root' )
    out_getpftparams%r_ntoc_root = 1.0 / out_getpftparams%r_cton_root

    ! y-axis intersection in the relationship of non-metabolic versus metabolic N per leaf area
    out_getpftparams%ncw_min = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'ncw_min' )

    ! slope in the relationship of non-metabolic versus metabolic N per leaf area
    out_getpftparams%r_n_cw_v = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_n_cw_v' )

    ! constant ratio of C to structural N
    out_getpftparams%r_ctostructn_leaf = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_ctostructn_leaf' )

  end function getpftparams


  subroutine initglobal_plant( plant, ngridcells )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of all _pools on all gridcells at the beginning
    !  of the simulation.
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_params_core, only: npft

    ! argument
    type( plant_type ), dimension(npft,ngridcells), intent(inout) :: plant
    integer, intent(in) :: ngridcells

    ! local variables
    integer :: pft
    integer :: jpngr

    !-----------------------------------------------------------------------------
    ! derive which PFTs are present from fpc_grid (which is prescribed)
    !-----------------------------------------------------------------------------
    do jpngr=1,ngridcells
      do pft=1,npft
        call initpft( plant(pft,jpngr) )
        plant(pft,jpngr)%pftno = pft
      end do
    end do

  end subroutine initglobal_plant


  subroutine initpft( plant )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of specified PFT on specified gridcell
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    ! argument
    type( plant_type ), intent(inout) :: plant

    plant%fpc_grid  = 0.0
    plant%lai_ind   = 0.0
    plant%fapar_ind = 0.0
    plant%acrown    = 0.0

    ! canpopy state variables
    plant%narea            = 0.0
    plant%narea_metabolic  = 0.0
    plant%narea_structural = 0.0
    plant%lma              = 0.0
    plant%sla              = 0.0
    plant%nmass            = 0.0
    plant%r_cton_leaf      = 0.0
    plant%r_ntoc_leaf      = 0.0

  end subroutine initpft


  subroutine initdaily_plant( plant_fluxes )
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    !----------------------------------------------------------------
    ! arguments
    type( plant_fluxes_type ), dimension(npft), intent(inout) :: plant_fluxes

    plant_fluxes(:)%dgpp = 0.0
    plant_fluxes(:)%dnpp = carbon(0.0)
    plant_fluxes(:)%drd  = 0.0

  end subroutine initdaily_plant


  subroutine initoutput_plant( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_interface, only: interface

    ! arguments
    integer, intent(in) :: ngridcells
    
    ! annual output variables
    if (interface%params_siml%loutplant) then

      if (interface%steering%init) then
        allocate( outanarea_mb(npft,ngridcells) )
        allocate( outanarea_cw(npft,ngridcells) )
        allocate( outalai(npft,ngridcells) )
        allocate( outalma(npft,ngridcells) )
        allocate( outacton_lm(npft,ngridcells) )
      end if
      
      outanarea_mb(:,:) = 0.0
      outanarea_cw(:,:) = 0.0
      outalai     (:,:) = 0.0
      outalma     (:,:) = 0.0
      outacton_lm (:,:) = 0.0

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


    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------
    if (interface%params_siml%loutplant) then

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

    ! local variables
    integer :: pft
    integer :: it

    !----------------------------------------------------------------
    ! DAILY FOR HIGH FREQUENCY OUTPUT
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    ! it = floor( real( doy - 1 ) / real( interface%params_siml%outdt ) ) + 1

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    ! if (interface%params_siml%loutplant) then
    !   ! nothing yet
    ! end if

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
    integer :: it, moy, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    ! if ( .not. interface%steering%spinup &
    !      .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
    !      .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

    !   ! Write daily output only during transient simulation
    !   do it=1,interface%params_siml%outnt

    !     ! Define 'itime' as a decimal number corresponding to day in the year + year
    !     itime = real(interface%steering%outyear) + real( it - 1 ) * interface%params_siml%outdt / real( ndayyear )
        
    !   end do
    ! end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutplant) then

      itime = real(interface%steering%outyear)

      write(319,999) itime, sum(outanarea_mb(:,jpngr))
      write(320,999) itime, sum(outanarea_cw(:,jpngr))
      write(321,999) itime, sum(outacton_lm(:,jpngr))
      write(322,999) itime, sum(outalma(:,jpngr))

    end if

    return

    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_plant

end module md_plant
