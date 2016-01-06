module _params_modl
  !////////////////////////////////////////////////////////////////
  !  Module contains (constrainable) model parameters.
  !  Model parameters adopted here are from LPX C3 grass PFT
  !  Litter and soil turnover parameters are divided by 365 to 
  !  convert from [1/yr] to [1/d].
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core

  implicit none

  ! NON PFT-DEPENDENT PARAMETERS
  real :: kbeer             ! canopy light extinction coefficient
  real :: F_NRETAIN         ! fraction of N retained at leaf abscission 
  real :: fpc_tree_max      ! maximum fractional plant coverage of trees
  real :: growtheff         ! growth respiration coefficient = yield factor [unitless]

  ! PFT-DEPENDENT PARAMETERS
  real, dimension(npft)        :: k_decay_leaf        ! leaf decay constant [year-1]
  real, dimension(npft)        :: k_decay_sapw        ! sapwood decay constant [year-1]
  real, dimension(npft)        :: k_decay_root        ! root decay constant [year-1]
  real, dimension(npft)        :: cton_soil           ! C:N ratio of soil organic matter [gC/gN]

  integer, dimension(npft)     :: lu_category      ! land use category associated with PFT
  integer, dimension(npft)     :: pftcode          ! code for identification of PFT (1=C3 grass)
  logical, dimension(npft,nlu) :: islu             ! islu(ipft,ilu) is true if ipft belongs to ilu
  logical, dimension(npft)     :: grass            ! boolean for growth form 'grass'
  logical, dimension(npft)     :: tree             ! boolean for growth form 'tree'
  logical, dimension(npft)     :: nfixer              ! whether plant is capable of symbiotically fixing N

  real, dimension(npft)        :: reinickerp          ! Reinicker-p for geometry
  real, dimension(npft)        :: wooddens            ! wood density (gC * m-3)
  real, dimension(npft)        :: allom1, allom2      ! allometry parameters
  real, dimension(npft)        :: allom3              ! allometry parameter
  real, dimension(npft)        :: crownarea_max       ! maximum crown area
  real, dimension(npft)        :: latosa              ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
  real, dimension(npft)        :: woodtosapw          ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
  real, dimension(npft)        :: lmtorm0             ! leaf to root ratio under non-water stressed conditions

  real, dimension(npft)        :: r_cton_root         ! C:N ratio in roots
  real, dimension(npft)        :: r_ntoc_root         ! N:C ratio in roots (inverse of 'r_cton_root')

  type outtype_getpftpar
    real    :: k_decay_leaf        ! leaf decay constant [year-1]
    real    :: k_decay_sapw        ! sapwood decay constant [year-1]
    real    :: k_decay_root        ! root decay constant [year-1]
    real    :: cton_soil           ! C:N ratio of soil organic matter [gC/gN]
    integer :: lu_category         ! land use category associated with PFT
    integer :: pftcode             ! code for identification of PFT (1=C3 grass)
    logical, dimension(nlu) :: islu! islu(ipft,ilu) is true if ipft belongs to ilu
    logical :: grass               ! boolean for growth form 'grass'
    logical :: tree                ! boolean for growth form 'tree'
    logical :: nfixer              ! whether plant is capable of symbiotically fixing N
    real    :: reinickerp          ! Reinicker-p for geometry
    real    :: wooddens            ! wood density (gC * m-3)
    real    :: allom1, allom2      ! allometry parameters
    real    :: allom3              ! allometry parameter
    real    :: crownarea_max       ! maximum crown area
    real    :: latosa              ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
    real    :: woodtosapw          ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    real    :: lmtorm0             ! leaf to root ratio under non-water stressed conditions
    real    :: r_cton_root         ! C:N ratio in roots
    real    :: r_ntoc_root         ! N:C ratio in roots (inverse of 'r_cton_root')
  end type outtype_getpftpar

contains

  subroutine getpar_modl()
    !////////////////////////////////////////////////////////////////
    !  Subroutine reads model parameters from input file.
    !  It was necessary to separate this SR from module _params_modl
    !  because this SR uses module _waterbal, which also uses
    !  _params_modl.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use _sofunutils, only: getparreal
    use _params_site, only: lTeBS, lGrC3

    ! local variables
    integer :: pft
    integer :: npft_site
    type(outtype_getpftpar), dimension(npft) :: out_getpftpar

    !----------------------------------------------------------------
    ! NON-PFT DEPENDENT PARAMETERS
    !----------------------------------------------------------------
    ! canopy light extinction coefficient for Beer's Law
    kbeer = getparreal( 'params/params_modl.dat', 'kbeer' )

    ! fraction of N retained at leaf abscission 
    F_NRETAIN = getparreal( 'params/params_modl.dat', 'F_NRETAIN' )
    
    ! maximum fractional plant coverage of trees (sum of all tree PFTs)
    fpc_tree_max = getparreal( 'params/params_modl.dat', 'fpc_tree_max' )

    ! growth efficiency = yield factor, central value: 0.6, range: 0.5-0.7; Zhang et al. (2009), see Li et al., 2014
    growtheff = getparreal( 'params/params_modl.dat', 'growtheff' )


    !----------------------------------------------------------------
    ! PFT DEPENDENT PARAMETERS
    ! read parameter input file and store values in single array
    ! important: Keep this order of reading PFT parameters fixed.
    !----------------------------------------------------------------
    pft = 0
    if ( lTeBS ) then
      pft = pft + 1
      out_getpftpar(pft) = getpftpar( 'TeBS' )
    end if 
    if ( lGrC3 ) then
      pft = pft + 1
      out_getpftpar(pft) = getpftpar( 'GrC3' )
    end if
    npft_site = pft

    do pft=1,npft_site
      k_decay_leaf(pft)  = out_getpftpar(pft)%k_decay_leaf
      k_decay_sapw(pft)  = out_getpftpar(pft)%k_decay_sapw
      k_decay_root(pft)  = out_getpftpar(pft)%k_decay_root
      cton_soil(pft)     = out_getpftpar(pft)%cton_soil
      lu_category(pft)   = out_getpftpar(pft)%lu_category
      pftcode(pft)       = out_getpftpar(pft)%pftcode
      islu(pft,:)        = out_getpftpar(pft)%islu(:)
      grass(pft)         = out_getpftpar(pft)%grass
      tree(pft)          = out_getpftpar(pft)%tree
      nfixer(pft)        = out_getpftpar(pft)%nfixer
      reinickerp(pft)    = out_getpftpar(pft)%reinickerp
      wooddens(pft)      = out_getpftpar(pft)%wooddens
      allom1(pft)        = out_getpftpar(pft)%allom1
      allom2(pft)        = out_getpftpar(pft)%allom2
      allom3(pft)        = out_getpftpar(pft)%allom3
      crownarea_max(pft) = out_getpftpar(pft)%crownarea_max
      latosa(pft)        = out_getpftpar(pft)%latosa
      woodtosapw(pft)    = out_getpftpar(pft)%woodtosapw
      lmtorm0(pft)       = out_getpftpar(pft)%lmtorm0
      r_cton_root(pft)   = out_getpftpar(pft)%r_cton_root
      r_ntoc_root(pft)   = out_getpftpar(pft)%r_ntoc_root
    end do

  end subroutine getpar_modl


  function getpftpar( pftname ) result( out_getpftpar )
    !----------------------------------------------------------------
    ! Read PFT parameters from respective file, given the PFT name
    !----------------------------------------------------------------
    use _sofunutils, only: getparreal

    ! arguments
    character(len=*) :: pftname

    ! local variables
    real :: lu_category_prov    ! land use category associated with PFT (provisional)
    real :: code_growthform
    real :: code_nfixer

    ! function return variable
    type(outtype_getpftpar) out_getpftpar

    ! leaf decay constant, read in as [years-1], central value: 0.0 yr-1 for deciduous plants
    out_getpftpar%k_decay_leaf = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'k_decay_leaf' ) / ndayyear 

    ! sapwood decay constant [days], read in as [years-1], central value: xxx
    out_getpftpar%k_decay_sapw = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'k_decay_sapw' ) / ndayyear 

    ! root decay constant [days], read in as [years-1], central value: 1.04 (Shan et al., 1993; see Li et al., 2014)
    out_getpftpar%k_decay_root = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'k_decay_root' ) / ndayyear 

    ! C:N ratio of soil organic matter [1]
    out_getpftpar%cton_soil = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'cton_soil' )

    ! land use category associated with PFT (provisional)
    lu_category_prov = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'lu_category_prov' )
    if (lu_category_prov==1.0) then
      out_getpftpar%lu_category = lunat
      out_getpftpar%islu(lunat) = .true.
    else
      out_getpftpar%islu(lunat) = .false.
    end if

    ! PFT identification code (same as in LPX)
    ! 3 = Temperate Needle Evergreen        TeNE
    ! 5 = Temperate Broadleaved Summergreen TeBS
    ! 8 = C3 grass                          TeH
    out_getpftpar%pftcode = int( getparreal( trim('params/params_pft_'//pftname//'.dat'), 'pftcode' ) ) 

    ! leaf type: broadleaved (1), needleleaved (2), grass (3) or moss (4)
    code_growthform = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'code_growthform' )
    if ( code_growthform==3.0) then
      out_getpftpar%grass = .true.
      out_getpftpar%tree  = .false.
    else
      out_getpftpar%tree = .true.
      out_getpftpar%grass = .false.
    end if     

    ! N-fixing plant? (0=false, 1=true)
    code_nfixer = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'code_nfixer' )
    if (code_nfixer==0.0) then
      out_getpftpar%nfixer = .false.
    else
      out_getpftpar%nfixer = .true.
    endif 

    ! Reinicker-p for geometry (not PFT-dependent in LPX)
    out_getpftpar%reinickerp = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'reinickerp' ) ! put into pft-loop to have all geometry parameters in one place 

    ! wood density (gC * m-3)
    out_getpftpar%wooddens = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'wooddens' )

    ! allometry parameters
    out_getpftpar%allom1 = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'allom1' )
    out_getpftpar%allom2 = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'allom2' )
    out_getpftpar%allom3 = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'allom3' )      ! put into pft-loop to have all geometry parameters in one place 

    ! maximum crown area
    out_getpftpar%crownarea_max = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'crownarea_max' ) ! put into pft-loop to have all geometry parameters in one place 

    ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
    out_getpftpar%latosa = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'latosa' )

    ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    out_getpftpar%woodtosapw = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'woodtosapw' )

    ! leaf to root ratio under non-water stressed conditionss
    out_getpftpar%lmtorm0 = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'lmtorm0' )

    ! C:N in root biomass
    out_getpftpar%r_cton_root = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'r_cton_root' )
    out_getpftpar%r_ntoc_root = 1.0 / out_getpftpar%r_cton_root

  end function getpftpar

end module _params_modl
