module md_params_modl
  !////////////////////////////////////////////////////////////////
  !  Module contains (constrainable) model parameters.
  !  Model parameters adopted here are from LPX C3 grass PFT
  !  Litter and soil turnover parameters are divided by 365 to 
  !  convert from [1/yr] to [1/d].
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core

  implicit none

  private
  public getpar_modl, params_glob, params_pft

  ! NON PFT-DEPENDENT PARAMETERS
  type paramstype
    real :: kbeer             ! canopy light extinction coefficient
    real :: f_nretain         ! fraction of N retained at leaf abscission 
    real :: fpc_tree_max      ! maximum fractional plant coverage of trees
    real :: growtheff         ! growth respiration coefficient = yield factor [unitless]
  end type paramstype

  type( paramstype ) :: params_glob

  ! PFT-DEPENDENT PARAMETERS
  type pftparamstype
    integer :: lu_category         ! land use category associated with PFT
    integer :: pftcode             ! code for identification of PFT (1=C3 grass)
    logical, dimension(nlu) :: islu! islu(ipft,ilu) is true if ipft belongs to ilu
    logical :: grass               ! boolean for growth form 'grass'
    logical :: tree                ! boolean for growth form 'tree'
    logical :: nfixer              ! whether plant is capable of symbiotically fixing N
    logical :: c3grass             ! whether plant follows C3 photosynthesis
    logical :: c4grass             ! whether plant follows C4 photosynthesis
    real    :: k_decay_leaf        ! leaf decay constant [year-1]
    real    :: k_decay_sapw        ! sapwood decay constant [year-1]
    real    :: k_decay_root        ! root decay constant [year-1]
    real    :: cton_soil           ! C:N ratio of soil organic matter [gC/gN]
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
  end type pftparamstype

  type( pftparamstype ), dimension(npft) :: params_pft

contains

  subroutine getpar_modl()
    !////////////////////////////////////////////////////////////////
    !  Subroutine reads model parameters from input file.
    !  It was necessary to separate this SR from module md_params_modl
    !  because this SR uses module md_waterbal, which also uses
    !  _params_modl.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal
    use md_params_site, only: lTeBS, lGrC3, lGrC4

    ! local variables
    integer :: pft
    integer :: npft_site

    !----------------------------------------------------------------
    ! NON-PFT DEPENDENT PARAMETERS
    !----------------------------------------------------------------
    ! canopy light extinction coefficient for Beer's Law
    params_glob%kbeer = getparreal( 'params/params_modl.dat', 'kbeer' )

    ! fraction of N retained at leaf abscission 
    params_glob%f_nretain = getparreal( 'params/params_modl.dat', 'f_nretain' )
    
    ! maximum fractional plant coverage of trees (sum of all tree PFTs)
    params_glob%fpc_tree_max = getparreal( 'params/params_modl.dat', 'fpc_tree_max' )

    ! growth efficiency = yield factor, central value: 0.6, range: 0.5-0.7; Zhang et al. (2009), see Li et al., 2014
    params_glob%growtheff = getparreal( 'params/params_modl.dat', 'growtheff' )


    !----------------------------------------------------------------
    ! PFT DEPENDENT PARAMETERS
    ! read parameter input file and store values in single array
    ! important: Keep this order of reading PFT parameters fixed.
    !----------------------------------------------------------------
    pft = 0
    if ( lTeBS ) then
      pft = pft + 1
      params_pft(pft) = getpftparams( 'TeBS' )
    end if 
    if ( lGrC3 ) then
      pft = pft + 1
      params_pft(pft) = getpftparams( 'GrC3' )
    end if
    if ( lGrC4 ) then
      pft = pft + 1
      params_pft(pft) = getpftparams( 'GrC4' )
    end if
    npft_site = pft

  end subroutine getpar_modl


  function getpftparams( pftname ) result( out_getpftparams )
    !----------------------------------------------------------------
    ! Read PFT parameters from respective file, given the PFT name
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    character(len=*) :: pftname

    ! local variables
    real :: lu_category_prov    ! land use category associated with PFT (provisional)
    real :: code_growthform
    real :: code_nfixer
    real :: pftcode

    ! function return variable
    type( pftparamstype ) :: out_getpftparams

    ! leaf decay constant, read in as [years-1], central value: 0.0 yr-1 for deciduous plants
    out_getpftparams%k_decay_leaf = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'k_decay_leaf' ) / ndayyear 

    ! sapwood decay constant [days], read in as [years-1], central value: xxx
    out_getpftparams%k_decay_sapw = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'k_decay_sapw' ) / ndayyear 

    ! root decay constant [days], read in as [years-1], central value: 1.04 (Shan et al., 1993; see Li et al., 2014)
    out_getpftparams%k_decay_root = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'k_decay_root' ) / ndayyear 

    ! C:N ratio of soil organic matter [1]
    out_getpftparams%cton_soil = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'cton_soil' )

    ! land use category associated with PFT (provisional)
    lu_category_prov = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'lu_category_prov' )
    if (lu_category_prov==1.0) then
      out_getpftparams%lu_category = lunat
      out_getpftparams%islu(lunat) = .true.
    else
      out_getpftparams%islu(lunat) = .false.
    end if

    ! PFT identification code (same as in LPX)
    ! 3 = Temperate Needle Evergreen        TeNE
    ! 5 = Temperate Broadleaved Summergreen TeBS
    ! 8 = C3 grass                          TeH
    ! 9 = C4 grass                          TrH
    out_getpftparams%c3grass = .false.
    out_getpftparams%c4grass = .false.
    pftcode = int( getparreal( trim('params/params_pft_'//pftname//'.dat'), 'pftcode' ) ) 
    if (pftcode==8.0) then
      out_getpftparams%c3grass = .true.
    else if (pftcode==9.0) then
      out_getpftparams%c4grass = .true.
    end if      

    ! leaf type: broadleaved (1), needleleaved (2), grass (3) or moss (4)
    code_growthform = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'code_growthform' )
    if ( code_growthform==3.0) then
      out_getpftparams%grass = .true.
      out_getpftparams%tree  = .false.
    else
      out_getpftparams%tree = .true.
      out_getpftparams%grass = .false.
    end if     

    ! N-fixing plant? (0=false, 1=true)
    code_nfixer = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'code_nfixer' )
    if (code_nfixer==0.0) then
      out_getpftparams%nfixer = .false.
    else
      out_getpftparams%nfixer = .true.
    endif 

    ! Reinicker-p for geometry (not PFT-dependent in LPX)
    out_getpftparams%reinickerp = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'reinickerp' ) ! put into pft-loop to have all geometry parameters in one place 

    ! wood density (gC * m-3)
    out_getpftparams%wooddens = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'wooddens' )

    ! allometry parameters
    out_getpftparams%allom1 = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'allom1' )
    out_getpftparams%allom2 = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'allom2' )
    out_getpftparams%allom3 = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'allom3' )      ! put into pft-loop to have all geometry parameters in one place 

    ! maximum crown area
    out_getpftparams%crownarea_max = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'crownarea_max' ) ! put into pft-loop to have all geometry parameters in one place 

    ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
    out_getpftparams%latosa = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'latosa' )

    ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    out_getpftparams%woodtosapw = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'woodtosapw' )

    ! leaf to root ratio under non-water stressed conditionss
    out_getpftparams%lmtorm0 = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'lmtorm0' )

    ! C:N in root biomass
    out_getpftparams%r_cton_root = getparreal( trim('params/params_pft_'//pftname//'.dat'), 'r_cton_root' )
    out_getpftparams%r_ntoc_root = 1.0 / out_getpftparams%r_cton_root

  end function getpftparams

end module md_params_modl
