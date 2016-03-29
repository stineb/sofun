subroutine getpar_modl
  !////////////////////////////////////////////////////////////////
  



  !             XXXXXXXXX DISCONTINUED XXXXXXXXX




  !  Subroutine reads model parameters from input file.
  !  It was necessary to separate this SR from module md_params_modl
  !  because this SR uses module md_waterbal, which also uses
  !  _params_modl.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_modl
  use md_waterbal, only: getpar_modl_waterbal
  use md_nuptake, only: getpar_modl_nuptake
  use md_npp, only: getpar_modl_npp
  use md_phenology, only: getpar_phenology
  use md_ntransform, only: getpar_modl_ntransform
  use md_littersom, only: getpar_modl_littersom
  use md_sofunutils, only: getparreal

  implicit none

  ! local variables
  integer :: pft

  integer, parameter :: npar_pft = 18
  real, dimension(npar_pft) :: pft_params_array

  real, dimension(npft) :: lu_category_prov    ! land use category associated with PFT (provisional)

  !----------------------------------------------------------------
  ! NON-PFT DEPENDENT PARAMETERS
  !----------------------------------------------------------------
  ! canopy light extinction coefficient for Beer's Law
  kbeer = getparreal( 'params_modl.dat', 'kbeer' )

  ! fraction of N retained at leaf abscission 
  F_NRETAIN = getparreal( 'params_modl.dat', 'F_NRETAIN' )
  
  ! maximum fractional plant coverage of trees (sum of all tree PFTs)
  fpc_tree_max = getparreal( 'params_modl.dat', 'fpc_tree_max' )

  ! growth efficiency = yield factor, central value: 0.6, range: 0.5-0.7; Zhang et al. (2009), see Li et al., 2014
  growtheff = getparreal( 'params_modl.dat', 'growtheff' )

  
  !----------------------------------------------------------------
  ! MODULE-SPECIFIC PARAMETERS
  !----------------------------------------------------------------
  call getpar_modl_waterbal()
  call getpar_modl_npp()
  call getpar_modl_nuptake()
  call getpar_modl_ntransform()
  call getpar_modl_littersom()

  !----------------------------------------------------------------
  ! PFT DEPENDENT PARAMETERS
  ! read parameter input file and store values in single array
  !----------------------------------------------------------------
  open(unit=02,file='params_pft.dat',status='OLD')      
  read(02,*) pft_params_array
  close(02)

  do pft=1,npft 

    ! leaf decay constant, read in as [years-1], central value: 0.0 yr-1 for deciduous plants
    k_decay_leaf(pft) = pft_params_array(1) / ndayyear 

    ! sapwood decay constant [days], read in as [years-1], central value: xxx
    k_decay_sapw(pft) = pft_params_array(2) / ndayyear 

    ! root decay constant [days], read in as [years-1], central value: 1.04 (Shan et al., 1993; see Li et al., 2014)
    k_decay_root(pft) = pft_params_array(3) / ndayyear 

    ! C:N ratio of soil organic matter [1]
    cton_soil(pft) = pft_params_array(4) 

    ! land use category associated with PFT (provisional)
    lu_category_prov(pft) = pft_params_array(5)
    if (lu_category_prov(pft)==1.0) then
      lu_category(pft) = lunat
      islu(pft,lunat) = .true.
    else
      islu(pft,lunat) = .false.
    end if

    ! PFT identification code (same as in LPX)
    ! 3 = Temperate Needle Evergreen        TeNE
    ! 5 = Temperate Broadleaved Summergreen TeBS
    ! 8 = C3 grass                          TeH
    pftcode(pft) = int(pft_params_array(6)) 

    ! leaf type: broadleaved (1), needleleaved (2), grass (3) or moss (4)
    if (pft_params_array(7)==3.0) then
      grass(pft) = .true.
      tree(pft)  = .false.
    else
      tree(pft) = .true.
      grass(pft) = .false.
    end if     

    ! N-fixing plant? (0=false, 1=true)
    if (pft_params_array(8)==0.0) then
      nfixer(pft) = .false.
    else
      nfixer(pft) = .true.
    endif 

    ! Reinicker-p for geometry (not PFT-dependent in LPX)
    reinickerp = pft_params_array(9) ! put into pft-loop to have all geometry parameters in one place 

    ! wood density (gC * m-3)
    wooddens(pft) = pft_params_array(10) 

    ! allometry parameters
    allom1(pft) = pft_params_array(11)
    allom2(pft) = pft_params_array(12)
    allom3 = pft_params_array(13)      ! put into pft-loop to have all geometry parameters in one place 

    ! maximum crown area
    crownarea_max(pft) = pft_params_array(14) ! put into pft-loop to have all geometry parameters in one place 

    ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
    latosa(pft) = pft_params_array(15)

    ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    woodtosapw(pft) = pft_params_array(16)

    ! leaf to root ratio under non-water stressed conditionss
    lmtorm0(pft) = pft_params_array(17)

    ! C:N in root biomass
    r_cton_root(pft) = pft_params_array(18)
    r_ntoc_root(pft) = 1.0 / r_cton_root(pft)

  end do

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC PARAMETERS
  !----------------------------------------------------------------
  call getpar_phenology()

  return
  888 write(0,*) 'GETPAR_MODL: error opening file params_modl.dat. Abort. '
  stop


end subroutine getpar_modl
