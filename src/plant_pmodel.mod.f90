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

  implicit none

  private
  public getpar_modl_plant, params_pft_plant, canopy, ispresent

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  ! Canopy state variables
  type canopy_type
    real :: fapar_ind
  end type canopy_type

  type( canopy_type ), dimension(npft) :: canopy


  ! This is not used interactively, therefore always true
  logical, dimension(npft,maxgrid), parameter :: ispresent = .true.

  !-----------------------------------------------------------------------
  ! Parameters. Runtime read-in
  !-----------------------------------------------------------------------
  type params_pft_plant_type
    character(len=4) :: pftname    ! standard PFT name with 4 characters length
    integer :: lu_category         ! land use category associated with PFT
    logical, dimension(nlu) :: islu! islu(ipft,ilu) is true if ipft belongs to ilu
    logical :: grass               ! boolean for growth form 'grass'
    logical :: tree                ! boolean for growth form 'tree'
    logical :: nfixer              ! whether plant is capable of symbiotically fixing N
    logical :: c3grass             ! whether plant follows C3 photosynthesis
    logical :: c4grass             ! whether plant follows C4 photosynthesis
  end type params_pft_plant_type

  type( params_pft_plant_type ), dimension(npft) :: params_pft_plant

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! MODULE-SPECIFIC, PRIVATE VARIABLES
  !----------------------------------------------------------------

contains

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
    ! PFT DEPENDENT PARAMETERS
    ! read parameter input file and store values in single array
    ! important: Keep this order of reading PFT parameters fixed.
    !----------------------------------------------------------------
    pft = 0
    if ( interface%params_siml%lTeBS ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TeBS' )
    end if 
    if ( interface%params_siml%lGrC3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GrC3' )
    end if
    if ( interface%params_siml%lGrC4 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GrC4' )
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
    real :: code_growthform
    real :: code_nfixer

    ! function return variable
    type( params_pft_plant_type ) :: out_getpftparams

    ! standard PFT name
    out_getpftparams%pftname = pftname

    ! PFT names
    ! GrC3 : C3 grass                          
    ! GrC4 : C4 grass     
    if (trim(pftname)=='GrC3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3grass = .true.
      out_getpftparams%c4grass = .false.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='GNC3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3grass = .true.
      out_getpftparams%c4grass = .false.
      out_getpftparams%nfixer  = .true.
    else if (trim(pftname)=='GrC4') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3grass = .false.
      out_getpftparams%c4grass = .true.
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
    out_getpftparams%k_decay_leaf = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_leaf' ) / ndayyear 

    ! sapwood decay constant [days], read in as [years-1], central value: xxx
    out_getpftparams%k_decay_sapw = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_sapw' ) / ndayyear 

    ! root decay constant [days], read in as [years-1], central value: 1.04 (Shan et al., 1993; see Li et al., 2014)
    out_getpftparams%k_decay_root = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'k_decay_root' ) / ndayyear 

    ! root C:N and N:C ratio (gC/gN and gN/gC)
    out_getpftparams%r_cton_root = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_cton_root' )
    out_getpftparams%r_ntoc_root = 1.0 / out_getpftparams%r_cton_root

  end function getpftparams

end module md_plant
