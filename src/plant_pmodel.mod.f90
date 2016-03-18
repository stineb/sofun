module _plant
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

  private
  public getpar_modl_plant, params_pft_plant, fapar_ind

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  real, dimension(npft,maxgrid) :: fapar_ind

  !-----------------------------------------------------------------------
  ! Parameters. Runtime read-in
  !-----------------------------------------------------------------------
  ! type params_plant_type
  ! end type params_plant_type

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
    !  It was necessary to separate this SR from module _plant
    !  because this SR uses module _waterbal, which also uses
    !  _plant.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------    
    use _sofunutils, only: getparreal
    use _params_site, only: lTeBS, lGrC3, lGrC4

    ! local variables
    integer :: pft
    integer :: npft_site

    !----------------------------------------------------------------
    ! PFT DEPENDENT PARAMETERS
    ! read parameter input file and store values in single array
    ! important: Keep this order of reading PFT parameters fixed.
    !----------------------------------------------------------------
    pft = 0
    if ( lTeBS ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TeBS' )
    end if 
    if ( lGrC3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GrC3' )
    end if
    if ( lGrC4 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GrC4' )
    end if
    npft_site = pft

  end subroutine getpar_modl_plant


  function getpftparams( pftname ) result( out_getpftparams )
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
    type( params_pft_plant_type ) :: out_getpftparams

    ! standard PFT name
    out_getpftparams%pftname = pftname

    ! PFT names
    ! GrC3 : C3 grass                          
    ! GrC4 : C4 grass     
    if (trim(pftname)=='GrC3') then
      out_getpftparams%grass = .true.
      out_getpftparams%tree  = .false.
      out_getpftparams%c3grass = .true.
      out_getpftparams%c4grass = .false.
      out_getpftparams%nfixer = .false.
    else if (trim(pftname)=='GrC4') then
      out_getpftparams%grass = .true.
      out_getpftparams%tree  = .false.
      out_getpftparams%c3grass = .false.
      out_getpftparams%c4grass = .true.
      out_getpftparams%nfixer = .false.
    end if      

    ! land use category associated with PFT (provisional)
    lu_category_prov = getparreal( trim('params/params_plant_'//trim(pftname)//'.dat'), 'lu_category_prov' )
    if (lu_category_prov==1.0) then
      out_getpftparams%lu_category = lunat
      out_getpftparams%islu(lunat) = .true.
    else
      out_getpftparams%islu(lunat) = .false.
    end if


  end function getpftparams

end module _plant
