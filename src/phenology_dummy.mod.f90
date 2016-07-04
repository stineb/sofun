module md_phenology
  !////////////////////////////////////////////////////////////////
  ! TEMPERATURE-DRIVEN PHENOLOGY 
  ! Not used (yet). Therefore, vegetation is always active (not 
  ! phenology-limited).
  ! Contains the "main" subroutine 'gettempphenology and phenology' and all 
  ! necessary subroutines for handling input/output. 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------  
  use md_params_core, only: ndayyear, npft

  implicit none

  private

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  real, dimension(ndayyear,npft)    :: dtphen       ! daily temperature-driven phenology (=dphen_t in LPX)
  logical, dimension(ndayyear,npft) :: sprout       ! boolean whether PFT is present
  logical, dimension(ndayyear,npft) :: shedleaves   ! boolean whether PFT is present

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! (none)

  !----------------------------------------------------------------
  ! Parameters
  !----------------------------------------------------------------
  ! (none)

contains

  subroutine gettempphenology( jpngr, dtemp )
    !//////////////////////////////////////////////////////////
    ! Defines dtphen, the temperature-driven phenology
    !----------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    real, dimension(ndayyear), intent(in) :: dtemp

    dtphen(:,:)     = 1.0
    sprout(:,:)     = .false.
    shedleaves(:,:) = .false.

  end subroutine gettempphenology


  subroutine getpar_modl_phenology()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------

    ! (empty)
  
  end subroutine getpar_modl_phenology


end module md_phenology





