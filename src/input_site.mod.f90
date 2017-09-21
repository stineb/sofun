module input_site
  !////////////////////////////////////////////////////////////////
  ! This module contains subroutines and variables needed when pre-
  ! scribing certain variables instead of simulating them online.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: npft, maxgrid

  implicit none

  real, dimension(npft,maxgrid) :: fpc_grid_data

contains

  subroutine getinput_site
    !////////////////////////////////////////////////////////////////
    !  Subroutine reads to prescribe variables that are otherwise 
    !  (per default) calculated online
    ! This is not appropriate here as its only for one year
    !----------------------------------------------------------------
    use md_params_siml, only: nyeartrend, sitename
    use md_forcing_siterun, only: read1year_daily  ! xxx put all reading functions/subroutines into a single module
    use md_sofunutils, only: getparreal  ! xxx put all reading functions/subroutines into a single module
    use md_params_core, only: ndayyear

    ! local variables
    integer :: day, mo, dm, pft, yr

    ! Get prescribed PFT selection 
    ! xxx try: make npft a definable variable depending on this selection 
    ! Get prescribed fractional plant cover (FPC) for each PFT
    fpc_grid_data(1,1) = getparreal( sitename//".parameter", 'in_fpc_grid_1' )
    fpc_grid_data(2,1) = getparreal( sitename//".parameter", 'in_fpc_grid_2' )
    fpc_grid_data(3,1) = getparreal( sitename//".parameter", 'in_fpc_grid_3' )
    fpc_grid_data(4,1) = getparreal( sitename//".parameter", 'in_fpc_grid_4' )
    fpc_grid_data(5,1) = getparreal( sitename//".parameter", 'in_fpc_grid_5' )
    fpc_grid_data(6,1) = getparreal( sitename//".parameter", 'in_fpc_grid_6' )
    fpc_grid_data(7,1) = getparreal( sitename//".parameter", 'in_fpc_grid_7' )
    fpc_grid_data(8,1) = getparreal( sitename//".parameter", 'in_fpc_grid_8' )
    fpc_grid_data(9,1) = getparreal( sitename//".parameter", 'in_fpc_grid_9' )

    return

  end subroutine getinput_site

end module input_site
