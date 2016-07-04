module md_vegdynamics
  !////////////////////////////////////////////////////////////////
  ! VEGETATION DYNAMICS
  ! Not used (yet): no change in number of individuals. Is always 1.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------  
  implicit none

  private
  public vegdynamics

contains

  subroutine vegdynamics( jpngr, doy )
    !//////////////////////////////////////////////////////////////////
    ! Updates canopy and stand variables and calls 'estab_daily' to 
    ! simulate establishment of new individuals
    !------------------------------------------------------------------

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! (empty)


  end subroutine vegdynamics


end module md_vegdynamics
