subroutine initoutput
  !////////////////////////////////////////////////////////////////
  ! Initialises all daily variables with zero.
  ! Called at the beginning of each year by 'biosphere'.
  !----------------------------------------------------------------
  use _outvars
  use _waterbal, only: initoutput_waterbal
  use _nuptake, only: initoutput_nuptake
  use _ntransform, only: initoutput_ntransform
  use _littersom, only: initoutput_littersom
  
  implicit none

  ! daily output variables
  outdgpp(:,:,:)      = 0.0
  outdnpp(:,:,:)      = 0.0
  outdnup(:,:,:)      = 0.0
  outdCleaf(:,:,:)    = 0.0
  outdninorg(:,:,:)   = 0.0

  ! annual output variables
  outagpp(:,:)        = 0.0
  outanpp(:,:)        = 0.0
  outanup(:,:)        = 0.0
  outaCveg(:,:)       = 0.0
  outaCalloc(:,:)     = 0.0
  outaNalloc(:,:)     = 0.0
  outaCveg2lit(:,:)   = 0.0
  outaNveg2lit(:,:)   = 0.0
  outaninorg(:,:)     = 0.0

  ! initialise module-specific output variables
  call initoutput_waterbal
  call initoutput_nuptake
  call initoutput_ntransform  ! nothing yet
  call initoutput_littersom

end subroutine initoutput
