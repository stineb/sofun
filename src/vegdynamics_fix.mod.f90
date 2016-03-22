module _vegdynamics
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk

  use _params_core

  implicit none

  private
  public vegdynamics

contains

  subroutine vegdynamics( jpngr, doy )
    !//////////////////////////////////////////////////////////////////
    ! Updates canopy and stand variables and calls 'estab_daily' to 
    ! simulate establishment of new individuals
    !------------------------------------------------------------------
    use _params_core, only: npft
    use _phenology, only: dtphen, sprout, params_pft_pheno

    ! xxx debug
    use _plant, only: ispresent

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    integer :: pft


    do pft=1,npft

      ! print*, 'xxx try: do nothing in vegdynamics'

      ! if (params_pft_pheno(pft)%summergreen) then
      !   !----------------------------------------------------------
      !   ! GRASSES, summergreen
      !   !----------------------------------------------------------

      !   if ( sprout(doy,pft) ) then
      !     !----------------------------------------------------------
      !     ! beginning of season
      !     !----------------------------------------------------------
      !     print*, 'starting to grow on day ',doy
      !     call estab_daily( pft, jpngr, doy )

      ! else

      !   stop 'estab_daily not implemented for non-summergreen'

      ! end if

      ! ! ! Update canopy state variables
      ! ! canopy = get_canopy( lai_ind(pft,jpngr) )

    end do

  end subroutine vegdynamics


  subroutine estab_daily( pft, jpngr, doy )
    !//////////////////////////////////////////////////////////////////
    ! Calculates leaf-level metabolic N content per unit leaf area as a
    ! function of Vcmax25.
    !------------------------------------------------------------------
    use _plant, only: initpft, params_pft_plant, ispresent, nind

    ! arguments
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! initialise all pools of this PFT with zero
    call initpft( pft, jpngr )

    ! add C (and N) to labile pool (available for allocation)
    call add_seed( pft, jpngr )
    
    ! set other state variables: 'ispresent' and 'nind'
    ispresent(pft,jpngr) = .true.

    if (params_pft_plant(pft)%grass) then
      nind(pft,jpngr) = 1.0
    else
      stop 'estab_daily not implemented for trees'
    end if


  end subroutine estab_daily


  subroutine add_seed( pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! To initialise plant pools, add "sapling" mass
    !------------------------------------------------------------------
    use _classdefs
    use _plant, only: plabl, seed, dnpp

    ! arguments
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    plabl(pft,jpngr) = seed

  end subroutine add_seed


end module _vegdynamics
