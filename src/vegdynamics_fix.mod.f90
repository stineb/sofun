module md_vegdynamics
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk

  use md_params_core

  implicit none

  private
  public vegdynamics

contains

  subroutine vegdynamics( jpngr, doy )
    !//////////////////////////////////////////////////////////////////
    ! Updates canopy and stand variables and calls 'estab_daily' to 
    ! simulate establishment of new individuals
    !------------------------------------------------------------------
    use md_params_core, only: npft
    use md_phenology, only: sprout, params_pft_pheno
    use md_plant, only: params_pft_plant

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    integer :: pft

    do pft=1,npft

      if (params_pft_plant(pft)%grass) then
        !----------------------------------------------------------
        ! GRASSES, summergreen
        !----------------------------------------------------------

        if ( sprout(doy,pft) ) then
          !----------------------------------------------------------
          ! beginning of season
          !----------------------------------------------------------
          print*, 'starting to grow on day ',doy
          call estab_daily( pft, jpngr, doy )

        end if

      else

        stop 'estab_daily not implemented for non-summergreen'

      end if

    end do

  end subroutine vegdynamics


  subroutine estab_daily( pft, jpngr, doy )
    !//////////////////////////////////////////////////////////////////
    ! Calculates leaf-level metabolic N content per unit leaf area as a
    ! function of Vcmax25.
    !------------------------------------------------------------------
    use md_plant, only: initpft, params_pft_plant, nind, frac_leaf
    use md_interface

    ! arguments
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! ! initialise all pools of this PFT with zero
    ! call initpft( pft, jpngr )

    ! add C (and N) to labile pool (available for allocation)
    call add_seed( pft, jpngr )
    if ( .not. interface%steering%dofree_alloc ) frac_leaf(pft) = 0.5
    
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
    use md_classdefs
    use md_plant, only: plabl, seed, dnpp

    ! arguments
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    plabl(pft,jpngr) = orgplus( plabl(pft,jpngr), seed )

  end subroutine add_seed


end module md_vegdynamics
