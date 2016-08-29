module md_turnover
  !////////////////////////////////////////////////////////////////
  ! NPP_LPJ MODULE
  ! Contains the "main" subroutine 'npp' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'npp' must contain this list 
  ! of subroutines (names that way).
  !   - npp
  !   - getpar_modl_npp
  !   - initio_npp
  !   - initoutput_npp
  !   - getout_daily_npp
  !   - getout_monthly_npp
  !   - writeout_ascii_npp
  ! Required module-independent model state variables (necessarily 
  ! updated by 'waterbal') are:
  !   - daily NPP ('dnpp')
  !   - soil temperature ('xxx')
  !   - inorganic N _pools ('no3', 'nh4')
  !   - xxx 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_classdefs
    
  implicit none

  private
  public turnover

contains

  subroutine turnover( plant, tile, jpngr )
    !////////////////////////////////////////////////////////////////
    !  Annual vegetation biomass turnover, called at the end of the
    !  year.
    !----------------------------------------------------------------
    use md_classdefs
    use md_params_core, only: npft, nlu, eps
    use md_plant, only: plant_type, plitt_af, plitt_bg, params_plant, params_pft_plant, &
      outacveg2lit, outanveg2lit
    use md_tile, only: tile_type

    ! arguments
    type( plant_type ), dimension(npft), intent(inout) :: plant
    type( tile_type ), dimension(nlu), intent(in)      :: tile
    integer, intent(in)                                :: jpngr

    ! local variables
    integer :: pft
    integer :: lu
    type( orgpool ) :: dleaf
    type( orgpool ) :: droot
    type( orgpool ) :: dturn
    type( orgpool ) :: dlabl

    do pft=1,npft

      lu = params_pft_plant(pft)%lu_category

      if (plant(pft)%plabl%c%c12 < -1.0*eps) stop 'before turnover labile C is neg.'
      if (plant(pft)%plabl%n%n14 < -1.0*eps) stop 'before turnover labile N is neg.'

      !--------------------------------------------------------------
      ! Get turnover fractions
      ! Turnover-rates are reciprocals of tissue longevity
      ! dleaf=1.0/long_leaf(pft)
      ! assuming no continuous leaf turnover
      !--------------------------------------------------------------
      if (params_pft_plant(pft)%tree) then

        ! determine amount of decaying leaf and root 
        dleaf = orgfrac( params_pft_plant(pft)%k_decay_leaf_base, plant(pft)%pleaf )
        droot = orgfrac( params_pft_plant(pft)%k_decay_root, plant(pft)%proot )

        ! some N is retained from leaves, therefore not all is to be "paid" for again
        dleaf%n = nfrac( ( 1.0 - params_plant%f_nretain ), dleaf%n )

        ! total of what is to be "paid" out of the labile pool
        dturn = orgplus( dleaf, droot )

        ! inflate C "cost" by inverse of growth efficiency
        dlabl%c = cfrac( 1.0 / params_plant%growtheff, dturn%c )
        dlabl%n = dturn%n

        ! reduce labile pool instead of actual plant pool assuming that the decaying plant components
        ! are continuously being renewed
        call orgsub( dlabl, plant(pft)%plabl )

        ! add decaying amounts to litter
        call orgmvRec( droot, droot, plitt_bg(pft,jpngr), outacveg2lit(pft,jpngr), outanveg2lit(pft,jpngr), scale=tile(lu)%nind(pft) )
        call orgmvRec( dleaf, dleaf, plitt_af(pft,jpngr), outacveg2lit(pft,jpngr), outanveg2lit(pft,jpngr), scale=tile(lu)%nind(pft) )

      end if

    enddo                     !pft

  end subroutine turnover


end module md_turnover
