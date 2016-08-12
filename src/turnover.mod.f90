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
  use md_plant
    
  implicit none

  private
  public turnover, turnover_leaf, turnover_root

contains

  subroutine turnover( jpngr, doy )
    !////////////////////////////////////////////////////////////////
    !  Annual vegetation biomass turnover, called at the end of the
    !  year.
    !----------------------------------------------------------------
    use md_classdefs
    use md_params_core, only: npft

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu
    real :: dleaf
    real :: droot

    ! xxx verbose
    logical, parameter :: verbose = .false.
    type( orgpool ) :: orgtmp, orgtmp2

    do pft=1,npft

      if (plabl(pft,jpngr)%c%c12<0.0) stop 'before turnover labile C is neg.'
      if (plabl(pft,jpngr)%n%n14<0.0) stop 'before turnover labile N is neg.'

      ! if (ispresent(pft,jpngr)) then
        !--------------------------------------------------------------
        ! Get turnover fractions
        ! Turnover-rates are reciprocals of tissue longevity
        ! dleaf=1.0/long_leaf(pft)
        ! assuming no continuous leaf turnover
        !--------------------------------------------------------------
        if (params_pft_plant(pft)%grass) then

          ! grasses have continuous turnover
          ! dleaf = 2.5 / 365.0
          droot = params_pft_plant(pft)%k_decay_root

          ! Alternative turnover function: increase turnover rate towards high LAI
          dleaf = (lai_ind(pft,jpngr)*params_pft_plant(pft)%k_decay_leaf_width)**8 + params_pft_plant(pft)%k_decay_leaf_base

        else

          stop 'turnover not implemented for non-grasses'

        endif

        !--------------------------------------------------------------
        ! Calculate leaf turnover in this day 
        !--------------------------------------------------------------
        if (verbose) print*, 'calling turnover_leaf() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              plitt = ', plitt_af(pft,jpngr)
        if (verbose) orgtmp  =  pleaf(pft,jpngr)
        if (verbose) orgtmp2 =  plitt_af(pft,jpngr)
        !--------------------------------------------------------------
        if ( dleaf>0.0 )                 call turnover_leaf( dleaf, pft, jpngr )
        !--------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              plitt = ', plitt_af(pft,jpngr)
        if (verbose) print*, '              --- balance: '
        if (verbose) print*, '                  dlitt - dleaf                = ',  orgminus( &
                                                                                      orgminus( &
                                                                                        plitt_af(pft,jpngr), &
                                                                                        orgtmp2 &
                                                                                        ), &
                                                                                      orgminus( &
                                                                                        orgtmp, &
                                                                                        pleaf(pft,jpngr) &
                                                                                        ) &
                                                                                      )

        !--------------------------------------------------------------
        ! Calculate root turnover in this day 
        !--------------------------------------------------------------
        if (verbose) print*, 'calling turnover_root() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              pleaf = ', proot(:,jpngr)
        if (verbose) print*, '              plitt = ', plitt_bg(pft,jpngr)
        if (verbose) orgtmp  =  proot(pft,jpngr)
        if (verbose) orgtmp2 =  plitt_bg(pft,jpngr)
        !--------------------------------------------------------------
        if ( droot>0.0 )                 call turnover_root( droot, pft, jpngr )
        !--------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              proot = ', proot(:,jpngr)
        if (verbose) print*, '              plitt = ', plitt_bg(pft,jpngr)
        if (verbose) print*, '              --- balance: '
        if (verbose) print*, '                  dlitt - droot                = ',  orgminus( &
                                                                                      orgminus( &
                                                                                        plitt_bg(pft,jpngr), &
                                                                                        orgtmp2 &
                                                                                        ), &
                                                                                      orgminus( &
                                                                                        orgtmp, &
                                                                                        proot(pft,jpngr) &
                                                                                        ) &
                                                                                      )

      ! endif                   !present
    enddo                     !pft

  end subroutine turnover


  subroutine turnover_leaf( dleaf, pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction dleaf for leaf pool
    !------------------------------------------------------------------
    ! arguments
    real, intent(in)    :: dleaf
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    ! local variables
    type(orgpool) :: lm_turn

    ! print*,'dleaf ', dleaf
    ! stop

    ! determine absolute turnover
    lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) ) ! leaf turnover

    ! reduce leaf mass and root mass
    call orgsub( lm_turn, pleaf(pft,jpngr) )

    ! add all organic (fixed) C to litter
    ! call cmvRec( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))
    call cmv( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, scale=nind(pft,jpngr) )

    ! retain fraction of N
    call nmv( nfrac( params_plant%f_nretain, lm_turn%n ), lm_turn%n, plabl(pft,jpngr)%n )

    ! rest goes to litter
    ! call nmvRec( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    call nmv( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, scale=nind(pft,jpngr) )

  end subroutine turnover_leaf


  subroutine turnover_root( droot, pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction droot for root pool
    !------------------------------------------------------------------
    ! arguments
    real, intent(in)    :: droot
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    ! local variables
    type(orgpool) :: rm_turn

    ! determine absolute turnover
    rm_turn = orgfrac( droot, proot(pft,jpngr) ) ! root turnover

    ! reduce leaf mass and root mass
    call orgsub( rm_turn, proot(pft,jpngr) )

    ! add all organic (fixed) C to litter
    ! call cmvRec( rm_turn%c, rm_turn%c, plitt_bg(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))
    call cmv( rm_turn%c, rm_turn%c, plitt_bg(pft,jpngr)%c, scale=nind(pft,jpngr))

    ! retain fraction of N
    call nmv( nfrac( params_plant%f_nretain, rm_turn%n ), rm_turn%n, plabl(pft,jpngr)%n )

    ! rest goes to litter
    ! call nmvRec( rm_turn%n, rm_turn%n, plitt_bg(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    call nmv( rm_turn%n, rm_turn%n, plitt_bg(pft,jpngr)%n, scale=nind(pft,jpngr) )

  end subroutine turnover_root


end module md_turnover
