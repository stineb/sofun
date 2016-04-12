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
  public turnover

contains

  subroutine turnover( jpngr, doy )
    !////////////////////////////////////////////////////////////////
    !  Annual vegetation biomass turnover, called at the end of the
    !  year.
    !----------------------------------------------------------------
    use md_classdefs
    use md_phenology, only: shedleaves, params_pft_pheno
    use md_params_core, only: npft

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu
    real :: dleaf
    real :: dsapw
    real :: droot
    real :: dlabl

    ! xxx verbose
    logical, parameter :: verbose = .false.
    type( orgpool ) :: orgtmp, orgtmp2

    do pft=1,npft

      if (plabl(pft,jpngr)%c%c12<0.0) stop 'before turnover labile C is neg.'
      if (plabl(pft,jpngr)%n%n14<0.0) stop 'before turnover labile N is neg.'

      if (ispresent(pft,jpngr)) then
    
        ! print*,'to A plabl(pft,jpngr): ', plabl(pft,jpngr)
        ! print*,'params_pft_plant(pft)%grass ', params_pft_plant(pft)%grass
        ! print*,'params_pft_pheno(pft)%summergreen', params_pft_pheno(pft)%summergreen
        ! print*,'shedleaves(doy,pft)', shedleaves(doy,pft)
        !--------------------------------------------------------------
        ! Get turnover fractions
        ! Turnover-rates are reciprocals of tissue longevity
        ! dleaf=1.0/long_leaf(pft)
        ! assuming no continuous leaf turnover
        !--------------------------------------------------------------
        if (params_pft_plant(pft)%grass) then

          ! grasses have continuous turnover
          dleaf = params_pft_plant(pft)%k_decay_leaf
          droot = params_pft_plant(pft)%k_decay_root
          dlabl = 0.0
          dsapw = 0.0

        else

          ! if (params_pft_pheno(pft)%summergreen) then
          !   if (shedleaves(doy,pft)) then
          !     ! print*, 'end of season: labile C, labile N', plabl(pft,jpngr)%c%c12, plabl(pft,jpngr)%n%n14
          !     ! stop 'do beni'
          !     ! dleaf = 1.0
          !     ! droot = 1.0
          !     dlabl = 1.0
          !     dsapw = 1.0
          !     ! print*, 'complete die-off on day', doy
          !   else
          !     dleaf = params_pft_plant(pft)%k_decay_leaf
          !     droot = params_pft_plant(pft)%k_decay_root
          !     dlabl = 0.0
          !     dsapw = 0.0
          !   end if
          ! else
          !   dleaf = params_pft_plant(pft)%k_decay_leaf
          !   droot = params_pft_plant(pft)%k_decay_root
          !   dlabl = 0.0
          !   dsapw = 0.0
          ! end if   

        endif

        ! print*,'dleaf ', dleaf
        ! print*,'droot ', droot
        ! print*,'dsapw ', dsapw
        ! print*,'dlabl ', dlabl

        !--------------------------------------------------------------
        ! Calculate biomass turnover in this year 
        ! xxx todo         
        ! Make consistent with turnover accounted for in allocation!
        !--------------------------------------------------------------
        if (verbose) print*, 'calling turnover_leaf() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              pleaf = ', pleaf(:,jpngr)
        if (verbose) print*, '              plitt = ', plitt_af(pft,jpngr)
        if (verbose) orgtmp  =  pleaf(pft,jpngr)
        if (verbose) orgtmp2 =  plitt_af(pft,jpngr)
        if ( dleaf>0.0 )                 call turnover_leaf( dleaf, pft, jpngr )
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
        if (verbose) print*, '... done'

        if ( droot>0.0 )                 call turnover_root( droot, pft, jpngr )
        if ( params_pft_plant(pft)%tree .and. dsapw>0.0 ) call turnover_sapw( dsapw, pft, jpngr )
        if ( dlabl>0.0 )                 call turnover_labl( dlabl, pft, jpngr )

        ! ! add labile C and N to litter as well
        ! if (dlabl>0.0) then
        !   print*, 'TURNOVER: WARNING LABILE C AND N ARE SET TO ZERO WITHOUT MASS CONSERVATION'
        !   lb_turn%c%c12 = 0.0
        !   lb_turn%n%n14 = 0.0
        !   ! call cmvRec( lb_turn%c, lb_turn%c, plitt_bg(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
        !   ! call nmvRec( lb_turn%n, lb_turn%n, plitt_bg(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
        !   ! print*, 'end of growing season-plabl:', plabl
        !   ! print*, 'moving C and N to exudates and ninorg:', lb_turn
        ! end if

      ! print*,'to B plabl(pft,jpngr): ', plabl(pft,jpngr)
      
        if (plabl(pft,jpngr)%c%c12<0.0) stop 'after turnover labile C is neg.'
        if (plabl(pft,jpngr)%n%n14<0.0) stop 'after turnover labile N is neg.'

      endif                   !present
    enddo                     !pft

  end subroutine turnover


  subroutine turnover_leaf( dleaf, pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction dleaf for leaf pool
    !------------------------------------------------------------------
    ! use md_plant, only: update_foliage_vars

    ! arguments
    real, intent(in)    :: dleaf
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    ! local variables
    type(orgpool) :: lm_turn

    ! xxx verbose
    logical, parameter :: verbose = .false.
    type( orgpool ) :: orgtmp, orgtmp2

    ! new version: >>>>>>>
    ! determine absolute turnover
    lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) ) ! leaf turnover

    ! ! xxx test
    ! orgtmp  = pleaf(pft,jpngr)
    ! orgtmp2 = plitt_af(pft,jpngr)

    ! reduce leaf mass and root mass
    call orgmv( lm_turn, pleaf(pft,jpngr), plitt_af(pft,jpngr) )
    ! print*, 'balance in turnover_leaf:'
    ! print*, 'lm_turn ', lm_turn
    ! print*, 'dleaf   ', orgminus(orgtmp, pleaf(pft,jpngr))
    ! print*, 'dlitt   ', orgminus(plitt_af(pft,jpngr), orgtmp2)
    ! print*, 'balance=', orgminus( orgminus(orgtmp, pleaf(pft,jpngr)), orgminus(plitt_af(pft,jpngr), orgtmp2) )
    ! <<<<<<<<<<<

    ! ! old version didn't conserve mass: >>>>>>>
    ! ! determine absolute turnover
    ! lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) ) ! leaf turnover

    ! ! reduce leaf mass and root mass
    ! call orgsub( lm_turn, pleaf(pft,jpngr) )

    ! ! add all organic (fixed) C to litter
    ! call cmvRec( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))

    ! ! retain fraction of N
    ! call nmv( nfrac( params_plant%f_nretain, lm_turn%n ), lm_turn%n, plabl(pft,jpngr)%n )

    ! ! rest goes to litter
    ! call nmvRec( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    ! ! <<<<<<<<<<<

  end subroutine turnover_leaf


  subroutine turnover_root( droot, pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction droot for root pool
    !------------------------------------------------------------------
    use md_classdefs

    ! arguments
    real, intent(in)    :: droot
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    ! local variables
    type(orgpool) :: rm_turn

    ! new version: >>>>>>>
    ! determine absolute turnover
    rm_turn = orgfrac( droot, proot(pft,jpngr) ) ! root turnover

    ! reduce root mass and root mass
    call orgmv( rm_turn, proot(pft,jpngr), plitt_bg(pft,jpngr) )
    ! <<<<<<<<<<<

    ! ! old version didn't conserve mass: >>>>>>>
    ! ! determine absolute turnover
    ! rm_turn = orgfrac( droot, proot(pft,jpngr) ) ! root turnover

    ! ! reduce leaf mass and root mass
    ! call orgsub( rm_turn, proot(pft,jpngr) )

    ! ! add all organic (fixed) C to litter
    ! call cmvRec( rm_turn%c, rm_turn%c, plitt_bg(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))

    ! ! retain fraction of N
    ! call nmv( nfrac( params_plant%f_nretain, rm_turn%n ), rm_turn%n, plabl(pft,jpngr)%n )

    ! ! rest goes to litter
    ! call nmvRec( rm_turn%n, rm_turn%n, plitt_bg(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    ! ! <<<<<<<<<<<

  end subroutine turnover_root


  subroutine turnover_labl( dlabl, pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction dlabl for labl pool
    !------------------------------------------------------------------
    ! arguments
    real, intent(in)    :: dlabl
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    ! local variables
    type(orgpool) :: lb_turn

    ! new version: >>>>>>>
    ! determine absolute turnover
    lb_turn = orgfrac( dlabl, plabl(pft,jpngr) ) ! labl turnover

    ! reduce labl mass and labl mass
    call orgmv( lb_turn, plabl(pft,jpngr), plitt_af(pft,jpngr) )
    ! <<<<<<<<<<<

    ! ! old version didn't conserve mass: >>>>>>>
    ! ! determine absolute turnover
    ! lb_turn = orgfrac( dlabl, plabl(pft,jpngr) ) ! labl turnover

    ! ! reduce labl mass and root mass
    ! call orgsub( lb_turn, plabl(pft,jpngr) )

    ! ! add all organic (fixed) C to litter
    ! call cmvRec( lb_turn%c, lb_turn%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))

    ! ! rest goes to litter
    ! call nmvRec( lb_turn%n, lb_turn%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    ! ! <<<<<<<<<<<

  end subroutine turnover_labl

  subroutine turnover_sapw( dsapw, pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction dsapw for sapw pool
    !------------------------------------------------------------------
    ! arguments
    real, intent(in)    :: dsapw
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    ! ! local variables
    ! type(orgpool) :: lm_turn

    ! ! determine absolute turnover
    ! lm_turn = orgfrac( dsapw, psapw(pft,jpngr) ) ! sapw turnover

    ! ! reduce sapw mass and root mass
    ! call orgsub( lm_turn, psapw(pft,jpngr) )

    ! ! add all organic (fixed) C to litter
    ! call cmvRec( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))

    ! ! retain fraction of N
    ! ! xxx try
    ! F_NRETAIN = 0.0
    ! call nmv( nfrac( F_NRETAIN, lm_turn%n ), lm_turn%n, plabl(pft,jpngr)%n )

    ! ! rest goes to litter
    ! call nmvRec( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )

    ! !--------------------------------------------------------------
    ! ! Update foliage-related state variables (lai_ind, fpc_grid, and fapar_ind)
    ! !--------------------------------------------------------------
    ! call update_foliage_vars( pft, jpngr )

  end subroutine turnover_sapw

end module md_turnover
