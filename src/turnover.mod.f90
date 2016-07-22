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
  public turnover, turnover_root, turnover_leaf, turnover_labl

contains

  subroutine turnover( jpngr, doy )
    !////////////////////////////////////////////////////////////////
    !  Annual vegetation biomass turnover, called at the end of the
    !  year.
    !----------------------------------------------------------------
    use md_classdefs
    use md_params_core, only: npft
    use md_phenology, only: shedleaves

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu
    real :: dlabl
    real :: dleaf
    real :: droot

    ! xxx verbose
    logical, parameter :: verbose = .false.
    type( orgpool ) :: orgtmp, orgtmp2

    do pft=1,npft

      if (plabl(pft,jpngr)%c%c12<0.0) stop 'before turnover labile C is neg.'
      if (plabl(pft,jpngr)%n%n14<0.0) stop 'before turnover labile N is neg.'

      !--------------------------------------------------------------
      ! Get turnover fractions
      ! Turnover-rates are reciprocals of tissue longevity
      ! dleaf=1.0/long_leaf(pft)
      ! assuming no continuous leaf turnover
      !--------------------------------------------------------------
      if (params_pft_plant(pft)%grass) then

        if (shedleaves(doy,pft)) then

          droot = 1.0
          dleaf = 1.0
          dlabl = 1.0
          
        else

          ! dlabl = 0.01
          dlabl = 0.0

          ! Alternative turnover function: increase turnover rate towards high LAI
          ! dleaf = (lai_ind(pft,jpngr)*params_pft_plant(pft)%k_decay_leaf_width)**8 + params_pft_plant(pft)%k_decay_leaf_base

          droot =  (lai_ind(pft,jpngr)*params_pft_plant(pft)%k_decay_leaf_width)**8 + params_pft_plant(pft)%k_decay_leaf_base
          dleaf =  (lai_ind(pft,jpngr)*params_pft_plant(pft)%k_decay_leaf_width)**8 + params_pft_plant(pft)%k_decay_leaf_base

          ! xxx try
          ! dleaf = 1.5 / 365.0
          ! droot = 1.5 / 365.0

        end if

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

      !--------------------------------------------------------------
      ! Calculate labile turnover in this day 
      !--------------------------------------------------------------
      if (verbose) print*, 'calling turnover_root() ... '
      if (verbose) print*, '              with state variables:'
      if (verbose) print*, '              pleaf = ', plabl(:,jpngr)
      if (verbose) print*, '              plitt = ', plitt_af(pft,jpngr)
      if (verbose) orgtmp  =  plabl(pft,jpngr)
      if (verbose) orgtmp2 =  plitt_af(pft,jpngr)
      !--------------------------------------------------------------
      if ( dlabl>0.0 )                 call turnover_labl( dlabl, pft, jpngr )
      !--------------------------------------------------------------
      if (verbose) print*, '              ==> returned: '
      if (verbose) print*, '              plabl = ', plabl(:,jpngr)
      if (verbose) print*, '              plitt = ', plitt_af(pft,jpngr)
      if (verbose) print*, '              --- balance: '
      if (verbose) print*, '                  dlitt - dlabl                = ',  orgminus( &
                                                                                    orgminus( &
                                                                                      plitt_af(pft,jpngr), &
                                                                                      orgtmp2 &
                                                                                      ), &
                                                                                    orgminus( &
                                                                                      orgtmp, &
                                                                                      proot(pft,jpngr) &
                                                                                      ) &
                                                                                    )
    enddo                     !pft

  end subroutine turnover


  subroutine turnover_leaf( dleaf, pft, jpngr )
    !//////////////////////////////////////////////////////////////////
    ! Leaf turnover governed by LAI turnover, leaf C and N follows from
    ! new LAI.
    !------------------------------------------------------------------
    use md_plant, only: lai_ind, leaftraits, canopy, get_canopy, get_leaftraits, get_lai
    use md_waterbal, only: solar
    use md_gpp, only: mactnv_unitiabs
    use md_params_core, only: eps

    ! arguments
    real, intent(in)    :: dleaf
    integer, intent(in) :: pft
    integer, intent(in) :: jpngr

    ! local variables
    type(orgpool) :: lm_init
    type(orgpool) :: lm_turn
    real :: nleaf
    real :: cleaf
    real :: dlai
    real :: lai_new
    real :: diff, tmp1, tmp2

    integer, save :: invoc = 0
    real, save :: iabs_times_nv

    invoc = invoc + 1

    ! store leaf C and N before turnover
    lm_init = pleaf(pft,jpngr)

    ! reduce LAI
    lai_new = (1.0 - dleaf) * lai_ind(pft,jpngr)

    ! update canopy state (only variable fAPAR so far implemented)
    canopy(pft) = get_canopy( lai_new )

    ! re-calculate metabolic and structural N, given new LAI and fAPAR
    leaftraits(pft) = get_leaftraits( pft, lai_new, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )

    nleaf = leaftraits(pft)%narea_canopy
    cleaf = leaftraits(pft)%leafc_canopy

    ! when light conditions change (annual maximum of solar%meanmppfd(:) * mactnv_unitiabs(pft,:)), 
    ! more N and C may be needed in spite of LAI-turnover > 0.0. 
    ! In this case, calculate LAI as a function of leaf C 
    if ( nleaf>pleaf(pft,jpngr)%n%n14 .or. cleaf>pleaf(pft,jpngr)%c%c12 ) then

      if ( invoc > 1 .and. maxval( solar%meanmppfd(:) * mactnv_unitiabs(pft,:) ) > iabs_times_nv ) then
        ! Light conditions change so that required C and N at current LAI are insufficient. Reduce LAI.
        diff = maxval( solar%meanmppfd(:) * mactnv_unitiabs(pft,:) ) / iabs_times_nv

        ! print*,'diff ', diff
        ! tmp1 = lai_ind(pft,jpngr) * params_plant%kbeer
        ! tmp2 = exp( tmp2 )
        ! dlai = 1.0 / tmp1 * log( diff * tmp2 / (diff * tmp2 - tmp2 + 1.0 ) )

        ! try        
        dlai = (1.0 / diff)**2
        lai_new = dlai * lai_ind(pft,jpngr)

        ! update canopy state (only variable fAPAR so far implemented)
        canopy(pft) = get_canopy( lai_new )

        ! re-calculate metabolic and structural N, given new LAI and fAPAR
        leaftraits(pft) = get_leaftraits( pft, lai_new, solar%meanmppfd(:), mactnv_unitiabs(pft,:) )

        nleaf = leaftraits(pft)%narea_canopy
        cleaf = leaftraits(pft)%leafc_canopy

        ! when light conditions change (annual maximum of solar%meanmppfd(:) * mactnv_unitiabs(pft,:)), 
        ! more N and C may be needed in spite of LAI-turnover > 0.0. 
        ! In this case, calculate LAI as a function of leaf C 
        if ( nleaf>pleaf(pft,jpngr)%n%n14 .or. cleaf>pleaf(pft,jpngr)%c%c12 ) then
          print*,'nleaf ', nleaf
          print*,'pleaf(pft,jpngr)%n%n14 ', pleaf(pft,jpngr)%n%n14
          stop 'shit'
        endif

      else

        if ( (cleaf - pleaf(pft,jpngr)%c%c12) > eps ) then
          print*,'dleaf ', dleaf
          print*,'cleaf ', cleaf
          print*,'pleaf(pft,jpngr)%c%c12 ', pleaf(pft,jpngr)%c%c12
          stop 'unknown reason'
        else
          pleaf(pft,jpngr)%c%c12 = cleaf
        end if

        if ( (nleaf - pleaf(pft,jpngr)%n%n14) > eps ) then
          print*,'dleaf ', dleaf
          print*,'nleaf ', nleaf
          print*,'pleaf(pft,jpngr)%n%n14 ', pleaf(pft,jpngr)%n%n14
          stop 'unknown reason'
        else
          pleaf(pft,jpngr)%n%n14 = nleaf
        end if

      end if

    end if

    ! update 
    lai_ind(pft,jpngr) = lai_new
    pleaf(pft,jpngr)%c%c12 = cleaf
    pleaf(pft,jpngr)%n%n14 = nleaf

    ! determine C and N turned over
    lm_turn = orgminus( lm_init, pleaf(pft,jpngr) )

    if ( lm_turn%c%c12 < -1.0*eps ) then
      stop 'negative turnover C'
    else if ( lm_turn%c%c12 < 0.0 ) then
       lm_turn%c%c12 = 0.0
    end if
    if ( lm_turn%n%n14 < -1.0*eps ) then
      stop 'negative turnover N'
    else if ( lm_turn%n%n14 < 0.0 ) then
       lm_turn%n%n14 = 0.0
    end if

    ! add all organic (fixed) C to litter
    call cmvRec( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))
    ! call cmv( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, scale=nind(pft,jpngr) )

    ! retain fraction of N
    call nmv( nfrac( params_plant%f_nretain, lm_turn%n ), lm_turn%n, plabl(pft,jpngr)%n )

    ! rest goes to litter
    call nmvRec( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    ! call nmv( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, scale=nind(pft,jpngr) )

    ! save for next time step
    iabs_times_nv = maxval( solar%meanmppfd(:) * mactnv_unitiabs(pft,:) )

  end subroutine turnover_leaf


  ! subroutine turnover_leaf( dleaf, pft, jpngr )
  !   !//////////////////////////////////////////////////////////////////
  !   ! Execute turnover of fraction dleaf for leaf pool
  !   !------------------------------------------------------------------
  !   ! arguments
  !   real, intent(in)    :: dleaf
  !   integer, intent(in) :: pft
  !   integer, intent(in) :: jpngr

  !   ! local variables
  !   type(orgpool) :: lm_turn

  !   ! print*,'dleaf ', dleaf
  !   ! stop

  !   ! determine absolute turnover
  !   lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) ) ! leaf turnover

  !   ! reduce leaf mass and root mass
  !   call orgsub( lm_turn, pleaf(pft,jpngr) )

  !   ! add all organic (fixed) C to litter
  !   call cmvRec( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))
  !   ! call cmv( lm_turn%c, lm_turn%c, plitt_af(pft,jpngr)%c, scale=nind(pft,jpngr) )

  !   ! retain fraction of N
  !   call nmv( nfrac( params_plant%f_nretain, lm_turn%n ), lm_turn%n, plabl(pft,jpngr)%n )

  !   ! rest goes to litter
  !   call nmvRec( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
  !   ! call nmv( lm_turn%n, lm_turn%n, plitt_af(pft,jpngr)%n, scale=nind(pft,jpngr) )

  ! end subroutine turnover_leaf


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
    call cmvRec( rm_turn%c, rm_turn%c, plitt_bg(pft,jpngr)%c, outaCveg2lit(pft,jpngr), scale=nind(pft,jpngr))
    ! call cmv( rm_turn%c, rm_turn%c, plitt_bg(pft,jpngr)%c, scale=nind(pft,jpngr))

    ! retain fraction of N
    call nmv( nfrac( params_plant%f_nretain, rm_turn%n ), rm_turn%n, plabl(pft,jpngr)%n )

    ! rest goes to litter
    call nmvRec( rm_turn%n, rm_turn%n, plitt_bg(pft,jpngr)%n, outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr) )
    ! call nmv( rm_turn%n, rm_turn%n, plitt_bg(pft,jpngr)%n, scale=nind(pft,jpngr) )

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

    if (dlabl>0.0) then

      ! detelbine absolute turnover
      lb_turn = orgfrac( dlabl, plabl(pft,jpngr) ) ! labl turnover

      ! reduce leaf mass and labl mass
      call orgsub( lb_turn, plabl(pft,jpngr) )

      call orgmvRec( lb_turn, lb_turn, plitt_af(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr), scale=nind(pft,jpngr))

    end if

  end subroutine turnover_labl


end module md_turnover
