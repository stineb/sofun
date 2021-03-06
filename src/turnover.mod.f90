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
  public turnover, turnover_root, turnover_leaf, turnover_labl, initoutput_turnover

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  real, dimension(:,:), allocatable :: outaCveg2lit
  real, dimension(:,:), allocatable :: outaNveg2lit

contains

  subroutine turnover( plant, out_pmodel, solar, temppheno )
    !////////////////////////////////////////////////////////////////
    !  Annual vegetation biomass turnover, called at the end of the
    !  year.
    !----------------------------------------------------------------
    use md_classdefs
    use md_params_core, only: npft, eps, nmonth
    use md_plant, only: plant_type
    use md_waterbal, only: solartype
    use md_gpp, only: outtype_pmodel
    use md_phenology, only: temppheno_type

    ! arguments
    type( plant_type ), dimension(npft), intent(inout)         :: plant ! npft counts over PFTs in all land units (tiles)
    type( outtype_pmodel ), dimension(npft,nmonth), intent(in) :: out_pmodel
    type( solartype ), intent(in)                              :: solar
    type(temppheno_type), dimension(npft), intent(in)          :: temppheno

    ! local variables
    integer :: pft
    integer :: lu
    real :: dlabl
    real :: dleaf
    real :: droot
    ! real :: balance

    ! xxx verbose
    logical, parameter :: verbose = .false.
    type( orgpool ) :: orgtmp, orgtmp2

    do pft=1,npft

      if (plant(pft)%plabl%c%c12 < -1.0*eps) stop 'before turnover labile C is neg.'
      if (plant(pft)%plabl%n%n14 < -1.0*eps) stop 'before turnover labile N is neg.'

      !--------------------------------------------------------------
      ! Get turnover fractions
      ! Turnover-rates are reciprocals of tissue longevity
      ! dleaf=1.0/long_leaf(pft)
      ! assuming no continuous leaf turnover
      !--------------------------------------------------------------
      if (params_pft_plant(pft)%grass) then

        ! balance = plant_fluxes(pft)%dnpp%c12 - plant_fluxes(pft)%dcex

        if (temppheno(pft)%shedleaves) then

          droot = 1.0
          dleaf = 1.0
          dlabl = 1.0

          stop 'shedding the fucking leaves'
          
        else

          ! Increase turnover rate towards high LAI ( when using non-zero value for k_decay_leaf_width, e.g. 0.08 )
          dleaf =  ( plant(pft)%lai_ind * params_pft_plant(pft)%k_decay_leaf_width )**8 + params_pft_plant(pft)%k_decay_leaf_base

          ! constant turnover rate
          droot = params_pft_plant(pft)%k_decay_root
          dlabl = 0.0 !    params_pft_plant(pft)%k_decay_labl

        end if

      else

        stop 'turnover not implemented for non-grasses'

      endif

      !--------------------------------------------------------------
      ! Calculate leaf turnover in this day 
      !--------------------------------------------------------------
      if (verbose) print*, 'calling turnover_leaf() ... '
      if (verbose) print*, '              with state variables:'
      if (verbose) print*, '              pleaf = ', plant(pft)%pleaf
      if (verbose) print*, '              plitt = ', plant(pft)%plitt_af
      if (verbose) orgtmp  =  plant(pft)%pleaf
      if (verbose) orgtmp2 =  plant(pft)%plitt_af
      !--------------------------------------------------------------
      if ( dleaf>0.0 ) call turnover_leaf( dleaf, plant(pft), out_pmodel(pft,:), solar, pft ) !, jpngr
      !--------------------------------------------------------------
      if (verbose) print*, '              ==> returned: '
      if (verbose) print*, '              pleaf = ', plant(pft)%pleaf
      if (verbose) print*, '              plitt = ', plant(pft)%plitt_af
      if (verbose) print*, '              --- balance: '
      if (verbose) print*, '                  dlitt - dleaf                = ',  orgminus( &
                                                                                    orgminus( &
                                                                                      plant(pft)%plitt_af, &
                                                                                      orgtmp2 &
                                                                                      ), &
                                                                                    orgminus( &
                                                                                      orgtmp, &
                                                                                      plant(pft)%pleaf &
                                                                                      ) &
                                                                                    )

      !--------------------------------------------------------------
      ! Calculate root turnover in this day 
      !--------------------------------------------------------------
      if (verbose) print*, 'calling turnover_root() ... '
      if (verbose) print*, '              with state variables:'
      if (verbose) print*, '              pleaf = ', plant(pft)%proot
      if (verbose) print*, '              plitt = ', plant(pft)%plitt_bg
      if (verbose) orgtmp  =  plant(pft)%proot
      if (verbose) orgtmp2 =  plant(pft)%plitt_bg
      !--------------------------------------------------------------
      if ( droot>0.0 ) call turnover_root( droot, plant(pft) )  ! pft, jpngr
      !--------------------------------------------------------------
      if (verbose) print*, '              ==> returned: '
      if (verbose) print*, '              proot = ', plant(pft)%proot
      if (verbose) print*, '              plitt = ', plant(pft)%plitt_bg
      if (verbose) print*, '              --- balance: '
      if (verbose) print*, '                  dlitt - droot                = ',  orgminus( &
                                                                                    orgminus( &
                                                                                      plant(pft)%plitt_bg, &
                                                                                      orgtmp2 &
                                                                                      ), &
                                                                                    orgminus( &
                                                                                      orgtmp, &
                                                                                      plant(pft)%proot &
                                                                                      ) &
                                                                                    )

      !--------------------------------------------------------------
      ! Calculate labile turnover in this day 
      !--------------------------------------------------------------
      if (verbose) print*, 'calling turnover_root() ... '
      if (verbose) print*, '              with state variables:'
      if (verbose) print*, '              pleaf = ', plant(:)%plabl
      if (verbose) print*, '              plitt = ', plant(pft)%plitt_af
      if (verbose) orgtmp  =  plant(pft)%plabl
      if (verbose) orgtmp2 =  plant(pft)%plitt_af
      !--------------------------------------------------------------
      if ( dlabl>0.0 ) call turnover_labl( dlabl, plant(pft) )  !  pft, jpngr
      !--------------------------------------------------------------
      if (verbose) print*, '              ==> returned: '
      if (verbose) print*, '              plabl = ', plant(:)%plabl
      if (verbose) print*, '              plitt = ', plant(pft)%plitt_af
      if (verbose) print*, '              --- balance: '
      if (verbose) print*, '                  dlitt - dlabl                = ',  orgminus( &
                                                                                    orgminus( &
                                                                                      plant(pft)%plitt_af, &
                                                                                      orgtmp2 &
                                                                                      ), &
                                                                                    orgminus( &
                                                                                      orgtmp, &
                                                                                      plant(pft)%proot &
                                                                                      ) &
                                                                                    )
    enddo                     !pft

  end subroutine turnover


  subroutine turnover_leaf( dleaf, plant, out_pmodel, solar, pft )  ! pft, jpngr
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction dleaf for leaf pool
    !------------------------------------------------------------------
    use md_params_core, only: npft, eps, nmonth
    use md_waterbal, only: solartype
    use md_gpp, only: outtype_pmodel
    use md_plant, only: plant_type, get_fapar

    ! arguments
    real, intent(in)    :: dleaf
    type( plant_type ), intent(inout)  :: plant ! npft counts over PFTs in all land units (tiles)
    type( outtype_pmodel ), dimension(nmonth), intent(in) :: out_pmodel
    type( solartype ), intent(in)      :: solar
    integer, intent(in) :: pft

    ! local variables
    type(orgpool) :: lm_turn
    type(orgpool) :: lm_init

    real :: nleaf
    real :: cleaf
    real :: dlai
    real :: lai_new
    real :: diff
    integer :: nitr

    ! number of iterations to match leaf C given leaf N
    nitr = 0

    ! store leaf C and N before turnover
    lm_init = plant%pleaf


    ! reduce leaf C (given by turnover rate)
    cleaf = ( 1.0 - dleaf ) *  plant%pleaf%c%c12

    ! get new LAI based on cleaf
    ! print*,'IN TURNOVER: out_pmodel(:)%actnv_unitiabs:'
    ! print*,out_pmodel(:)%actnv_unitiabs
    lai_new = get_lai( pft, cleaf, solar%meanmppfd(:), out_pmodel(:)%actnv_unitiabs )

    ! update canopy state (only variable fAPAR so far implemented)
    plant%fapar_ind = get_fapar( lai_new )

    ! re-calculate metabolic and structural N, given new LAI and fAPAR
    call update_leaftraits( plant, pft, lai_new, solar%meanmppfd(:), out_pmodel(:)%actnv_unitiabs )
    ! leaftraits(pft) = get_leaftraits( pft, lai_new, solar%meanmppfd(:), out_pmodel%actnv_unitiabs )

    ! get updated leaf N
    nleaf = plant%narea_canopy

    do while ( nleaf > lm_init%n%n14 )

      nitr = nitr + 1

      ! reduce leaf C a bit more
      cleaf = cleaf * lm_init%n%n14 / nleaf

      ! get new LAI based on cleaf
      lai_new = get_lai( pft, cleaf, solar%meanmppfd(:), out_pmodel(:)%actnv_unitiabs )

      ! update canopy state (only variable fAPAR so far implemented)
      plant%fapar_ind = get_fapar( lai_new )

      ! re-calculate metabolic and structural N, given new LAI and fAPAR
      call update_leaftraits( plant, pft, lai_new, solar%meanmppfd(:), out_pmodel(:)%actnv_unitiabs )
      ! leaftraits(pft) = get_leaftraits( pft, lai_new, solar%meanmppfd(:), out_pmodel(:)%actnv_unitiabs )

      ! get updated leaf N
      nleaf = plant%narea_canopy

      if (nitr>30) exit

    end do

    ! if (nitr>0) print*,'no. of iterations ', nitr
    ! if (nitr>0) print*,'final reduction of leaf C ', cleaf / lm_init%c%c12
    ! if (nitr>0) print*,'final reduction of leaf N ', nleaf / lm_init%n%n14

    ! update 
    plant%lai_ind = lai_new
    plant%pleaf%c%c12 = cleaf
    plant%pleaf%n%n14 = nleaf

    ! determine C and N turned over
    lm_turn = orgminus( lm_init, plant%pleaf )

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
    ! call cmvRec( lm_turn%c, lm_turn%c, plant%plitt_af%c, outaCveg2lit(pft,jpngr), scale=plant%nind)
    call cmv( lm_turn%c, lm_turn%c, plant%plitt_af%c, scale=plant%nind )

    ! retain fraction of N
    call nmv( nfrac( params_plant%f_nretain, lm_turn%n ), lm_turn%n, plant%plabl%n )

    ! rest goes to litter
    ! call nmvRec( lm_turn%n, lm_turn%n, plant%plitt_af%n, outaNveg2lit(pft,jpngr), scale=plant%nind )
    call nmv( lm_turn%n, lm_turn%n, plant%plitt_af%n, scale=plant%nind )

  end subroutine turnover_leaf


  subroutine turnover_root( droot, plant ) ! pft, jpngr
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction droot for root pool
    !------------------------------------------------------------------
    use md_plant, only: plant_type

    ! arguments
    real, intent(in)    :: droot
    type( plant_type ), intent(inout) :: plant ! npft counts over PFTs in all land units (tiles)
    ! integer, intent(in) :: pft
    ! integer, intent(in) :: jpngr

    ! local variables
    type(orgpool) :: rm_turn

    ! determine absolute turnover
    rm_turn = orgfrac( droot, plant%proot ) ! root turnover

    ! reduce leaf mass and root mass
    call orgsub( rm_turn, plant%proot )

    ! add all organic (fixed) C to litter
    ! call cmvRec( rm_turn%c, rm_turn%c, plant%plitt_bg%c, outaCveg2lit(pft,jpngr), scale=plant%nind)
    call cmv( rm_turn%c, rm_turn%c, plant%plitt_bg%c, scale=plant%nind)

    ! retain fraction of N
    call nmv( nfrac( params_plant%f_nretain, rm_turn%n ), rm_turn%n, plant%plabl%n )

    ! rest goes to litter
    ! call nmvRec( rm_turn%n, rm_turn%n, plant%plitt_bg%n, outaNveg2lit(pft,jpngr), scale=plant%nind )
    call nmv( rm_turn%n, rm_turn%n, plant%plitt_bg%n, scale=plant%nind )

  end subroutine turnover_root


  subroutine turnover_labl( dlabl, plant ) ! pft, jpngr
    !//////////////////////////////////////////////////////////////////
    ! Execute turnover of fraction dlabl for labl pool
    !------------------------------------------------------------------
    use md_plant, only: plant_type

    ! arguments
    real, intent(in)    :: dlabl
    type( plant_type ), intent(inout) :: plant ! npft counts over PFTs in all land units (tiles)
    ! integer, intent(in) :: pft
    ! integer, intent(in) :: jpngr

    ! local variables
    type(orgpool) :: lb_turn

    ! detelbine absolute turnover
    lb_turn = orgfrac( dlabl, plant%plabl ) ! labl turnover

    !! xxx think of something more plausible to put the labile C and N to

    ! reduce leaf mass and labl mass
    call orgsub( lb_turn, plant%plabl )

    ! call orgmvRec( lb_turn, lb_turn, plant%plitt_af, outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr), scale=plant%nind )
    call orgmv( lb_turn, lb_turn, plant%plitt_af, scale=plant%nind )

  end subroutine turnover_labl


  subroutine initoutput_turnover( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_interface, only: interface
    use md_params_core, only: npft

    ! arguments
    integer, intent(in) :: ngridcells
    
    ! annual output variables
    if (interface%params_siml%loutturnover) then

      if (interface%steering%init) then
        allocate( outaCveg2lit(npft,ngridcells) )
        allocate( outaCveg2lit(npft,ngridcells) )
      end if
      
      outaCveg2lit(:,:) = 0.0
      outaCveg2lit(:,:) = 0.0

    end if

  end subroutine initoutput_turnover

end module md_turnover
