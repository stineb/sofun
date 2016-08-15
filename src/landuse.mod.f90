module md_landuse
  !////////////////////////////////////////////////////////////////
  !----------------------------------------------------------------
  use md_classdefs
  use md_plant
  use md_params_core, only: npft, maxgrid
    
  implicit none

  private
  public mharv, grharvest, initoutput_landuse, initio_landuse, &
    getout_annual_landuse, writeout_ascii_landuse, initglobal_landuse, &
    init_mharv

  !----------------------------------------------------------------
  ! Private, module-specific state variables
  !----------------------------------------------------------------
  type(orgpool), dimension(npft,maxgrid) :: mharv  ! harvested biomass [gC/ind.] (=lm_ind)

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, dimension(npft,maxgrid) :: outacharv
  real, dimension(npft,maxgrid) :: outanharv


contains


  subroutine grharvest( jpngr, doy )
    !////////////////////////////////////////////////////////////////
    !  Annual vegetation biomass turnover, called at the end of the
    !  year.
    !----------------------------------------------------------------
    use md_classdefs
    use md_params_core, only: npft, eps
    use md_interface
    use md_plant, only: pleaf, plabl
    use md_waterbal, only: solar
    use md_gpp, only: out_pmodel

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    type( orgpool ) :: lm_init
    type( orgpool ) :: lm_turn
    type( orgpool ) :: lb_turn
    real :: cleaf
    real :: nleaf
    real :: dleaf
    real :: dlabl
    real :: lai_new

    integer :: nitr
    integer :: pft
    integer :: lu

    ! xxx try
    real, parameter :: min_cleaf_left = 30.0

    do pft=1,npft

      if ( params_pft_plant(pft)%grass .and. interface%landuse(jpngr)%do_grharvest(doy) ) then

        if (pleaf(pft,jpngr)%c%c12>min_cleaf_left) then

          ! number of iterations to match leaf C given leaf N
          nitr = 0

          ! store leaf C and N before turnover
          lm_init = pleaf(pft,jpngr)

          ! reduce leaf C (given by turnover rate)
          cleaf = min_cleaf_left

          ! get new LAI based on cleaf
          lai_new = get_lai( pft, cleaf, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )

          ! update canopy state (only variable fAPAR so far implemented)
          canopy(pft) = get_canopy( lai_new )

          ! re-calculate metabolic and structural N, given new LAI and fAPAR
          leaftraits(pft) = get_leaftraits( pft, lai_new, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )

          ! get updated leaf N
          nleaf = leaftraits(pft)%narea_canopy

          do while ( nleaf > lm_init%n%n14 )

            print*,'initial leaf C reduction: ', cleaf / lm_init%c%c12
            stop

            nitr = nitr + 1

            ! reduce leaf C a bit more
            cleaf = cleaf * lm_init%n%n14 / nleaf

            ! get new LAI based on cleaf
            lai_new = get_lai( pft, cleaf, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )

            ! update canopy state (only variable fAPAR so far implemented)
            canopy(pft) = get_canopy( lai_new )

            ! re-calculate metabolic and structural N, given new LAI and fAPAR
            leaftraits(pft) = get_leaftraits( pft, lai_new, solar%meanmppfd(:), out_pmodel(pft,:)%actnv_unitiabs )

            ! get updated leaf N
            nleaf = leaftraits(pft)%narea_canopy

            print*,'leaf C reduced by ', cleaf / lm_init%c%c12
            print*,'leaf N reduced by ', nleaf / lm_init%n%n14

            stop 'should never occurr, no?'

          end do

          if (nitr>0) print*,'no. of iterations ', nitr
          if (nitr>0) print*,'final reduction of leaf C ', cleaf / lm_init%c%c12
          if (nitr>0) print*,'final reduction of leaf N ', nleaf / lm_init%n%n14

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

          ! reduce leaf mass and root mass
          call orgmv( lm_turn, lm_turn, mharv(pft,jpngr) )

          ! labile pool harvest
          dleaf = ( 1.0 - cleaf / lm_init%c%c12 )
          ! xxx try: 
          dlabl = min( 1.0, max( 0.0, dleaf ) )
          lb_turn = orgfrac( dlabl, plabl(pft,jpngr) ) ! leaf turnover
          call orgmv( lb_turn, plabl(pft,jpngr), mharv(pft,jpngr) )

        end if

      end if

    end do

  end subroutine grharvest


  ! subroutine grharvest( jpngr, doy )
  !   !////////////////////////////////////////////////////////////////
  !   !  Annual vegetation biomass turnover, called at the end of the
  !   !  year.
  !   !----------------------------------------------------------------
  !   use md_classdefs
  !   use md_params_core, only: npft
  !   use md_interface
  !   use md_plant, only: pleaf, plabl

  !   ! arguments
  !   integer, intent(in) :: jpngr
  !   integer, intent(in) :: doy

  !   ! local variables
  !   integer :: pft
  !   integer :: lu

  !   ! real, parameter :: dleaf = 0.5
  !   ! real, parameter :: dlabl = 0.25

  !   ! xxx try
  !   real, parameter :: min_cleaf_left = 28.0
  !   real :: dleaf
  !   real :: dlabl

  !   type( orgpool ) :: lm_turn
  !   type( orgpool ) :: lb_turn


  !   do pft=1,npft

  !     if (plabl(pft,jpngr)%c%c12<0.0) stop 'labile C is neg.'
  !     if (plabl(pft,jpngr)%n%n14<0.0) stop 'labile N is neg.'

  !     if ( params_pft_plant(pft)%grass .and. interface%landuse(jpngr)%do_grharvest(doy) ) then
  !       ! grasses are harvested

  !       ! print*,'harvest on day ', doy

  !       if (pleaf(pft,jpngr)%c%c12>min_cleaf_left) then
  !         dleaf = ( 1.0 - min_cleaf_left / pleaf(pft,jpngr)%c%c12 )
  !         ! dlabl = 0.0
  !         dlabl = min( 1.0, max( 0.0, dleaf / 0.5 ) )
  !       end if

  !       ! determine absolute turnover
  !       lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) ) ! leaf turnover
  !       lb_turn = orgfrac( dlabl, plabl(pft,jpngr) ) ! leaf turnover

  !       ! reduce leaf mass and root mass
  !       call orgmv( lm_turn, pleaf(pft,jpngr), mharv(pft,jpngr) )
  !       call orgmv( lb_turn, plabl(pft,jpngr), mharv(pft,jpngr) )

  !       ! call orgsub( lm_turn, pleaf(pft,jpngr) )
  !       ! call orgsub( lb_turn, plabl(pft,jpngr) )
        
  !       ! ! add harvested biomass to harvest pool (off site decay, 100%/yr)
  !       ! call orgmvRec( lm_turn, lm_turn, mharv(pft,jpngr), outacharv(pft,jpngr), outanharv(pft,jpngr), scale=nind(pft,jpngr) )
  !       ! call orgmvRec( lb_turn, lb_turn, mharv(pft,jpngr), outacharv(pft,jpngr), outanharv(pft,jpngr), scale=nind(pft,jpngr) )

  !     endif

  !   end do

  ! end subroutine grharvest

  subroutine initglobal_landuse()
    !////////////////////////////////////////////////////////////////
    !  Initialisation of all pools on all gridcells at the beginning
    !  of the simulation.
    !----------------------------------------------------------------
    mharv(:,:) = orgpool(carbon(0.0),nitrogen(0.0))  

  end subroutine initglobal_landuse


  subroutine init_mharv(jpngr)
    !////////////////////////////////////////////////////////////////
    !  Initialisation harvest pool on first day of the year
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr

    mharv(:,jpngr) = orgpool(carbon(0.0),nitrogen(0.0))  

  end subroutine init_mharv


  subroutine initoutput_landuse()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    ! annual output variables
    outacharv(:,:) = 0.0
    outanharv(:,:) = 0.0

  end subroutine initoutput_landuse


  subroutine initio_landuse()
    !////////////////////////////////////////////////////////////////
    ! Opens input/output files.
    !----------------------------------------------------------------
    use md_interface
    ! use md_params_siml, only: runname

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------

    ! HARVESTED C
    filnam=trim(prefix)//'.a.charv.out'
    open(360,file=filnam,err=999,status='unknown')

    ! HARVESTED N
    filnam=trim(prefix)//'.a.nharv.out'
    open(361,file=filnam,err=999,status='unknown')

    return

    999  stop 'INITIO_LANDUSE: error opening output files'

  end subroutine initio_landuse


  subroutine getout_annual_landuse( jpngr )
    !////////////////////////////////////////////////////////////////
    !  SR called once a year to gather annual output variables.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft

    ! arguments
    integer, intent(in) :: jpngr

    ! Output annual value at day of peak LAI
    outacharv(:,:) = mharv(:,:)%c%c12
    outanharv(:,:) = mharv(:,:)%n%n14

  end subroutine getout_annual_landuse


  subroutine writeout_ascii_landuse( year )
    !/////////////////////////////////////////////////////////////////////////
    ! Write daily ASCII output
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !-------------------------------------------------------------------------
    ! use md_params_siml, only: spinup, interface%params_siml%daily_out_startyr, &
    !   interface%params_siml%daily_out_endyr, outyear
    use md_params_core, only: ndayyear, nlu
    use md_interface

    ! arguments
    integer, intent(in) :: year       ! simulation year

    ! local variables
    real :: itime
    integer :: day, moy, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1

    if (nlu>1) stop 'Output only for one LU category implemented.'

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    itime = real(interface%steering%outyear)

    write(360,999) itime, sum(outacharv(:,jpngr))
    write(361,999) itime, sum(outanharv(:,jpngr))

    return

  999     format (F20.8,F20.8)

  end subroutine writeout_ascii_landuse

end module md_landuse