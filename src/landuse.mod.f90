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
    use md_params_core, only: npft
    use md_interface
    use md_plant, only: pleaf, plabl

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: doy

    ! local variables
    integer :: pft
    integer :: lu

    ! real, parameter :: dleaf = 0.5
    ! real, parameter :: dlabl = 0.25

    ! xxx try
    real, parameter :: min_cleaf_left = 28.0
    real :: dleaf
    real :: dlabl

    type( orgpool ) :: lm_turn
    type( orgpool ) :: lb_turn


    do pft=1,npft

      if (plabl(pft,jpngr)%c%c12<0.0) stop 'labile C is neg.'
      if (plabl(pft,jpngr)%n%n14<0.0) stop 'labile N is neg.'

      if ( params_pft_plant(pft)%grass .and. interface%landuse(jpngr)%do_grharvest(doy) ) then
        ! grasses are harvested

        ! print*,'harvest on day ', doy

        if (pleaf(pft,jpngr)%c%c12>min_cleaf_left) then
          dleaf = ( 1.0 - min_cleaf_left / pleaf(pft,jpngr)%c%c12 )
          ! dlabl = 0.0
          dlabl = min( 1.0, max( 0.0, dleaf / 0.5 ) )
        end if

        ! determine absolute turnover
        lm_turn = orgfrac( dleaf, pleaf(pft,jpngr) ) ! leaf turnover
        lb_turn = orgfrac( dlabl, plabl(pft,jpngr) ) ! leaf turnover

        ! reduce leaf mass and root mass
        call orgmv( lm_turn, pleaf(pft,jpngr), mharv(pft,jpngr) )
        call orgmv( lb_turn, plabl(pft,jpngr), mharv(pft,jpngr) )

        ! call orgsub( lm_turn, pleaf(pft,jpngr) )
        ! call orgsub( lb_turn, plabl(pft,jpngr) )
        
        ! ! add harvested biomass to harvest pool (off site decay, 100%/yr)
        ! call orgmvRec( lm_turn, lm_turn, mharv(pft,jpngr), outacharv(pft,jpngr), outanharv(pft,jpngr), scale=nind(pft,jpngr) )
        ! call orgmvRec( lb_turn, lb_turn, mharv(pft,jpngr), outacharv(pft,jpngr), outanharv(pft,jpngr), scale=nind(pft,jpngr) )

      endif

    end do

  end subroutine grharvest

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
    use md_interface

    if (interface%params_siml%loutlanduse) then

      ! annual output variables   
      outacharv(:,:) = 0.0
      outanharv(:,:) = 0.0

    end if

  end subroutine initoutput_landuse


  subroutine initio_landuse()
    !////////////////////////////////////////////////////////////////
    ! Opens input/output files.
    !----------------------------------------------------------------
    use md_interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    if (interface%params_siml%loutlanduse) then
      !////////////////////////////////////////////////////////////////
      ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
      !----------------------------------------------------------------

      ! HARVESTED C
      filnam=trim(prefix)//'.a.charv.out'
      open(360,file=filnam,err=999,status='unknown')

      ! HARVESTED N
      filnam=trim(prefix)//'.a.nharv.out'
      open(361,file=filnam,err=999,status='unknown')

    end if

    return

    999  stop 'INITIO_LANDUSE: error opening output files'

  end subroutine initio_landuse


  subroutine getout_annual_landuse( jpngr )
    !////////////////////////////////////////////////////////////////
    !  SR called once a year to gather annual output variables.
    !----------------------------------------------------------------
    use md_interface
    use md_params_core, only: ndayyear, npft

    ! arguments
    integer, intent(in) :: jpngr

    ! Output annual value at day of peak LAI
    if (interface%params_siml%loutlanduse) then
      outacharv(:,:) = mharv(:,:)%c%c12
      outanharv(:,:) = mharv(:,:)%n%n14
    end if

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

    if (interface%params_siml%loutlanduse) then
      !-------------------------------------------------------------------------
      ! ANNUAL OUTPUT
      ! Write annual value, summed over all PFTs / LUs
      ! xxx implement taking sum over PFTs (and gridcells) in this land use category
      !-------------------------------------------------------------------------
      itime = real(interface%steering%outyear)

      write(360,999) itime, sum(outacharv(:,jpngr))
      write(361,999) itime, sum(outanharv(:,jpngr))

    end if

    return

  999     format (F20.8,F20.8)

  end subroutine writeout_ascii_landuse

end module md_landuse