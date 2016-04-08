module md_landuse
  !////////////////////////////////////////////////////////////////
  !----------------------------------------------------------------
  use md_classdefs
  use md_plant
  use md_params_core, only: npft, maxgrid
    
  implicit none

  private
  public grharvest

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
    real, parameter :: min_cleaf_left = 50.0
    real :: dleaf
    real :: dlabl

    type( orgpool ) :: lm_turn
    type( orgpool ) :: lb_turn

    ! set to zero on the first day of the year
    if (doy==1) then
      mharv(:,jpngr) = orgpool( carbon(0.0), nitrogen(0.0))
    end if

    if (plabl(pft,jpngr)%c%c12<0.0) stop 'labile C is neg.'
    if (plabl(pft,jpngr)%n%n14<0.0) stop 'labile N is neg.'

    do pft=1,npft

      if ( params_pft_plant(pft)%grass .and. interface%landuse(jpngr)%do_grharvest(doy) ) then
        ! grasses are harvested

        print*,'harvest on day ', doy

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

end module md_landuse