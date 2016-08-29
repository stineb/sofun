module md_allocation
  !////////////////////////////////////////////////////////////////
  ! ALLOCATION MODULE
  ! Binary allocation formulation: either to leaves or to roots.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: npft, nlu, maxgrid
  use md_plant, only: plant_type

  implicit none

  private 
  public allocation_annual, update_tree, initio_allocation, &
    initoutput_allocation, getout_daily_allocation,         &
    writeout_ascii_allocation

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  real, dimension(npft) :: dcleaf
  real, dimension(npft) :: dcroot
  real, dimension(npft) :: dcsapw
  real, dimension(npft) :: dcwood
  real, dimension(npft) :: dnleaf
  real, dimension(npft) :: dnroot
  real, dimension(npft) :: dnsapw
  real, dimension(npft) :: dnwood

  !----------------------------------------------------------------
  ! Module-specific parameters (PFT dependent)
  !----------------------------------------------------------------
  type params_tmodel_type 
    real :: a_par          ! initial slope of the relationship between height and diameter
    real :: cr_par         ! initial ratio of crown area to stem cross - sectional area ('c')
    real :: maxheight      ! asymptotic maximum height
    real :: rho            ! density of wood
    real :: z_par          ! Ratio of fine-root mass to foliage area (0.17 kg C m-2; White et al. (2000))
  end type params_tmodel_type
  type( params_tmodel_type ), dimension(npft) :: params_tmodel

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! output variables
  real, dimension(npft,maxgrid) :: outaCalclm
  real, dimension(npft,maxgrid) :: outaNalclm
  real, dimension(npft,maxgrid) :: outaCalcrm
  real, dimension(npft,maxgrid) :: outaNalcrm

contains

  subroutine allocation_annual( tree )
    !//////////////////////////////////////////////////////////////////
    ! Finds optimal shoot:root growth ratio to balance C:N stoichiometry
    ! of a grass (no wood allocation).
    !------------------------------------------------------------------
    use md_classdefs
    use md_plant, only: params_plant, params_pft_plant, plant_type, &
      drgrow
    use md_params_core, only: pi, eps

    ! arguments
    type( plant_type ), dimension(npft), intent(inout) :: tree 

    ! local variables
    type( plant_type ) :: before    ! tree before allocation
    real               :: diam_inc  ! diameter increment (m)
    real               :: dWs_per_dD
    real               :: dHP_per_dD
    real               :: calloc    ! C to be allocated to new growth (gC)
    real               :: dclabl
    integer            :: lu
    integer            :: pft

    do pft=1,npft

      lu = params_pft_plant(pft)%lu_category

      ! save tree instance before allocation
      before = tree(pft)

      !------------------------------------------------------------------------
      ! Get allocatable C
      !------------------------------------------------------------------------
      calloc = params_plant%growtheff * tree(pft)%plabl%c%c12 * 1e-3   ! conversion from gC to kgC 

      !------------------------------------------------------------------------
      ! Get geometric relationships given current state variables
      !------------------------------------------------------------------------
      dWs_per_dD = get_dWs_per_dD( tree(pft) )
      dHP_per_dD = get_dHP_per_dD( tree(pft) )

      !------------------------------------------------------------------------
      ! Calculate change in diameter, stem mass and height
      !------------------------------------------------------------------------
      ! change in diameter
      diam_inc = calloc / ( dHP_per_dD + dWs_per_dD )

      !------------------------------------------------------------------------
      ! Update state variables
      !------------------------------------------------------------------------
      call update_tree( tree(pft), diam_inc )

      !------------------------------------------------------------------------
      ! Calculate change in C and N pool sizes
      !------------------------------------------------------------------------
      dcleaf(pft) = tree(pft)%pleaf%c%c12 - before%pleaf%c%c12
      dcroot(pft) = tree(pft)%proot%c%c12 - before%proot%c%c12
      dcsapw(pft) = tree(pft)%psapw%c%c12 - before%psapw%c%c12
      dcwood(pft) = tree(pft)%pwood%c%c12 - before%pwood%c%c12

      dnleaf(pft) = tree(pft)%pleaf%n%n14 - before%pleaf%n%n14
      dnroot(pft) = tree(pft)%proot%n%n14 - before%proot%n%n14
      dnsapw(pft) = tree(pft)%psapw%n%n14 - before%psapw%n%n14
      dnwood(pft) = tree(pft)%pwood%n%n14 - before%pwood%n%n14

      ! test if it makes sense
      if ( abs((dWs_per_dD * diam_inc) - dcwood(pft)) > eps ) stop 'stem mass increment not equal to change in wood C'
      if ( abs((dHP_per_dD * diam_inc) - (dcleaf(pft) + dcroot(pft))) > eps ) stop 'actual root and foliage growth not equal to change in sum of pools'

      !------------------------------------------------------------------------
      ! Growth respiration
      !------------------------------------------------------------------------
      drgrow(pft) = ( 1.0 / params_plant%growtheff - 1.0 ) * calloc

      !------------------------------------------------------------------------
      ! Depletion of labile pool (ignoring N balance)
      !------------------------------------------------------------------------
      dclabl = drgrow(pft) + dcleaf(pft) + dcroot(pft) + dcwood(pft)

      ! labile pool should be depleted to zero
      if ( abs( tree(pft)%plabl%c%c12 - dclabl ) > eps ) stop 'did not deplete labile N pool'

      tree(pft)%plabl%c%c12 = 0.0

      ! xxx try: crude - N is supplied by god.
      tree(pft)%plabl%n%n14 = 0.0

    end do

  end subroutine allocation_annual


  subroutine update_tree( tree, diam_inc )
    !/////////////////////////////////////////////////////////////////////////
    ! Updates tree pools, given diameter increment following geometry rules by
    ! T-model.
    !-------------------------------------------------------------------------
    use md_plant, only: params_pft_plant, plant_type
    use md_params_core, only: pi

    ! arguments
    type( plant_type ), intent(inout) :: tree
    real, intent(in)                  :: diam_inc

    !-------------------------------------------------------------------------
    ! Geometric relationships from T-model
    !-------------------------------------------------------------------------
    ! Change in diameter, determined from previous year; is zero in first year
    tree%diam = tree%diam + diam_inc

    ! update height (Eq. 4)
    tree%height = params_tmodel(tree%pftno)%maxheight * ( 1.0 - exp( -1.0 * params_tmodel(tree%pftno)%a_par * tree%diam / params_tmodel(tree%pftno)%maxheight ) )

    ! update crown fraction ("thickness of crown as a fraction of total height", Eq. 11)
    tree%fcrown = tree%height / ( params_tmodel(tree%pftno)%a_par * tree%diam ) 

    ! update projected crown area (Eq. 8)
    tree%acrown = ( ( pi * params_tmodel(tree%pftno)%cr_par ) / ( 4.0 * params_tmodel(tree%pftno)%a_par) ) * tree%diam * tree%height

    ! update total stem mass (Eq. 6)
    tree%pwood%c%c12 = ( pi / 8.0 ) * ( tree%diam**2.0 ) * tree%height * params_tmodel(tree%pftno)%rho

    ! update foliage mass (ss is specific leaf area)
    tree%pleaf%c%c12 = tree%acrown * params_pft_plant(tree%pftno)%lai_ind / params_pft_plant(tree%pftno)%lma

    ! update sapwood mass (tree%acrown/cr = As; As:= sapwood cross-sectional area; Eq. 14; L*cr*v=1; Hf=H(1-tree%fcrown/2))
    tree%psapw%c%c12 = tree%acrown * params_tmodel(tree%pftno)%rho * tree%height * ( 1.0 - tree%fcrown / 2.0 ) / params_tmodel(tree%pftno)%cr_par

    ! update root mass
    tree%proot%c%c12 = params_tmodel(tree%pftno)%z_par * params_pft_plant(tree%pftno)%lma * tree%pleaf%c%c12

    !-------------------------------------------------------------------------
    ! Update N in pools with fixed or given stoichiometry
    !-------------------------------------------------------------------------
    tree%pleaf%n%n14 = tree%r_ntoc_leaf
    tree%proot%n%n14 = params_pft_plant(tree%pftno)%r_ntoc_root
    tree%psapw%n%n14 = params_pft_plant(tree%pftno)%r_ntoc_wood
    tree%pwood%n%n14 = params_pft_plant(tree%pftno)%r_ntoc_wood


  end subroutine update_tree


  function get_dWs_per_dD( tree ) result( dWs_per_dD )
    !///////////////////////////////////////////////////////////////////
    ! Returns stem mass increment proportionality to diameter increment 
    ! Based on Eq. 7, using Eq. 5 in Li et al., 2014
    !-------------------------------------------------------------------
    use md_params_core, only: pi

    ! arguments
    type( plant_type ), intent(in) :: tree

    ! function return variable
    real :: dWs_per_dD

    dWs_per_dD = pi / 8.0 * params_tmodel(tree%pftno)%rho * tree%diam &
      * ( params_tmodel(tree%pftno)%a_par * tree%diam * (1.0 - ( tree%height / params_tmodel(tree%pftno)%maxheight) ) + 2.0 * tree%height )

  end function get_dWs_per_dD
    

  function get_dHP_per_dD( tree ) result( dHP_per_dD )
    !///////////////////////////////////////////////////////////////////
    ! Returns "root and foliage allocation" per unit diameter increment.
    ! Corresponds to second term in Eq. 16 without dD/dt in Li et al., 2014
    !-------------------------------------------------------------------
    use md_params_core, only: pi

    ! arguments
    type( plant_type ), intent(in) :: tree

    ! function return variable
    real :: dHP_per_dD

    dHP_per_dD = params_pft_plant(tree%pftno)%lai_ind * ( ( pi * params_tmodel(tree%pftno)%cr_par) / ( 4.0 * params_tmodel(tree%pftno)%a_par) ) &
      * ( params_tmodel(tree%pftno)%a_par * tree%diam * ( 1.0 - ( tree%height / params_tmodel(tree%pftno)%maxheight ) + tree%height ) )      &
      * ( 1.0 / params_tmodel(tree%pftno)%lma + params_tmodel(tree%pftno)%z_par ) 
  
  end function get_dHP_per_dD


  subroutine getpar_tree()
    !////////////////////////////////////////////////////////////////
    !  Subroutine reads model parameters from input file.
    !  It was necessary to separate this SR from module md_plant
    !  because this SR uses module md_waterbal, which also uses
    !  _plant.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------    
    use md_sofunutils, only: getparreal
    use md_interface

    ! local variables
    integer :: pft
    integer :: npft_site

    !----------------------------------------------------------------
    ! PFT DEPENDENT PARAMETERS
    ! read parameter input file and store values in single array
    ! important: Keep this order of reading PFT parameters fixed.
    !----------------------------------------------------------------
    do pft=1,npft
      params_tmodel(pft) = getpftparams( params_pft_plant(pft)%pftname )
    end do

  end subroutine getpar_tree

  function getpftparams( pftname ) result( out_getpar )
    !----------------------------------------------------------------
    ! Read PFT parameters from respective file, given the PFT name
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    character(len=*), intent(in) :: pftname

    ! function return variable
    type( params_tmodel_type ) :: out_getpar

    out_getpar%a_par     = getparreal( trim('params/params_tmodel_'//pftname//'.dat'), 'a_par    ' ) ! initial slope of the relationship between height and diameter
    out_getpar%cr_par    = getparreal( trim('params/params_tmodel_'//pftname//'.dat'), 'cr_par   ' ) ! initial ratio of crown area to stem cross - sectional area ('c')
    out_getpar%maxheight = getparreal( trim('params/params_tmodel_'//pftname//'.dat'), 'maxheight' ) ! asymptotic maximum height
    out_getpar%rho       = getparreal( trim('params/params_tmodel_'//pftname//'.dat'), 'rho      ' ) ! density of wood
    out_getpar%z_par     = getparreal( trim('params/params_tmodel_'//pftname//'.dat'), 'z_par    ' ) ! Ratio of fine-root mass to foliage area (0.17 kg C m-2; White et al. (2000))

  end function getpftparams


  subroutine initio_allocation()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use md_interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------
    if (interface%params_siml%loutalloc) then

      ! C ALLOCATED TO LEAF GROWTH 
      filnam=trim(prefix)//'.a.calclm.out'
      open(350,file=filnam,err=999,status='unknown')

      ! N ALLOCATED TO LEAF GROWTH 
      filnam=trim(prefix)//'.a.nalclm.out'
      open(351,file=filnam,err=999,status='unknown')

      ! C ALLOCATED TO ROOT GROWTH 
      filnam=trim(prefix)//'.a.calcrm.out'
      open(352,file=filnam,err=999,status='unknown')

      ! N ALLOCATED TO ROOT GROWTH 
      filnam=trim(prefix)//'.a.nalcrm.out'
      open(353,file=filnam,err=999,status='unknown')

    end if

    return

    999  stop 'INITIO_ALLOCATION: error opening output files'

  end subroutine initio_allocation


  subroutine initoutput_allocation()
    !////////////////////////////////////////////////////////////////
    !  Initialises nuptake-specific output variables
    !----------------------------------------------------------------
    use md_interface

    ! xxx remove their day-dimension
    if (interface%params_siml%loutalloc) then
      outaCalclm(:,:) = 0.0
      outaNalclm(:,:) = 0.0
      outaCalcrm(:,:) = 0.0
      outaNalcrm(:,:) = 0.0
    end if
    
    ! print*, 'initialising outaCalloc',outaCalloc

  end subroutine initoutput_allocation


  subroutine getout_daily_allocation( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    if (interface%params_siml%loutalloc) then
      outaCalclm(:,jpngr) = outaCalclm(:,jpngr) + dcleaf(:) 
      outaNalclm(:,jpngr) = outaNalclm(:,jpngr) + dnleaf(:)
      outaCalcrm(:,jpngr) = outaCalcrm(:,jpngr) + dcroot(:) 
      outaNalcrm(:,jpngr) = outaNalcrm(:,jpngr) + dnroot(:)
    end if

    ! print*, 'collecting outaCalloc',outaCalloc

  end subroutine getout_daily_allocation


  subroutine writeout_ascii_allocation( year )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: year       ! simulation year

    ! local variables
    real    :: itime
    integer :: jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutalloc) then

      itime = real(interface%steering%outyear)

      write(350,999) itime, sum(outaCalclm(:,jpngr))
      write(351,999) itime, sum(outaNalclm(:,jpngr))
      write(352,999) itime, sum(outaCalcrm(:,jpngr))
      write(353,999) itime, sum(outaNalcrm(:,jpngr))
    end if

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_allocation

end module md_allocation