module md_plant
  !////////////////////////////////////////////////////////////////
  !  Module contains (constrainable) model parameters.
  !  Model parameters adopted here are from LPX C3 grass PFT
  !  Litter and soil turnover parameters are divided by 365 to 
  !  convert from [1/yr] to [1/d].
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core

  implicit none

  private
  public dgpp, lai_ind, canopy, getpar_modl_plant, params_pft_plant,         &
    initdaily_plant, initoutput_plant, initio_plant, getout_daily_plant,     &
    getout_annual_plant, writeout_ascii_plant, maxdoy

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  ! fluxes
  real, dimension(npft) :: dgpp             ! daily gross primary production [gC/m2/d]

  ! Canopy state variables
  real, dimension(npft,maxgrid) :: lai_ind

  type canopy_type
    real :: fapar_ind
  end type canopy_type

  type( canopy_type ), dimension(npft) :: canopy

  !-----------------------------------------------------------------------
  ! Parameters. Runtime read-in
  !-----------------------------------------------------------------------
  type params_pft_plant_type
    character(len=4) :: pftname    ! standard PFT name with 4 characters length
    integer :: lu_category         ! land use category associated with PFT
    logical, dimension(nlu) :: islu! islu(ipft,ilu) is true if ipft belongs to ilu
    logical :: grass               ! boolean for growth form 'grass'
    logical :: tree                ! boolean for growth form 'tree'
    logical :: nfixer              ! whether plant is capable of symbiotically fixing N
    logical :: c3                  ! whether plant follows C3 photosynthesis
    logical :: c4                  ! whether plant follows C4 photosynthesis
  end type params_pft_plant_type

  type( params_pft_plant_type ), dimension(npft) :: params_pft_plant

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, allocatable, dimension(:,:,:) :: outdgpp    ! daily gross primary production [gC/m2/d]
  real, allocatable, dimension(:,:,:) :: outdlai

  ! annual
  real, dimension(npft,maxgrid) :: outagpp

  ! required for outputting leaf trait variables in other modules
  integer, dimension(npft) :: maxdoy  ! DOY of maximum LAI

contains

  subroutine getpar_modl_plant()
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
    pft = 0
    if ( interface%params_siml%lTeBS ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TeBS' )

    else if ( interface%params_siml%lGrC3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GrC3' )

    else if ( interface%params_siml%lGNC3 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GNC3' )

    else if ( interface%params_siml%lGrC4 ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'GrC4' )

    else
      stop 'PLANT:GETPAR_MODL_PLANT: PFT name not valid. See run/<simulationname>.sofun.parameter'
    end if
    npft_site = pft

  end subroutine getpar_modl_plant


  function getpftparams( pftname ) result( out_getpftparams )
    !----------------------------------------------------------------
    ! Read PFT parameters from respective file, given the PFT name
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    character(len=*), intent(in) :: pftname

    ! local variables
    real :: lu_category_prov    ! land use category associated with PFT (provisional)
    real :: code_growthform
    real :: code_nfixer

    ! function return variable
    type( params_pft_plant_type ) :: out_getpftparams

    ! standard PFT name
    out_getpftparams%pftname = pftname

    ! PFT names
    ! GrC3 : C3 grass                          
    ! GrC4 : C4 grass     
    if (trim(pftname)=='GrC3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .false.
    else if (trim(pftname)=='GNC3') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .true.
      out_getpftparams%c4      = .false.
      out_getpftparams%nfixer  = .true.
    else if (trim(pftname)=='GrC4') then
      out_getpftparams%grass   = .true.
      out_getpftparams%tree    = .false.
      out_getpftparams%c3      = .false.
      out_getpftparams%c4      = .true.
      out_getpftparams%nfixer  = .false.
    end if      

    ! land use category associated with PFT (provisional) 
    lu_category_prov = getparreal( trim('params/params_plant_'//trim(pftname)//'.dat'), 'lu_category_prov' )
    if (lu_category_prov==1.0) then
      out_getpftparams%lu_category = lunat
      out_getpftparams%islu(lunat) = .true.
    else
      out_getpftparams%islu(lunat) = .false.
    end if

  end function getpftparams


  subroutine initdaily_plant()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    !----------------------------------------------------------------
    dgpp(:) = 0.0

  end subroutine initdaily_plant


  subroutine initoutput_plant()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_interface

    if (interface%steering%init .and. interface%params_siml%loutdgpp ) allocate( outdgpp(npft,ndayyear,maxgrid) )
    outdgpp  (:,:,:) = 0.0
    
    ! this is needed also for other (annual) output variables
    allocate( outdlai(npft,ndayyear,maxgrid) )
    outdlai  (:,:,:) = 0.0
 
    ! annual output variables
    if (interface%params_siml%loutplant) then
      outagpp(:,:) = 0.0
    end if

  end subroutine initoutput_plant


  subroutine initio_plant()
    !////////////////////////////////////////////////////////////////
    ! Opens input/output files.
    !----------------------------------------------------------------
    use md_interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !////////////////////////////////////////////////////////////////
    ! DAILY OUTPUT: OPEN ASCII OUTPUT FILES 
    !----------------------------------------------------------------
    ! GPP
    if (interface%params_siml%loutdgpp) then
      filnam=trim(prefix)//'.d.gpp.out'
      open(101,file=filnam,err=999,status='unknown')
    end if 

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------
    if (interface%params_siml%loutplant) then

      ! GPP 
      filnam=trim(prefix)//'.a.gpp.out'
      open(310,file=filnam,err=999,status='unknown')

    end if

    return

    999  stop 'INITIO: error opening output files'

  end subroutine initio_plant


  subroutine getout_daily_plant( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    ! LOCAL VARIABLES
    integer :: pft

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    if (interface%params_siml%loutdgpp ) outdgpp(:,doy,jpngr) = dgpp(:)

    ! this is needed also for other (annual) output variables
    outdlai(:,doy,jpngr) = lai_ind(:,jpngr)
      
    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    if (interface%params_siml%loutplant) then
      outagpp(:,jpngr)     = outagpp(:,jpngr) + dgpp(:)
    end if

  end subroutine getout_daily_plant


  subroutine getout_annual_plant( jpngr )
    !////////////////////////////////////////////////////////////////
    !  SR called once a year to gather annual output variables.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr

    ! local variables
    integer :: pft
    integer :: doy
    integer, dimension(npft) :: maxlai

    maxlai(:) = 0.0
    maxdoy(:) = 1

    ! Output annual value at day of peak LAI
    do pft=1,npft
      maxdoy(pft) = 1
      maxlai(pft) = outdlai(pft,1,jpngr)
      do doy=2,ndayyear
        if ( outdlai(pft,doy,jpngr) > maxlai(pft) ) then
          maxlai(pft) = outdlai(pft,doy,jpngr)
          maxdoy(pft) = doy
        end if
      end do
    end do

  end subroutine getout_annual_plant


  subroutine writeout_ascii_plant( year )
    !/////////////////////////////////////////////////////////////////////////
    ! Write daily ASCII output
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !-------------------------------------------------------------------------
    ! use md_params_siml, only: spinup, interface%params_siml%daily_out_startyr, &
    !   interface%params_siml%daily_out_endyr, outyear
    use md_params_core, only: ndayyear
    use md_interface

    ! arguments
    integer, intent(in) :: year       ! simulation year

    ! local variables
    real :: itime
    integer :: day, moy, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    ! Write daily value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if ( .not. interface%steering%spinup &
      .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
      .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      ! Write daily output only during transient simulation
      do day=1,ndayyear

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(interface%steering%outyear) + real(day-1)/real(ndayyear)
        
        if (interface%params_siml%loutdgpp  ) write(101,999) itime, sum(outdgpp(:,day,jpngr))

      end do
    end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutplant) then

      itime = real(interface%steering%outyear)

      write(310,999) itime, sum(outagpp(:,jpngr))

    end if

    return

    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_plant


end module md_plant
