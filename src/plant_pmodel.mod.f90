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
  public plant_type, dgpp, getpar_modl_plant, params_pft_plant,              &
    initdaily_plant, initoutput_plant, initio_plant, getout_daily_plant,     &
    writeout_ascii_plant, maxdoy, initglobal_plant, get_leaftraits, initpft, &
    getout_annual_plant, initio_nc_plant, writeout_nc_plant

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  type plant_type

    ! PFT index that goes along with this instance of 'plant'
    integer :: pftno

    ! canopy
    real :: lai_ind             ! fraction of absorbed photosynthetically active radiation
    real :: fapar_ind           ! fraction of absorbed photosynthetically active radiation
    real :: acrown              ! crown area

    ! leaf traits
    real :: narea               ! total leaf N per unit leaf area (gN m-2)
    real :: narea_metabolic     ! metabolic leaf N per unit leaf area (gN m-2)
    real :: narea_structural    ! structural leaf N per unit leaf area (gN m-2)
    real :: lma                 ! leaf mass per area (gC m-2)
    real :: sla                 ! specific leaf area (m2 gC-1)
    real :: nmass               ! leaf N per unit leaf mass, g N / g-dry mass
    real :: r_cton_leaf         ! leaf C:N ratio [gC/gN] 
    real :: r_ntoc_leaf         ! leaf N:C ratio [gN/gC]

  end type plant_type

  ! fluxes
  real, dimension(npft) :: dgpp             ! daily gross primary production [gC/m2/d]

  !-----------------------------------------------------------------------
  ! Parameters. Runtime read-in
  !-----------------------------------------------------------------------
  ! NON PFT-DEPENDENT PARAMETERS
  type params_plant_type
    real :: kbeer             ! canopy light extinction coefficient
  end type params_plant_type

  type( params_plant_type ) :: params_plant

  ! PFT-DEPENDENT PARAMETERS
  type params_pft_plant_type
    character(len=4) :: pftname    ! standard PFT name with 4 characters length
    integer :: lu_category         ! land use category associated with PFT
    logical, dimension(nlu) :: islu! islu(ipft,ilu) is true if ipft belongs to ilu
    logical :: grass               ! boolean for growth form 'grass'
    logical :: tree                ! boolean for growth form 'tree'
    logical :: nfixer              ! whether plant is capable of symbiotically fixing N
    logical :: c3                  ! whether plant follows C3 photosynthesis
    logical :: c4                  ! whether plant follows C4 photosynthesis
    real    :: sla                 ! specific leaf area (m2 gC-1)
    real    :: lma                 ! leaf mass per area (gC m-2)
    real    :: r_ntolma            ! constant ratio of structural N to C (LMA) (gN/gC)
  end type params_pft_plant_type

  type( params_pft_plant_type ), dimension(npft) :: params_pft_plant

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, allocatable, dimension(:,:,:) :: outdgpp    ! daily gross primary production [gC/m2/d]

  ! annual
  real, dimension(:,:), allocatable :: outagpp
  real, dimension(:,:), allocatable :: outanarea_mb
  real, dimension(:,:), allocatable :: outanarea_cw
  real, dimension(:,:), allocatable :: outalai
  real, dimension(:,:), allocatable :: outalma
  real, dimension(:,:), allocatable :: outacton_lm

  ! required for outputting leaf trait variables in other modules
  integer, dimension(npft) :: maxdoy  ! DOY of maximum LAI

  !----------------------------------------------------------------
  ! Module-specific NetCDF output file and variable names
  !----------------------------------------------------------------
  character(len=256) :: ncoutfilnam_gpp
  character(len=*), parameter :: GPP_NAME="gpp"

contains

  function get_fapar( lai ) result( fapar )
    !////////////////////////////////////////////////////////////////
    ! FOLIAGE PROJECTIVE COVER 
    ! = Fraction of Absorbed Photosynthetically Active Radiation
    ! Function returns fractional plant cover an individual
    ! Eq. 7 in Sitch et al., 2003
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: lai

    ! function return variable
    real :: fapar

    fapar = ( 1.0 - exp( -1.0 * params_plant%kbeer * lai) )

  end function get_fapar


  function get_leaf_n_metabolic_canopy( mylai, meanmppfd, nv, myfapar ) result( mynleaf_metabolic )
    !////////////////////////////////////////////////////////////////
    ! Calculates metabolic leaf N at canopy-level, determined by 
    ! light conditions (meanmppfd) and the Rubisco-N per unit absorbed
    ! light.
    !----------------------------------------------------------------
    use md_params_core, only: nmonth

    ! arguments
    real, intent(in)                    :: mylai
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv
    real, intent(in), optional          :: myfapar

    ! function return variable
    real :: mynleaf_metabolic  ! mol N m-2-ground

    ! local variables
    real :: maxnv

    ! Metabolic N is predicted and is optimised at a monthly time scale. 
    ! Leaf traits are calculated based on metabolic N => cellwall N => cellwall C / LMA
    ! Leaves get thinner at the bottom of the canopy => increasing LAI through the season comes at a declining C and N cost
    ! Monthly variations in metabolic N, determined by variations in meanmppfd and nv should not result in variations in leaf traits. 
    ! In order to prevent this, assume annual maximum metabolic N, part of which is deactivated during months with lower insolation (and Rd reduced.)
    maxnv = maxval( meanmppfd(:) * nv(:) )

    if (present(myfapar)) then
      mynleaf_metabolic = maxnv * myfapar
    else
      mynleaf_metabolic = maxnv * get_fapar( mylai )
    end if

  end function get_leaf_n_metabolic_canopy


  subroutine get_leaftraits( plant, meanmppfd, nv )
    !////////////////////////////////////////////////////////////////
    ! Calculates leaf traits based on (predicted) metabolic Narea and
    ! (prescribed) parameters that relate structural to metabolic
    ! Narea and Carea to structural Narea:
    ! Narea_metabolic  = predicted
    ! Narea_structural = rN:C_struct * LMA
    !----------------------------------------------------------------
    use md_params_core, only: c_content_of_biomass, nmonth, n_molmass, c_molmass

    ! arguments
    type( plant_type ), intent(inout)   :: plant
    real, dimension(nmonth), intent(in) :: meanmppfd
    real, dimension(nmonth), intent(in) :: nv

    ! local variables
    real :: narea_metabolic_canopy   ! g N m-2-ground

    ! canopy-level, in units of gN / m2-ground 
    narea_metabolic_canopy  = n_molmass * get_leaf_n_metabolic_canopy(  -9999.9, meanmppfd(:), nv(:), plant%fapar_ind )

    ! leaf-level, in units of gN / m2-leaf 
    ! assume narea_metabolic is representative of the outer canopy, therefore divide by 1.0 (or just leave)
    plant%narea_metabolic  = narea_metabolic_canopy / 1.0
    plant%narea_structural = params_pft_plant(plant%pftno)%r_ntolma * params_pft_plant(plant%pftno)%lma
    plant%narea            = plant%narea_metabolic + plant%narea_structural
    plant%lma              = params_pft_plant(plant%pftno)%lma

    ! additional traits
    plant%nmass            = plant%narea / ( plant%lma / c_content_of_biomass )
    plant%r_cton_leaf      = params_pft_plant(plant%pftno)%lma / plant%narea
    plant%r_ntoc_leaf      = 1.0 / plant%r_cton_leaf

  end subroutine get_leaftraits


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
    ! NON-PFT DEPENDENT PARAMETERS
    !----------------------------------------------------------------
    ! canopy light extinction coefficient for Beer's Law
    params_plant%kbeer = getparreal( 'params/params_plant.dat', 'kbeer' )

    !----------------------------------------------------------------
    ! PFT DEPENDENT PARAMETERS
    ! read parameter input file and store values in single array
    ! important: Keep this order of reading PFT parameters fixed.
    !----------------------------------------------------------------
    pft = 0
    if ( interface%params_siml%lTeBE ) then
      pft = pft + 1
      params_pft_plant(pft) = getpftparams( 'TeBE' )

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

    ! leaf mass per area (gC m-2)
    out_getpftparams%lma = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'lma' )
    out_getpftparams%sla = 1.0 / out_getpftparams%lma

    ! constant ratio of leaf structural N to LMA
    out_getpftparams%r_ntolma = getparreal( trim('params/params_plant_'//pftname//'.dat'), 'r_ntolma' )

  end function getpftparams


  subroutine initglobal_plant( plant, ngridcells )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of all _pools on all gridcells at the beginning
    !  of the simulation.
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_params_core, only: npft

    ! argument
    type( plant_type ), dimension(npft,ngridcells), intent(inout) :: plant
    integer, intent(in) :: ngridcells

    ! local variables
    integer :: pft
    integer :: jpngr

    !-----------------------------------------------------------------------------
    ! derive which PFTs are present from fpc_grid (which is prescribed)
    !-----------------------------------------------------------------------------
    do jpngr=1,ngridcells
      do pft=1,npft
        call initpft( plant(pft,jpngr) )
        plant(pft,jpngr)%pftno = pft
      end do
    end do

  end subroutine initglobal_plant


  subroutine initpft( plant )
    !////////////////////////////////////////////////////////////////
    !  Initialisation of specified PFT on specified gridcell
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    ! argument
    type( plant_type ), intent(inout) :: plant

    plant%acrown = 1.0 ! dummy

    ! canpopy state variables
    plant%narea            = 0.0
    plant%narea_metabolic  = 0.0
    plant%narea_structural = 0.0
    plant%lma              = 0.0
    plant%sla              = 0.0
    plant%nmass            = 0.0
    plant%r_cton_leaf      = 0.0
    plant%r_ntoc_leaf      = 0.0

  end subroutine initpft


  subroutine initdaily_plant()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    !----------------------------------------------------------------
    dgpp(:) = 0.0

  end subroutine initdaily_plant


  subroutine initoutput_plant( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    use md_interface, only: interface

    ! arguments
    integer, intent(in) :: ngridcells

    if ( interface%steering%init .and. interface%params_siml%loutdgpp ) allocate( outdgpp(npft,ndayyear,ngridcells) )
    if ( interface%params_siml%loutdgpp ) outdgpp  (:,:,:) = 0.0
    
    ! annual output variables
    if (interface%params_siml%loutplant) then

      if (interface%steering%init) then
        allocate( outagpp(npft,ngridcells) )
        allocate( outanarea_mb(npft,ngridcells) )
        allocate( outanarea_cw(npft,ngridcells) )
        allocate( outalai(npft,ngridcells) )
        allocate( outalma(npft,ngridcells) )
        allocate( outacton_lm(npft,ngridcells) )

        outagpp(:,:)      = 0.0
        outanarea_mb(:,:) = 0.0
        outanarea_cw(:,:) = 0.0
        outalai     (:,:) = 0.0
        outalma     (:,:) = 0.0
        outacton_lm (:,:) = 0.0
      end if

    end if

  end subroutine initoutput_plant


  subroutine initio_plant()
    !////////////////////////////////////////////////////////////////
    ! Opens input/output files.
    !----------------------------------------------------------------
    use md_interface, only: interface

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
      print*,'filnam ', filnam
      open(101,file=filnam,err=999,status='unknown')
    end if 

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------
    if (interface%params_siml%loutplant) then

      ! GPP 
      filnam=trim(prefix)//'.a.gpp.out'
      open(310,file=filnam,err=999,status='unknown')

      ! METABOLIC NAREA (AT ANNUAL LAI MAXIMUM)
      filnam=trim(prefix)//'.a.narea_mb.out'
      open(319,file=filnam,err=999,status='unknown')

      ! CELL WALL NAREA (AT ANNUAL LAI MAXIMUM)
      filnam=trim(prefix)//'.a.narea_cw.out'
      open(320,file=filnam,err=999,status='unknown')

      ! LEAF C:N RATIO (AT ANNUAL LAI MAXIMUM)
      filnam=trim(prefix)//'.a.cton_lm.out'
      open(321,file=filnam,err=999,status='unknown')

      ! LMA (AT ANNUAL LAI MAXIMUM)
      filnam=trim(prefix)//'.a.lma.out'
      open(322,file=filnam,err=999,status='unknown')

    end if

    return

    999  stop 'INITIO: error opening output files'

  end subroutine initio_plant


  subroutine initio_nc_plant()
    !////////////////////////////////////////////////////////////////
    ! Opens NetCDF output files.
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: init_nc_3D, check
    use md_interface, only: interface

    ! local variables
    character(len=256) :: prefix

    character(len=*), parameter :: TITLE = "SOFUN GP-model output, module md_plant"
    character(len=4) :: year_char

    integer :: jpngr, doy

    write(year_char,999) interface%steering%outyear

    prefix = "./output_nc/"//trim(interface%params_siml%runname)

    if (interface%params_siml%lncoutdgpp) then

      ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
      ! overwrite this file, if it already exists.
      print*,'initialising gpp NetCDF file ...'
      ncoutfilnam_gpp = trim(prefix)//'.'//year_char//".d.gpp.nc"
      call init_nc_3D( filnam  = ncoutfilnam_gpp, &
                      nlon     = interface%domaininfo%nlon, &
                      nlat     = interface%domaininfo%nlat, &
                      lon      = interface%domaininfo%lon, &
                      lat      = interface%domaininfo%lat, &
                      outyear  = interface%steering%outyear, &
                      varnam   = GPP_NAME, &
                      varunits = "gC m-2 d-1", &
                      longnam  = "daily gross primary productivivty", &
                      title    = TITLE &
                      )

    end if

    999  format (I4.4)
    
  end subroutine initio_nc_plant


  subroutine getout_daily_plant( plant, jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft
    use md_interface, only: interface

    ! arguments
    type( plant_type ), dimension(npft), intent(in) :: plant
    integer, intent(in)                             :: jpngr
    integer, intent(in)                             :: moy
    integer, intent(in) :: doy

    ! LOCAL VARIABLES
    integer :: pft

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    if ( interface%params_siml%loutdgpp ) outdgpp(:,doy,jpngr) = dgpp(:)

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    if (interface%params_siml%loutplant) then
      outagpp(:,jpngr)     = outagpp(:,jpngr) + dgpp(:)
    end if

  end subroutine getout_daily_plant


  subroutine getout_annual_plant( plant, jpngr )
    !////////////////////////////////////////////////////////////////
    !  SR called once a year to gather annual output variables.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft
    use md_interface, only: interface

    ! arguments
    type( plant_type ), dimension(npft), intent(in) :: plant
    integer, intent(in)                             :: jpngr

    ! local variables
    integer :: pft
    integer :: doy

    ! Output annual value at day of peak LAI
    if (interface%params_siml%loutplant) then
      outanarea_mb(:,jpngr) = plant(:)%narea_metabolic 
      outanarea_cw(:,jpngr) = plant(:)%narea_structural
      outalai     (:,jpngr) = plant(:)%lai_ind
      outalma     (:,jpngr) = plant(:)%lma
      outacton_lm (:,jpngr) = plant(:)%r_cton_leaf
    end if

  end subroutine getout_annual_plant


  subroutine writeout_ascii_plant()
    !/////////////////////////////////////////////////////////////////////////
    ! Write daily ASCII output
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !-------------------------------------------------------------------------
    use md_params_core, only: ndayyear
    use md_interface, only: interface

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
      write(319,999) itime, sum(outanarea_mb(:,jpngr))
      write(320,999) itime, sum(outanarea_cw(:,jpngr))
      write(321,999) itime, sum(outacton_lm(:,jpngr))
      write(322,999) itime, sum(outalma(:,jpngr))

    end if

    return

    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_plant


  subroutine writeout_nc_plant()
    !/////////////////////////////////////////////////////////////////////////
    ! Write NetCDF output
    !-------------------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: write_nc_3D, check
    use md_interface, only: interface

    ! local variables
    integer :: doy, jpngr
    integer :: ncid
    integer :: varid_temp

    real, dimension(:,:,:,:), allocatable :: outarr

    if (npft>1) stop 'writeout_nc_plant(): npft > 1. Think of something clever.'

    ! if ( .not. interface%steering%spinup &
    !       .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
    !       .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

      !-------------------------------------------------------------------------
      ! fapar
      !-------------------------------------------------------------------------
      print*,'writing gpp NetCDF file ...'
      if (interface%params_siml%lncoutdgpp) call write_nc_3D( ncoutfilnam_gpp, &
                                                              GPP_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              ndayyear, &
                                                              outdgpp(1,:,:) &
                                                              )

    ! end if

  end subroutine writeout_nc_plant


end module md_plant
