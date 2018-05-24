program main
  !////////////////////////////////////////////////////////////////
  !  Main program for site scale simulations, here used for 
  !  SOFUN (Seasonal Optimisation of Fixation and Uptake of 
  !  Nitrogen)
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_interface, only: interface, outtype_biosphere
  use md_params_siml, only: getpar_siml, getsteering
  use md_params_domain, only: getpar_domain, type_params_domain
  use md_grid, only: get_domaininfo, getgrid
  use md_params_soil, only: getsoil
  use md_forcing, only: get_fpc_grid, getclimate, getninput, ninput_type, gettot_ninput, getfapar, getlanduse, getco2
  use md_params_core, only: dummy, maxgrid, ndayyear, npft
  use md_biosphere, only: biosphere_annual
  use md_sofunutils, only: median

  implicit none

  ! local variables
  integer :: yr                                  ! simulation year
  character(len=245) :: runname
  integer, parameter :: maxlen_runname = 50      ! maximum length of runname (arbitrary)
  type( ninput_type ), dimension(maxgrid) :: nfert_field, ndep_field 
  type( type_params_domain ) :: params_domain
  logical, parameter :: verbose = .true.

  type(outtype_biosphere) :: out_biosphere       ! holds all the output used for calculating the cost or maximum likelihood function 
  real :: cost_annual = 0.0                      ! annual cost (model-observation fit after Choler et al., 2010 BG)
  real, allocatable, dimension(:) :: cost_bysite ! cost over multiple years for each site in simsuite
  real :: cost                                   ! overall cost as median across sites
  real :: k_decay_tissue                         ! parameter read from standard input
  character(len=20) :: simsuite

  ! For runname file reading
  integer :: ios
  integer, parameter :: read_unit = 99
  character(len=6), allocatable :: runname_list(:)
  character(len=6) :: line
  integer :: nruns, irun

  !----------------------------------------------------------------
  ! READ SIMULATION SUITE NAME AND PARAMETERS FROM STANDARD INPUT
  !----------------------------------------------------------------
  ! read (*,100) k_decay_tissue
  read (*,*) simsuite, k_decay_tissue
  print*,'PARAMETERS FOR CALIBRATION:'
  print*,'k_decay_tissue: ', k_decay_tissue

  ! translate parameters from standard input to appropriate derived type 
  interface%params_calib%k_decay_tissue = k_decay_tissue

  ! set parameter to define that this is a calibration run (no output, overriding parameter values)
  interface%params_siml%is_calib = .true.

  !----------------------------------------------------------------
  ! READ RUNNAMES FOR THIS SIMULATION SUITE
  !----------------------------------------------------------------
  open(unit=read_unit, file='run/runnames_'//trim(simsuite)//'.txt', iostat=ios)
  if ( ios /= 0 ) stop "Error opening file run/runnames_<simsuite>.txt"

  nruns = 0
  do
    read(read_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    nruns = nruns + 1
  end do

  print*, "Simsuite "//trim(simsuite)//" contains ", nruns, "simulations"
  allocate(runname_list(nruns))
  rewind(read_unit)

  do irun = 1, nruns
    read(read_unit, '(A)') runname_list(irun)
  end do

  close(read_unit)

  ! allocate variables that have site-dimension 
  allocate( cost_bysite(nruns) )

  !----------------------------------------------------------------
  ! LOOP THROUGH RUNS
  !----------------------------------------------------------------
  runloop: do irun = 1, nruns

    runname = runname_list(irun)
    
    ! make sure runname length is smaller/equal than maxlen_runname
    if (len_trim(runname)>=maxlen_runname) then
      stop 'runname too long'
    endif

    ! write simulation name to standard output (screen)
    print*, '------------SOFUN : '//trim(runname)//'-------------'

    !----------------------------------------------------------------
    ! GET SIMULATION PARAMETERS FROM FILE <runname>.sofun.parameter
    ! SR getpar_siml is defined in _params_siml.mod.F
    !----------------------------------------------------------------
    interface%params_siml = getpar_siml( trim(runname) )

    !----------------------------------------------------------------
    ! GET SITE PARAMETERS AND INPUT DATA
    ! site location (lon,lat), soil type, vegetation type
    ! SR getpar_site is defined in _params_site.mod.F. 
    ! 'sitename' is global variable
    !----------------------------------------------------------------
    params_domain = getpar_domain( trim(interface%params_siml%sitename) )

    !----------------------------------------------------------------
    ! GET GRID INFORMATION
    ! longitude, latitude, elevation
    !----------------------------------------------------------------
    ! temporarily get full grid information (full lon-lat-array)
    interface%domaininfo = get_domaininfo( params_domain )

    ! allocate variable size arrays
    if (.not.allocated(interface%grid))         allocate( interface%grid(         interface%domaininfo%maxgrid ) )
    if (.not.allocated(interface%climate))      allocate( interface%climate(      interface%domaininfo%maxgrid ) )
    if (.not.allocated(interface%ninput_field)) allocate( interface%ninput_field( interface%domaininfo%maxgrid ) )
    if (.not.allocated(interface%landuse))      allocate( interface%landuse(      interface%domaininfo%maxgrid ) )
    if (.not.allocated(interface%soilparams))   allocate( interface%soilparams(   interface%domaininfo%maxgrid ) )
    if (.not.allocated(interface%dfapar_field)) allocate( interface%dfapar_field( ndayyear, interface%domaininfo%maxgrid ) )
    if (.not.allocated(interface%fpc_grid))     allocate( interface%fpc_grid( npft, interface%domaininfo%maxgrid ) )

    ! vectorise 2D array, keeping only land gridcells
    interface%grid(:) = getgrid( interface%domaininfo, params_domain )

    ! Obtain land unit dependent parameters, define decomposition _rates
    !call luparameters

    !----------------------------------------------------------------
    ! GET SOIL PARAMETERS
    !----------------------------------------------------------------
    interface%soilparams(:) = getsoil( interface%domaininfo, interface%grid(:) )

    !----------------------------------------------------------------
    ! GET VEGETATION COVER (fractional projective cover by PFT)
    !----------------------------------------------------------------
    interface%fpc_grid(:,:) = get_fpc_grid( interface%domaininfo, interface%grid(:), interface%params_siml )


    ! LOOP THROUGH YEARS
    print*, '-------------------START OF SIMULATION--------------------'
    yearloop: do yr=1,interface%params_siml%runyears

      !----------------------------------------------------------------
      ! Define simulations "steering" variables (forcingyear, etc.)
      !----------------------------------------------------------------
      ! print*,'getting steering'
      interface%steering = getsteering( yr, interface%params_siml )

      if (yr == interface%params_siml%spinupyears+1 ) then
        print*, '------------------TRANSIENT SIMULATION--------------------'
      endif

      !----------------------------------------------------------------
      ! Get external (environmental) forcing
      !----------------------------------------------------------------
      ! Climate
      interface%climate(:) = getclimate( &
                                        interface%domaininfo, &
                                        interface%grid, &
                                        interface%steering%init, &
                                        interface%steering%climateyear, &
                                        interface%params_siml%in_ppfd,  &
                                        interface%params_siml%in_netrad &
                                        )

      ! CO2
      interface%pco2 = getco2( &
                              trim(runname), &
                              interface%domaininfo, &
                              interface%steering%forcingyear, &
                              interface%params_siml%const_co2_year, &
                              interface%params_siml%firstyeartrend,&
                              interface%params_siml%co2_forcing_file&
                              )

      ! Atmospheric N deposition
      ndep_field(:) = getninput( &
                                "ndep", &
                                trim(runname), &
                                interface%domaininfo, &
                                interface%steering%forcingyear, &
                                interface%params_siml%firstyeartrend, &
                                interface%params_siml%const_ndep_year, &
                                interface%params_siml%ndep_noy_forcing_file, &
                                interface%params_siml%ndep_nhx_forcing_file, &
                                interface%climate(:)&
                                )
      
      ! N fertiliser input
      nfert_field(:) = getninput( &
                                "nfert", &
                                trim(runname), &
                                interface%domaininfo, &
                                interface%steering%forcingyear, &
                                interface%params_siml%firstyeartrend, &
                                interface%params_siml%const_nfert_year, &
                                interface%params_siml%nfert_noy_forcing_file, &
                                interface%params_siml%nfert_nhx_forcing_file, &
                                interface%climate(:)&
                                )

      ! Interface holds only total reactive N input (N deposition + N fertiliser)                             
      interface%ninput_field(:) = gettot_ninput( nfert_field(:), ndep_field(:) )
                                   
      ! 'SOFUN: holding harvesting regime constant at 1993 level.'
      interface%landuse(:) = getlanduse( &
                                        trim(runname), &
                                        interface%domaininfo, &
                                        interface%steering%forcingyear, &
                                        interface%params_siml%do_grharvest_forcing_file, &
                                        interface%params_siml%const_lu_year, &
                                        interface%params_siml%firstyeartrend &
                                        )

      !----------------------------------------------------------------
      ! Get prescribed fAPAR if required (otherwise set to dummy value)
      !----------------------------------------------------------------
      interface%dfapar_field(:,:) = getfapar( &
                                            interface%domaininfo, &
                                            interface%grid, &
                                            interface%steering%forcingyear, &
                                            interface%params_siml%fapar_forcing_source &
                                            )    

      !----------------------------------------------------------------
      ! Call SR biosphere at an annual time step but with vectors 
      ! containing data for each day of this year.
      !----------------------------------------------------------------
      print*,'--------------------------------------------------------'
      print*,'Simulation year: ', interface%steering%year, ' - Real year: ', interface%steering%outyear
      print*,'--------------------------------------------------------'
      
      !----------------------------------------------------------------
      ! Call biosphere (wrapper for all modules, contains gridcell loop)
      !----------------------------------------------------------------
      out_biosphere = biosphere_annual() 
      !----------------------------------------------------------------

      ! calculate cost of mod-obs fit (sum annual costs)
      if (yr > interface%params_siml%spinupyears ) &
        cost_annual = cost_annual &
          + ( sum( abs( out_biosphere%fapar(:) - interface%dfapar_field(:,1) ) ) / ndayyear ) / ( sum( interface%dfapar_field(:,1) ) / ndayyear )

    enddo yearloop

    print*, '--------------END OF SIMULATION---------------'

    ! get final cost
    cost_bysite(irun) = cost_annual / interface%params_siml%nyeartrend

  end do runloop

  ! Get overall cost as median across sites
  cost = median( cost_bysite(:), nruns ) 
  print*,'cost_bysite: ', cost_bysite
  print*, cost

100  format(F10.7)

end program main

