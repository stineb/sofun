program main
  !////////////////////////////////////////////////////////////////
  !  Main program for site scale simulations, here used for 
  !  SOFUN (Seasonal Optimisation of Fixation and Uptake of 
  !  Nitrogen)
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_interface, only: interface
  use md_params_siml, only: getpar_siml, getsteering
  use md_params_domain, only: getpar_domain, type_params_domain
  use md_grid, only: get_domaininfo, getgrid
  use md_params_soil, only: getsoil
  use md_forcing, only: getclimate, getninput, ninput_type, gettot_ninput, &
    getfapar, getlanduse, getco2
  use md_params_core, only: dummy, maxgrid, ndayyear
  use md_biosphere, only: biosphere_annual

  implicit none

  ! local variables
  integer :: yr           ! simulation year
  real    :: c_uptake     ! annual net global C uptake by biosphere
  character(len=245) :: runname
  integer, parameter :: maxlen_runname = 50      ! maximum length of runname (arbitrary)
  type( ninput_type ), dimension(maxgrid) :: nfert_field, ndep_field 
  type( type_params_domain ) :: params_domain
  logical, parameter :: verbose = .true.

  !----------------------------------------------------------------
  ! READ RUNNAME FROM STANDARD INPUT
  !----------------------------------------------------------------
  read (*,*) runname
  ! make sure runname length is smaller/equal than maxlen_runname
  if (len_trim(runname)>=maxlen_runname) then
    stop 'runname too long'
  endif

  ! write simulation name to standard output (screen)
  write(0,*) '------------SOFUN : '//trim(runname)//'-------------'

  !----------------------------------------------------------------
  ! GET SIMULATION PARAMETERS FROM FILE <runname>.sofun.parameter
  ! SR getpar_siml is defined in _params_siml.mod.F
  !----------------------------------------------------------------
  print*,'starting'
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
  allocate( interface%grid(         interface%domaininfo%maxgrid ) )
  allocate( interface%climate(      interface%domaininfo%maxgrid ) )
  allocate( interface%ninput_field( interface%domaininfo%maxgrid ) )
  allocate( interface%landuse(      interface%domaininfo%maxgrid ) )
  allocate( interface%soilparams(   interface%domaininfo%maxgrid ) )
  allocate( interface%dfapar_field( ndayyear, interface%domaininfo%maxgrid ) )

  ! vectorise 2D array, keeping only land gridcells
  interface%grid(:) = getgrid( interface%domaininfo, params_domain )

  ! Obtain land unit dependent parameters, define decomposition _rates
  !call luparameters

  !----------------------------------------------------------------
  ! GET SOIL PARAMETERS
  !----------------------------------------------------------------
  interface%soilparams(:) = getsoil( interface%domaininfo, interface%grid(:) )

  ! LOOP THROUGH YEARS
  write(0,*) '-------------------START OF SIMULATION--------------------'

  do yr=1,interface%params_siml%runyears

    !----------------------------------------------------------------
    ! Define simulations "steering" variables (forcingyear, etc.)
    !----------------------------------------------------------------
    ! print*,'getting steering'
    interface%steering = getsteering( yr, interface%params_siml )

    if (yr == interface%params_siml%spinupyears+1 ) then
      write(0,*) '------------------TRANSIENT SIMULATION--------------------'
    endif

    !----------------------------------------------------------------
    ! Get external (environmental) forcing
    !----------------------------------------------------------------
    ! Climate
    if (verbose) print*,'getting WFDEI climate ...'
    interface%climate(:) = getclimate( &
                                      interface%domaininfo, &
                                      interface%grid, &
                                      interface%steering%init, &
                                      interface%steering%climateyear, &
                                      interface%params_siml%in_ppfd,  &
                                      interface%params_siml%in_netrad &
                                      )
    if (verbose) print*,'... done.'

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
                                 
    ! write(0,*) 'SOFUN: holding harvesting regime constant at 1993 level.'
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
    if (verbose) print*,'getting fAPAR ...'
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
    c_uptake = biosphere_annual() 
    !----------------------------------------------------------------

  enddo

  write(0,*) '--------------END OF SIMULATION---------------'

100  format(A,I6,I6,F8.2)
! 888  write(0,*) 'error opening file'
777  format (F20.8,F20.8)
999  format (I4.4)

end program main
