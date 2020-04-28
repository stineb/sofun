program main
  !////////////////////////////////////////////////////////////////
  ! Main program for SOFUN, used for single-site simulations and 
  ! "lonlat" (=multiple points on a spatial grid) simulations.
  ! 
  ! - reads run name from standard input
  ! - invokes functions to get simulation parameters, grid definition, etc.
  ! - loops over simulation years, including spinup
  ! - invokes functions to get annually varying forcing
  ! - call function biosphere()
  !
  ! Example (setup 'pmodel'):
  ! echo RUNNAME ./runpmodel
  ! 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  !----------------------------------------------------------------
  use md_interface, only: interface, outtype_biosphere
  use md_params_siml, only: getpar_siml, getsteering
  use md_params_domain, only: getpar_domain, type_params_domain
  use md_grid, only: get_domaininfo, getgrid
  use md_params_soil, only: getsoil
  use md_forcing, only: get_fpc_grid, getclimate, getninput, ninput_type, gettot_ninput, &
    getfapar, getlanduse, getco2
  use md_params_core, only: dummy, maxgrid, ndayyear, npft
  use md_biosphere, only: biosphere_annual

  implicit none

  ! local variables
  integer :: yr                                  ! simulation year
  character(len=245) :: runname                  ! run name
  integer, parameter :: maxlen_runname = 50      ! maximum length of runname (arbitrary)
  type( ninput_type ), dimension(maxgrid) :: nfert_field, ndep_field 
  type( type_params_domain ) :: params_domain
  logical, parameter :: verbose = .true.         ! set this true for additional messages

  type(outtype_biosphere) :: out_biosphere       ! derived type containing certain quantities calculated within biosphere(), not used here.

  !----------------------------------------------------------------
  ! READ RUNNAME FROM STANDARD INPUT
  !----------------------------------------------------------------
  read (*,*) runname
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

  ! set parameter to define that this is not a calibration run (otherwise sofun.f90 would not have been compiled, but sofun_simsuite.f90)
  interface%params_siml%is_calib = .false.

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
  allocate( interface%grid( interface%domaininfo%maxgrid ) )
  allocate( interface%climate( ndayyear, interface%domaininfo%maxgrid ) )
  ! allocate( interface%ninput_field(   interface%domaininfo%maxgrid ) )
  ! allocate( interface%landuse(        interface%domaininfo%maxgrid ) )
  allocate( interface%soilparams(     interface%domaininfo%maxgrid ) )
  allocate( interface%vegcover( ndayyear, interface%domaininfo%maxgrid ) )
  allocate( interface%fpc_grid( npft, interface%domaininfo%maxgrid ) )

  ! vectorise 2D array, keeping only land gridcells
  interface%grid(:) = getgrid( interface%domaininfo, params_domain )

  !----------------------------------------------------------------
  ! GET SOIL PARAMETERS
  !----------------------------------------------------------------
  interface%soilparams(:) = getsoil( interface%domaininfo )

  !----------------------------------------------------------------
  ! GET VEGETATION COVER (fractional projective cover by PFT)
  !----------------------------------------------------------------
  interface%fpc_grid(:,:) = get_fpc_grid( interface%domaininfo, interface%grid(:), interface%params_siml )


  ! LOOP THROUGH YEARS
  write(0,*) 'SOFUN site-level run: ', runname
  print*, '-------------------START OF SIMULATION--------------------'

  do yr=1,interface%params_siml%runyears

    !----------------------------------------------------------------
    ! Define simulations "steering" variables (forcingyear, etc.)
    !----------------------------------------------------------------
    ! print*,'getting steering'
    interface%steering = getsteering( yr, interface%params_siml )

    if (yr == interface%params_siml%spinupyears+1 ) then
      print*, '------------------TRANSIENT SIMULATION--------------------'
    endif

    print*,'--------------------------------------------------------'
    print*,'Simulation year: ', interface%steering%year, ' - Real year: ', interface%steering%outyear
    print*,'--------------------------------------------------------'
    
    !----------------------------------------------------------------
    ! Get prescribed fAPAR if required (otherwise set to dummy value)
    !----------------------------------------------------------------
    if (verbose) print*,'getting fAPAR ...'
    interface%vegcover(:,:) = getfapar( &
                                        interface%domaininfo, &
                                        interface%grid, &
                                        interface%steering%forcingyear, &
                                        interface%params_siml%fapar_forcing_source &
                                        )    


    interface%fpc_grid(:,:) = get_fpc_grid( &
                                            interface%domaininfo, &
                                            interface%grid, &
                                            interface%params_siml &
                                            )

    !----------------------------------------------------------------
    ! Get external (environmental) forcing
    !----------------------------------------------------------------
    ! Get climate variables for this year (full fields and 365 daily values for each variable)
    if (verbose) print*,'getting climate ...'
    interface%climate(:,:) = getclimate( &
                                          interface%domaininfo, &
                                          interface%grid, &
                                          interface%steering%init, &
                                          interface%steering%climateyear, &
                                          interface%params_siml%in_ppfd,  &
                                          interface%params_siml%in_netrad &
                                          )
    if (verbose) print*,'... done.'

    ! Get annual, gobally uniform CO2
    if (verbose) print*,'getting CO2 ...'
    interface%pco2 = getco2( &
                            trim(runname), &
                            interface%domaininfo, &
                            interface%steering%forcingyear, &
                            interface%params_siml%const_co2_year, &
                            interface%params_siml%firstyeartrend,&
                            interface%params_siml%co2_forcing_file&
                            )
    if (verbose) print*,'... done.'

    ! ! Atmospheric N deposition (note that actual data is not read in all SOFUN setups)
    ! ndep_field(:) = getninput( &
    !                           "ndep", &
    !                           trim(runname), &
    !                           interface%domaininfo, &
    !                           interface%steering%forcingyear, &
    !                           interface%params_siml%firstyeartrend, &
    !                           interface%params_siml%const_ndep_year, &
    !                           interface%params_siml%ndep_noy_forcing_file, &
    !                           interface%params_siml%ndep_nhx_forcing_file, &
    !                           interface%climate(:)&
    !                           )
    
    ! ! N fertiliser input (note that actual data is not read in all SOFUN setups)
    ! nfert_field(:) = getninput( &
    !                           "nfert", &
    !                           trim(runname), &
    !                           interface%domaininfo, &
    !                           interface%steering%forcingyear, &
    !                           interface%params_siml%firstyeartrend, &
    !                           interface%params_siml%const_nfert_year, &
    !                           interface%params_siml%nfert_noy_forcing_file, &
    !                           interface%params_siml%nfert_nhx_forcing_file, &
    !                           interface%climate(:)&
    !                           )

    ! ! Interface holds only total reactive N input (N deposition + N fertiliser)                             
    ! interface%ninput_field(:) = gettot_ninput( nfert_field(:), ndep_field(:) )
                                 
    ! ! Get land use information (note that actual data is not read in all SOFUN setups)
    ! interface%landuse(:) = getlanduse( &
    !                                   trim(runname), &
    !                                   interface%domaininfo, &
    !                                   interface%steering%forcingyear, &
    !                                   interface%params_siml%do_grharvest_forcing_file, &
    !                                   interface%params_siml%const_lu_year, &
    !                                   interface%params_siml%firstyeartrend &
    !                                   )

    !----------------------------------------------------------------
    ! Call SR biosphere at an annual time step but with vectors 
    ! containing data for each day of this year.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Call biosphere (wrapper for all modules, contains gridcell loop)
    !----------------------------------------------------------------
    if (verbose) print*,'calling biosphere ...'
    out_biosphere = biosphere_annual() 
    if (verbose) print*,'... done.'
    !----------------------------------------------------------------

  enddo

  print*, '--------------END OF SIMULATION---------------'

100  format (A,I6,I6,F8.2)
777  format (F20.8,F20.8)
999  format (I4.4)

end program main
