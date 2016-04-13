program main
  !////////////////////////////////////////////////////////////////
  !  Main program for site scale simulations, here used for 
  !  SOFUN (Seasonal Optimisation of Fixation and Uptake of 
  !  Nitrogen)
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  use md_interface
  use md_params_siml, only: getpar_siml, getsteering
  use md_params_site, only: getpar_site
  use md_grid, only: getgrid
  use md_params_soil, only: getsoil_field
  use md_forcing_siterun, only: getclimate_site, getndep, getfapar, getlanduse, getco2

  implicit none

  integer :: yr
  real    :: c_uptake     ! annual net global C uptake by biosphere
  character(len=245) :: runname

  integer, parameter :: maxlen_runname = 50      ! maximum length of runname (arbitrary)

  ! xxx try
  integer :: day

  !----------------------------------------------------------------
  ! READ RUNNAME FROM STANDARD INPUT
  !----------------------------------------------------------------
  read (*,*) runname
  ! make sure runname length is smaller/equal than maxlen_runname
  if (len_trim(runname)>=maxlen_runname) then
    stop'runname too long'
  endif

  ! write simulation name to standard output (screen)
  write(0,*) '------------SOFUN : '//trim(runname)//'-------------'

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
  interface%params_site = getpar_site( trim(interface%params_siml%sitename) )

  !----------------------------------------------------------------
  ! GET GRID INFORMATION
  ! longitude, latitude, elevation
  !----------------------------------------------------------------
  interface%grid(:) = getgrid( trim(interface%params_siml%sitename) )

  ! Get soil parameters (if not defined in <sitename>.parameter)
  !call getsoilpar

  ! Obtain land unit dependent parameters, define decomposition _rates
  !call luparameters

  !----------------------------------------------------------------
  ! GET SOIL PARAMETERS
  !----------------------------------------------------------------
  interface%soilparams(:) = getsoil_field( interface%params_site%soilcode )

  ! LOOP THROUGH YEARS
  write(0,*) '------------START OF SIMULATION-------------'

  do yr=1,interface%params_siml%runyears

    !----------------------------------------------------------------
    ! Define simulations "steering" variables (forcingyear, etc.)
    !----------------------------------------------------------------
    ! print*,'getting steering'
    interface%steering = getsteering( yr, interface%params_siml )

    if (yr == interface%params_siml%spinupyears+1 ) then
      write(0,*) '-----------TRANSIENT SIMULATION-------------'
    endif

    !----------------------------------------------------------------
    ! Get external (environmental) forcing
    !----------------------------------------------------------------
    interface%climate(:) = getclimate_site( &
                                          trim(runname), &
                                          trim(interface%params_siml%sitename), &
                                          ! 1993 &
                                          interface%steering%climateyear &
                                          )
    interface%pco2 = getco2( &
                            trim(runname), &
                            trim(interface%params_siml%sitename), &
                            ! 1993, &
                            interface%steering%forcingyear, &
                            interface%params_siml%const_co2, &
                            interface%params_siml%firstyeartrend,&
                            interface%params_siml%co2_forcing_file&
                            )
    interface%ndep_field(:) = getndep( &
                                      trim(runname), &
                                      trim(interface%params_siml%sitename), &
                                      ! 1993, &
                                      interface%steering%forcingyear, &
                                      interface%params_siml%firstyeartrend, &
                                      interface%params_siml%const_ndep, &
                                      interface%params_siml%ndep_noy_forcing_file, &
                                      interface%params_siml%ndep_nhx_forcing_file, &
                                      interface%climate(:)&
                                      )
    write(0,*) 'SOFUN: holding harvesting regime constant at 1993 level.'
    interface%landuse(:) = getlanduse( &
                                      trim(runname), &
                                      trim(interface%params_siml%sitename), &
                                      1993, &
                                      ! interface%steering%forcingyear, &
                                      interface%params_siml%do_grharvest_forcing_file &
                                      )
    ! pft_field(:)     = getpft( trim(interface%params_siml%sitename), interface%steering%forcingyear )

    !----------------------------------------------------------------
    ! Get prescribed fAPAR
    !----------------------------------------------------------------
    interface%mfapar_field(:,:) = getfapar( &
                                          trim(runname), &
                                          trim(interface%params_siml%sitename), &
                                          interface%steering%forcingyear, &
                                          interface%params_siml%prescr_monthly_fapar &
                                          )

    !----------------------------------------------------------------
    ! Call SR biosphere at an annual time step but with vectors 
    ! containing data for each day of this year.
    !----------------------------------------------------------------
    write(0,100) 'sim. year, year AD, pco2', yr, interface%steering%forcingyear, interface%pco2

    !================================================================
    !================================================================
    call biosphere( c_uptake ) 
    !================================================================
    !================================================================

  enddo

  write(0,*) '--------------END OF SIMULATION---------------'

100  format(A,I6,I6,F8.2)
! 888  write(0,*) 'error opening file'
777  format (F20.8,F20.8)

end program
