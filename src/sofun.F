program main
  !////////////////////////////////////////////////////////////////
  !  Main program for site scale simulations, here used for 
  !  SOFUN (Seasonal Optimisation of Fixation and Uptake of 
  !  Nitrogen)
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  use _params_core
  use _params_siml
  use _params_site
  use _forcing_siterun
  use _gridvars
  use _params_soil

  implicit none

  integer :: year

  real :: c_uptake     ! annual net global C uptake by biosphere

  ! local variables
  real                                     :: pco2
  real, dimension(ndayyear,maxgrid)        :: dndep_field
  real, dimension(nlu,maxgrid)             :: lu_area
  real, dimension(nmonth,maxgrid)          :: mfapar_field
  type(paramtype_soil), dimension(maxgrid) :: params_soil_field

  ! xxx try
  integer :: day

  type( outtype_steering ) :: out_steering

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
  call getpar_siml(trim(runname))

  !----------------------------------------------------------------
  ! GET SITE PARAMETERS AND INPUT DATA
  ! site location (lon,lat), soil type, vegetation type
  ! SR getpar_site is defined in _params_site.mod.F. 
  ! 'sitename' is global variable
  !----------------------------------------------------------------
  call getpar_site

  !----------------------------------------------------------------
  ! GET GRID INFORMATION
  ! longitude, latitude, elevation
  !----------------------------------------------------------------
  call getgrid

  ! Get soil parameters (if not defined in <sitename>.parameter)
  !call getsoilpar

  ! Obtain land unit dependent parameters, define decomposition _rates
  !call luparameters

  !----------------------------------------------------------------
  ! GET SOIL PARAMETERS
  !----------------------------------------------------------------
  params_soil_field(:) = getsoil_field( soilcode_field(:) )

  ! LOOP THROUGH YEARS
  write(0,*) '------------START OF SIMULATION-------------'

  do year=1,runyears

    !----------------------------------------------------------------
    ! Define simulations "steering" variables (forcingyear, etc.)
    !----------------------------------------------------------------
    ! print*,'getting steering'
    out_steering = getsteering( year )

    if (year == spinupyears+1 ) then
      write(0,*) '-----------TRANSIENT SIMULATION-------------'
    endif

    !----------------------------------------------------------------
    ! Get external (environmental) forcing for all simulation years
    ! xxx move this into annual loop
    !----------------------------------------------------------------
    ! out_climate = getclimate_site( trim(runname), trim(sitename), out_steering%climateyear )
    call getclimate_site( trim(runname), trim(sitename), out_steering%climateyear )

    !----------------------------------------------------------------
    ! Get external (environmental) forcing
    !----------------------------------------------------------------
    pco2             = getco2( trim(runname), trim(sitename), out_steering%forcingyear )
    dndep_field(:,:) = getndep( trim(runname), trim(sitename), out_steering%forcingyear )
    lu_area(:,:)     = getlanduse( trim(runname), out_steering%forcingyear )

    !----------------------------------------------------------------
    ! Get prescribed fAPAR
    !----------------------------------------------------------------
    mfapar_field(:,:) = getfapar( trim(runname), trim(sitename), out_steering%forcingyear )

    !----------------------------------------------------------------
    ! Call SR biosphere at an annual time step but with vectors 
    ! containing data for each day of this year.
    !----------------------------------------------------------------
    write(0,100) 'sim. year, year AD, pco2', year, out_steering%forcingyear, pco2

    ! write(0,*) 'lon, lat, elv', lon, lat, elv

    ! write(0,*) 'params_soil_field', params_soil_field(:)
    ! write(0,*) 'andep', sum(dndep_field(:,:))
    ! write(0,*) 'aprec', sum(dprec_field(:,:))
    ! write(0,*) 'dtemp ave', sum(dtemp_field(:,:))/ndayyear
    ! write(0,*) 'dfsun ave', sum(dfsun_field(:,:))/ndayyear
    ! write(0,*) 'dvpd ave', sum(dvpd_field(:,:))/ndayyear
    ! write(0,*) 'mfapar_field', mfapar_field
    ! write(0,*) 'calling elvis ...'
    ! stop 'before calling biosphere()'

    call biosphere( &
      year, lon, lat, elv &
      , params_soil_field, lu_area, pco2 &
      , dtemp_field, dprec_field &
      , dfsun_field, dvpd_field, dndep_field &
      , c_uptake &
      , mfapar_field &
      ) 

  enddo

  write(0,*) '--------------END OF SIMULATION---------------'

100  format(A,I6,I6,F8.2)
! 888  write(0,*) 'error opening file'
777  format (F20.8,F20.8)

end program
