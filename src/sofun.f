# 1 "sofun.f90"
program main
!////////////////////////////////////////////////////////////////
!  Main program for site scale simulations, here used for
!  SOFUN (Seasonal Optimisation of Fixation and Uptake of
!  Nitrogen)
! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
! contact: b.stocker@imperial.ac.uk
!----------------------------------------------------------------
# 1 "./sofun_module_control.inc" 1 
! //////////////////////////////////////////////////////////////////////
! SELECT MODULES FOR DIFFERENT PARAMETRISATIONS
! ----------------------------------------------------------------------
! This is used only to define setup (prescribing/simulating model parts)
! on the top level (program sofun). Within 'biosphere', select respective
! modules instead.

! ----------------------------------------------------------------------
! MAINTENANCE RESPIRATION
! Default: like in LPX
! Option:   After Mäkelä et al., 2008, New Phyt.

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! PHOTOSYNTHESIS
! Default: xxx
! Option:  After Mäkelä et al., 2008, New Phyt.

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! N UPTAKE
! Default: FUN type after Fisher et al., 2011
! Option:  After Mäkelä et al., 2008, New Phyt.

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ALLOCATION
! Default: annually called, like in LPX for grass PFT
! Option:  Pipe model as described in Mäkelä et al., 2008, New Phyt.
! Note that _allocation_pipe does not account for the separation between
! sapwood and heartwood. The pool 'psapw' is taken to represent
! 'live wood'. This affects how turnover is formulated (see turnover.F)

! ----------------------------------------------------------------------

! //////////////////////////////////////////////////////////////////////
! SELECT MODEL INTEGRATION
! defining which model parts are to be calculated online and which to
! be prescribed. Default is always to calculate it online (not prescr.)
! In general, if selected, daily values for the course of the whole
! transient simulation are read in when initialising the model. This
! requires input files formatted like ''
! This is only useful for site scale simulations.
! ----------------------------------------------------------------------

! Prescribe GPP per day


! ----------------------------------------------------------------------
! FIXED/DYNAMIC VEGETATION
! Default: fixed, reading in FPC for each PFT from site parameter file)
! Option:  dynamic, like in LPX

! ----------------------------------------------------------------------


! //////////////////////////////////////////////////////////////////////
! SANITY
! Compile additional code for online sanity checks. This is recommended
! for simulations during development stage. De-activate later to reduce
! computational costs.
! ----------------------------------------------------------------------

# 10 "sofun.f90" 2 
  use md_params_core
  use md_params_siml
  use md_params_site
  use md_forcing_siterun
  use md_gridvars
  use md_params_soil

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
  call getpar_siml( trim(runname) )

!----------------------------------------------------------------
! GET SITE PARAMETERS AND INPUT DATA
! site location (lon,lat), soil type, vegetation type
! SR getpar_site is defined in _params_site.mod.F.
! 'sitename' is global variable
!----------------------------------------------------------------
  call getpar_site()

!----------------------------------------------------------------
! GET GRID INFORMATION
! longitude, latitude, elevation
!----------------------------------------------------------------
  call getgrid()

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
! Get external (environmental) forcing
!----------------------------------------------------------------
    climate_field = getclimate_site( trim(runname), trim(sitename), out_steering%climateyear )

!----------------------------------------------------------------
! Get external (environmental) forcing
!----------------------------------------------------------------
    pco2             = getco2( trim(runname), trim(sitename), out_steering%forcingyear )
    dndep_field(:,:) = getndep( trim(runname), trim(sitename), out_steering%forcingyear )
    lu_area(:,:)     = getlanduse( trim(runname), out_steering%forcingyear )
! pft_field(:)     = getpft( trim(sitename), out_steering%forcingyear )

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
      year, lon(:), lat(:), elv(:) &
      , params_soil_field(:), lu_area(:,:), pco2 &
      , climate_field%dtemp(:,:), climate_field%dprec(:,:) &
      , climate_field%dfsun(:,:), climate_field%dvpd(:,:) &
      , dndep_field(:,:) &
      , c_uptake &
      , mfapar_field &
      ) 

  enddo

  write(0,*) '--------------END OF SIMULATION---------------'

100  format(A,I6,I6,F8.2)
! 888  write(0,*) 'error opening file'
777  format (F20.8,F20.8)

end program
