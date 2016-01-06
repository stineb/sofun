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
  use input_site
  use gridvars_siterun
  ! use _vars_core

  implicit none

  integer :: year

  real :: c_uptake     ! annual net global C uptake by biosphere

  ! local variables
  real                              :: pco2
  real, dimension(ndayyear,maxgrid) :: dndep_field
  real, dimension(nlu,maxgrid)      :: lu_area

  type( outtype_steering ) :: out_steering
  ! type( outtype_climate )  :: out_climate

  ! This is necessary because of keyword-specification of arguments to 'biosphere'
  interface
    subroutine biosphere( year, lon, lat, elv, lu_area, pco2 &
      , dtemp_field, dprec_field, dfsun_field, dvpd_field, dndep_field &
      , c_uptake &
      , dgpp_data &
      , fpc_grid_data &
      )

    ! required arguments
    integer, intent(in)                   :: year       ! simulation year
    real, intent(in), dimension(maxgrid)  :: lon        ! longitude vector/field (degrees E)              
    real, intent(in), dimension(maxgrid)  :: lat        ! latitude vector/field (degrees N)             
    real, intent(in), dimension(maxgrid)  :: elv        ! elevation (altitude) vector/field (m above sea level)                  
    real, dimension(3)                    :: lu_area    ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside 'getlanduse'
    real, intent(in)                      :: pco2
    real, intent(in), dimension(ndayyear,maxgrid) :: dtemp_field
    real, intent(in), dimension(ndayyear,maxgrid) :: dprec_field
    real, intent(in), dimension(ndayyear,maxgrid) :: dfsun_field
    real, intent(in), dimension(ndayyear,maxgrid) :: dvpd_field
    real, intent(in), dimension(ndayyear,maxgrid) :: dndep_field
    real, intent(out)                     :: c_uptake   ! annual net global C uptake by biosphere
    logical                               :: spinup

    ! optional arguments
    real, intent(in), dimension(ndayyear,npft), optional :: dgpp_data       ! for runs with prescribed GPP
    real, intent(in), dimension(nlu,maxgrid), optional   :: fpc_grid_data   ! for runs with prescribed GPP

    end subroutine biosphere
  end interface

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
  !       write (*,200)
  !200    format ('SIMULATION NAME: '$)
  !      write (*,*) runname


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
  call getgrid

  ! Get soil parameters (if not defined in <sitename>.parameter)
  !call getsoilpar

  ! Obtain land unit dependent parameters, define decomposition _rates
  !call luparameters

  !----------------------------------------------------------------
  ! GET MODEL INPUT DATA
  ! Prescribe variables instead of calculating these parts online
  ! This reads daily values for each day of the transient simulation
  !----------------------------------------------------------------
  call getinput_site


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
    pco2        = getco2( trim(runname), out_steering%forcingyear )
    dndep_field = getndep( trim(runname), trim(sitename), out_steering%forcingyear )
    lu_area     = getlanduse( trim(runname), out_steering%forcingyear )

    ! ! xxx test
    ! dndep_field(:) = dndep_field(:) * 0.4
    ! dprec_field(100:250,1) = 0.0

    !----------------------------------------------------------------
    ! Call SR biosphere at an annual time step but with vectors 
    ! containing data for each day of this year.
    !----------------------------------------------------------------
    write(0,100) 'sim. year, year AD, pco2', year, out_steering%forcingyear, pco2
    call biosphere( year, lon, lat, elv, lu_area, pco2 &
      , dtemp_field, dprec_field &
      , dfsun_field, dvpd_field, dndep_field &
      , c_uptake &
#if _prescr_gpp_site
      , dgpp_data=dgpp_data(transyear,:,:) &  ! optional argument
#endif
      , fpc_grid_data=fpc_grid_data &
      ) 

  enddo

  write(0,*) '--------------END OF SIMULATION---------------'

100     format(A,I6,I6,F8.2)

end program
