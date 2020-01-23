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
  use md_params_siml, only: getsteering
  ! use md_grid, only: get_domaininfo, getgrid, type_params_domain
  ! use md_params_soil, only: getsoil
  use md_forcing, only: getclimate, getco2, forcingData, climate_type
  use md_interface, only: interfacetype_biosphere, outtype_biosphere, myinterface
  use md_params_core, only: n_dim_soil_types, MSPECIES, MAX_INIT_COHORTS, ntstepsyear
  use md_biosphere, only: biosphere_annual

  use datatypes

  implicit none

  !----------------------------------------------------------------
  ! LOCAL VARIABLES READ FROM NAMELIST
  !----------------------------------------------------------------
  ! Simulation parameters
  logical :: spinup
  integer :: spinupyears
  integer :: recycle
  integer :: firstyeartrend
  integer :: nyeartrend

  ! ! integer :: model_run_years
  ! integer :: equi_days
  ! logical :: outputhourly
  ! logical :: outputdaily
  ! logical :: do_U_shaped_mortality
  ! logical :: update_annaulLAImax
  ! logical :: do_closedN_run

  ! site information
  real :: longitude
  real :: latitude
  real :: altitude

  ! ! Tile parameters (are public in datatypes)
  ! integer :: soiltype
  ! real :: FLDCAP
  ! real :: WILTPT
  ! real :: K1
  ! real :: K2
  ! real :: K_nitrogen
  ! real :: etaN
  ! real :: MLmixRatio
  ! real :: l_fract
  ! real :: retransN
  ! real :: fNSNmax
  ! real :: f_N_add
  ! real :: f_initialBSW

  ! naked arrays
  real, dimension(0:15,17) :: params_species
  real, dimension(9,8)     :: params_soil
  real, dimension(10,5)    :: init_cohort

  ! ! initial soil pool size
  ! real :: init_fast_soil_C
  ! real :: init_slow_soil_C
  ! real :: init_Nmineral
  ! real :: N_input

  !----------------------------------------------------------------
  ! LOCAL VARIABLES READ FROM/TO FILE/OUTPUT
  !----------------------------------------------------------------
  ! integer, parameter :: nt = 11 * ntstepsyear
  ! real, dimension(nt,13) :: forcing
  ! real :: output

  ! local variables
  integer :: datalines
  integer :: yr_data   ! Years of the forcing data
  integer :: totyears
  integer :: totdays
  integer :: days_data ! days of the forcing data
  real    :: timestep  ! hour, Time step of forcing data, usually hourly (1.0)
  type(outtype_biosphere) :: out_biosphere  ! holds all the output used for calculating the cost or maximum likelihood function 
  integer :: yr
  logical, parameter :: verbose = .false.

  character(len=100) :: namelistfile = '/Users/lmarques/sofun/params/parameters_Allocation.nml' !'parameters_WC_biodiversity.nml' ! 'parameters_CN.nml'

  !----------------------------------------------------------------
  ! DECLARATIONS TO READ FROM NAMELIST FILE
  !----------------------------------------------------------------
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  integer :: nml_unit

  namelist /vegn_parameters_nml/  &
  soiltype, FLDCAP, WILTPT, &
  pt, phenotype, lifeform, &
  Vmax, Vannual,wet_leaf_dreg,   &
  gamma_L, gamma_LN, gamma_SW, gamma_FR,  &
  rho_FR, root_r, root_zeta,Kw_root, &
  !rho_N_up0, N_roots0, &
  leaf_size, leafLS, LAImax, LAI_light,   &
  LMA, LNbase, CNleafsupport, c_LLS,      &
  K1,K2, K_nitrogen, etaN, MLmixRatio,    &
  LMAmin, fsc_fine, fsc_wood, &
  GR_factor, l_fract, retransN,f_N_add,   &
  f_initialBSW,f_LFR_max,  &
  gdd_crit,tc_crit, tc_crit_on, &
  alphaHT, thetaHT, alphaCA, thetaCA, alphaBM, thetaBM, &
  maturalage, v_seed, seedlingsize, prob_g,prob_e,      &
  mortrate_d_c, mortrate_d_u, A_mort, B_mort,DBHtp,     &
  phiRL, phiCSA, rho_wood, taperfactor, &
  tauNSC, fNSNmax, understory_lai_factor, &
  CNleaf0,CNsw0,CNwood0,CNroot0,CNseed0, &
  NfixRate0, NfixCost0,  &
  internal_gap_frac

  namelist /soil_data_nml/ &
  GMD, GSD, vwc_sat,k_sat_ref, psi_sat_ref, &
  chb, alphaSoil,heat_capacity_dry

  namelist /initial_state_nml/ &
  init_n_cohorts, init_cohort_species, init_cohort_nindivs, &
  init_cohort_bl, init_cohort_br, init_cohort_bsw, &
  init_cohort_bHW, init_cohort_seedC, init_cohort_nsc, &
  init_fast_soil_C, init_slow_soil_C,    & 
  init_Nmineral, N_input,  &
  filepath_in,climfile, &
  equi_days, &
  outputhourly, &
  outputdaily, &
  do_U_shaped_mortality, &
  update_annualLAImax, &
  do_closedN_run, &
  spinup, spinupyears, recycle, firstyeartrend, nyeartrend, &
  longitude, latitude, altitude

  !----------------------------------------------------------------
  ! READ FROM NAMELIST FILE
  !----------------------------------------------------------------
  ! Vegetation parameters: tile and species
  if(read_from_parameter_file)then
    nml_unit = 999
    open(nml_unit, file=trim(namelistfile), form='formatted', action='read', status='old')
    read (nml_unit, nml=vegn_parameters_nml, iostat=io, end=10)
  10    close (nml_unit)
  endif
  write(*,nml=vegn_parameters_nml)

  ! Soil parameters
  if(read_from_parameter_file)then
    nml_unit = 999
    open(nml_unit, file=trim(namelistfile), form='formatted', action='read', status='old')
    read (nml_unit, nml=soil_data_nml, iostat=io, end=20)
  20   close (nml_unit)
    write (*, nml=soil_data_nml)
  endif

  ! Initial soil and cohort states
  ! xxx todo: take info from myinterface instead of from namelist
  ! --- Generate cohorts according to "initial_state_nml" ---
  nml_unit = 999
  open(nml_unit, file=trim(namelistfile), form='formatted', action='read', status='old')
  read (nml_unit, nml=initial_state_nml, iostat=io, end=30)
  30    close (nml_unit)
  write(*,nml=initial_state_nml)

  !----------------------------------------------------------------
  ! POPULATE MYINTERFACE WITH ARGUMENTS FROM R
  !----------------------------------------------------------------
  myinterface%params_siml%do_spinup        = spinup
  myinterface%params_siml%spinupyears      = spinupyears
  myinterface%params_siml%recycle          = recycle
  myinterface%params_siml%firstyeartrend   = firstyeartrend
  myinterface%params_siml%nyeartrend       = nyeartrend

  if (myinterface%params_siml%do_spinup) then
    myinterface%params_siml%runyears = myinterface%params_siml%nyeartrend + myinterface%params_siml%spinupyears
  else
    myinterface%params_siml%runyears = myinterface%params_siml%nyeartrend
    myinterface%params_siml%spinupyears = 0
  endif

  ! Simulation parameters
  ! myinterface%params_siml%model_run_years       = myinterface%params_siml%runyears    ! xxx delete model_run_years from arguments
  !myinterface%params_siml%equi_days             = 0   ! to always write output; xxx todo: remove once output is passed back to R
  myinterface%params_siml%outputhourly          = outputhourly
  myinterface%params_siml%outputdaily           = outputdaily
  myinterface%params_siml%do_U_shaped_mortality = do_U_shaped_mortality
  ! myinterface%params_siml%update_annaulLAImax   = update_annaulLAImax      ! xxx actually not used at all. todo: remove from arguments etc.
  ! myinterface%params_siml%do_closedN_run        = do_closedN_run           ! xxx actually not used at all. todo: remove from arguments etc.

  ! Tile parameters
  myinterface%params_tile%soiltype     = soiltype
  myinterface%params_tile%FLDCAP       = FLDCAP
  myinterface%params_tile%WILTPT       = WILTPT
  myinterface%params_tile%K1           = K1
  myinterface%params_tile%K2           = K2
  myinterface%params_tile%K_nitrogen   = K_nitrogen
  myinterface%params_tile%etaN         = etaN
  myinterface%params_tile%MLmixRatio   = MLmixRatio
  myinterface%params_tile%l_fract      = l_fract
  myinterface%params_tile%retransN     = retransN
  myinterface%params_tile%f_initialBSW = f_initialBSW

  ! Species parameters
  myinterface%params_species%lifeform(:)     = lifeform(:)              !  params_species(:,1)
  myinterface%params_species%phenotype(:)    = phenotype(:)              !  params_species(:,2)
  myinterface%params_species%pt(:)           = pt(:)              !  params_species(:,3)
  myinterface%params_species%seedlingsize(:) = seedlingsize(:)              !  params_species(:,4)
  myinterface%params_species%LMA(:)          = LMA(:)              !  params_species(:,5)
  myinterface%params_species%phiRL(:)        = phiRL(:)              !  params_species(:,6)
  myinterface%params_species%LNbase(:)       = LNbase(:)              !  params_species(:,7)
  myinterface%params_species%laimax(:)       = laimax(:)              !  params_species(:,8)
  myinterface%params_species%LAI_light(:)    = LAI_light(:)              !  params_species(:,9)
  myinterface%params_species%Nfixrate0(:)    = Nfixrate0(:)              !  params_species(:,10)
  myinterface%params_species%NfixCost0(:)    = NfixCost0(:)              !  params_species(:,11)
  myinterface%params_species%phiCSA(:)       = phiCSA(:)              !  params_species(:,12)
  myinterface%params_species%mortrate_d_c(:) = mortrate_d_c(:)              !  params_species(:,13)
  myinterface%params_species%mortrate_d_u(:) = mortrate_d_u(:)              !  params_species(:,14)
  myinterface%params_species%maturalage(:)   = maturalage(:)              !  params_species(:,15)
  myinterface%params_species%fNSNmax(:)      = fNSNmax(:)              !  params_species(:,16)
  myinterface%params_species%f_N_add         = f_N_add                !  params_species(:,17)

  ! Initial cohort sizes
  myinterface%init_cohort%init_cohort_species(:) = init_cohort_species(:)      ! init_cohort(:,1)
  myinterface%init_cohort%init_cohort_nindivs(:) = init_cohort_nindivs(:)      ! init_cohort(:,2)
  myinterface%init_cohort%init_cohort_bsw(:)     = init_cohort_bsw(:)      ! init_cohort(:,3)
  myinterface%init_cohort%init_cohort_bHW(:)     = init_cohort_bHW(:)      ! init_cohort(:,4)
  myinterface%init_cohort%init_cohort_nsc(:)     = init_cohort_nsc(:)      ! init_cohort(:,5)

  ! Initial soil pools
  myinterface%init_soil%init_fast_soil_C = init_fast_soil_C
  myinterface%init_soil%init_slow_soil_C = init_slow_soil_C
  myinterface%init_soil%init_Nmineral    = init_Nmineral
  myinterface%init_soil%N_input          = N_input

  !----------------------------------------------------------------
  ! GET GRID INFORMATION
  !----------------------------------------------------------------
  myinterface%grid%lon = longitude    ! real( longitude )
  myinterface%grid%lat = latitude    ! real( latitude )
  myinterface%grid%elv = altitude    ! real( altitude )
  myinterface%grid%dogridcell = .true.  ! xxx todo remove all dogridcell statemetns
  myinterface%grid%landfrac = 1.0
  myinterface%grid%area     = 1.0

  !----------------------------------------------------------------
  ! GET SOIL PARAMETERS
  !----------------------------------------------------------------
  ! myinterface%params_soil = getsoil( params_soil )
  myinterface%params_soil%GMD               = GMD
  myinterface%params_soil%GSD               = GSD
  myinterface%params_soil%vwc_sat           = vwc_sat
  myinterface%params_soil%chb               = chb
  myinterface%params_soil%psi_sat_ref       = psi_sat_ref
  myinterface%params_soil%k_sat_ref         = k_sat_ref
  myinterface%params_soil%alphaSoil         = alphaSoil
  myinterface%params_soil%heat_capacity_dry = heat_capacity_dry

  !----------------------------------------------------------------
  ! READ FORCING FILE
  !----------------------------------------------------------------
  ! call read_FACEforcing( forcingData, datalines, days_data, yr_data, timestep )
  call read_NACPforcing( forcingData, datalines, days_data, yr_data, timestep )
  myinterface%steps_per_day = int(24.0/timestep)
  myinterface%dt_fast_yr = 1.0/(365.0 * myinterface%steps_per_day)
  myinterface%step_seconds = 24.0*3600.0/myinterface%steps_per_day ! seconds_per_year * dt_fast_yr
  ntstepsyear = myinterface%steps_per_day * 365
  !print*,'ntstepsyear ', ntstepsyear
  write(*,*) myinterface%steps_per_day, myinterface%dt_fast_yr, myinterface%step_seconds
  
  totyears = myinterface%params_siml%runyears
  totdays  = int(totyears/yr_data+1)*days_data
  myinterface%params_siml%equi_days = totdays - days_data

  !print*,'myinterface%params_siml%equi_days ', myinterface%params_siml%equi_days
  !print*,'days_data                         ', days_data
  !print*,'myinterface%params_siml%runyears  ', myinterface%params_siml%runyears
   !stop  

  ! record some variables that are determined by the SR that reads the forcing
  myinterface%datalines = datalines

  allocate(myinterface%climate(ntstepsyear))
  allocate(myinterface%pco2(ntstepsyear))


  !----------------------------------------------------------------
  ! GET CALIBRATABLE MODEL PARAMETERS (so far a small list)
  !----------------------------------------------------------------
  ! XXX warning: converting from double to single may cause a problem
  ! when calibrating and parameters are varied in their Nth digit after
  ! the comma!  
  ! myinterface%params_calib%kphio = real(par(1))

  ! !----------------------------------------------------------------
  ! ! GET VEGETATION COVER (fractional projective cover by PFT)
  ! !----------------------------------------------------------------
  ! myinterface%fpc_grid(:) = get_fpc_grid( myinterface%domaininfo, myinterface%params_siml )


  ! LOOP THROUGH YEARS
  ! print*, '-------------------START OF SIMULATION--------------------'


  do yr=1,myinterface%params_siml%runyears

    !----------------------------------------------------------------
    ! Define simulations "steering" variables (forcingyear, etc.)
    !----------------------------------------------------------------
    myinterface%steering = getsteering( yr, myinterface%params_siml )

    ! if (yr == myinterface%params_siml%spinupyears+1 ) then
    !   print*, '------------------TRANSIENT SIMULATION--------------------'
    ! endif

    !----------------------------------------------------------------
    ! Get external (environmental) forcing
    !----------------------------------------------------------------
    ! Get climate variables for this year (full fields and 365 daily values for each variable)
    myinterface%climate(:) = getclimate( &
                                        datalines, & ! nt, &
                                        forcingData, & ! forcing, &
                                        myinterface%steering%climateyear_idx, &
                                        myinterface%steering%climateyear &
                                        )

    ! print*,'sofun(): myinterface%climate(:)'
    ! print*, myinterface%climate(8000)%year
    ! print*, myinterface%climate(8000)%doy
    ! print*, myinterface%climate(8000)%hod
    ! print*, myinterface%climate(8000)%PAR
    ! print*, myinterface%climate(8000)%radiation
    ! print*, myinterface%climate(8000)%Tair
    ! print*, myinterface%climate(8000)%Tsoil
    ! print*, myinterface%climate(8000)%rain
    ! print*, myinterface%climate(8000)%windU
    ! print*, myinterface%climate(8000)%P_air
    ! print*, myinterface%climate(8000)%RH
    ! print*, myinterface%climate(8000)%CO2
    ! print*, myinterface%climate(8000)%soilwater


    ! Get annual, gobally uniform CO2
    myinterface%pco2(:) = getco2(  &
                                    datalines, & ! nt, &
                                    forcingData, & ! forcing, &
                                    myinterface%steering%climateyear_idx, &  ! to make it equivalent to BiomeE
                                    myinterface%steering%climateyear &
                                    ! myinterface%steering%forcingyear_idx, &
                                    ! myinterface%steering%forcingyear &
                                    )

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

  contains

  !========================================================================
  ! read in forcing data (Users need to write their own data input procedure)
  subroutine read_FACEforcing(forcingData,datalines,days_data,yr_data,timestep)
    type(climate_type),pointer,intent(inout) :: forcingData(:)
    integer,intent(inout) :: datalines,days_data,yr_data
    real, intent(inout)   :: timestep
    !------------local var -------------------
    type(climate_type), pointer :: climateData(:)
    character(len=80)  commts
    integer, parameter :: niterms=9       ! MDK data for Oak Ridge input
    integer, parameter :: ilines=22*366*24 ! the maxmum records of Oak Ridge FACE, 1999~2007
    integer,dimension(ilines) :: year_data
    real,   dimension(ilines) :: doy_data,hour_data
    real input_data(niterms,ilines)
    real inputstep
    integer :: istat1,istat2,istat3
    integer :: doy,idays
    integer :: i,j,k
    integer :: m,n

    ! xxx temporary
    character(len=80) :: filepath_in = '/Users/lmarques/BiomeE-Allocation/model/input/'
    character(len=80) :: climfile    = 'ORNL_forcing.txt'

    climfile=trim(filepath_in)//trim(climfile)

    ! open forcing data
    open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
    write(*,*)istat2
    ! skip 2 lines of input met data file
    read(11,'(a160)') commts
    ! read(11,'(a160)') commts ! MDK data only has one line comments
    m       = 0  ! to record the lines in a file
    idays   = 1  ! the total days in a data file
    yr_data = 0 ! to record years of a dataset
    do    ! read forcing files
      m=m+1
      read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
      (input_data(n,m),n=1,niterms)
      if(istat3<0)exit
      if(m == 1) then
        doy = doy_data(m)
      else
        doy = doy_data(m-1)
      endif
      if(doy /= doy_data(m)) idays = idays + 1
      !write(*,*)year_data(m),doy_data(m),hour_data(m)
    enddo ! end of reading the forcing file

    timestep = hour_data(2) - hour_data(1)
    write(*,*)"forcing", datalines, yr_data, timestep, myinterface%dt_fast_yr
    if (timestep==1.0)then
      write(*,*)"the data freqency is hourly"
    elseif(timestep==0.5)then
      write(*,*)"the data freqency is half hourly"
    else
      write(*,*)"Please check time step!"
      stop
    endif
    close(11)    ! close forcing file
    ! Put the data into forcing 
    datalines = m - 1
    days_data = idays
    yr_data  = year_data(datalines-1) - year_data(1) + 1

    ! allocate(climateData(datalines))

    ! xxx try
    allocate(climateData(datalines - 144))  !72
    days_data = days_data - 3

    do i=1,datalines
     if (.not. (doy_data(i)==60 .and. (year_data(i)==2000 .or. year_data(i)==2004 .or. year_data(i)==2008))) then
       climateData(i)%year      = year_data(i)          ! Year
       climateData(i)%doy       = doy_data(i)           ! day of the year
       climateData(i)%hod       = hour_data(i)          ! hour of the day
       climateData(i)%PAR       = input_data(1,i)       ! umol/m2/s
       climateData(i)%radiation = input_data(2,i)       ! W/m2
       climateData(i)%Tair      = input_data(3,i) + 273.16  ! air temperature, K
       climateData(i)%Tsoil     = input_data(4,i) + 273.16  ! soil temperature, K
       climateData(i)%RH        = input_data(5,i) * 0.01    ! relative humidity (0.xx)
       climateData(i)%rain      = input_data(6,i)/(timestep * 3600)! ! kgH2O m-2 s-1
       climateData(i)%windU     = input_data(7,i)        ! wind velocity (m s-1)
       climateData(i)%P_air     = input_data(8,i)        ! pa
       climateData(i)%CO2       = input_data(9,i) * 1.0e-6       ! mol/mol
       climateData(i)%soilwater = 0.8    ! soil moisture, vol/vol
     end if
   enddo
   forcingData => climateData

   ! xxx try
   datalines = datalines - 144  !144 !72

   write(*,*)"forcing", datalines,days_data,yr_data
 end subroutine read_FACEforcing


  !=============================================================
  ! for reading in NACP site synthesis forcing
  subroutine read_NACPforcing(forcingData,datalines,days_data,yr_data,timestep)
    type(climate_type),pointer,intent(inout) :: forcingData(:)
    integer,intent(inout) :: datalines,days_data,yr_data
    real, intent(inout)   :: timestep
    !------------local var -------------------
    type(climate_type), pointer :: climateData(:)
    character(len=80)  commts
    integer, parameter :: niterms=15       ! NACP site forcing
    integer, parameter :: ilines=22*366*48 ! the maxmum records
    integer,dimension(ilines) :: year_data
    real,   dimension(ilines) :: doy_data,hour_data
    real input_data(niterms,ilines)
    real inputstep
    integer :: istat1,istat2,istat3
    integer :: doy,idays
    integer :: i,j,k
    integer :: m,n

    ! xxx try
    integer :: idx_climatedata

    ! xxx temporary
    character(len=80) :: filepath_in = '/Users/lmarques/BiomeE-Allocation/model/input/'
    character(len=80) :: climfile    = 'US-WCrforcing.txt'

    climfile=trim(filepath_in)//trim(climfile)
    write(*,*)'inputfile: ',climfile
    ! open forcing data
    open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
    write(*,*)istat2
    ! skip 2 lines of input met data file
    read(11,'(a160)') commts
    read(11,'(a160)') commts
    m       = 0  ! to record the lines in a file
    idays   = 1  ! the total days in a data file
    yr_data = 0 ! to record years of a dataset
    do    ! read forcing files
      m=m+1
      read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
      (input_data(n,m),n=1,niterms)

      if(istat3<0)exit
      if(m == 1) then
        doy = doy_data(m)
      else
        doy = doy_data(m-1)
      endif
      if(doy /= doy_data(m)) idays = idays + 1
      !write(*,*)year_data(m),doy_data(m),hour_data(m)
      ! discard one line
      !read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
      !                        (input_data(n,m),n=1,niterms)
    enddo ! end of reading the forcing file

    timestep = hour_data(2) - hour_data(1)
    write(*,*)"forcing",datalines,yr_data,timestep,myinterface%dt_fast_yr
    if (timestep==1.0)then
      write(*,*)"the data freqency is hourly"
    elseif(timestep==0.5)then
      write(*,*)"the data freqency is half hourly"
    else
      write(*,*)"Please check time step!"
      stop
    endif
    close(11)    ! close forcing file
    ! Put the data into forcing 
    datalines = m - 1
    days_data = idays
    yr_data  = year_data(datalines-1) - year_data(1) + 1

    ! allocate(climateData(datalines))

    ! xxx try
    allocate(climateData(datalines- 96))
    days_data = days_data - 2

    idx_climatedata = 0
    do i=1,datalines

      if (.not. (doy_data(i)==60 .and. (year_data(i)==2000 .or. year_data(i)==2004))) then

        idx_climatedata = idx_climatedata + 1

        climateData(idx_climatedata)%year      = year_data(i)          ! Year
        climateData(idx_climatedata)%doy       = doy_data(i)           ! day of the year
        climateData(idx_climatedata)%hod       = hour_data(i)          ! hour of the day
        climateData(idx_climatedata)%PAR       = input_data(11,i)*2.0  ! umol/m2/s
        climateData(idx_climatedata)%radiation = input_data(11,i)      ! W/m2
        climateData(idx_climatedata)%Tair      = input_data(1,i)       ! air temperature, K
        climateData(idx_climatedata)%Tsoil     = input_data(1,i)       ! soil temperature, K
        climateData(idx_climatedata)%rain      = input_data(7,i)       ! kgH2O m-2 s-1
        climateData(idx_climatedata)%windU     = input_data(5,i)        ! wind velocity (m s-1)
        climateData(idx_climatedata)%P_air     = input_data(9,i)        ! pa
        climateData(idx_climatedata)%RH        = input_data(3,i)/mol_h2o*mol_air* & ! relative humidity (0.xx)
        climateData(idx_climatedata)%P_air/esat(climateData(idx_climatedata)%Tair-273.16)
        climateData(idx_climatedata)%CO2       = input_data(15,i) * 1.0e-6       ! mol/mol
        climateData(idx_climatedata)%soilwater = 0.8    ! soil moisture, vol/vol

        ! if (abs(climateData(idx_climatedata)%hod-12.0)<0.001) then
        !   print*,'A year, doy, Tair ', climateData(idx_climatedata)%year, &
        !    climateData(idx_climatedata)%doy, &
        !    climateData(idx_climatedata)%Tair
        ! end if 

      end if

    enddo
    forcingData => climateData

    ! xxx try
    datalines = datalines - 96

    write(*,*)"forcing", datalines,days_data,yr_data

  end subroutine read_NACPforcing


end program main
