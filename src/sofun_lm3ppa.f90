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
  use md_forcing, only: getclimate, getco2, forcingData, climate_type
  use md_interface, only: interfacetype_biosphere, outtype_biosphere, myinterface
  use md_params_core, only: n_dim_soil_types, MSPECIES, MAX_INIT_COHORTS, ntstepsyear, out_max_cohorts, &
    ndayyear, nvars_daily_tile, nvars_hourly_tile, nvars_daily_cohorts, nvars_annual_cohorts, nvars_annual_tile
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
  integer :: runyears

  ! site information
  real :: longitude
  real :: latitude
  real :: altitude

  integer :: idx_hourly_start, idx_hourly_end, idx_daily_start, idx_daily_end ! year index for which climate is read in.

  ! naked arrays
  real, dimension(0:15,17) :: params_species
  real, dimension(9,8)     :: params_soil
  real, dimension(10,5)    :: init_cohort

  !----------------------------------------------------------------
  ! LOCAL VARIABLES READ FROM/TO FILE/OUTPUT
  !----------------------------------------------------------------
  ! local variables
  integer :: datalines
  integer :: yr_data   ! Years of the forcing data
  integer :: totyears
  integer :: totdays
  integer :: days_data ! days of the forcing data
  real    :: timestep  ! hour, Time step of forcing data, usually hourly (1.0)
  integer :: ntstepsyear_forcing                 ! 365*48 when half-hourly inputs, 365*24 when hourly inputs
  type(outtype_biosphere) :: out_biosphere  ! holds all the output used for calculating the cost or maximum likelihood function 
  integer :: yr
  logical, parameter :: verbose = .false.
  integer :: iday

  character(len=100) :: namelistfile = '/Users/benjaminstocker/sofun/params/parameters_Allocation.nml' !'parameters_WC_biodiversity.nml' ! 'parameters_CN.nml'

  ! output arrays (naked) to be passed back to C/R
  real, dimension(:,:), allocatable  :: out_hourly_tile 
  real, dimension(:,:), allocatable  :: out_daily_tile       !fno4
  real, dimension(:,:,:), allocatable:: out_daily_cohorts    !fno3
  real, dimension(:,:), allocatable  :: out_annual_tile      !fno5
  real, dimension(:,:,:), allocatable:: out_annual_cohorts   !fno2

  ! whether fast time step processes are simulated. If .false., then C, N, and W balance is simulated daily.
  logical, parameter :: daily = .false.

  !----------------------------------------------------------------
  ! DECLARATIONS TO READ FROM NAMELIST FILE
  !----------------------------------------------------------------
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  integer :: nml_unit

  integer :: j

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
  call read_FACEforcing( forcingData, datalines, days_data, yr_data, timestep ) !! ORNL
  ! call read_NACPforcing( forcingData, datalines, days_data, yr_data, timestep ) !!US-WCrforcing
  
  ! record some useful variables that are determined by the SR that reads the forcing
  if (daily) timestep = 24.0
  myinterface%steps_per_day = int(24.0/timestep)
  myinterface%dt_fast_yr = 1.0/(365.0 * myinterface%steps_per_day)
  myinterface%step_seconds = 24.0 * 3600.0 / myinterface%steps_per_day ! seconds_per_year * dt_fast_yr
  ntstepsyear = myinterface%steps_per_day * 365
  if (daily) then 
    ntstepsyear_forcing = ntstepsyear * timestep
  else
    ntstepsyear_forcing = ntstepsyear
  end if
  totyears = myinterface%params_siml%runyears
  totdays  = int(totyears/yr_data+1) * days_data
  myinterface%params_siml%equi_days = totdays - days_data
  myinterface%datalines = datalines

  print*, myinterface%steps_per_day, myinterface%dt_fast_yr, myinterface%step_seconds

  ! allocate memory
  allocate(myinterface%climate(ntstepsyear))
  allocate(myinterface%pco2(ntstepsyear))

  allocate(out_biosphere%hourly_tile(ntstepsyear))

  allocate(out_hourly_tile(     ntstepsyear * myinterface%params_siml%nyeartrend,  nvars_hourly_tile                        ))
  allocate(out_daily_cohorts(   ndayyear * myinterface%params_siml%nyeartrend,     out_max_cohorts,    nvars_daily_cohorts  ))
  allocate(out_daily_tile(      ndayyear * myinterface%params_siml%nyeartrend,     nvars_daily_tile                         ))
  allocate(out_annual_cohorts(  myinterface%params_siml%runyears,                  out_max_cohorts,    nvars_annual_cohorts ))
  allocate(out_annual_tile(     myinterface%params_siml%runyears,                  nvars_annual_tile                        ))
  
  ! LOOP THROUGH YEARS
  print*, '--------------START OF SIMULATION---------------'

  do yr=1, myinterface%params_siml%runyears

    !----------------------------------------------------------------
    ! Define simulations "steering" variables (forcingyear, etc.)
    !----------------------------------------------------------------
    myinterface%steering = getsteering( yr, myinterface%params_siml )

    if (yr == myinterface%params_siml%spinupyears+1 ) then
      print*, '--------------TRANSIENT SIMULATION---------------'

    endif

    !----------------------------------------------------------------
    ! Get external (environmental) forcing
    !----------------------------------------------------------------
    ! Get climate variables for this year (full fields and 365 daily values for each variable)
    myinterface%climate(:) = getclimate( &
                                        datalines, &
                                        ntstepsyear, &
                                        ntstepsyear_forcing, &
                                        daily, &
                                        forcingData, &
                                        myinterface%steering%climateyear_idx, &
                                        myinterface%steering%climateyear, &
                                        myinterface%grid%elv &
                                        )

    ! ! Get annual, gobally uniform CO2
    ! myinterface%pco2(:) = getco2( &
    !                               datalines, &
    !                               ntstepsyear, &
    !                               ntstepsyear_forcing, &
    !                               daily, &
    !                               forcingData, &
    !                               myinterface%steering%climateyear_idx, &  ! to make it equivalent to BiomeE
    !                               myinterface%steering%climateyear &
    !                               )

    !----------------------------------------------------------------
    ! Call SR biosphere at an annual time step but with vectors 
    ! containing data for each day of this year.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Call biosphere (wrapper for all modules, contains gridcell loop)
    !----------------------------------------------------------------
    if (verbose) print*,'calling biosphere ...'
    call biosphere_annual( out_biosphere ) 
    if (verbose) print*,'... done.'
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Populate big output arrays
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Print out_hourly_tile
    !----------------------------------------------------------------
    if (.not. myinterface%steering%spinup) then

      idx_hourly_start = (yr - myinterface%params_siml%spinupyears - 1) * ntstepsyear + 1          ! To exclude the spinup years and include only the transient years
      idx_hourly_end   = idx_hourly_start + ntstepsyear - 1
      call populate_outarray_hourly_tile( out_biosphere%hourly_tile(:), out_hourly_tile(idx_hourly_start:idx_hourly_end, :) )
    
    end if

    !----------------------------------------------------------------
    ! Print out_daily_tile
    !----------------------------------------------------------------
    ! Output only for transient years

    if (.not. myinterface%steering%spinup) then  

      idx_daily_start = (yr - myinterface%params_siml%spinupyears - 1) * ndayyear + 1 
      idx_daily_end   = idx_daily_start + ndayyear - 1

      call populate_outarray_daily_tile( out_biosphere%daily_tile(:), out_daily_tile(idx_daily_start:idx_daily_end, :) )
    
    end if

    !----------------------------------------------------------------
    ! Print out_daily_cohorts
    !----------------------------------------------------------------
    ! Output only for transient years

    if (.not. myinterface%steering%spinup) then  
      call populate_outarray_daily_cohorts( out_biosphere%daily_cohorts(:,:), out_daily_cohorts(idx_daily_start:idx_daily_end,:,:) )
    end if

    !----------------------------------------------------------------
    ! Print out_annual_tile
    !----------------------------------------------------------------
    call populate_outarray_annual_tile( out_biosphere%annual_tile, out_annual_tile(yr,:) )

    !----------------------------------------------------------------
    ! Print out_annual_cohorts
    !----------------------------------------------------------------
    call populate_outarray_annual_cohorts( out_biosphere%annual_cohorts(:), out_annual_cohorts(yr,:,:) )

  enddo

  deallocate(myinterface%climate)
  deallocate(myinterface%pco2)
  deallocate(out_biosphere%hourly_tile)
  deallocate(out_hourly_tile)
  deallocate(out_daily_cohorts)
  deallocate(out_daily_tile)
  deallocate(out_annual_cohorts)
  deallocate(out_annual_tile)

  print*, '--------------END OF SIMULATION---------------'

  100  format (A,I6,I6,F8.2)
  777  format (F20.8,F20.8)
  999  format (I4.4)

contains

  subroutine populate_outarray_hourly_tile( hourly_tile, out_hourly_tile ) !, idx_daily_start, idx_daily_end

    use md_interface, only: outtype_hourly_tile

    ! arguments
    type(outtype_hourly_tile), dimension(ntstepsyear), intent(in) :: hourly_tile    ! dimension(ntstepsyear)
    real, dimension(ntstepsyear, nvars_hourly_tile), intent(inout) :: out_hourly_tile

    out_hourly_tile(:, 1)  = hourly_tile(:)%year
    out_hourly_tile(:, 2)  = hourly_tile(:)%doy
    out_hourly_tile(:, 3)  = hourly_tile(:)%hour
    out_hourly_tile(:, 4)  = hourly_tile(:)%rad
    out_hourly_tile(:, 5)  = hourly_tile(:)%Tair
    out_hourly_tile(:, 6)  = hourly_tile(:)%Prcp
    out_hourly_tile(:, 7)  = hourly_tile(:)%GPP
    out_hourly_tile(:, 8)  = hourly_tile(:)%Resp
    out_hourly_tile(:, 9)  = hourly_tile(:)%Transp
    out_hourly_tile(:, 10) = hourly_tile(:)%Evap
    out_hourly_tile(:, 11) = hourly_tile(:)%Runoff
    out_hourly_tile(:, 12) = hourly_tile(:)%Soilwater
    out_hourly_tile(:, 13) = hourly_tile(:)%wcl
    out_hourly_tile(:, 14) = hourly_tile(:)%FLDCAP
    out_hourly_tile(:, 15) = hourly_tile(:)%WILTPT

  end subroutine populate_outarray_hourly_tile


  subroutine populate_outarray_daily_tile( daily_tile, out_daily_tile ) !, idx_daily_start, idx_daily_end

    use md_interface, only: outtype_daily_tile

    ! arguments
    type(outtype_daily_tile), dimension(ndayyear), intent(in) :: daily_tile
    real, dimension(ndayyear, nvars_daily_tile), intent(inout) :: out_daily_tile

    out_daily_tile(:, 1)  = daily_tile(:)%year 
    out_daily_tile(:, 2)  = daily_tile(:)%doy
    out_daily_tile(:, 3)  = daily_tile(:)%Tc
    out_daily_tile(:, 4)  = daily_tile(:)%Prcp
    out_daily_tile(:, 5)  = daily_tile(:)%totWs
    out_daily_tile(:, 6)  = daily_tile(:)%Trsp
    out_daily_tile(:, 7)  = daily_tile(:)%Evap
    out_daily_tile(:, 8)  = daily_tile(:)%Runoff
    out_daily_tile(:, 9)  = daily_tile(:)%ws1
    out_daily_tile(:, 10) = daily_tile(:)%ws2
    out_daily_tile(:, 11) = daily_tile(:)%ws3
    out_daily_tile(:, 12) = daily_tile(:)%LAI
    out_daily_tile(:, 13) = daily_tile(:)%GPP
    out_daily_tile(:, 14) = daily_tile(:)%Rauto
    out_daily_tile(:, 15) = daily_tile(:)%Rh
    out_daily_tile(:, 16) = daily_tile(:)%NSC
    out_daily_tile(:, 17) = daily_tile(:)%seedC
    out_daily_tile(:, 18) = daily_tile(:)%leafC
    out_daily_tile(:, 19) = daily_tile(:)%rootC
    out_daily_tile(:, 20) = daily_tile(:)%SW_C
    out_daily_tile(:, 21) = daily_tile(:)%HW_C
    out_daily_tile(:, 22) = daily_tile(:)%NSN
    out_daily_tile(:, 23) = daily_tile(:)%seedN
    out_daily_tile(:, 24) = daily_tile(:)%leafN
    out_daily_tile(:, 25) = daily_tile(:)%rootN
    out_daily_tile(:, 26) = daily_tile(:)%SW_N
    out_daily_tile(:, 27) = daily_tile(:)%HW_N
    out_daily_tile(:, 28) = daily_tile(:)%McrbC
    out_daily_tile(:, 29) = daily_tile(:)%fastSOM
    out_daily_tile(:, 30) = daily_tile(:)%slowSOM
    out_daily_tile(:, 31) = daily_tile(:)%McrbN
    out_daily_tile(:, 32) = daily_tile(:)%fastSoilN
    out_daily_tile(:, 33) = daily_tile(:)%slowSoilN
    out_daily_tile(:, 34) = daily_tile(:)%mineralN
    out_daily_tile(:, 35) = daily_tile(:)%N_uptk

  end subroutine populate_outarray_daily_tile


  subroutine populate_outarray_daily_cohorts( daily_cohorts, out_daily_cohorts ) 

    use md_interface, only: outtype_daily_cohorts

    ! arguments
    type(outtype_daily_cohorts), dimension(ndayyear, out_max_cohorts), intent(in) :: daily_cohorts
    real, dimension(ndayyear, out_max_cohorts,nvars_daily_cohorts), intent(inout) :: out_daily_cohorts

    out_daily_cohorts(:,:, 1)  = daily_cohorts(:,:)%year
    out_daily_cohorts(:,:, 2)  = daily_cohorts(:,:)%doy
    out_daily_cohorts(:,:, 3)  = daily_cohorts(:,:)%hour
    out_daily_cohorts(:,:, 4)  = daily_cohorts(:,:)%cID
    out_daily_cohorts(:,:, 5)  = daily_cohorts(:,:)%PFT
    out_daily_cohorts(:,:, 6)  = daily_cohorts(:,:)%layer
    out_daily_cohorts(:,:, 7)  = daily_cohorts(:,:)%density
    out_daily_cohorts(:,:, 8)  = daily_cohorts(:,:)%f_layer
    out_daily_cohorts(:,:, 9)  = daily_cohorts(:,:)%LAI
    out_daily_cohorts(:,:, 10) = daily_cohorts(:,:)%gpp
    out_daily_cohorts(:,:, 11) = daily_cohorts(:,:)%resp
    out_daily_cohorts(:,:, 12) = daily_cohorts(:,:)%transp
    out_daily_cohorts(:,:, 13) = daily_cohorts(:,:)%NPPleaf
    out_daily_cohorts(:,:, 14) = daily_cohorts(:,:)%NPProot
    out_daily_cohorts(:,:, 15) = daily_cohorts(:,:)%NPPwood    
    out_daily_cohorts(:,:, 16) = daily_cohorts(:,:)%NSC
    out_daily_cohorts(:,:, 17) = daily_cohorts(:,:)%seedC
    out_daily_cohorts(:,:, 18) = daily_cohorts(:,:)%leafC
    out_daily_cohorts(:,:, 19) = daily_cohorts(:,:)%rootC
    out_daily_cohorts(:,:, 20) = daily_cohorts(:,:)%SW_C
    out_daily_cohorts(:,:, 21) = daily_cohorts(:,:)%HW_C
    out_daily_cohorts(:,:, 22) = daily_cohorts(:,:)%NSN
    out_daily_cohorts(:,:, 23) = daily_cohorts(:,:)%seedN
    out_daily_cohorts(:,:, 24) = daily_cohorts(:,:)%leafN
    out_daily_cohorts(:,:, 25) = daily_cohorts(:,:)%rootN
    out_daily_cohorts(:,:, 26) = daily_cohorts(:,:)%SW_N
    out_daily_cohorts(:,:, 27) = daily_cohorts(:,:)%HW_N

  end subroutine populate_outarray_daily_cohorts


  subroutine populate_outarray_annual_tile( annual_tile, out_annual_tile )

    use md_interface, only: outtype_annual_tile

    ! arguments
    type(outtype_annual_tile), intent(in) :: annual_tile
    real, dimension(nvars_annual_tile), intent(inout) :: out_annual_tile

    out_annual_tile(1)  = annual_tile%year
    out_annual_tile(2)  = annual_tile%CAI
    out_annual_tile(3)  = annual_tile%LAI
    out_annual_tile(4)  = annual_tile%GPP
    out_annual_tile(5)  = annual_tile%Rauto
    out_annual_tile(6)  = annual_tile%Rh
    out_annual_tile(7)  = annual_tile%rain
    out_annual_tile(8)  = annual_tile%SoilWater
    out_annual_tile(9)  = annual_tile%Transp
    out_annual_tile(10) = annual_tile%Evap
    out_annual_tile(11) = annual_tile%Runoff
    out_annual_tile(12) = annual_tile%plantC
    out_annual_tile(13) = annual_tile%soilC
    out_annual_tile(14) = annual_tile%plantN
    out_annual_tile(15) = annual_tile%soilN
    out_annual_tile(16) = annual_tile%totN
    out_annual_tile(17) = annual_tile%NSC
    out_annual_tile(18) = annual_tile%SeedC
    out_annual_tile(19) = annual_tile%leafC
    out_annual_tile(20) = annual_tile%rootC
    out_annual_tile(21) = annual_tile%SapwoodC
    out_annual_tile(22) = annual_tile%WoodC
    out_annual_tile(23) = annual_tile%NSN
    out_annual_tile(24) = annual_tile%SeedN
    out_annual_tile(25) = annual_tile%leafN
    out_annual_tile(26) = annual_tile%rootN
    out_annual_tile(27) = annual_tile%SapwoodN
    out_annual_tile(28) = annual_tile%WoodN
    out_annual_tile(29) = annual_tile%McrbC
    out_annual_tile(30) = annual_tile%fastSOM
    out_annual_tile(31) = annual_tile%SlowSOM
    out_annual_tile(32) = annual_tile%McrbN
    out_annual_tile(33) = annual_tile%fastSoilN
    out_annual_tile(34) = annual_tile%slowSoilN
    out_annual_tile(35) = annual_tile%mineralN
    out_annual_tile(36) = annual_tile%N_fxed
    out_annual_tile(37) = annual_tile%N_uptk
    out_annual_tile(38) = annual_tile%N_yrMin
    out_annual_tile(39) = annual_tile%N_P2S
    out_annual_tile(40) = annual_tile%N_loss
    out_annual_tile(41) = annual_tile%totseedC
    out_annual_tile(42) = annual_tile%totseedN
    out_annual_tile(43) = annual_tile%Seedling_C
    out_annual_tile(44) = annual_tile%Seedling_N

  end subroutine populate_outarray_annual_tile


subroutine populate_outarray_annual_cohorts( annual_cohorts, out_annual_cohorts ) 

    use md_interface, only: outtype_annual_cohorts

    ! arguments
    type(outtype_annual_cohorts), dimension(out_max_cohorts), intent(in) :: annual_cohorts
    real, dimension(out_max_cohorts,nvars_annual_cohorts), intent(inout) :: out_annual_cohorts    

    out_annual_cohorts(:, 1)  = annual_cohorts(:)%year
    out_annual_cohorts(:, 2)  = annual_cohorts(:)%cID
    out_annual_cohorts(:, 3)  = annual_cohorts(:)%PFT
    out_annual_cohorts(:, 4)  = annual_cohorts(:)%layer
    out_annual_cohorts(:, 5)  = annual_cohorts(:)%density
    out_annual_cohorts(:, 6)  = annual_cohorts(:)%f_layer
    out_annual_cohorts(:, 7)  = annual_cohorts(:)%dDBH
    out_annual_cohorts(:, 8)  = annual_cohorts(:)%dbh
    out_annual_cohorts(:, 9)  = annual_cohorts(:)%height
    out_annual_cohorts(:, 10)  = annual_cohorts(:)%Acrown
    out_annual_cohorts(:, 11) = annual_cohorts(:)%wood
    out_annual_cohorts(:, 12) = annual_cohorts(:)%nsc
    out_annual_cohorts(:, 13) = annual_cohorts(:)%NSN
    out_annual_cohorts(:, 14) = annual_cohorts(:)%NPPtr
    out_annual_cohorts(:, 15) = annual_cohorts(:)%seed
    out_annual_cohorts(:, 16) = annual_cohorts(:)%NPPL
    out_annual_cohorts(:, 17) = annual_cohorts(:)%NPPR
    out_annual_cohorts(:, 18) = annual_cohorts(:)%NPPW
    out_annual_cohorts(:, 19) = annual_cohorts(:)%GPP
    out_annual_cohorts(:, 20) = annual_cohorts(:)%NPP
    out_annual_cohorts(:, 21) = annual_cohorts(:)%N_uptk
    out_annual_cohorts(:, 22) = annual_cohorts(:)%N_fix
    out_annual_cohorts(:, 23) = annual_cohorts(:)%maxLAI

  end subroutine populate_outarray_annual_cohorts

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
    integer :: idx_climatedata

    character(len=80) :: filepath_in = '/Users/benjaminstocker/sofun/input/'
    character(len=80) :: climfile    = 'ORNL_forcing.txt'

    climfile=trim(filepath_in)//trim(climfile)
    write(*,*)'inputfile: ',climfile

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
      write(*,*)"the data frequency is hourly"
    elseif(timestep==0.5)then
      write(*,*)"the data frequency is half hourly"
    else
      write(*,*)"Please check time step!"
      stop
    endif
    close(11)    ! close forcing file

    ! print*,'1'
    ! print*,'input par', input_data(1,10)

    ! Put the data into forcing 
    datalines = m - 1
    days_data = idays
    yr_data  = year_data(datalines-1) - year_data(1) + 1

    ! print*,'2'  
    allocate(climateData(datalines - 72))  !3*24
    days_data = days_data - 3

    ! print*,'3'
    ! print*,'datalines', datalines
    ! print*,'size(input_data)', shape(input_data)
    ! print*,'length(climateData)', size(climateData)
    idx_climatedata = 0
    do i=1,datalines
      ! print*,'i, doy_data(i), year_data(i)', i, doy_data(i), year_data(i)
     if (.not. (doy_data(i)==60 .and. (year_data(i)==2000 .or. year_data(i)==2004 .or. year_data(i)==2008))) then
       
       idx_climatedata = idx_climatedata + 1 !xxx debug 

       climateData(idx_climatedata)%year      = year_data(i)          ! Year
       climateData(idx_climatedata)%doy       = doy_data(i)           ! day of the year
       climateData(idx_climatedata)%hod       = hour_data(i)          ! hour of the day
       climateData(idx_climatedata)%PAR       = input_data(1,i)       ! umol/m2/s
       climateData(idx_climatedata)%radiation = input_data(2,i)       ! W/m2
       climateData(idx_climatedata)%Tair      = input_data(3,i) + 273.16  ! air temperature, K
       climateData(idx_climatedata)%Tsoil     = input_data(4,i) + 273.16  ! soil temperature, K
       climateData(idx_climatedata)%RH        = input_data(5,i) * 0.01    ! relative humidity (0.xx)
       climateData(idx_climatedata)%rain      = input_data(6,i)/(timestep * 3600)! ! kgH2O m-2 s-1
       climateData(idx_climatedata)%windU     = input_data(7,i)        ! wind velocity (m s-1)
       climateData(idx_climatedata)%P_air     = input_data(8,i)        ! pa
       climateData(idx_climatedata)%CO2       = input_data(9,i) * 1.0e-6       ! mol/mol
       climateData(idx_climatedata)%soilwater = 0.8    ! soil moisture, vol/vol
      else
     end if
   enddo
   forcingData => climateData

   ! print*,'4'
   datalines = datalines - 72  !3*24

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
    integer :: idx_climatedata

    character(len=80) :: filepath_in = '/Users/benjaminstocker/sofun/input/'
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

      end if

    enddo
    forcingData => climateData

    datalines = datalines - 96

    write(*,*)"forcing", datalines,days_data,yr_data

  end subroutine read_NACPforcing

end program main