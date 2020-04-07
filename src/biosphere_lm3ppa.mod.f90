module md_biosphere

  use datatypes
  use md_vegetation
  use md_soil
  use md_params_core

  implicit none

  private
  public biosphere_annual

  type(vegn_tile_type),  pointer :: vegn
  type(soil_tile_type),  pointer :: soil
  type(cohort_type),     pointer :: cx, cc

contains

  subroutine biosphere_annual(out_biosphere)
    !////////////////////////////////////////////////////////////////
    ! function BIOSPHERE_annual calculates net ecosystem exchange (nee)
    ! in response to environmental boundary conditions (atmospheric 
    ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
    ! LPJ, also formulated as subroutine.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_interface, only: myinterface, outtype_biosphere
    use md_params_core, only: ntstepsyear !diff
  
    ! return variable
    type(outtype_biosphere), intent(inout) :: out_biosphere

    ! ! local variables
    integer :: dm, moy, jpngr, doy
    ! logical, save           :: init_daily = .true.   ! is true only on the first day of the simulation 
    logical, parameter :: verbose = .false.       ! change by hand for debugging etc.

    !----------------------------------------------------------------
    ! Biome-E stuff
    !----------------------------------------------------------------
    integer, parameter :: rand_seed = 86456
    integer, parameter :: totalyears = 10
    integer, parameter :: nCohorts = 1
    ! integer :: steps_per_day ! 24 or 48
    real    :: tsoil, soil_theta
    real    :: NPPtree,fseed, fleaf, froot, fwood ! for output
    real    :: dDBH ! yearly growth of DBH, mm
    real    :: plantC,plantN, soilC, soilN
    real    :: dSlowSOM  ! for multiple tests only
    character(len=150) :: plantcohorts, plantCNpools, soilCNpools, allpools, faststepfluxes  ! output file names
    logical :: new_annual_cycle = .False.
    logical :: switch = .True.
    integer :: istat1, istat2, istat3
    integer :: year0,  year1
    integer :: fno1, fno2, fno3, fno4, fno5 ! output files
    integer :: totyears
    integer :: i, j, k
    integer :: idata
    integer, save :: simu_steps 
    integer, save :: iyears
    integer, save :: idays
    integer, save :: idoy
    character(len=50) :: filepath_out, filesuffix
    character(len=50) :: parameterfile(10), chaSOM(10)

    ! benjaminstocker debug
    integer :: idx, forcingyear
    !integer :: ntstepsyear          ! number of days in a year


    !------------------------------------------------------------------------
    ! Create output files
    ! XXX add this to output instead
    !------------------------------------------------------------------------
    filepath_out   = '/Users/benjaminstocker/sofun/output/'
    filesuffix     = '_test.csv' ! tag for simulation experiments
    plantcohorts   = trim(filepath_out)//'Annual_cohorts'//trim(filesuffix)  ! has 22 columns
    plantCNpools   = trim(filepath_out)//'Daily_cohorts'//trim(filesuffix)  ! daily has 27 columns
    soilCNpools    = trim(filepath_out)//'Daily_tile'//trim(filesuffix)
    allpools       = trim(filepath_out)//'Annual_tile'//trim(filesuffix)
    faststepfluxes = trim(filepath_out)//'Hourly_tile'//trim(filesuffix) ! hourly, has 15 columns and 

    fno1=91; fno2=101; fno3=102; fno4=103; fno5=104
    open(fno1, file=trim(faststepfluxes), ACTION='write', IOSTAT=istat1)
    open(fno2, file=trim(plantcohorts),   ACTION='write', IOSTAT=istat1)
    open(fno3, file=trim(plantCNpools),   ACTION='write', IOSTAT=istat2)
    open(fno4, file=trim(soilCNpools),    ACTION='write', IOSTAT=istat3)
    open(fno5, file=trim(allpools),       ACTION='write', IOSTAT=istat3)

    !----------------------------------------------------------------
    ! INITIALISATIONS
    !----------------------------------------------------------------
    if (myinterface%steering%init) then

      !------------------------------------------------------------------------
      ! Translation to LM3-PPA variables
      !------------------------------------------------------------------------

          ! hourly tile
          write(fno1,'(5(a8,","),25(a12,","))')                    &
          'year','doy','hour','rad',                               &
          'Tair','Prcp', 'GPP', 'Resp',                            &
          'Transp','Evap','Runoff','Soilwater',                    &
          'wcl','FLDCAP','WILTPT'

          ! annual cohorts
          write(fno2,'(3(a5,","),25(a9,","))') 'year',             &
          'cID','PFT','layer','density', 'f_layer',                &
          'dDBH','dbh','height','Acrown',                          &
          'wood','nsc', 'NSN','NPPtr','seed',                      &
          'NPPL','NPPR','NPPW','GPP_yr','NPP_yr',                  &
          'N_uptk','N_fix','maxLAI'
          
          ! daily cohorts
          write(fno3,'(5(a5,","),25(a8,","))')                     &
          'year','doy','hour','cID','PFT',                         &
          'layer','density', 'f_layer', 'LAI',                     &
          'gpp','resp','transp',                                   &
          'NPPleaf','NPProot','NPPwood',                           &
          'NSC','seedC','leafC','rootC','SW_C','HW_C',             &
          'NSN','seedN','leafN','rootN','SW_N','HW_N'
          
          ! daily tile
          write(fno4,'(2(a5,","),55(a10,","))')  'year','doy',     &
          'Tc','Prcp', 'totWs',  'Trsp', 'Evap','Runoff',          &
          'ws1','ws2','ws3', 'LAI','GPP', 'Rauto', 'Rh',           &
          'NSC','seedC','leafC','rootC','SW_C','HW_C',             &
          'NSN','seedN','leafN','rootN','SW_N','HW_N',             &
          'McrbC', 'fastSOM',   'slowSOM',                         &
          'McrbN', 'fastSoilN', 'slowSoilN',                       &
          'mineralN', 'N_uptk'
          
          ! annual tile
          write(fno5,'(1(a5,","),80(a12,","))')  'year',           &
          'CAI','LAI','GPP', 'Rauto',   'Rh',                      &
          'rain','SoilWater','Transp','Evap','Runoff',             &
          'plantC','soilC',    'plantN', 'soilN','totN',           &
          'NSC', 'SeedC', 'leafC', 'rootC', 'SapwoodC', 'WoodC',   &
          'NSN', 'SeedN', 'leafN', 'rootN', 'SapwoodN', 'WoodN',   &
          'McrbC','fastSOM',   'SlowSOM',                          &
          'McrbN','fastSoilN', 'slowSoilN',                        &
          'mineralN', 'N_fxed','N_uptk','N_yrMin','N_P2S','N_loss',&
          'totseedC','totseedN','Seedling_C','Seedling_N'

      !------------------------------------------------------------------------
      ! Initialisations
      !------------------------------------------------------------------------
      ! Parameter initialization: Initialize PFT parameters
      ! print*,'1: vegn%cohorts ', vegn%n_cohorts
      call initialize_PFT_data()

      ! Initialize vegetation tile and plant cohorts
      allocate(vegn)
      
      ! print*,'2: vegn%cohorts , vegn%CAI, vegn%LAI', vegn%n_cohorts, vegn%CAI, vegn%LAI
      call initialize_vegn_tile(vegn, nCohorts)
      
      ! Sort and relayer cohorts
      ! print*,'3: vegn%cohorts , vegn%CAI, vegn%LAI', vegn%n_cohorts, vegn%CAI, vegn%LAI
      call relayer_cohorts(vegn)

      ! print*,'4: vegn%cohorts , vegn%CAI, vegn%LAI', vegn%n_cohorts, vegn%CAI, vegn%LAI
      call Zero_diagnostics(vegn)

      !------------------------------------------------------------------------
      ! Read in forcing data
      ! This reads it all into memory and then extracts from the huge array in
      ! daily/hourly steps
      !------------------------------------------------------------------------
      ! !call read_FACEforcing(forcingData,datalines,days_data,yr_data,timestep)
      ! call read_NACPforcing(forcingData,datalines,days_data,yr_data,timestep)
      ! steps_per_day = int(24.0/timestep)
      ! dt_fast_yr = 1.0/(365.0 * steps_per_day)
      ! step_seconds = 24.0*3600.0/steps_per_day ! seconds_per_year * dt_fast_yr
      ! write(*,*)steps_per_day,dt_fast_yr,step_seconds

      year0 = myinterface%climate(1)%year  ! forcingData(1)%year
      iyears = 1
      idoy   = 0
      idays  = 0

    endif 

    ! print*, 'year', myinterface%climate(1)%year
    ! print*, 'doy',myinterface%climate(1)%doy
    ! print*, 'hour',myinterface%climate(1)%hod
    ! print*, 'PAR',myinterface%climate(1)%PAR
    ! print*, 'radiation',myinterface%climate(1)%radiation
    ! print*, 'Tair', myinterface%climate(1)%Tair
    ! print*, 'Tsoil', myinterface%climate(1)%Tsoil
    ! print*, 'windU', myinterface%climate(1)%windU
    ! print*, 'P_air', myinterface%climate(1)%P_air
    ! print*, 'CO2', myinterface%climate(1)%CO2
    ! print*, 'soilwater', myinterface%climate(1)%soilwater


    ! every year, start reading from beginning of myinterface%climate (annual chunks are passed!)
    simu_steps = 0

    !----------------------------------------------------------------
    ! LOOP THROUGH MONTHS
    !----------------------------------------------------------------
    doy = 0
    monthloop: do moy=1,nmonth

      !----------------------------------------------------------------
      ! LOOP THROUGH DAYS
      !----------------------------------------------------------------
        ! print*,'5: vegn%cohorts , vegn%CAI, vegn%LAI', vegn%n_cohorts, vegn%CAI, vegn%LAI
      dayloop: do dm=1,ndaymonth(moy)

        doy = doy + 1
        idays = idays + 1

        if (verbose) print*,'----------------------'
        if (verbose) print*,'YEAR, DOY ', myinterface%steering%year, doy
        if (verbose) print*,'----------------------'

        idoy = idoy + 1

        !----------------------------------------------------------------
        ! Get daily mean air and soil temperature
        !----------------------------------------------------------------
        ! vegn%Tc_daily = myinterface%climate(jpngr)%dtemp(doy)
        ! tsoil         =                                       ! soil temp. is prescribed
        ! soil_theta    = 

        ! ForestESS:
        ! get daily mean temperature from hourly/half-hourly data
        vegn%Tc_daily = 0.0
        tsoil         = 0.0

        !----------------------------------------------------------------
        ! FAST LOOP
        !----------------------------------------------------------------
        fastloop: do i=1,myinterface%steps_per_day

          ! idata = MOD(simu_steps, myinterface%datalines)+1
          idata = simu_steps + 1
          year0 =  myinterface%climate(idata)%year  ! Current year
          vegn%Tc_daily = vegn%Tc_daily +  myinterface%climate(idata)%Tair
          tsoil         = myinterface%climate(idata)%tsoil
          simu_steps    = simu_steps + 1

          ! print*,'5.0.1'
          call vegn_CNW_budget_fast(vegn, myinterface%climate(idata))

          ! ! benjaminstocker consistency check
          ! myinterface%climate(idata)%tsoil = myinterface%climate(idata)%Tair
          
          ! diagnostics
          ! print*,'5.0.2'
          call hourly_diagnostics(vegn, myinterface%climate(idata), iyears, idoy, i, idays, fno1, out_biosphere%hourly_tile(idata) )

          ! print*,'5.0.3'

        enddo fastloop

        !if (doy == 100 .and. iyears == 3) stop 'consistency?'

        vegn%Tc_daily = vegn%Tc_daily/myinterface%steps_per_day
        tsoil         = tsoil/myinterface%steps_per_day
        soil_theta    = vegn%thetaS

        !-------------------------------------------------
        ! Daily calls
        !-------------------------------------------------
        ! print*,'5: vegn%cohorts , vegn%CAI, vegn%LAI', vegn%n_cohorts, vegn%CAI, vegn%LAI
        call daily_diagnostics(vegn, myinterface%climate(idata), iyears, idoy, idays, fno3, fno4, out_biosphere%daily_cohorts(doy,:), out_biosphere%daily_tile(doy) )

        ! Determine start and end of season and maximum leaf (root) mass
        ! print*,'6: vegn%cohorts , vegn%CAI, vegn%LAI', vegn%n_cohorts, vegn%CAI, vegn%LAI
        call vegn_phenology(vegn, j)


        ! Kill all individuals of a cohort if NSC falls below threshold
        !call vegn_starvation(vegn)

        ! Produce new biomass from 'carbon_gain' (is zero afterwards)
        ! print*,'7: vegn%cohorts , vegn%CAI, vegn%LAI', vegn%n_cohorts, vegn%CAI, vegn%LAI
        call vegn_growth_EW(vegn)

      end do dayloop

    end do monthloop

    !----------------------------------------------------------------
    ! Annual calls
    !----------------------------------------------------------------
    idoy = 0

    print*,'sim. year  ', iyears
    print*,'real year: ', year0

    if (update_annualLAImax) call vegn_annualLAImax_update(vegn)
      
    ! print*,'8: vegn%cohorts , vegn%CAI, vegn%LAI', vegn%n_cohorts, vegn%CAI, vegn%LAI
    call annual_diagnostics(vegn, iyears, fno2, fno5, out_biosphere%annual_cohorts(:), out_biosphere%annual_tile)

    !---------------------------------------------
    ! Reproduction and mortality
    !---------------------------------------------        
    ! Kill all individuals in a cohort if NSC falls below critical point
    call vegn_annual_starvation(vegn)

    ! Natural mortality (reducing number of individuals 'nindivs')
    ! (~Eq. 2 in Weng et al., 2015 BG)
    call vegn_nat_mortality(vegn, real( seconds_per_year ))
    ! print*,'6: ', vegn%n_cohorts

    ! seed C and germination probability (~Eq. 1 in Weng et al., 2015 BG)
    call vegn_reproduction(vegn)
    ! print*,'7: ', vegn%n_cohorts

     !---------------------------------------------
    ! Re-organize cohorts
    !---------------------------------------------
    call kill_lowdensity_cohorts(vegn)
    ! print*,'8: ', vegn%n_cohorts
    call relayer_cohorts(vegn)
    ! print*,'9: ', vegn%n_cohorts
    call vegn_mergecohorts(vegn)
    ! print*,'10: ', vegn%n_cohorts

    ! !---------------------------------------------
    ! ! Set annual variables zero
    ! !---------------------------------------------
    call Zero_diagnostics(vegn)
    ! print*,'11: ', vegn%n_cohorts

    ! update the years of model run
    iyears = iyears + 1

    if (myinterface%steering%finalize) then
      !----------------------------------------------------------------
      ! Finazlize run: deallocating memory
      !----------------------------------------------------------------
      close(91)
      close(101)
      close(102)
      close(103)
      close(104)
      deallocate(vegn%cohorts)
      !deallocate(myinterface%climate) !differ
      !deallocate(out_biosphere%hourly_tile)
      ! stop 'actually finalizing'

    end if

    if (verbose) print*,'Done with biosphere for this year. Guete Rutsch!'

  end subroutine biosphere_annual

end module md_biosphere
