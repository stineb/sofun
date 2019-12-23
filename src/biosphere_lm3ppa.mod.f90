module md_biosphere

  ! use md_params_core
  ! use md_classdefs
  ! use md_plant, only: plant_type, plant_fluxes_type, initdaily_plant, initglobal_plant, getout_daily_plant, getout_annual_plant, getpar_modl_plant, initoutput_plant
  ! use md_params_soil, only: paramtype_soil
  ! use md_waterbal, only: solartype, waterbal, get_solar, getout_daily_waterbal, initoutput_waterbal, getpar_modl_waterbal, initio_nc_waterbal, writeout_nc_waterbal, init_rlm_waterbal, get_rlm_waterbal, getrlm_daily_waterbal
  ! use md_gpp, only: outtype_pmodel, getpar_modl_gpp, initoutput_gpp, gpp, getout_daily_gpp, getout_annual_gpp, initio_nc_gpp, writeout_nc_gpp
  ! use md_vegdynamics, only: vegdynamics
  ! use md_tile, only: tile_type, tile_fluxes_type, initglobal_tile, initdaily_tile
  ! use md_interface, only: getout_daily_forcing, initoutput_forcing, initio_nc_forcing, writeout_nc_forcing
  ! use md_soiltemp, only: getout_daily_soiltemp, soiltemp, initoutput_soiltemp
  ! use md_sofunutils, only: calc_patm

  use esdvm_mod
  use md_params_core

  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! ForestESS stuff
  !----------------------------------------------------------------
  type(tile_type),   pointer :: vegn
  type(cohort_type), pointer :: cp,cc

  ! !----------------------------------------------------------------
  ! ! Module-specific (private) variables
  ! !----------------------------------------------------------------
  ! ! derived types from L1 modules
  ! type( tile_type ),         allocatable, dimension(:,:) :: tile
  ! type( tile_fluxes_type ),  allocatable, dimension(:)   :: tile_fluxes
  ! type( plant_type ),        allocatable, dimension(:,:) :: plant
  ! type( plant_fluxes_type ), allocatable, dimension(:)   :: plant_fluxes

  ! ! derived types from L2 modules
  ! type( solartype )                              :: solar
  ! type( outtype_pmodel ), dimension(npft,nmonth) :: out_pmodel ! P-model output variables for each month and PFT determined beforehand (per unit fAPAR and PPFD only)

contains

  function biosphere_annual() result( out_biosphere )
    !////////////////////////////////////////////////////////////////
    ! function BIOSPHERE_annual calculates net ecosystem exchange (nee)
    ! in response to environmental boundary conditions (atmospheric 
    ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
    ! LPJ, also formulated as subroutine.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_interface, only: interface, outtype_biosphere
  
    ! return variable
    type(outtype_biosphere) :: out_biosphere

    ! ! local variables
    integer :: dm, moy, jpngr, doy
    ! logical, save           :: init_daily = .true.   ! is true only on the first day of the simulation 
    logical, parameter      :: verbose = .false.       ! change by hand for debugging etc.

    !----------------------------------------------------------------
    ! ForestESS stuff
    !----------------------------------------------------------------
    integer,parameter :: rand_seed = 86456
    integer,parameter :: totalyears = 10
    integer,parameter :: nCohorts = 1
    integer           :: datalines                                          ! the total lines in forcing data file
    integer           :: yr_data                                            ! Years of the forcing data
    integer           :: days_data                                          ! days of the forcing data
    integer           :: steps_per_day                                      ! 24 or 48
    real              :: timestep                                           ! hour, Time step of forcing data, usually hourly (1.0)
    real              :: tsoil, soil_theta
    real              :: NPPtree, fseed, fleaf, froot, fwood                 ! for output
    real              :: dDBH                                               ! yearly growth of DBH, mm
    character(len=50) :: plantcohorts, plantCNpools, soilCNpools, allpools  ! output file names
    logical           :: new_annual_cycle = .False.
    integer           :: istat1, istat2, istat3
    integer           :: year0, year1, iyears
    integer           :: totyears, totdays
    integer           :: i, j, k, idays, idoy
    integer           :: idata
    integer, save     :: simu_steps

    ! xxx debug
    integer :: idx, forcingyear

    !------------------------------------------------------------------------
    ! Create output files
    ! XXX add this to output instead
    !------------------------------------------------------------------------
    plantcohorts = 'Annual_cohorts.txt'
    plantCNpools = 'Plant_C_N_pools.csv'  ! daily
    soilCNpools  = 'Soil_C_N_pools.csv'
    allpools     = 'AnnualEcosystemDynamics.csv'

    open(101,file=plantcohorts, ACTION='write', IOSTAT=istat1)
    open(102,file=plantCNpools, ACTION='write', IOSTAT=istat2)
    open(103,file=soilCNpools,  ACTION='write', IOSTAT=istat3)
    open(104,file=allpools,     ACTION='write', IOSTAT=istat3)
    
    ! write header
    write(101,'(3(a5,","),25(a9,","))')               &
    'cID','PFT','layer','density', 'f_layer',    &
    'dDBH','dbh','height','Acrown','totwood',    &
    'nsc', 'NSN',                                &
    'NPPtr','f_seed','f_leaf','f_root','f_wood', &
    'GPP-yr','NPP-yr','Ra_yr','N_uptk','maxLAI'


    write(102,'(5(a5,","),25(a8,","))')               &
    'year','doy','cID','PFT',                    &
    'layer','density', 'f_layer', 'LAI',         &
    'NSC','seedC','leafC','rootC','SW-C','HW-C', &
    'NSN','seedN','leafN','rootN','SW-N','HW-N'

    write(103,'(2(a5,","),25(a8,","))')  'year','doy',         &
    'GPP', 'NPP', 'Rh',   &
    'McrbC', 'fineL', 'struL', 'McrbN', 'fineN', 'struN', &
    'mineralN', 'N_uptk'

    write(104,'(1(a5,","),25(a8,","))')  'year',               &
    'GPP', 'NPP', 'Rh',                                   &
    'NSC','SeedC','leafC','rootC','SW-C', 'HW-C','LAI',   &
    'McrbC', 'fineL', 'struL', 'McrbN', 'fineN', 'struN', &
    'mineralN', 'N_uptk'


    !------------------------------------------------------------------------
    ! Read in forcing data
    ! Requires data frame with one row for each time step (hourly) and columns:
    !
    !   YEAR     : year
    !   DOY      : day of the year  
    !   HOUR     : hour of the day
    !   PAR      : umol/m2/s
    !   Swdown   : W/m2         
    !   TEMP     : air temperature, deg C       
    !   STEMP    : soil temperature, deg C  
    !   RH       : relative humidity, %     
    !   RAIN     : kgH2O m-2 s-1
    !   WIND     : wind velocity (m s-1)       
    !   PRESSURE : Pa         
    !   aCO2_AW  : ???
    !   amb_co2  : ???       
    !
    ! To define array 'forcingData'
    ! This is how it's read from file input/Temperate_forcing.txt:
    !
    ! forcingData(i)%year      = year_data(i)                     ! Year
    ! forcingData(i)%doy       = doy_data(i)                      ! day of the year
    ! forcingData(i)%hod       = hour_data(i)                     ! hour of the day
    ! forcingData(i)%PAR       = input_data(1,i)                  ! umol/m2/s
    ! forcingData(i)%radiation = input_data(2,i)                  ! W/m2
    ! forcingData(i)%Tair      = input_data(3,i) + 273.16         ! air temperature, K
    ! forcingData(i)%Tsoil     = input_data(4,i) + 273.16         ! soil temperature, K
    ! forcingData(i)%RH        = input_data(5,i)                  ! relative humidity
    ! forcingData(i)%rain      = input_data(6,i)/(timestep * 3600)! kgH2O m-2 s-1
    ! forcingData(i)%windU     = input_data(7,i)                  ! wind velocity (m s-1)
    ! forcingData(i)%pressure  = input_data(8,i)                  ! pa
    ! forcingData(i)%soilwater = 0.8                              ! soil moisture, vol/vol
    !------------------------------------------------------------------------
    forcingyear = modulo(interface%steering%year - 1, interface%params_siml%recycle) + 1998
    call read_forcingdata( forcingData, datalines, days_data, yr_data, timestep, forcingyear)  ! interface%steering%climateyear
    steps_per_day = int( 24.0 / timestep )

    !----------------------------------------------------------------
    ! INITIALISATIONS
    !----------------------------------------------------------------
    if (interface%steering%init) then

      !------------------------------------------------------------------------
      ! Initialisations
      !------------------------------------------------------------------------
      ! Parameter initialization: Initialize PFT parameters
      call initialize_PFT_data()

      ! print*,'1'

      ! Initialize vegetation tile and plant cohorts
      allocate( vegn )
      call initialize_vegn_tile( vegn, nCohorts )

      ! total years of model run
      totyears = model_run_years
      totdays  = INT( totyears / yr_data + 1) * days_data

      ! Sort and relayer cohorts
      call relayer_cohorts(vegn)

      ! ----- model run ----------
      year0      = forcingData(1)%year
      iyears     = 1
      idoy       = 0
      simu_steps = 0
      call vegn_annual_diagnostics_zero(vegn)

      ! output initial cohorts
      write(*,*)"Initial cohorts:", vegn%n_cohorts
      write(*,'(3(I6,","))')0, vegn%n_cohorts
      do i=1,vegn%n_cohorts
        cc => vegn%cohorts(i)
        write(*,*)i,vegn%cohorts(i)%species,vegn%cohorts(i)%height,vegn%cohorts(i)%crownarea
        dDBH = 0.0  ! mm in diameter
        write(*,'(3(I5,","),1(F9.1,","),25(F9.3,","))')          &
              cc%ccID,cc%species,cc%layer,                       &
              cc%nindivs*10000, cc%layerfrac, dDBH,              &
              cc%dbh,cc%height,cc%crownarea,                     &
              cc%bsw+cc%bHW,cc%nsc,cc%NSN,                       &
              NPPtree,fseed, fleaf, froot, fwood,                &
              cc%annualGPP,cc%annualNPP, cc%annualResp,          &
              cc%N_up_yr*1000,spdata(cc%species)%laimax
      enddo

      !   !----------------------------------------------------------------
      !   ! GET MODEL PARAMETERS
      !   ! read model parameters that may be varied for optimisation
      !   !----------------------------------------------------------------
      !   if (verbose) print*, 'getpar_modl() ...'
      !   call getpar_modl_plant()
      !   call getpar_modl_waterbal()
      !   call getpar_modl_gpp()
      !   if (verbose) print*, '... done'

      !   !----------------------------------------------------------------
      !   ! Initialise pool variables and/or read from restart file (not implemented)
      !   !----------------------------------------------------------------
      !   if (verbose) print*, 'initglobal_() ...'
      !   allocate( tile(  nlu,  size(interface%grid) ) )
      !   allocate( tile_fluxes(  nlu ) )
      !   allocate( plant( npft, size(interface%grid) ) )
      !   allocate( plant_fluxes( npft ) )

      !   call initglobal_plant( plant(:,:), size(interface%grid) )
      !   call initglobal_tile(  tile(:,:),  size(interface%grid) )
      !   if (verbose) print*, '... done'

    endif 

    ! !----------------------------------------------------------------
    ! ! Open NetCDF output files (one for each year)
    ! !----------------------------------------------------------------
    ! if (.not.interface%params_siml%is_calib) then
    !   if (verbose) print*, 'initio_nc_() ...'
    !   call initio_nc_forcing()
    !   call initio_nc_gpp()
    !   call initio_nc_waterbal()
    !   if (verbose) print*, '... done'
    ! end if
    
    ! !----------------------------------------------------------------
    ! ! Initialise output variables for this year
    ! !----------------------------------------------------------------
    ! if (.not.interface%params_siml%is_calib) then
    !   if (verbose) print*, 'initoutput_() ...'
    !   call initoutput_waterbal( size(interface%grid) )
    !   call initoutput_gpp(      size(interface%grid) )
    !   call initoutput_plant(    size(interface%grid) )
    !   call initoutput_forcing(  size(interface%grid) )
    !   call initoutput_soiltemp( size(interface%grid) )
    !   if (verbose) print*, '... done'
    ! end if

    ! ! additional initialisation for rolling annual mean calculations (also needed in calibration mode)
    ! call init_rlm_waterbal( size(interface%grid) )

    !----------------------------------------------------------------
    ! LOOP THROUGH GRIDCELLS
    !----------------------------------------------------------------
    if (verbose) print*,'looping through gridcells ...'
    gridcellloop: do jpngr=1,size(interface%grid)

      if (interface%grid(jpngr)%dogridcell) then

        ! if (verbose) print*,'----------------------'
        ! if (verbose) print*,'JPNGR: ', jpngr
        ! if (verbose) print*,'----------------------'

        ! !----------------------------------------------------------------
        ! ! Get radiation based on daily temperature, sunshine fraction, and 
        ! ! elevation.
        ! ! This is not compatible with a daily biosphere-climate coupling. I.e., 
        ! ! there is a daily loop within 'get_solar'!
        ! !----------------------------------------------------------------
        ! if (verbose) print*,'calling get_solar() ... '
        ! if (verbose) print*,'    with argument lat = ', interface%grid(jpngr)%lat
        ! if (verbose) print*,'    with argument elv = ', interface%grid(jpngr)%elv
        ! if (verbose) print*,'    with argument dfsun (ann. mean) = ', sum( interface%climate(jpngr)%dfsun(:) / ndayyear )
        ! if (verbose) print*,'    with argument dppfd (ann. mean) = ', sum( interface%climate(jpngr)%dppfd(:) / ndayyear )
        ! solar = get_solar( &
        !                   interface%grid(jpngr)%lat, & 
        !                   interface%grid(jpngr)%elv, & 
        !                   interface%climate(jpngr)%dfsun(:), & 
        !                   interface%climate(jpngr)%dppfd(:)  &
        !                   )
        ! if (verbose) print*,'... done'

        ! !----------------------------------------------------------------
        ! ! calculate constant atmospheric pressure as a function of elevation
        ! !----------------------------------------------------------------
        ! interface%climate(jpngr)%dpatm(:) = calc_patm(interface%grid(jpngr)%elv)

        !----------------------------------------------------------------
        ! LOOP THROUGH MONTHS
        !----------------------------------------------------------------
        doy=0
        monthloop: do moy=1,nmonth

          !----------------------------------------------------------------
          ! LOOP THROUGH DAYS
          !----------------------------------------------------------------
          dayloop: do dm=1,ndaymonth(moy)
            doy=doy+1

            if (verbose) print*,'----------------------'
            if (verbose) print*,'YEAR, Doy ', interface%steering%year, doy
            if (verbose) print*,'----------------------'

            idoy = idoy + 1

            !----------------------------------------------------------------
            ! Get daily mean air and soil temperature
            !----------------------------------------------------------------
            ! vegn%Tc_daily = interface%climate(jpngr)%dtemp(doy)
            ! tsoil         =                                       ! soil temp. is prescribed
            ! soil_theta    = 

            ! ForestESS:
            ! get daily mean temperature from hourly/half-hourly data
            vegn%Tc_daily = 0.0
            tsoil         = 0.0

            do i=1,steps_per_day

               idata         = MOD(simu_steps, datalines)+1
               year0         = forcingData(idata)%year  ! Current year
               vegn%Tc_daily = vegn%Tc_daily + forcingData(idata)%Tair
               tsoil         = tsoil + forcingData(idata)%tsoil
               simu_steps    = simu_steps + 1

               ! fast-step calls
               ! ***** none ******
               ! leaf photosynthesis
               ! transpiration
               ! plant respiration
               ! soil water dynamics

            enddo

            vegn%Tc_daily = vegn%Tc_daily / steps_per_day
            tsoil         = tsoil / steps_per_day
            soil_theta    = vegn%soil_theta

            ! if (simu_steps==17520) then
            !   print*,'vegn%n_cohorts', vegn%n_cohorts
            !   do idx = 1,vegn%n_cohorts
            !     print*,'idx, nindivs', idx, vegn%cohorts(idx)
            !   end do
            !   print*,'tsoil ', tsoil
            !   stop 'beni'
            ! end if


            !-------------------------------------------------
            ! Daily calls
            !-------------------------------------------------
            ! Determine start and end of season and maximum leaf (root) mass
            call vegn_phenology(vegn,j)

            ! Determine 'carbon_gain' available for growth from NSC
            call vegn_C_N_budget(vegn, tsoil, soil_theta)

            ! Kill all individuals of a cohort if NSC falls below threshold
            call vegn_starvation(vegn)

            ! Produce new biomass from 'carbon_gain' (is zero afterwards)
            call vegn_growth_EW(vegn)

            ! daily output
            do i = 1, vegn%n_cohorts

              cc => vegn%cohorts(i)
              ! cohorts
              write(102,'(5(I5,","),1(F8.1,","),6(F8.3,","),2(F8.2,","),25(F8.2,","))')  &
                  iyears,idoy,cc%ccID,cc%species,cc%layer,   &
                  cc%nindivs*10000, cc%layerfrac, cc%LAI, &
                  cc%NSC, cc%seedC, cc%bl, cc%br, cc%bsw, cc%bHW, &
                  cc%NSN*1000, cc%seedN*1000, cc%leafN*1000, &
                  cc%rootN*1000,cc%sapwdN*1000,cc%woodN*1000

            enddo

            ! Tile level, daily
            write(103,'(2(I5,","),15(F8.4,","))') iyears, idoy,  &
                  vegn%GPP, vegn%NPP, vegn%Rh, &
                  vegn%MicrobialC, vegn%metabolicL, vegn%structuralL, &
                  vegn%MicrobialN*1000, vegn%metabolicN*1000, vegn%structuralN*1000, &
                  vegn%mineralN*1000,   vegn%N_uptake*1000

            ! !----------------------------------------------------------------
            ! ! initialise daily updated variables 
            ! !----------------------------------------------------------------
            ! if (verbose) print*,'calling initdaily_() ...'
            ! call initdaily_plant( plant_fluxes(:) )
            ! call initdaily_tile( tile_fluxes(:) )
            ! if (verbose) print*,'... done.'


            ! !----------------------------------------------------------------
            ! ! update canopy and tile variables and simulate daily 
            ! ! establishment / sprouting
            ! !----------------------------------------------------------------
            ! if (verbose) print*,'calling vegdynamics() ... '
            ! call vegdynamics( tile(:,jpngr), &
            !                   plant(:,jpngr), &
            !                   solar, &
            !                   out_pmodel(:,:), &
            !                   interface%vegcover(jpngr)%dfapar(doy), &
            !                   interface%fpc_grid(:,jpngr) &
            !                   )
            ! if (verbose) print*,'... done'


            ! !----------------------------------------------------------------
            ! ! calculate GPP
            ! !----------------------------------------------------------------
            ! if (verbose) print*,'calling gpp() ... '
            ! call gpp( &
            !           interface%pco2, & 
            !           interface%climate(jpngr)%dtemp(doy), & 
            !           interface%climate(jpngr)%dvpd(doy), & 
            !           interface%climate(jpngr)%dpatm(doy), &
            !           interface%climate(jpngr)%dppfd(doy), &
            !           interface%vegcover(jpngr)%dfapar(doy), &
            !           plant(:,jpngr)%fpc_grid, &
            !           solar%dayl(doy), &
            !           solar%meanmppfd(moy), &
            !           tile(:,jpngr)%soil%phy%wscal, &
            !           tile(:,jpngr)%soil%phy%rlmalpha, &
            !           interface%params_siml%soilmstress, &
            !           interface%params_siml%tempstress, &
            !           plant_fluxes(:)%dgpp, &
            !           plant_fluxes(:)%drd, &
            !           plant_fluxes(:)%dtransp, &
            !           init_daily &
            !           )
            ! if (verbose) print*,'... done'

            ! !----------------------------------------------------------------
            ! ! get soil moisture, and runoff
            ! !----------------------------------------------------------------
            ! if (verbose) print*,'calling waterbal() ... '
            ! ! print*,'lon,lat,ilon,ilat,jpngr', interface%grid(jpngr)%lon, interface%grid(jpngr)%lat, interface%grid(jpngr)%ilon, interface%grid(jpngr)%ilat, jpngr
            ! call waterbal( &
            !               tile(:,jpngr)%soil, &
            !               tile_fluxes(:), &
            !               plant_fluxes(:), &
            !               doy, &
            !               jpngr, & 
            !               interface%grid(jpngr)%lat,             & 
            !               interface%grid(jpngr)%elv,             & 
            !               interface%climate(jpngr)%dprec(doy),   & 
            !               interface%climate(jpngr)%dsnow(doy),   & 
            !               interface%climate(jpngr)%dtemp(doy),   & 
            !               interface%climate(jpngr)%dfsun(doy),   &
            !               interface%climate(jpngr)%dnetrad(doy), &
            !               interface%climate(jpngr)%dvpd(doy),    &
            !               interface%vegcover(jpngr)%dfapar(doy)  &
            !               )
            ! if (verbose) print*,'... done'

            ! ! !----------------------------------------------------------------
            ! ! ! calculate soil temperature
            ! ! !----------------------------------------------------------------
            ! ! if (verbose) print*, 'calling soiltemp() ... '
            ! ! call soiltemp(&
            ! !               tile(:,jpngr)%soil, &
            ! !               interface%climate(jpngr)%dtemp(:), &
            ! !               size(interface%grid), &
            ! !               interface%steering%init, &
            ! !               jpngr, & 
            ! !               moy, & 
            ! !               doy & 
            ! !               )
            ! ! if (verbose) print*, '... done'

            ! !----------------------------------------------------------------
            ! ! collect from daily updated state variables for annual variables
            ! !----------------------------------------------------------------
            ! if (.not.interface%params_siml%is_calib) then
            !   if (verbose) print*,'calling getout_daily() ... '
            !   call getout_daily_waterbal( jpngr, moy, doy, solar, tile(:,jpngr)%soil%phy, tile_fluxes(:) )
            !   call getout_daily_gpp( out_pmodel(:,moy), plant_fluxes(:), jpngr, doy )
            !   call getout_daily_plant( plant(:,jpngr), plant_fluxes(:), jpngr, moy, doy )
            !   call getout_daily_forcing( jpngr, moy, doy )
            !   call getout_daily_soiltemp( jpngr, moy, doy, tile(:,jpngr)%soil%phy )
            !   if (verbose) print*,'... done'
            ! end if

            ! call getrlm_daily_waterbal( jpngr, doy )

            ! !----------------------------------------------------------------
            ! ! populate function return variable
            ! !----------------------------------------------------------------
            ! !if (npft>1) stop 'think about npft > 1'
            ! out_biosphere%fapar(doy)   = plant(1,jpngr)%fapar_ind
            ! out_biosphere%gpp(doy)     = plant_fluxes(1)%dgpp
            ! out_biosphere%transp(doy)  = plant_fluxes(1)%dtransp
            ! out_biosphere%latenth(doy) = plant_fluxes(1)%dlatenth

            ! init_daily = .false.

          end do dayloop

        end do monthloop

        !----------------------------------------------------------------
        ! collect annual output
        !----------------------------------------------------------------

        ! if (.not.interface%params_siml%is_calib) then
        !   if (verbose) print*,'calling getout_annual_() ... '
        !   call getout_annual_plant( plant(:,jpngr), jpngr )
        !   call getout_annual_gpp( jpngr )
        !   if (verbose) print*,'... done'
        ! end if

        ! ! annual calls
        ! idata = MOD(simu_steps+1, datalines)+1
        ! year1 = forcingData(idata)%year         ! Check if it is the last day of a year
        ! new_annual_cycle = ((year0 /= year1).OR. &
        !                  (idata == steps_per_day .and. simu_steps > datalines))

        idoy = 0

        write(101,'(2(I6,","),1(F9.2,","))') iyears, vegn%n_cohorts,vegn%annualN*1000
        write(*,  '(2(I6,","),1(F9.2,","))') iyears, vegn%n_cohorts,vegn%annualN*1000

        ! output yearly variables
        write(*,'(3(a5,","),25(a9,","))') &
            'chtID','PFT','layer','density', 'f_layer',  &
            'dDBH','dbh','height','Acrown', &
            'wood','nsc', 'NSN','NPPtr',     &
            'NPPL','NPPR','NPPW','GPP-yr','NPP-yr','N_uptk','spLAI'

        do i = 1, vegn%n_cohorts

          cc => vegn%cohorts(i)
          NPPtree = cc%seedC + cc%NPPleaf + cc%NPProot + cc%NPPwood

          fseed = cc%seedC/NPPtree
          fleaf = cc%NPPleaf/NPPtree
          froot = cc%NPProot/NPPtree
          fwood = cc%NPPwood/NPPtree

          dDBH = (cc%DBH - cc%DBH_ys) * 1000.0

          write(101,'(3(I5,","),1(F9.1,","),5(F9.3,","),1(F9.2,","),20(F9.3,","))')       &
              cc%ccID,cc%species,cc%layer,                       &
              cc%nindivs*10000, cc%layerfrac, dDBH,              &
              cc%dbh,cc%height,cc%crownarea,                     &
              cc%bsw+cc%bHW,cc%nsc,cc%NSN,                       &
              NPPtree,fseed, fleaf, froot, fwood,                &
              cc%annualGPP,cc%annualNPP, cc%annualResp,          &
              cc%N_up_yr*1000,spdata(cc%species)%laimax

          ! Screen output
          write(*,'(3(I5,","),1(F9.1,","),5(F9.3,","),1(F9.2,","),20(F9.3,","))')       &
              cc%ccID,cc%species,cc%layer,                        &
              cc%nindivs*10000, cc%layerfrac,dDBH,                &
              cc%dbh,cc%height,cc%crownarea,                      &
              cc%bsw+cc%bHW,cc%nsc,cc%NSN,                        &
              fseed, fleaf, froot, fwood,                         &
              cc%annualGPP,cc%annualNPP,                          &
              cc%N_up_yr*1000,      &
              spdata(cc%species)%laimax

          ! Vegn pools:
          vegn%maxNSC     = vegn%maxNSC   + cc%NSC      * cc%nindivs
          vegn%maxSeedC   = vegn%maxSeedC + cc%seedC    * cc%nindivs
          vegn%maxleafC   = vegn%maxleafC + cc%bl       * cc%nindivs
          vegn%maxrootC   = vegn%maxrootC + cc%br       * cc%nindivs
          vegn%SapwoodC   = vegn%SapwoodC + cc%bsw      * cc%nindivs
          vegn%WoodC      = vegn%WoodC    + cc%bHW      * cc%nindivs
          vegn%maxLAI     = vegn%maxLAI   + cc%leafarea * cc%nindivs
              
        enddo

        ! tile pools
        write(104,'(1(I5,","),13(F8.4,","),6(F8.2,","),2(F8.3,","))') &
              iyears,       &
              vegn%annualGPP, vegn%annualNPP, vegn%annualRh, &
              vegn%maxNSC, vegn%maxSeedC, vegn%maxleafC, vegn%maxrootC,  &
              vegn%SapwoodC, vegn%WoodC, vegn%maxLAI,                    &
              vegn%MicrobialC, vegn%metabolicL, vegn%structuralL, &
              vegn%MicrobialN*1000, vegn%metabolicN*1000, vegn%structuralN*1000, &
              vegn%mineralN*1000,   vegn%accu_Nup*1000

        ! 'year',                 &
        !'NSC','SeedC','leafC', 'rootC', 'SW-C',  'HW-C', 'LAI', &
        !'McrbC', 'fineL', 'struL', 'McrbN', 'fineN', 'struN',   &
        !'mineralN', 'N_uptake'

        ! ---------- annual call -------------
        ! update the LAImax of each PFT according to available N for next year
        call vegn_annualLAImax_update(vegn)

        !---------------------------------------------
        ! Reproduction and mortality: initilising a new cohort based on 
        !---------------------------------------------
        ! seed C and germination probability
        call vegn_reproduction(vegn)

        ! Natural mortality (reducing number of individuals 'nindivs')
        call vegn_nat_mortality(vegn, real(seconds_per_year))

        ! Kill all individuals in a cohort if NSC falls below critical point
        call vegn_starvation(vegn)

        !---------------------------------------------
        ! Re-organize cohorts
        !---------------------------------------------
        call relayer_cohorts(vegn)
        call vegn_mergecohorts(vegn)
        call kill_lowdensity_cohorts(vegn)

        !---------------------------------------------
        ! set annual variables zero
        !---------------------------------------------
        call vegn_annual_diagnostics_zero(vegn)
        iyears = iyears + 1

      end if
    end do gridcellloop

    ! !----------------------------------------------------------------
    ! ! Get rolling multi-year averages (needs to store entire arrays)
    ! !----------------------------------------------------------------
    ! call get_rlm_waterbal( tile(:,:)%soil%phy, interface%steering%init )

    !----------------------------------------------------------------
    ! Write to NetCDF output
    !----------------------------------------------------------------
    ! if (.not.interface%params_siml%is_calib) then
    !   if (verbose) print*,'calling writeout_nc_() ... '
    !   call writeout_nc_forcing()
    !   call writeout_nc_gpp()
    !   call writeout_nc_waterbal()
    !   if (verbose) print*,'... done'
    ! end if

    if (interface%steering%finalize) then
      !----------------------------------------------------------------
      ! Finazlize run: deallocating memory
      !----------------------------------------------------------------
      ! deallocate( tile )
      ! deallocate( tile_fluxes )
      ! deallocate( plant )
      ! deallocate( plant_fluxes )

      ! deallocate(cc)
      close(101)
      close(102)
      close(103)
      close(104)
      deallocate(vegn%cohorts)
      deallocate(forcingData)

    end if

    if (verbose) print*,'Done with biosphere for this year. Guete Rutsch!'

  end function biosphere_annual


  subroutine read_forcingdata(forcingData,datalines,days_data,yr_data,timestep,readyear)

    ! read in forcing data (Users need to write their own data input procedure)
    type(climate_data_type),pointer,intent(inout) :: forcingData(:)
    integer,intent(inout) :: datalines,days_data,yr_data
    integer, intent(in)   :: readyear
    real, intent(inout)   :: timestep

    !------------local var -------------------
    type(climate_data_type), pointer :: climateData(:)
    character(len=50)  filepath_in
    character(len=50)  climfile
    character(len=50)  parafile       ! parameter file
    character(len=80)  commts
    character(len=4)  readyear_char
    integer, parameter :: iiterms=9        ! MDK data for Oak Ridge input
    integer, parameter :: ilines=12*366*24 ! the maxmum records of Oak Ridge FACE, 1999~2007
    integer,dimension(ilines) :: year_data
    real,   dimension(ilines) :: doy_data,hour_data
    real input_data(iiterms,ilines)
    real inputstep
    integer :: istat1,istat2,istat3
    integer :: doy,idays
    integer :: i,j,k
    integer :: m,n

    write(readyear_char,999) readyear

    filepath_in = 'input/'
    climfile    = 'Temperate_forcing_'//trim(readyear_char)//'.txt'
    climfile    = trim(filepath_in)//trim(climfile)

    print*,'Reading forcing from ', trim(climfile)

    ! open forcing data
    open(11,file=climfile,status='old',ACTION='read',IOSTAT=istat2)
    write(*,*)istat2
    ! skip 2 lines of input met data file
    read(11,'(a160)') commts
    ! read(11,'(a160)') commts ! MDK data only has one line comments
    m       = 0  ! to record the lines in a file
    idays   = 1  ! the total days in a data file
    yr_data = 0  ! to record years of a dataset

    ! read forcing files
    do
      m=m+1
      read(11,*,IOSTAT=istat3)year_data(m),doy_data(m),hour_data(m),   &
                              (input_data(n,m),n=1,iiterms)
      if (istat3<0) exit
      if (m == 1) then
          doy = doy_data(m)
      else
          doy = doy_data(m-1)
      endif
      if (doy /= doy_data(m)) idays = idays + 1

    enddo ! end of reading the forcing file

    timestep = hour_data(2) - hour_data(1)
    write(*,*)"forcing",datalines,yr_data,timestep
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

    allocate(climateData(datalines))
    do i=1,datalines
       climateData(i)%year      = year_data(i)                     ! Year
       climateData(i)%doy       = doy_data(i)                      ! day of the year
       climateData(i)%hod       = hour_data(i)                     ! hour of the day
       climateData(i)%PAR       = input_data(1,i)                  ! umol/m2/s
       climateData(i)%radiation = input_data(2,i)                  ! W/m2
       climateData(i)%Tair      = input_data(3,i) + 273.16         ! air temperature, K
       climateData(i)%Tsoil     = input_data(4,i) + 273.16         ! soil temperature, K
       climateData(i)%RH        = input_data(5,i)                  ! relative humidity
       climateData(i)%rain      = input_data(6,i)/(timestep * 3600)! kgH2O m-2 s-1
       climateData(i)%windU     = input_data(7,i)                  ! wind velocity (m s-1)
       climateData(i)%pressure  = input_data(8,i)                  ! pa
       climateData(i)%soilwater = 0.8                              ! soil moisture, vol/vol
    enddo
    forcingData => climateData
    write(*,*)"forcing", datalines,days_data,yr_data

    ! xxx debug
    ! print*,'forcingData(1)%year      ', forcingData(1)%year      
    ! print*,'forcingData(1)%doy       ', forcingData(1)%doy       
    ! print*,'forcingData(1)%hod       ', forcingData(1)%hod       
    ! print*,'forcingData(1)%PAR       ', forcingData(1)%PAR       
    ! print*,'forcingData(1)%radiation ', forcingData(1)%radiation 
    ! print*,'forcingData(1)%Tair      ', forcingData(1)%Tair      
    ! print*,'forcingData(:)%Tsoil     ', forcingData(:)%Tsoil     
    ! print*,'forcingData(1)%RH        ', forcingData(1)%RH        
    ! print*,'forcingData(1)%rain      ', forcingData(1)%rain      
    ! print*,'forcingData(1)%windU     ', forcingData(1)%windU     
    ! print*,'forcingData(1)%pressure  ', forcingData(1)%pressure  
    ! print*,'forcingData(1)%soilwater ', forcingData(1)%soilwater 
    ! stop

    999  format (I4.4)

  end subroutine read_forcingdata


end module md_biosphere
