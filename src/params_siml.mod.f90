module md_params_siml
  !////////////////////////////////////////////////////////////////
  ! Module for handling simulation parameters.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_sofunutils, only: getparint, getparstring, getparlogical

  implicit none

  private
  public paramstype_siml, getsteering, outtype_steering, getpar_siml

  !----------------------------------------------------------------
  ! Derived type for simulation parameters
  !----------------------------------------------------------------
  type paramstype_siml
    
    integer :: runyears        ! number of years of entire simulation (spinup+transient)
    integer :: spinupyears     ! number of spinup years
    integer :: nyeartrend      ! number of transient years
    integer :: firstyeartrend  ! year AD of first transient year
    integer :: recycle         ! length of standard recycling period
    integer :: daily_out_startyr! first year where daily output is written
    integer :: daily_out_endyr ! last year where daily output is written
    integer :: outdt           ! output periodicity
    integer :: outnt           ! number of output time steps per year
    integer :: secs_per_tstep  ! number of seconds per time step (now daily => 60 * 60 * 24)
    
    logical :: do_spinup            ! whether this simulation does spinup 
    logical :: is_calib             ! whether this simulation is a calibration simulation (overriding parameters and no output)

    integer :: const_co2_year       ! is true when using constant CO2, given by first transient year in 'co2_forcing_file'
    integer :: const_ndep_year      ! is true when using constant N deposition, given by first transient year in 'ndep_forcing_file'
    integer :: const_nfert_year     ! is true when using constant N fertilisation, given by first transient year in 'nfert_forcing_file'
    integer :: const_clim_year      ! is true when using constant climate, given by year 'firstyeartrend'
    integer :: const_lu_year        ! is true when using constant land use, given by year 'firstyeartrend'

    logical :: soilmstress          ! when true, an empirical soil moisture stress function is applied to GPP
    logical :: tempstress           ! when true, an empirical temperature stress function is applied to GPP
    
    character(len=256) :: runname
    character(len=256) :: sitename
    character(len=256) :: co2_forcing_file
    character(len=256) :: ndep_noy_forcing_file
    character(len=256) :: ndep_nhx_forcing_file
    character(len=256) :: nfert_noy_forcing_file
    character(len=256) :: nfert_nhx_forcing_file 
    character(len=256) :: do_grharvest_forcing_file
    character(len=256) :: fapar_forcing_source

    ! optionally prescribed variables (if false, then simulated internally)
    logical :: in_netrad    ! net radiation
    logical :: in_ppfd      ! photosynthetic photon flux density 

    ! activated PFTs
    logical :: lTrE        ! evergreen tree
    logical :: lTNE        ! evergreen tree, N-fixing
    logical :: lTrD        ! deciduous tree
    logical :: lTND        ! deciduous tree, N-fixing
    logical :: lGr3        ! grass, C3 photosynthetic pathway
    logical :: lGN3        ! grass, C3 photosynthetic pathway, N-fixing
    logical :: lGr4        ! grass, C4 photosynthetic pathway

    integer :: npft        ! number of activated PFTs

    ! Module-specific booleans defining whether a set of variables is written to annual output
    logical :: loutplant
    logical :: loutgpp
    logical :: loutwaterbal
    logical :: loutforcing

    ! Module-specific booleans whether a single variable is written to daily output
    logical :: loutdgpp
    logical :: loutdrd
    logical :: loutdtransp
    logical :: loutdwcont
    logical :: loutdaet
    logical :: loutdpet
    logical :: loutdnetrad
    logical :: loutdwbal
    logical :: loutdtemp
    logical :: loutdfapar
    logical :: loutdtemp_soil

    ! booleans defining whether variable is used as calibration target
    logical :: lcalibgpp
    logical :: lcalibfapar
    logical :: lcalibtransp
    logical :: lcaliblatenth

  end type paramstype_siml

  type outtype_steering
    integer :: year
    integer :: forcingyear     ! year AD for which forcings are read in (=firstyeartrend during spinup)
    integer :: climateyear     ! year AD for which climate is read in (recycling during spinup or when climate is held const.)
    integer :: outyear         ! year AD written to output
    logical :: spinup          ! is true during spinup
    logical :: init            ! is true in first simulation year
    logical :: finalize        ! is true in the last simulation year
    logical :: do_soilequil    ! true in year of analytical soil equilibration (during spinup)
    logical :: average_soil    ! true in years before analytical soil equilibration, when average in and out are taken
    logical :: project_nmin    ! true in all years before analytical soil equilibration, when projected soil N mineralisation is used
    logical :: dofree_alloc    ! true if allocation is not fixed by 'frac_leaf'
    logical :: add_ninorg      ! true in the first few years to get it started
  end type

contains

  function getsteering( year, params_siml ) result( out_steering )
    !////////////////////////////////////////////////////////////////
    ! Gets variables used for steering simulation for each 
    ! simulation year (setting booleans for opening files, doing   
    ! spinup etc.)
    !----------------------------------------------------------------
    use md_params_core, only: dummy

    ! arguments
    integer, intent(in) :: year
    type( paramstype_siml ), intent(in) :: params_siml

    ! function return variable
    type( outtype_steering ) :: out_steering

    ! local variables
    integer :: cycleyear

    integer, parameter :: spinupyr_soilequil_1 = 600   ! year of analytical soil equilibration, based on mean litter -> soil input flux
    integer, parameter :: spinupyr_soilequil_2 = 1200  ! year of analytical soil equilibration, based on mean litter -> soil input flux
    integer, parameter :: spinup_add_ninorg    = 100   ! year until which inorganic N is added to get it started

    out_steering%year = year

    if (params_siml%do_spinup) then

      if (year<=spinup_add_ninorg) then
        out_steering%add_ninorg = .true.
      else
        out_steering%add_ninorg = .false.
      end if

      if (year<=params_siml%spinupyears) then
        ! during spinup
        out_steering%spinup = .true.
        out_steering%forcingyear = params_siml%firstyeartrend
        cycleyear = get_cycleyear( year, params_siml%spinupyears, params_siml%recycle )
        out_steering%climateyear = cycleyear + params_siml%firstyeartrend - 1

      else  
        ! during transient simulation
        out_steering%spinup = .false.
        out_steering%forcingyear =  year - params_siml%spinupyears + params_siml%firstyeartrend - 1

        if (params_siml%const_clim_year/=int(dummy)) then
          ! constant climate year specified
          cycleyear = get_cycleyear( year, params_siml%spinupyears, params_siml%recycle )
          out_steering%climateyear = cycleyear + params_siml%const_clim_year - 1
        
        else
          ! constant climate year not specified
          out_steering%climateyear = out_steering%forcingyear
        
        end if

      endif
      out_steering%outyear = year + params_siml%firstyeartrend - params_siml%spinupyears - 1

      if ( year > 3 ) then
      ! if (year > (spinupyr_soilequil_1 + 1) ) then
      ! if (out_steering%forcingyear > 2003 ) then
        out_steering%dofree_alloc = .true.
      else
        out_steering%dofree_alloc = .false.
      end if

      if ( (year==spinupyr_soilequil_1 .or. year==spinupyr_soilequil_2 ) .and. year<=params_siml%spinupyears) then
        out_steering%do_soilequil = .true.
      else
        out_steering%do_soilequil = .false.
      end if

      if ( year<=params_siml%spinupyears .and. ( year > ( spinupyr_soilequil_1 - params_siml%recycle ) .and. year <= spinupyr_soilequil_1 &
              .or. year > ( spinupyr_soilequil_2 - params_siml%recycle ) .and. year <= spinupyr_soilequil_2 ) ) then
        out_steering%average_soil = .true.
      else
        out_steering%average_soil = .false.
      end if

      if ( year<=params_siml%spinupyears .and. year <= spinupyr_soilequil_1 ) then
        out_steering%project_nmin = .true.
      else
        out_steering%project_nmin = .false.
      end if

    else

      out_steering%dofree_alloc = .false.
      out_steering%do_soilequil = .false.
      out_steering%average_soil = .false.
      out_steering%project_nmin = .false.
      out_steering%forcingyear = year + params_siml%firstyeartrend - 1 
      out_steering%climateyear = out_steering%forcingyear
      out_steering%outyear = year + params_siml%firstyeartrend - 1

    endif

    if (year==1) then
      out_steering%init = .true.
    else
      out_steering%init = .false.
    endif 

    if (year==params_siml%runyears) then
      out_steering%finalize = .true.
    else
      out_steering%finalize = .false.
    end if

    ! print*, 'out_steering%climateyear'
    ! print*, out_steering%climateyear
    ! if (year>30) stop

  end function getsteering


  function get_cycleyear( year, spinupyears, recycle ) result( cycleyear )
    !////////////////////////////////////////////////////////////////
    ! Returns cycle year for climate recycling, given number of spinup
    ! years and recycle period length, so that in the last year of
    ! of the spinup, 'cycleyear' is equal to 'recycle'.
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: year
    integer, intent(in) :: spinupyears
    integer, intent(in) :: recycle

    ! local variables
    integer :: remainder
    integer :: nfits
    integer :: first_cycleyear

    ! function return variable
    integer :: cycleyear

    remainder = mod( spinupyears, recycle )
    nfits = (spinupyears - remainder) / recycle
    first_cycleyear = recycle - remainder + 1
    cycleyear = modulo( first_cycleyear + year - 1, recycle )  
    if (cycleyear==0) cycleyear = recycle

  end function get_cycleyear


  function getpar_siml( runname ) result( out_getpar_siml )
    !////////////////////////////////////////////////////////////////
    ! Reads simulation parameters from file 
    ! <runname>.sofun.parameter. Only once at start of simulation.
    !----------------------------------------------------------------
    use md_params_core, only: ndayyear, npft

    ! argument
    character(len=*), intent(in) :: runname

    ! function return variable
    type( paramstype_siml ) :: out_getpar_siml

    ! local variables
    integer :: npft_local, pft

    ! Read in main model parameters
    print*, 'reading parameter file ', runname//".sofun.parameter ..."

    out_getpar_siml%runname = runname

    !------------------------------------------------------------------
    ! read simulation parameters that may change between runs within an ensemble (simsuite)
    !------------------------------------------------------------------
    call getparstring( 'run/'//runname//'.sofun.parameter', 'sitename', out_getpar_siml%sitename )

    out_getpar_siml%firstyeartrend    = getparint( 'run/'//runname//'.sofun.parameter', 'firstyeartrend' )
    out_getpar_siml%nyeartrend        = getparint( 'run/'//runname//'.sofun.parameter', 'nyeartrend' )
    out_getpar_siml%daily_out_startyr = getparint( 'run/'//runname//'.sofun.parameter', 'daily_out_startyr' )
    out_getpar_siml%daily_out_endyr   = getparint( 'run/'//runname//'.sofun.parameter', 'daily_out_endyr' )

    out_getpar_siml%spinupyears       = getparint( 'run/'//runname//'.sofun.parameter', 'spinupyears' )
    out_getpar_siml%recycle           = getparint( 'run/'//runname//'.sofun.parameter', 'recycle' )
    out_getpar_siml%outdt             = getparint( 'run/'//runname//'.sofun.parameter', 'outdt' )
    out_getpar_siml%outnt = ceiling( real( ndayyear ) / real( out_getpar_siml%outdt ) )

    out_getpar_siml%do_spinup         = getparlogical( 'run/'//runname//'.sofun.parameter', 'spinup' )

    if (out_getpar_siml%do_spinup) then
      out_getpar_siml%runyears = out_getpar_siml%nyeartrend + out_getpar_siml%spinupyears
    else
      out_getpar_siml%runyears = out_getpar_siml%nyeartrend
      out_getpar_siml%spinupyears = 0
    endif

    ! time step in number of seconds
    out_getpar_siml%secs_per_tstep = getparint( 'run/'//runname//'.sofun.parameter', 'secs_per_tstep' )

    ! activated PFTs
    out_getpar_siml%lTrE = getparlogical( 'run/'//runname//'.sofun.parameter', 'lTrE' )
    out_getpar_siml%lTNE = getparlogical( 'run/'//runname//'.sofun.parameter', 'lTNE' )
    out_getpar_siml%lTrD = getparlogical( 'run/'//runname//'.sofun.parameter', 'lTrD' )
    out_getpar_siml%lTND = getparlogical( 'run/'//runname//'.sofun.parameter', 'lTND' )
    out_getpar_siml%lGr3 = getparlogical( 'run/'//runname//'.sofun.parameter', 'lGr3' )
    out_getpar_siml%lGN3 = getparlogical( 'run/'//runname//'.sofun.parameter', 'lGN3' )
    out_getpar_siml%lGr4 = getparlogical( 'run/'//runname//'.sofun.parameter', 'lGr4' )

    npft_local = 0
    if (out_getpar_siml%lTrE) npft_local = npft_local + 1
    if (out_getpar_siml%lTNE) npft_local = npft_local + 1
    if (out_getpar_siml%lTrD) npft_local = npft_local + 1
    if (out_getpar_siml%lTND) npft_local = npft_local + 1
    if (out_getpar_siml%lGr3) npft_local = npft_local + 1
    if (out_getpar_siml%lGr4) npft_local = npft_local + 1
    if (out_getpar_siml%lGN3) npft_local = npft_local + 1
    ! temporary solution to this
    print*,'found ', npft_local, ' activated PFTs.'
    if (npft/=npft_local) stop 'GETPAR_SIML: adjust number of activated PFTs by hand in params_core.'

    !------------------------------------------------------------------
    ! read simulation parameters that do not change change between runs within an ensemble (simsuite)
    !------------------------------------------------------------------
    call getparstring( 'run/'//runname//'.sofun.parameter', 'co2_forcing_file', out_getpar_siml%co2_forcing_file )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'ndep_noy_forcing_file', out_getpar_siml%ndep_noy_forcing_file )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'ndep_nhx_forcing_file', out_getpar_siml%ndep_nhx_forcing_file )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'nfert_noy_forcing_file', out_getpar_siml%nfert_noy_forcing_file )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'nfert_nhx_forcing_file', out_getpar_siml%nfert_nhx_forcing_file )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'do_grharvest_forcing_file', out_getpar_siml%do_grharvest_forcing_file )

    ! Reading from a different file!
    call getparstring( 'run/'//runname//'.sofun.parameter', 'fapar_forcing_source', out_getpar_siml%fapar_forcing_source )

    out_getpar_siml%in_netrad         = getparlogical( 'run/'//runname//'.sofun.parameter', 'in_netrad' )
    out_getpar_siml%in_ppfd           = getparlogical( 'run/'//runname//'.sofun.parameter', 'in_ppfd' )
        
    out_getpar_siml%const_co2_year    = getparint( 'run/'//runname//'.sofun.parameter', 'const_co2_year' )
    out_getpar_siml%const_ndep_year   = getparint( 'run/'//runname//'.sofun.parameter', 'const_ndep_year' )
    out_getpar_siml%const_nfert_year  = getparint( 'run/'//runname//'.sofun.parameter', 'const_nfert_year' )
    out_getpar_siml%const_clim_year   = getparint( 'run/'//runname//'.sofun.parameter', 'const_clim_year' )
    out_getpar_siml%const_lu_year     = getparint( 'run/'//runname//'.sofun.parameter', 'const_lu_year' )

    out_getpar_siml%soilmstress       = getparlogical( 'run/'//runname//'.sofun.parameter', 'soilmstress' )
    out_getpar_siml%tempstress        = getparlogical( 'run/'//runname//'.sofun.parameter', 'tempstress' )
    
    ! Module-specific booleans defining whether a set of variables is written to annual output
    out_getpar_siml%loutplant     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutplant' )
    out_getpar_siml%loutgpp       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutgpp' )
    out_getpar_siml%loutwaterbal  = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutwaterbal' )

    ! Module-specific booleans whether a single variable is written to daily output
    out_getpar_siml%loutdgpp       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdgpp' )
    out_getpar_siml%loutdrd        = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdrd' )
    out_getpar_siml%loutdtransp    = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdtransp' )
    out_getpar_siml%loutdwcont     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdwcont' )
    out_getpar_siml%loutdaet       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdaet' )
    out_getpar_siml%loutdpet       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdpet' )
    out_getpar_siml%loutdnetrad    = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdnetrad' )
    out_getpar_siml%loutdwbal      = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdwbal' )
    out_getpar_siml%loutdtemp      = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdtemp' )
    out_getpar_siml%loutdfapar     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdfapar' )
    out_getpar_siml%loutdtemp_soil = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdtemp_soil' )

    ! If NetCDF output writing is true and ascii output writing is false, then overwrite
    if (.not.out_getpar_siml%loutgpp      .and. out_getpar_siml%loutdgpp   ) out_getpar_siml%loutgpp      = .true.
    if (.not.out_getpar_siml%loutgpp      .and. out_getpar_siml%loutdrd    ) out_getpar_siml%loutgpp      = .true.
    if (.not.out_getpar_siml%loutgpp      .and. out_getpar_siml%loutdtransp) out_getpar_siml%loutgpp      = .true.
    if (.not.out_getpar_siml%loutwaterbal .and. out_getpar_siml%loutdpet   ) out_getpar_siml%loutwaterbal = .true.
    if (.not.out_getpar_siml%loutwaterbal .and. out_getpar_siml%loutdnetrad) out_getpar_siml%loutwaterbal = .true.
    if (.not.out_getpar_siml%loutwaterbal .and. out_getpar_siml%loutdwcont ) out_getpar_siml%loutwaterbal = .true.
    if (.not.out_getpar_siml%loutwaterbal .and. out_getpar_siml%loutdaet   ) out_getpar_siml%loutwaterbal = .true.
    if (.not.out_getpar_siml%loutwaterbal .and. out_getpar_siml%loutdwbal  ) out_getpar_siml%loutwaterbal = .true.
    if (.not.out_getpar_siml%loutforcing  .and. out_getpar_siml%loutdtemp  ) out_getpar_siml%loutforcing  = .true.
    if (.not.out_getpar_siml%loutforcing  .and. out_getpar_siml%loutdfapar ) out_getpar_siml%loutforcing  = .true.

    ! boolean to define which variables are used as calibration target
    out_getpar_siml%lcalibgpp     = getparlogical( 'run/'//runname//'.sofun.parameter', 'lcalibgpp' )
    out_getpar_siml%lcalibfapar   = getparlogical( 'run/'//runname//'.sofun.parameter', 'lcalibfapar' )
    out_getpar_siml%lcalibtransp  = getparlogical( 'run/'//runname//'.sofun.parameter', 'lcalibtransp' )
    out_getpar_siml%lcaliblatenth = getparlogical( 'run/'//runname//'.sofun.parameter', 'lcaliblatenth' )

    print*, "... done"

  end function getpar_siml

end module md_params_siml


