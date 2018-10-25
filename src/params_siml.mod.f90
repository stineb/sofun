module md_params_siml
  !////////////////////////////////////////////////////////////////
  !  Module contains simulation parameters read in by getpar_siml
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_sofunutils, only: getparint, getparstring, getparlogical

  implicit none

  private
  public paramstype_siml, getsteering, outtype_steering, getpar_siml

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

    ! booleans defining whether variable is written to ascii output
    logical :: loutdgpp       
    logical :: loutdrd
    logical :: loutdtransp    
    logical :: loutdnpp       
    logical :: loutdnup       
    logical :: loutdcex       
    logical :: loutdcleaf     
    logical :: loutdcroot     
    logical :: loutdclabl     
    logical :: loutdnlabl     
    logical :: loutdclitt     
    logical :: loutdnlitt     
    logical :: loutdlai       
    logical :: loutdfapar
    logical :: loutdtemp_soil 
    logical :: loutdtemp

    ! booleans defining whether module-specific output variables are to be written to output
    logical :: loutplant
    logical :: loutalloc
    logical :: loutgpp
    logical :: loutnpp
    logical :: loutntransform
    logical :: loutwaterbal  
    logical :: loutlittersom
    logical :: loutnuptake
    logical :: loutlanduse
    logical :: loutforcing
    logical :: loutturnover

    ! booleans defining whether variable is written to NetCDF output
    logical :: lncoutdtemp
    logical :: lncoutdgpp
    logical :: lncoutdfapar
    logical :: lncoutdwaterbal

    ! booleans defining whether variable is used as calibration target
    logical :: lcalibgpp
    logical :: lcalibfapar
    logical :: lcalibtransp

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
    !  SR defines variables used for steering simulation for each 
    !  simulation year (setting booleans for opening files, doing   
    !  spinup etc.)
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
    ! Returns cyce year for climate recycling, given number of spinup
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
    !  SR for reading and defining simulation parameters from file 
    !  <runname>.sofun.parameter. Only once at start of simulation.
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

    if (out_getpar_siml%do_spinup) then
      out_getpar_siml%runyears = out_getpar_siml%nyeartrend + out_getpar_siml%spinupyears
    else
      out_getpar_siml%runyears = out_getpar_siml%nyeartrend
      out_getpar_siml%spinupyears = 0
    endif

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
    call getparstring( 'input/dfapar_source.txt', 'fapar_forcing_source', out_getpar_siml%fapar_forcing_source )

    out_getpar_siml%in_netrad         = getparlogical( 'run/'//runname//'.sofun.parameter', 'in_netrad' )
    out_getpar_siml%in_ppfd           = getparlogical( 'run/'//runname//'.sofun.parameter', 'in_ppfd' )
    
    out_getpar_siml%do_spinup         = getparlogical( 'run/'//runname//'.sofun.parameter', 'spinup' )
    
    out_getpar_siml%const_co2_year    = getparint( 'run/'//runname//'.sofun.parameter', 'const_co2_year' )
    out_getpar_siml%const_ndep_year   = getparint( 'run/'//runname//'.sofun.parameter', 'const_ndep_year' )
    out_getpar_siml%const_nfert_year  = getparint( 'run/'//runname//'.sofun.parameter', 'const_nfert_year' )
    out_getpar_siml%const_clim_year   = getparint( 'run/'//runname//'.sofun.parameter', 'const_clim_year' )
    out_getpar_siml%const_lu_year     = getparint( 'run/'//runname//'.sofun.parameter', 'const_lu_year' )

    out_getpar_siml%soilmstress       = getparlogical( 'run/'//runname//'.sofun.parameter', 'soilmstress' )
    out_getpar_siml%tempstress        = getparlogical( 'run/'//runname//'.sofun.parameter', 'tempstress' )
    
    ! boolean for ascii output writing (core variables)
    out_getpar_siml%loutdgpp       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdgpp' )
    out_getpar_siml%loutdrd        = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdrd' )
    out_getpar_siml%loutdtransp    = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdtransp' )
    out_getpar_siml%loutdnpp       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdnpp' )
    out_getpar_siml%loutdnup       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdnup' )
    out_getpar_siml%loutdcex       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdcex' )
    out_getpar_siml%loutdcleaf     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdCleaf' )
    out_getpar_siml%loutdcroot     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdCroot' )
    out_getpar_siml%loutdclabl     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdClabl' )
    out_getpar_siml%loutdnlabl     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdNlabl' )
    out_getpar_siml%loutdclitt     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdClitt' )
    out_getpar_siml%loutdnlitt     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdNlitt' )
    out_getpar_siml%loutdlai       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdlai' )
    out_getpar_siml%loutdfapar     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdfapar' )
    out_getpar_siml%loutdtemp_soil = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdtemp_soil' )
    out_getpar_siml%loutdtemp      = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdtemp' )

    ! boolean for ascii output writing (module variables)
    out_getpar_siml%loutplant      = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutplant' )
    out_getpar_siml%loutalloc      = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutalloc' )
    out_getpar_siml%loutgpp        = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutgpp' )
    out_getpar_siml%loutnpp        = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutnpp' )
    out_getpar_siml%loutntransform = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutntransform') 
    out_getpar_siml%loutwaterbal   = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutwaterbal')
    out_getpar_siml%loutlittersom  = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutlittersom' )
    out_getpar_siml%loutnuptake    = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutnuptake' )
    out_getpar_siml%loutlanduse    = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutlanduse' )
    out_getpar_siml%loutforcing    = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutforcing' )
    out_getpar_siml%loutturnover   = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutturnover' )

    ! boolean for NetCDF output writing
    out_getpar_siml%lncoutdtemp    = getparlogical( 'run/'//runname//'.sofun.parameter', 'lncoutdtemp' )
    out_getpar_siml%lncoutdfapar   = getparlogical( 'run/'//runname//'.sofun.parameter', 'lncoutdfapar' )
    out_getpar_siml%lncoutdgpp     = getparlogical( 'run/'//runname//'.sofun.parameter', 'lncoutdgpp' )
    out_getpar_siml%lncoutdwaterbal= getparlogical( 'run/'//runname//'.sofun.parameter', 'lncoutdwaterbal' )

    ! If NetCDF output writing is true and ascii output writing is false, then overwrite
    if (.not.out_getpar_siml%loutdtemp    .and. out_getpar_siml%lncoutdtemp    )  out_getpar_siml%loutdtemp    = .true.
    if (.not.out_getpar_siml%loutdfapar   .and. out_getpar_siml%lncoutdfapar   )  out_getpar_siml%loutdfapar   = .true.
    if (.not.out_getpar_siml%loutdgpp     .and. out_getpar_siml%lncoutdgpp     )  out_getpar_siml%loutdgpp     = .true.
    if (.not.out_getpar_siml%loutwaterbal .and. out_getpar_siml%lncoutdwaterbal ) out_getpar_siml%loutwaterbal = .true.

    ! boolean to define which variables are used as calibration target
    out_getpar_siml%lcalibgpp    = getparlogical( 'run/'//runname//'.sofun.parameter', 'lcalibgpp' )
    out_getpar_siml%lcalibfapar  = getparlogical( 'run/'//runname//'.sofun.parameter', 'lcalibfapar' )
    out_getpar_siml%lcalibtransp = getparlogical( 'run/'//runname//'.sofun.parameter', 'lcalibtransp' )

    print*, "... done"

  end function getpar_siml

end module md_params_siml


