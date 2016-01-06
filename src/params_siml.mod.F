module _params_siml
  !////////////////////////////////////////////////////////////////
  !  Module contains simulation parameters read in by getpar_siml
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: recycle
  use _sofunutils, only: getparint, getparstring, getparlogical

  implicit none

  integer :: runyears        ! number of years of entire simulation (spinup+transient)
  integer :: spinupyears     ! number of spinup years
  integer :: nyeartrend      ! number of transient years
  integer :: firstyeartrend  ! year AD of first transient year
  logical :: do_spinup       ! whether this simulation does spinup 
  logical :: const_co2       ! is true when using constant CO2, given by first transient year in 'co2_forcing_file'
  logical :: const_ndep      ! is true when using constant N deposition, given by first transient year in 'ndep_forcing_file'
  character(len=256) :: runname
  character(len=6)   :: sitename
  character(len=256) :: input_dir
  character(len=256) :: co2_forcing_file
  character(len=256) :: ndep_forcing_file
  character(len=256) :: pftparfil

  logical :: spinup          ! is true during spinup
  logical :: init            ! is true in first simulation year

  type outtype_steering
    integer :: forcingyear     ! year AD for which forcings are read in (=firstyeartrend during spinup)
    integer :: climateyear     ! year AD for which climate is read in (recycling during spinup or when climate is held const.)
  end type

contains

  function getsteering( year ) result( out_steering )
    !////////////////////////////////////////////////////////////////
    !  SR defines variables used for steering simulation for each 
    !  simulation year (setting booleans for opening files, doing   
    !  spinup etc.)
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: year

    ! local variables
    integer :: first_cycleyear, cycleyear

    ! function return variable
    type( outtype_steering ) :: out_steering

    ! xxx try
    integer, parameter :: recycle_test = 5

    integer :: remainder, nfits

    if (do_spinup) then

      if (year<=spinupyears) then

        spinup = .true.
        out_steering%forcingyear = firstyeartrend

        remainder = mod( spinupyears, recycle_test )
        nfits = (spinupyears - remainder) / recycle_test
        first_cycleyear = recycle_test - remainder + 1
        cycleyear = modulo( first_cycleyear + year - 1, recycle_test )  
        if (cycleyear==0) cycleyear = recycle_test
        out_steering%climateyear = cycleyear + firstyeartrend - 1

      else  

        spinup = .false.
        out_steering%forcingyear =  year - spinupyears + firstyeartrend - 1 
        out_steering%climateyear = out_steering%forcingyear

      endif

    else

      out_steering%forcingyear = year + firstyeartrend - 1 
      out_steering%climateyear = out_steering%forcingyear

    endif

    ! write(0,*) 'recycle_test, spinupyears, firstyeartrend', recycle_test, spinupyears, firstyeartrend
    ! write(0,*) 'first_cycleyear, cycleyear', first_cycleyear, cycleyear
    ! write(0,*) 'year, forcingyear, climateyear', year, out_steering%forcingyear, out_steering%climateyear
    ! stop

    if (year==1) then
      init = .true.
    else
      init = .false.
    endif 

  end function getsteering


  subroutine getpar_siml( runname )
    !////////////////////////////////////////////////////////////////
    !  SR for reading and defining simulation parameters from file 
    !  <runname>.sofun.parameter. Only once at start of simulation.
    !----------------------------------------------------------------
    character(len=*), intent(in) :: runname

    ! Read in main model parameters
    write(0,*) 'reading parameter file ', runname//".sofun.parameter ..."

    ! sitename         = getparstring( runname//".sofun.parameter", 'sitename' )
    ! input_dir        = getparstring( runname//".sofun.parameter", 'input_dir' )
    ! co2_forcing_file = getparstring( runname//".sofun.parameter", 'co2_forcing_file' )
    ! pftparfil        = getparstring( runname//".sofun.parameter", 'pftparfil' )

    call getparstring( runname//".sofun.parameter", 'sitename', sitename )
    call getparstring( runname//".sofun.parameter", 'input_dir', input_dir )
    call getparstring( runname//".sofun.parameter", 'co2_forcing_file', co2_forcing_file )
    call getparstring( runname//".sofun.parameter", 'ndep_forcing_file', ndep_forcing_file )
    call getparstring( runname//".sofun.parameter", 'pftparfil', pftparfil )

    do_spinup        = getparlogical( runname//".sofun.parameter", 'spinup' )
    const_co2        = getparlogical( runname//".sofun.parameter", 'const_co2' )
    const_ndep       = getparlogical( runname//".sofun.parameter", 'const_ndep' )
    
    spinupyears      = getparint( runname//".sofun.parameter", 'spinupyears' )
    firstyeartrend   = getparint( runname//".sofun.parameter", 'firstyeartrend' )
    nyeartrend       = getparint( runname//".sofun.parameter", 'nyeartrend' )
    
    if (do_spinup) then
      runyears = nyeartrend + spinupyears
    else
      runyears = nyeartrend
      spinupyears = 0
    endif

  end subroutine getpar_siml

end module _params_siml

