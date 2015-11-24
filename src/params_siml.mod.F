module _params_siml
  !////////////////////////////////////////////////////////////////
  !  Module contains simulation parameters read in by getpar_siml
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: recycle
  use _sofunutils, only: getparint, getparstring, getparlogical

  implicit none
  integer :: runyears                    ! number of years of entire simulation (spinup+transient)
  integer :: spinupyears                 ! number of spinup years
  integer :: nyeartrend                  ! number of transient years
  integer :: firstyeartrend              ! year AD of first transient year
  logical :: do_spinup                   ! whether this simulation does spinup 
  logical :: init                        ! is true in first simulation year
  logical :: spinup                      ! is true during spinup
  logical :: const_co2                   ! is true when using constant CO2, given by first transient year in 'co2_forcing_file'
  logical :: const_ndep                  ! is true when using constant N deposition, given by first transient year in 'ndep_forcing_file'
  integer :: transyear                   ! year in transient simulation
  character(len=256) :: runname
  character(len=6)   :: sitename
  character(len=256) :: input_dir
  character(len=256) :: co2_forcing_file
  character(len=256) :: ndep_forcing_file
  character(len=256) :: pftparfil

contains

  subroutine getsteering( year, realyear)
    !////////////////////////////////////////////////////////////////
    !  SR defines variables used for steering simulation for each 
    !  simulation year (setting booleans for opening files, doing   
    !  spinup etc.)
    !----------------------------------------------------------------
    integer, intent(in) :: year
    integer, intent(out) :: realyear

    if (do_spinup) then
      realyear = max( firstyeartrend, year - spinupyears + firstyeartrend - 1 )
      if (year<=spinupyears) then
        spinup = .true.
      else  
        spinup = .false.
      endif
    else
      realyear = year + firstyeartrend
    endif

    if (year==1) then
      init = .true.
    else
      init = .false.
    endif 

    if (spinup) then
      transyear = mod(year,recycle)
      if (transyear==0) transyear = recycle
      ! use this if data starts from a given year that is not equal to the required start 
      ! year of the transient simulation. Here '1901' is start year in data set.
      !first_cycleyear = mod(recycle-mod(spinupyears+1901-firstyeartrend,recycle)+1,recycle)
      !cycleyear = modulo(first_cycleyear+year-1,recycle)
      !if (cycleyear.eq.0) cycleyear=31
    else
      transyear = realyear-firstyeartrend+1
    endif

  end subroutine getsteering


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

