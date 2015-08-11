module _params_siml
  !////////////////////////////////////////////////////////////////
  !  Module contains simulation parameters read in by getpar_siml
  ! Copyright (C) 2015  Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: recycle

  implicit none
  integer :: runyears                    ! number of years of entire simulation (spinup+transient)
  integer :: spinupyears                 ! number of spinup years
  integer :: nyeartrend                  ! number of transient years
  integer :: firstyeartrend              ! year AD of first transient year
  logical :: do_spinup                   ! whether this simulation does spinup 
  logical :: init                        ! is true in first simulation year
  logical :: spinup                      ! is true during spinup
  integer :: transyear                   ! year in transient simulation
  character(len=256) :: runname
  character(len=6)   :: sitename
  character(len=256) :: input_dir
  character(len=256) :: co2_forcing_file
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
    use _sofunutils

    character(len=*), intent(in) :: runname

    ! Read in main model parameters
    write(0,*) 'reading parameter file ', runname//".sofun.parameter ..."

    call getstring(runname//".sofun.parameter",'sitename',sitename)
    call getlogical(runname//".sofun.parameter",'spinup',do_spinup)
    call getint(runname//".sofun.parameter",'spinupyears',spinupyears)
    call getint(runname//".sofun.parameter",'firstyeartrend',firstyeartrend)
    call getint(runname//".sofun.parameter",'nyeartrend',nyeartrend)

    call getstring(runname//".sofun.parameter",'input_dir',input_dir)
    call getstring(runname//".sofun.parameter",'co2_forcing_file',co2_forcing_file)
    call getstring(runname//".sofun.parameter",'pftparfil',pftparfil)

    if (do_spinup) then
      runyears = nyeartrend + spinupyears
    else
      runyears = nyeartrend
      spinupyears = 0
    endif

  end subroutine getpar_siml

end module _params_siml

