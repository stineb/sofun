module _params_siml
  !////////////////////////////////////////////////////////////////
  !  Module contains simulation parameters read in by getpar_siml
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _sofunutils, only: getparint, getparstring, getparlogical

  implicit none

  integer :: runyears        ! number of years of entire simulation (spinup+transient)
  integer :: spinupyears     ! number of spinup years
  integer :: nyeartrend      ! number of transient years
  integer :: firstyeartrend  ! year AD of first transient year
  integer :: recycle         ! length of standard recycling period
  integer :: outyear         ! year AD written to output
  integer :: daily_out_startyr! first year where daily output is written
  integer :: daily_out_endyr ! last year where daily output is written
  logical :: do_spinup       ! whether this simulation does spinup 
  logical :: const_co2       ! is true when using constant CO2, given by first transient year in 'co2_forcing_file'
  logical :: const_ndep      ! is true when using constant N deposition, given by first transient year in 'ndep_forcing_file'
  character(len=256) :: runname
  character(len=256) :: sitename
  character(len=256) :: input_dir
  character(len=256) :: co2_forcing_file
  character(len=256) :: ndep_noy_forcing_file
  character(len=256) :: ndep_nhx_forcing_file
  character(len=256) :: pftparfil

  logical :: spinup          ! is true during spinup
  logical :: init            ! is true in first simulation year
  logical :: prescr_monthly_fapar

  ! booleans defining whether variable is written to output
  logical :: loutdgpp       
  logical :: loutdtransp    
  logical :: loutdnpp       
  logical :: loutdnup       
  logical :: loutdCleaf     
  logical :: loutdCroot     
  logical :: loutdClabl     
  logical :: loutdNlabl     
  logical :: loutdClitt     
  logical :: loutdNlitt     
  logical :: loutdCsoil     
  logical :: loutdNsoil     
  logical :: loutdlai       
  logical :: loutdninorg    
  logical :: loutdtemp_soil 

  logical :: loutntransform
  logical :: loutwaterbal  


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
    integer :: remainder, nfits

    ! function return variable
    type( outtype_steering ) :: out_steering
    

    if (do_spinup) then

      if (year<=spinupyears) then

        spinup = .true.
        out_steering%forcingyear = firstyeartrend

        remainder = mod( spinupyears, recycle )
        nfits = (spinupyears - remainder) / recycle
        first_cycleyear = recycle - remainder + 1
        cycleyear = modulo( first_cycleyear + year - 1, recycle )  
        if (cycleyear==0) cycleyear = recycle
        out_steering%climateyear = cycleyear + firstyeartrend - 1

      else  

        spinup = .false.
        out_steering%forcingyear =  year - spinupyears + firstyeartrend - 1 
        out_steering%climateyear = out_steering%forcingyear

      endif
      outyear = year + firstyeartrend - spinupyears - 1

    else

      out_steering%forcingyear = year + firstyeartrend - 1 
      out_steering%climateyear = out_steering%forcingyear
      outyear = year + firstyeartrend - 1

    endif

    ! write(0,*) 'recycle, spinupyears, firstyeartrend', recycle, spinupyears, firstyeartrend
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
    ! argument
    character(len=*), intent(in) :: runname

    ! Read in main model parameters
    write(0,*) 'reading parameter file ', runname//".sofun.parameter ..."

    ! sitename         = getparstring( 'run/'//runname//'.sofun.parameter', 'sitename' )
    ! input_dir        = getparstring( 'run/'//runname//'.sofun.parameter', 'input_dir' )
    ! co2_forcing_file = getparstring( 'run/'//runname//'.sofun.parameter', 'co2_forcing_file' )
    ! pftparfil        = getparstring( 'run/'//runname//'.sofun.parameter', 'pftparfil' )

    call getparstring( 'run/'//runname//'.sofun.parameter', 'sitename', sitename )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'input_dir', input_dir )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'co2_forcing_file', co2_forcing_file )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'ndep_noy_forcing_file', ndep_noy_forcing_file )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'ndep_nhx_forcing_file', ndep_nhx_forcing_file )
    call getparstring( 'run/'//runname//'.sofun.parameter', 'pftparfil', pftparfil )

    do_spinup            = getparlogical( 'run/'//runname//'.sofun.parameter', 'spinup' )
    const_co2            = getparlogical( 'run/'//runname//'.sofun.parameter', 'const_co2' )
    const_ndep           = getparlogical( 'run/'//runname//'.sofun.parameter', 'const_ndep' )
    
    spinupyears          = getparint( 'run/'//runname//'.sofun.parameter', 'spinupyears' )
    firstyeartrend       = getparint( 'run/'//runname//'.sofun.parameter', 'firstyeartrend' )
    nyeartrend           = getparint( 'run/'//runname//'.sofun.parameter', 'nyeartrend' )
    recycle              = getparint( 'run/'//runname//'.sofun.parameter', 'recycle' )
    
    prescr_monthly_fapar = getparlogical( 'run/'//runname//'.sofun.parameter', 'prescr_monthly_fapar' )
    
    daily_out_startyr    = getparint( 'run/'//runname//'.sofun.parameter', 'daily_out_startyr' )
    daily_out_endyr      = getparint( 'run/'//runname//'.sofun.parameter', 'daily_out_endyr' )

    if (do_spinup) then
      runyears = nyeartrend + spinupyears
    else
      runyears = nyeartrend
      spinupyears = 0
    endif

    loutdgpp       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdgpp' )
    loutdtransp    = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdtransp' )
    loutdnpp       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdnpp' )
    loutdnup       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdnup' )
    loutdCleaf     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdCleaf' )
    loutdCroot     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdCroot' )
    loutdClabl     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdClabl' )
    loutdNlabl     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdNlabl' )
    loutdClitt     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdClitt' )
    loutdNlitt     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdNlitt' )
    loutdCsoil     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdCsoil' )
    loutdNsoil     = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdNsoil' )
    loutdlai       = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdlai' )
    loutdninorg    = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdninorg' )
    loutdtemp_soil = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutdtemp_soil' )

    loutntransform = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutntransform') 
    loutwaterbal   = getparlogical( 'run/'//runname//'.sofun.parameter', 'loutwaterbal')

  end subroutine getpar_siml

end module _params_siml

