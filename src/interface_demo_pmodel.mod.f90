module md_interface
  !////////////////////////////////////////////////////////////////
  ! Module interface defines nested derived type that contains all
  ! the necessary simulation parameters, etc.
  !
  ! Author: Benjamin D. Stocker
  !----------------------------------------------------------------
  implicit none

  private
  public interface, get_interface

  type paramstype_siml

    integer :: daily_out_startyr! first year where daily output is written
    integer :: daily_out_endyr ! last year where daily output is written
    integer :: outdt           ! output periodicity
    integer :: outnt           ! number of output time steps per year
    
    logical :: is_calib             ! whether this simulation is a calibration simulation (overriding parameters and no output)

    logical :: soilmstress          ! when true, an empirical soil moisture stress function is applied to GPP
    logical :: tempstress           ! when true, an empirical temperature stress function is applied to GPP
    
    character(len=256) :: runname
    character(len=256) :: fapar_forcing_source

    ! booleans defining whether variable is written to ascii output
    logical :: loutdgpp       
    logical :: loutdrd
    logical :: loutdtransp    

    ! booleans defining whether module-specific output variables are to be written to output
    logical :: loutgpp

    ! booleans defining whether variable is written to NetCDF output
    logical :: lncoutdgpp

    ! booleans defining whether variable is used as calibration target
    logical :: lcalibgpp

  end type paramstype_siml


  type outtype_steering
    integer :: outyear         ! year AD written to output
    logical :: spinup          ! is true during spinup
    logical :: init            ! is true in first simulation year
  end type outtype_steering


  type domaininfo_type
    integer :: nlon
    integer :: nlat
    integer :: maxgrid
    real, dimension(:), allocatable :: lon
    real, dimension(:), allocatable :: lat
  end type domaininfo_type


  type gridtype
    integer :: ilon
    integer :: ilat
    logical :: dogridcell
  end type gridtype


  type paramstype_calib
    real :: kphio
    real :: temp_ramp_edge
    real :: soilm_par_a
    real :: soilm_par_b
  end type paramstype_calib  


  type interfacetype_biosphere
    type( gridtype ), dimension(:),         allocatable :: grid
    type( domaininfo_type )                             :: domaininfo
    type( outtype_steering )                            :: steering
    type( paramstype_siml )                             :: params_siml
    type( paramstype_calib )                            :: params_calib    ! calibratable parameters
  end type interfacetype_biosphere

  !----------------------------------------------------------------
  ! Interface instance is created here 
  ! (instead of locally defined and passed on as argument. Both are 
  ! ok but this has the advantage that unknown-size arguments are
  ! avoided).
  !----------------------------------------------------------------
  type( interfacetype_biosphere ) :: interface

contains

  subroutine get_interface()

    ! simulation parameters
    interface%params_siml%runname           = "demo"
    interface%params_siml%daily_out_startyr = 1
    interface%params_siml%daily_out_endyr   = 1
    interface%params_siml%outdt             = 1
    interface%params_siml%outnt             = 1
    interface%params_siml%is_calib          = .false.
    interface%params_siml%soilmstress       = .false.
    interface%params_siml%loutdgpp          = .false.
    interface%params_siml%loutdrd           = .false. 
    interface%params_siml%loutdtransp       = .false.     
    interface%params_siml%loutgpp           = .false. 
    interface%params_siml%lncoutdgpp        = .false. 
    interface%params_siml%lcalibgpp         = .false. 

    ! others don't need to be defined as they are not used in demo setup

  end subroutine get_interface

end module md_interface
