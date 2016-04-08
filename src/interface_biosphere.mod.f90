module md_interface

  use md_params_core, only: maxgrid, nmonth
  use md_grid, only: gridtype
  use md_forcing_siterun, only: landuse_type, climate_type, ndep_type
  use md_params_site, only: paramstype_site
  use md_params_soil, only: paramtype_soil
  use md_params_siml, only: outtype_steering, paramstype_siml

  implicit none

  private
  public interfacetype_biosphere, interface

  type interfacetype_biosphere
    integer                                           :: year
    real                                              :: pco2
    type( gridtype )      , dimension(maxgrid)        :: grid
    type( paramtype_soil ), dimension(maxgrid)        :: soilparams
    type( landuse_type)   , dimension(maxgrid)        :: landuse
    type( climate_type )  , dimension(maxgrid)        :: climate
    type( ndep_type)      , dimension(maxgrid)        :: ndep_field
    real                  , dimension(nmonth,maxgrid) :: mfapar_field
    type( paramstype_site )                           :: params_site
    type( outtype_steering )                          :: steering
    type( paramstype_siml )                           :: params_siml
  end type interfacetype_biosphere

  type( interfacetype_biosphere ) :: interface

end module md_interface
