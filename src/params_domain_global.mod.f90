module md_params_domain
  !////////////////////////////////////////////////////////////////
  !  Module contains site-specific parameters read by getpar_domain
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: maxgrid

  implicit none

  private 
  public getpar_domain, type_params_domain

  type type_params_domain

    ! For global simulations
    character(len=256) :: filnam_landmask
    character(len=256) :: filnam_topography
    character(len=256) :: filnam_soil

  end type type_params_domain

contains

  function getpar_domain( sitename ) result( out_params_domain )
    !////////////////////////////////////////////////////////////////
    !  SR reads site-specific parameters from <sitename>.parameter
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal, getparlogical, getparint, getparstring

    ! arguments
    character(len=*) :: sitename

    ! function return variable
    type( type_params_domain ) :: out_params_domain
  
    print*, 'reading site parameter file ', 'site_paramfils/'//trim(sitename)//'.parameter ...'

    call getparstring( 'site_paramfils/'//trim(sitename)//'.parameter', 'filnam_landmask',   out_params_domain%filnam_landmask )
    call getparstring( 'site_paramfils/'//trim(sitename)//'.parameter', 'filnam_topography', out_params_domain%filnam_topography )
    call getparstring( 'site_paramfils/'//trim(sitename)//'.parameter', 'filnam_soil',       out_params_domain%filnam_soil )

  end function getpar_domain

end module md_params_domain
