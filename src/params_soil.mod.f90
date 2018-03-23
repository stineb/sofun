module md_params_soil
  !////////////////////////////////////////////////////////////////
  ! Module containing soil parameters and functions to read them
  !----------------------------------------------------------------
  implicit none

  private
  public paramtype_soil, getsoil_field

  type paramtype_soil
    real :: whc
  end type

contains

  function getsoil_field( grid ) result( params_soil_field )
    !////////////////////////////////////////////////////////////////
    ! Function returns array containing all soil parameter values
    !----------------------------------------------------------------
    use md_params_core, only: maxgrid
    use md_grid, only: gridtype

    ! arguments
    type( gridtype ), intent(in), dimension(maxgrid) :: grid

    ! local variables
    integer :: jpngr

    ! function return variable
    type(paramtype_soil), dimension(maxgrid) :: params_soil_field

    if (maxgrid>1) stop 'in getsoil_field: think of something'

    do jpngr=1,maxgrid
      params_soil_field(jpngr) = getsoil( 1 )
    end do

  end function getsoil_field


  function getsoil( soilcode ) result( params_soil )
    !////////////////////////////////////////////////////////////////
    ! Function returns soil parameter values given soil code
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    integer, intent(in) :: soilcode

    ! local variables
    character(len=2) :: soilcode_char

    ! function return variable
    type(paramtype_soil) :: params_soil

    write( soilcode_char, 999 ) soilcode

    params_soil%perc_k1      = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'perc_k1' )
    params_soil%whc_eff      = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'whc_eff' )
    params_soil%thdiff_wp    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_wp' )
    params_soil%thdiff_whc15 = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_whc15' )
    params_soil%thdiff_fc    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_fc' )
    params_soil%forg         = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'forg' )
    params_soil%whc_blwp     = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'whc_blwp' )
    params_soil%por          = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'por' )
    params_soil%fsand        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsand' )
    params_soil%fclay        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fclay' )
    params_soil%fsilt        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsilt' )

    return
    999  format (I2.2)

  end function getsoil

end module md_params_soil
