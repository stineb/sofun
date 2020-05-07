module md_gpp
  !//////////////////////////////////////////////////////////////////////
  ! GPP MODULE
  ! Uses LM3-PPA structure to call the P-model photosynthesis routine
  !------------------------------------------------------------------------
  use datatypes

  implicit none

  private
  public gpp

  !-----------------------------------------------------------------------
  ! Uncertain (unknown) parameters. Hard-coded here (should be read runtime)
  !-----------------------------------------------------------------------
  type paramstype_gpp
    real :: beta = 146.0         ! Unit cost of carboxylation (dimensionless)
    real :: soilm_par_a = 0.0
    real :: soilm_par_b = 0.0
    real :: rd_to_vcmax = 0.014  ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
    real :: tau_acclim = 10.0    ! acclimation time scale of photosynthesis (d)
  end type paramstype_gpp

  type(paramstype_gpp) :: params_gpp
  ! params_gpp%beta = 146.0         ! Unit cost of carboxylation (dimensionless)
  ! params_gpp%soilm_par_a = 0.0
  ! params_gpp%soilm_par_b = 0.0
  ! params_gpp%rd_to_vcmax = 0.014  ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  ! params_gpp%tau_acclim = 10.0    ! acclimation time scale of photosynthesis (d)

  ! PFT-DEPENDENT PARAMETERS
  type pftparamstype_gpp
    real :: kphio = 0.24    ! hard-coded here, is a calibratable parameter in P-model, unrealistically high here to match ballpark of original model
  end type pftparamstype_gpp

  type(pftparamstype_gpp) :: params_pft_gpp

contains

  subroutine gpp( forcing, vegn, init)
    !//////////////////////////////////////////////////////////////////////
    ! GPP
    ! Calculates light availability and photosynthesis for each cohort 
    ! and sets the following cohort-level variables:
    ! - An_op   
    ! - An_cl   
    ! - w_scale 
    ! - transp  
    ! - gpp     
    !
    ! Calling the P-model for photosynthesis.
    !
    ! Adjustments to the forcing type have been made, requiring:
    ! %dtemp
    ! %dvpd
    ! %dpatm
    ! %dppfd
    ! 
    ! The module gpp_pmodel.mod.f90 cannot be used yet because the tile 
    ! structure has not been made consistent yet. Required structures are
    ! (SOFUN structures for LM3-PPA):
    ! - tile
    ! - tile_fluxes
    ! - climate (using the variables SOFUN uses)
    ! - params_pft_plant
    ! - time step length in seconds (base units should always be seconds)
    ! - remove 'd' from forcing component variables
    ! 
    !------------------------------------------------------------------------
    use md_forcing, only: climate_type
    use md_photosynth, only: pmodel, outtype_pmodel
    use md_params_core, only: kTkelvin, kfFEC
    use md_sofunutils, only: dampen_variability

    type(climate_type), intent(in):: forcing
    type(vegn_tile_type), intent(inout) :: vegn
    logical, intent(in) :: init   ! is true on the very first simulation day (first subroutine call of each gridcell)

    ! local variables
    type(cohort_type), pointer :: cc
    real   :: kappa  ! light extinction coefficient of corwn layers
    real   :: f_light(10) = 0.0, f_apar(10) = 0.0      ! incident light fraction, and aborbed light fraction of each layer
    real   :: LAIlayer(10), crownarea_layer(10), accuCAI, f_gap, fapar_tree, rad_top ! additional GPP for lower layer cohorts due to gaps
    integer:: i, layer

    real :: tk
    real, save :: co2_memory
    real, save :: vpd_memory
    real, save :: temp_memory
    real, save :: patm_memory

    type(outtype_pmodel) :: out_pmodel      ! list of P-model output variables

    logical, parameter :: verbose = .false.

    !----------------------------------------------------------------
    ! Calculate environmental conditions with memory, time scale 
    ! relevant for Rubisco turnover
    !----------------------------------------------------------------
    if (init) then
      co2_memory  = forcing%CO2 * 1.0e6
      temp_memory = (forcing%Tair - kTkelvin)
      vpd_memory  = forcing%vpd
      patm_memory = forcing%P_air
    end if 

    co2_memory  = dampen_variability( forcing%CO2 * 1.0e6,        params_gpp%tau_acclim, co2_memory )
    temp_memory = dampen_variability( (forcing%Tair - kTkelvin),  params_gpp%tau_acclim, temp_memory )
    vpd_memory  = dampen_variability( forcing%vpd,                params_gpp%tau_acclim, vpd_memory )
    patm_memory = dampen_variability( forcing%P_air,              params_gpp%tau_acclim, patm_memory )

    tk = forcing%Tair + kTkelvin

    !----------------------------------------------------------------
    ! Sum leaf area over cohorts in each crown layer -> LAIlayer(layer)
    !----------------------------------------------------------------
    f_gap = 0.1 ! 0.1
    accuCAI = 0.0
    LAIlayer(:) = 0.0
    crownarea_layer(:) = 0.0
    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      layer = max(1, min(cc%layer,9))
      LAIlayer(layer) = LAIlayer(layer) + cc%leafarea * cc%nindivs / (1.0 - f_gap)
      crownarea_layer(layer) = crownarea_layer(layer) + cc%crownarea
    end do

    
    !----------------------------------------------------------------
    ! Get light fraction received at each crown layer, relative to top-of-canopy -> f_light(layer) 
    !----------------------------------------------------------------
    ! Calculate kappa according to sun zenith angle 
    ! kappa = cc%extinct/max(cosz,0.01)

    ! Use constant light extinction coefficient
    kappa = cc%extinct
    f_light(:) = 0.0
    f_light(1) = 1.0
    do i=2,layer+1
      f_light(i) = f_light(i-1) * exp(0.0 - kappa * LAIlayer(i-1))
    end do

    ! lighten up the canopy (gaps)
    f_light(:) = f_light(:) * (1.0 - f_gap) + f_gap

    ! fraction of light absorbed by layer
    do i=1,layer
      f_apar(i) = f_light(i) - f_light(i+1)
      if (f_apar(i) < 0.0 ) stop 'negative fapar'
    end do

    !----------------------------------------------------------------
    ! Photosynthesis for each cohort
    !----------------------------------------------------------------
    accuCAI = 0.0

    do i = 1, vegn%n_cohorts

      cc => vegn%cohorts(i)
      associate ( sp => spdata(cc%species) )

      if (cc%status == LEAF_ON .and. cc%lai > 0.1 .and. temp_memory > -5.0) then

        !----------------------------------------------------------------
        ! Get light aborbed by cohort, dividing fAPAR up by crown areas
        !----------------------------------------------------------------
        layer = max(1, min(cc%layer,9))
        fapar_tree = 1.0 - exp(-kappa * cc%leafarea / cc%crownarea)   ! at individual-level: cc%leafarea represents leaf area index within the crown 

        !----------------------------------------------------------------
        ! P-model call for C3 plants to get a list of variables that are 
        ! acclimated to slowly varying conditions
        !----------------------------------------------------------------
        if (fapar_tree > 0.0 .and. forcing%PAR > 0.0) then

          ! !===============================
          ! ! XXX constant LUE hack:
          ! !===============================
          ! cc%An_op   = 1.0e-7 * fapar_tree * f_light(layer) * forcing%PAR / kfFEC  ! molC s-1 m-2 of leaves
          ! cc%An_cl   = 0.5e-9 * fapar_tree * f_light(layer) * forcing%PAR / kfFEC  ! molC s-1 m-2 of leaves
          ! cc%w_scale = 0.0
          ! cc%transp  = 0.0
          ! cc%resl    = cc%An_cl              * mol_C * myinterface%step_seconds ! kgC tree-1 step-1
          ! cc%gpp     = (cc%An_op + cc%An_cl) * mol_C * myinterface%step_seconds ! kgC step-1 tree-1
          ! !===============================
          if (verbose) print*,'calling pmodel...'
          if (verbose) print*,'      fapar'          , fapar_tree
          if (verbose) print*,'      ppfd'           , f_light(layer) * forcing%PAR * 1.0e-6    ! required in mol m-2 s-1
          if (verbose) print*,'      co2'            , co2_memory
          if (verbose) print*,'      tc'             , temp_memory
          if (verbose) print*,'      vpd'            , vpd_memory
          if (verbose) print*,'      patm'           , patm_memory
          if (verbose) print*,'      c4'             , .false.
          if (verbose) print*,'      method_optci'   , "prentice14"
          if (verbose) print*,'      method_jmaxlim' , "wang17"
          if (verbose) print*,'      kphio'          , params_pft_gpp%kphio
          if (verbose) print*,'      beta'           , params_gpp%beta
          if (verbose) print*,'      rd_to_vcmax'    , params_gpp%rd_to_vcmax

          out_pmodel = pmodel( &
                              fapar          = fapar_tree, &
                              ppfd           = f_light(layer) * forcing%PAR * 1.0e-6, &    ! required in mol m-2 s-1
                              co2            = co2_memory, &
                              tc             = temp_memory, &
                              vpd            = vpd_memory, &
                              patm           = patm_memory, &
                              c4             = .false., &
                              method_optci   = "prentice14", &
                              method_jmaxlim = "wang17", &
                              kphio          = params_pft_gpp%kphio, &
                              beta           = params_gpp%beta, &
                              rd_to_vcmax    = params_gpp%rd_to_vcmax &
                              )
          if (verbose) print*,'done.'

          ! irrelevant variables for this setup  
          cc%An_op   = 0.0
          cc%An_cl   = 0.0
          cc%transp  = 0.0
          cc%w_scale = -9999

          ! copy to cohort variables
          cc%resl    = out_pmodel%rd  * cc%crownarea * myinterface%step_seconds * mol_C     ! kgC step-1 tree-1
          cc%gpp     = out_pmodel%gpp * cc%crownarea * myinterface%step_seconds * 1.0e-3    ! kgC step-1 tree-1


        else

          cc%An_op   = 0.0
          cc%An_cl   = 0.0
          cc%transp  = 0.0
          cc%w_scale = -9999
          cc%resl    = 0.0
          cc%gpp     = 0.0

        end if

      else

        ! no leaves means no photosynthesis and no stomatal conductance either
        cc%An_op   = 0.0
        cc%An_cl   = 0.0
        cc%transp  = 0.0
        cc%w_scale = -9999
        cc%resl    = 0.0
        cc%gpp     = 0.0

      endif

      end associate

    end do

  end subroutine gpp


  function calc_soilmstress( soilm, meanalpha, isgrass ) result( outstress )
    !//////////////////////////////////////////////////////////////////
    ! Calculates empirically-derived stress (fractional reduction in light 
    ! use efficiency) as a function of soil moisture
    ! Input:  soilm (unitless, within [0,1]): daily varying soil moisture
    ! Output: outstress (unitless, within [0,1]): function of alpha to reduce GPP 
    !         in strongly water-stressed months
    !-----------------------------------------------------------------------
    ! argument
    real, intent(in) :: soilm                 ! soil water content (fraction)
    real, intent(in) :: meanalpha             ! mean annual AET/PET, average over multiple years (fraction)
    logical, intent(in), optional :: isgrass  ! vegetation cover information to distinguish sensitivity to low soil moisture

    real, parameter :: x0 = 0.0
    real, parameter :: x1 = 0.6

    real :: y0, beta

    ! function return variable
    real :: outstress

    if (soilm > x1) then
      outstress = 1.0
    else
      ! print*,'soilm_par_a, soilm_par_b, meanalpha', params_gpp%soilm_par_a, params_gpp%soilm_par_b, meanalpha

      y0 = (params_gpp%soilm_par_a + params_gpp%soilm_par_b * meanalpha)

      ! if (present(isgrass)) then
      !   if (isgrass) then
      !     y0 = apar_grass + bpar_grass * meanalpha
      !   else
      !     y0 = apar + bpar * meanalpha
      !   end if
      ! else
      !   y0 = apar + bpar * meanalpha
      ! end if

      beta = (1.0 - y0) / (x0 - x1)**2
      outstress = 1.0 - beta * ( soilm - x1 )**2
      outstress = max( 0.0, min( 1.0, outstress ) )
    end if

  end function calc_soilmstress


  function calc_ftemp_kphio( dtemp ) result( ftemp )
    !////////////////////////////////////////////////////////////////
    ! Calculates the instantaneous temperature response of the quantum
    ! yield efficiency based on Bernacchi et al., 2003 PCE (Equation
    ! and parameter values taken from Appendix B)
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: dtemp

    ! function return variable
    real :: ftemp

    ftemp = 0.352 + 0.022 * dtemp - 3.4e-4 * dtemp**2

  end function calc_ftemp_kphio


end module md_gpp
