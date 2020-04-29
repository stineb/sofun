module md_gpp
  !//////////////////////////////////////////////////////////////////////
  ! GPP MODULE
  ! Uses LM3-PPA structure to call the gs_Leuning() photosynthesis routine
  !------------------------------------------------------------------------
  use datatypes
  use md_soil, only: water_supply_layer

  implicit none

  private
  public gpp

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
    ! Subroutines from BiomeE-Allocation
    !------------------------------------------------------------------------
    use md_forcing, only: climate_type

    type(climate_type), intent(in):: forcing
    type(vegn_tile_type), intent(inout) :: vegn
    logical, intent(in) :: init   ! is true on the very first simulation day (first subroutine call of each gridcell)

    ! local variables
    type(cohort_type), pointer :: cc
    real   :: rad_top  ! downward radiation at the top of the canopy, W/m2
    real   :: rad_net  ! net radiation absorbed by the canopy, W/m2
    real   :: Tair, TairK     ! air temperature, degC and degK
    real   :: cana_q   ! specific humidity in canopy air space, kg/kg
    real   :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
    real   :: p_surf   ! surface pressure, Pa
    real   :: water_supply ! water supply per m2 of leaves
    real   :: fw, fs ! wet and snow-covered fraction of leaves
    real   :: psyn   ! net photosynthesis, mol C/(m2 of leaves s)
    real   :: resp   ! leaf respiration, mol C/(m2 of leaves s)
    real   :: tempLAI,w_scale2, transp ! mol H20 per m2 of leaf per second
    real   :: kappa  ! light extinction coefficient of corwn layers
    real   :: f_light(10)=0.0      ! light fraction of each layer
    real   :: LAIlayer(10),accuCAI,f_gap ! additional GPP for lower layer cohorts due to gaps
    integer:: i, layer

    !===========================================================
    ! Original BiomeE-Allocation
    !-----------------------------------------------------------
    ! Water supply for photosynthesis, Layers
    call water_supply_layer(forcing, vegn)

    ! Sum leaf area over cohorts in each crown layer -> LAIlayer(layer)
    f_gap = 0.1 ! 0.1
    accuCAI = 0.0
    LAIlayer = 0.0
    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      layer = Max(1, Min(cc%layer,9))
      LAIlayer(layer) = LAIlayer(layer) + cc%leafarea * cc%nindivs /(1.-f_gap)
    enddo

    ! ! Calculate kappa according to sun zenith angle 
    ! kappa = cc%extinct/max(cosz,0.01)
    
    ! Use constant light extinction coefficient
    kappa = cc%extinct

    ! Get light received at each crown layer as a fraction of top-of-canopy -> f_light(layer) 
    f_light(:) = 0.0
    f_light(1) = 1.0
    do i=2,layer
      f_light(i) = f_light(i-1) * (exp(0.0 - kappa * LAIlayer(i-1)) + f_gap)
    enddo

    ! Photosynthesis
    accuCAI = 0.0

    do i = 1, vegn%n_cohorts

      cc => vegn%cohorts(i)
      associate ( sp => spdata(cc%species) )

      if (cc%status == LEAF_ON .and. cc%lai > 0.1) then

        ! Convert forcing data
        layer    = Max (1, Min(cc%layer,9))
        rad_top  = f_light(layer) * forcing%radiation ! downward radiation at the top of the canopy, W/m2
        rad_net  = f_light(layer) * forcing%radiation * 0.9 ! net radiation absorbed by the canopy, W/m2
        p_surf   = forcing%P_air  ! Pa
        TairK    = forcing%Tair ! K
        Tair     = forcing%Tair - 273.16 ! degC
        cana_q   = (esat(Tair) * forcing%RH * mol_h2o) / (p_surf * mol_air)  ! air specific humidity, kg/kg
        cana_co2 = forcing%CO2 ! co2 concentration in canopy air space, mol CO2/mol dry air

        ! recalculate the water supply to mol H20 per m2 of leaf per second
        water_supply = cc%W_supply / (cc%leafarea * myinterface%step_seconds * mol_h2o) ! mol m-2 leafarea s-1

        !call get_vegn_wet_frac (cohort, fw=fw, fs=fs)
        fw = 0.0
        fs = 0.0

        call gs_Leuning(rad_top, rad_net, TairK, cana_q, cc%lai, &
          p_surf, water_supply, cc%species, sp%pt, &
          cana_co2, cc%extinct, fs+fw, cc%layer, &
          ! output:
          psyn, resp, w_scale2, transp )

        ! store the calculated photosynthesis, photorespiration, and transpiration for future use in growth
        cc%An_op   = psyn  ! molC s-1 m-2 of leaves
        cc%An_cl   = -resp  ! molC s-1 m-2 of leaves
        cc%w_scale = w_scale2
        cc%transp  = transp * mol_h2o * cc%leafarea * myinterface%step_seconds ! Transpiration (kgH2O/(tree step), Weng, 2017-10-16
        cc%resl    = -resp         * mol_C * cc%leafarea * myinterface%step_seconds ! fnsc*spdata(sp)%gamma_LN  * cc%leafN * tf * myinterface%dt_fast_yr  ! tree-1 step-1
        cc%gpp     = (psyn - resp) * mol_C * cc%leafarea * myinterface%step_seconds ! kgC step-1 tree-1

        ! if (rad_top > 0.0) print*,'psyn/rad_top, resp/rad_top ', psyn/rad_top, resp/rad_top

        ! ! psyn/rad_top is on the order of 1e-8; psyn/rad_top is on the order of -1e-9
        ! if (rad_top > 0.0) print*,'psyn/rad_top, resp/rad_top', psyn/rad_top, resp/rad_top

        if (isnan(cc%gpp)) stop '"gpp" is a NaN'

      else

        ! no leaves means no photosynthesis and no stomatal conductance either
        cc%An_op   = 0.0
        cc%An_cl   = 0.0
        cc%gpp     = 0.0
        cc%transp  = 0.0
        cc%w_scale = -9999

      endif

      end associate
    enddo

  end subroutine gpp


  subroutine gs_Leuning( rad_top, rad_net, tl, ea, lai, &
    p_surf, ws, pft, pt, ca, kappa, leaf_wet, layer, &
    apot, acl,w_scale2, transp )

    real,    intent(in)    :: rad_top ! PAR dn on top of the canopy, w/m2
    real,    intent(in)    :: rad_net ! PAR net on top of the canopy, w/m2
    real,    intent(in)    :: tl   ! leaf temperature, degK
    real,    intent(in)    :: ea   ! specific humidity in the canopy air, kg/kg
    real,    intent(in)    :: lai  ! leaf area index
    !real,    intent(in)    :: leaf_age ! age of leaf since budburst (deciduos), days
    real,    intent(in)    :: p_surf ! surface pressure, Pa
    real,    intent(in)    :: ws   ! water supply, mol H20/(m2 of leaf s)
    integer, intent(in)    :: pft  ! species
    integer, intent(in)    :: pt   ! physiology type (C3 or C4)
    real,    intent(in)    :: ca   ! concentartion of CO2 in the canopy air space, mol CO2/mol dry air
    real,    intent(in)    :: kappa! canopy extinction coefficient (move inside f(pft))
    real,    intent(in)    :: leaf_wet ! fraction of leaf that's wet or snow-covered
    integer, intent(in)    :: layer  ! the layer of this canopy
    ! note that the output is per area of leaf; to get the quantities per area of
    ! land, multiply them by LAI
    !real,    intent(out)   :: gs   ! stomatal conductance, m/s
    real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
    real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
    real,    intent(out)   :: w_scale2,transp  ! transpiration, mol H20/(m2 of leaf s)
    ! local variables     
    ! photosynthesis
    real :: vm
    real :: kc, ko ! Michaelis-Menten constants for CO2 and O2, respectively
    real :: ci
    real :: capgam ! CO2 compensation point
    real :: f2, f3
    real :: coef0, coef1
    real :: Resp
    ! conductance related
    real :: gs
    real :: b
    real :: ds  ! humidity deficit, kg/kg
    real :: hl  ! saturated specific humidity at the leaf temperature, kg/kg
    real :: do1
    ! misceleneous
    real :: dum2
    real, parameter :: light_crit = 0
    real, parameter :: gs_lim = 0.25
    real, parameter :: Rgas = 8.314 ! J mol-1 K-1, universal gas constant
    ! new average computations
    real :: lai_eq
    real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta 
    real :: light_top
    real :: par_net
    real :: Ag
    real :: An
    real :: Ag_l
    real :: Ag_rb
    real :: anbar
    real :: gsbar
    real :: w_scale
    real, parameter :: p_sea = 1.0e5 ! sea level pressure, Pa
    ! soil water stress
    real :: Ed, an_w, gs_w
    b = 0.01
    do1 = 0.09 ! kg/kg
    if (pft < 2) do1 = 0.15
    ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
    ! empirical relationship from McCree is light=rn*0.0000046
    light_top = rad_top*rad_phot;
    par_net   = rad_net*rad_phot;
    ! calculate humidity deficit, kg/kg
    call qscomp(tl, p_surf, hl)
    ds = max(hl - ea,0.0)
    !  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
    !  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
    !  vm=spdata(pft)%Vmax*exp(3000.0*(1.0/288.2-1.0/tl));
    ! corrected by Weng, 2013-01-17
    ! Weng, 2013-01-10
    ko=0.248    * exp(35948/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
    kc=0.000404 * exp(59356/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10
    vm=spdata(pft)%Vmax*exp(24920/Rgas*(1.0/298.2-1.0/tl)) ! / ((layer-1)*1.0+1.0) ! Ea = 33920
    !decrease Vmax due to aging of temperate deciduous leaves 
    !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
    !! Turned off by Weng, 2013-02-01, since we can't trace new leaves
    !  if (spdata(pft)%leaf_age_tau>0 .and. leaf_age>spdata(pft)%leaf_age_onset) then
    !     vm=vm*exp(-(leaf_age-spdata(pft)%leaf_age_onset)/spdata(pft)%leaf_age_tau)
    !  endif

    ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986
    capgam=0.5*kc/ko*0.21*0.209; ! Farquhar & Caemmerer 1982



    ! Find respiration for the whole canopy layer

    !  Resp=spdata(pft)%gamma_resp*vm*lai /((layer-1)*1.0+1.0)  ! Weng, 2013-01-17 add '/ ((layer-1)*1.0+1.0)'

    ! 2014-09-03, for Nitrogen model: resp = D*(A + B*LMA)
    ! (A+B*LMA) = LNA, D=Vmax/LNA = 25E-6/0.0012 = 0.02 for a standard deciduous species
    !! Leaf resp as a function of nitrogen
    !  Resp=spdata(pft)%gamma_resp*0.04*spdata(pft)%LNA  & ! basal rate, mol m-2 s-1
    !       * exp(24920/Rgas*(1.0/298.2-1.0/tl))         & ! temperature scaled
    !       * lai                                        & ! whole canopy
    !       /((layer-1)*1.0+1.0)                         !
    !! as a function of LMA
    !  Resp=(spdata(pft)%gamma_LNbase*spdata(pft)%LNbase+spdata(pft)%gamma_LMA*spdata(pft)%LMA)  & ! basal rate, mol m-2 s-1
    !  Resp=spdata(pft)%gamma_LNbase*(2.5*spdata(pft)%LNA-1.5*spdata(pft)%LNbase)     & ! basal rate, mol m-2 s-1
    Resp = spdata(pft)%gamma_LN/seconds_per_year & ! per seconds, ! basal rate, mol m-2 s-1
            * spdata(pft)%LNA * lai / mol_c    &     ! whole canopy, ! basal rate, mol m-2 s-1
            * exp(24920/Rgas*(1.0/298.2-1.0/tl))     ! temperature scaled
    !                                  &
    !       /((layer-1)*1.0+1.0)
    ! Temperature effects
    Resp=Resp/((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));


    ! ignore the difference in concentrations of CO2 near
    !  the leaf and in the canopy air, rb=0.
    Ag_l=0.;
    Ag_rb=0.;
    Ag=0.;
    anbar=-Resp/lai;
    gsbar=b;
    ! find the LAI level at which gross photosynthesis rates are equal
    ! only if PAR is positive
    if ( light_top > light_crit ) then

      if (pt==PT_C4) then ! C4 species

        coef0=(1+ds/do1)/spdata(pft)%m_cond
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0)

        if (ci>capgam) then
          f2=vm
          f3=18000.0*vm*ci ! 18000 or 1800?
          dum2=min(f2,f3)

          ! find LAI level at which rubisco limited rate is equal to light limited rate
          lai_eq = -log(dum2/(kappa*spdata(pft)%alpha_phot*light_top))/kappa
          lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

          ! gross photosynthesis for light-limited part of the canopy
          Ag_l   = spdata(pft)%alpha_phot * par_net &
                  * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))

          ! gross photosynthesis for rubisco-limited part of the canopy
          Ag_rb  = dum2*lai_eq

          Ag=(Ag_l+Ag_rb)/((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))))
          An=Ag-Resp
          anbar=An/lai

          if (anbar>0.0) then
            gsbar=anbar/(ci-capgam)/coef0
          endif

        endif ! ci>capgam

      else ! C3 species

        coef0=(1+ds/do1)/spdata(pft)%m_cond
        coef1=kc*(1.0+0.209/ko)
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0)
        f2=vm*(ci-capgam)/(ci+coef1)
        f3=vm/2.
        dum2=min(f2,f3)

        if (ci>capgam) then
          ! find LAI level at which rubisco limited rate is equal to light limited rate
          lai_eq=-log(dum2*(ci+2.*capgam)/(ci-capgam)/ &
           (spdata(pft)%alpha_phot*light_top*kappa))/kappa
          lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

          ! gross photosynthesis for light-limited part of the canopy
          Ag_l   = spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * par_net &
                   * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1.0-exp(-lai*kappa))

          ! gross photosynthesis for rubisco-limited part of the canopy
          Ag_rb  = dum2*lai_eq

          Ag=(Ag_l+Ag_rb) /((1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE))));
          An=Ag-Resp;
          anbar=An/lai
          !write(*,*)'An,Ag,AG_l,Ag_rb,lai',An,Ag, Ag_l, Ag_rb,lai

          if (anbar>0.0) then
            gsbar=anbar/(ci-capgam)/coef0;
          endif

        endif

      endif

    endif ! light is available for photosynthesis

    !write(898,'(1(I4,","),10(E10.4,","))') &
    !     layer, light_top, par_net, kappa, lai, lai_eq, ci, capgam, Ag_l, Ag_rb, Ag

    an_w=anbar

    if (an_w > 0.) then
      an_w=an_w*(1-spdata(pft)%wet_leaf_dreg*leaf_wet);
    endif
    gs_w = 1.56 * gsbar *(1-spdata(pft)%wet_leaf_dreg*leaf_wet); !Weng: 1.56 for H2O?

    if (gs_w > gs_lim) then
      if (an_w > 0.) an_w = an_w*gs_lim/gs_w
      gs_w = gs_lim
    endif

    ! find water availability diagnostic demand
    Ed = gs_w * ds * mol_air / mol_h2o ! ds*mol_air/mol_h2o is the humidity deficit in [mol_h2o/mol_air]

    ! the factor mol_air/mol_h2o makes units of gs_w and humidity deficit ds compatible:
    if (Ed>ws) then
      w_scale = ws/Ed
      gs_w = w_scale * gs_w
      if (an_w > 0.0) an_w = an_w * w_scale
      if (an_w < 0.0 .and. gs_w >b) gs_w = b
    endif

    gs=gs_w
    apot=an_w
    acl=-Resp/lai
    transp = min(ws,Ed) ! mol H20/(m2 of leaf s)
    ! just for reporting
    if (Ed>0.0) then
      w_scale2=min(1.0,ws/Ed)
    else
      w_scale2=1.0
    end if 

    ! finally, convert units of stomatal conductance to m/s from mol/(m2 s) by
    ! multiplying it by a volume of a mole of gas
    gs = gs * Rugas * Tl / p_surf
    !write(899, '(25(E12.4,","))') rad_net,par_net,apot*3600*12,acl*3600*12,Ed

  end subroutine gs_Leuning


  subroutine calc_solarzen(td, latdegrees, cosz, solarelev, solarzen)
    ! Calculate solar zenith angle **in radians**
    ! From Spitters, C. J. T. (1986), AgForMet 38: 231-242.
    implicit none
    real, intent(in) :: td             ! day(to minute fraction)
    real, intent(in) :: latdegrees     ! latitude in degrees
    real :: hour,latrad
    real :: delta    ! declination angle
    real :: pi, rad
    real, intent(out) :: cosz        ! cosz=cos(zen angle)=sin(elev angle)
    real, intent(out) :: solarelev    ! solar elevation angle (rad)
    real, intent(out) :: solarzen     ! solar zenith angle (rad)

    pi  = 3.1415926
    rad = pi / 180.0 ! Conversion from degrees to radians.
    hour = (td-floor(td))*24.0
    latrad = latdegrees*rad
    delta  = asin(-sin(rad*23.450)*cos(2.0*pi*(td+10.0)/365.0))
    cosz = sin(latrad)*sin(delta) + &
    cos(latrad)*cos(delta)*cos(rad* 15.0*(hour-12.0))
    cosz = max (cosz, 0.01)  ! Sun's angular is 0.01
    ! compute the solar elevation and zenth angles below
    solarelev = asin(cosz)/pi*180.0  !since asin(cos(zen))=pi/2-zen=elev
    solarzen = 90.0 - solarelev ! pi/2.d0 - solarelev

  end subroutine calc_solarzen


end module md_gpp
