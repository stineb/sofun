module md_ntransform
  !////////////////////////////////////////////////////////////////
  ! INORGANIC NITROGEN DYNAMICS MODULE AFTER XURI & PRENTICE 2008
  ! Contains the "main" subroutine 'ntransform' and all necessary 
  ! subroutines for handling input/output. 
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_classdefs
  use md_params_core, only: nlu, maxgrid, ndayyear

  implicit none

  private 
  public pno3, pnh4, ntransform, getpar_modl_ntransform, initglobal_ntransform, initdaily_ntransform, &
    initio_ntransform, initoutput_ntransform, getout_daily_ntransform, writeout_ascii_ntransform

  !----------------------------------------------------------------
  ! Public, module-specific state variables
  !----------------------------------------------------------------
  ! pools
  type( nitrogen ), dimension(nlu,maxgrid) :: pno3   ! soil nitrate pool [gN/m2]
  type( nitrogen ), dimension(nlu,maxgrid) :: pnh4   ! soil ammonium pool [gN/m2]

  !-----------------------------------------------------------------------
  ! Uncertain (unknown) parameters. Runtime read-in
  !-----------------------------------------------------------------------
  type params_ntransform_type
    real :: maxnitr                           ! maximum nitrification rate
    real :: non                               ! maximum NO from nitrification (day-1)
    real :: n2on                              ! maximum N2O from nitrification (day-1)
    real :: kn                                ! Michaelis-Menten coefficient [gN/m2]. Use this value if soil represents top 100 cm 
    real :: kdoc                              ! Michaelis-Menten coefficient [gC/m2]. Use this value if soil represents top 100 cm 
    real :: docmax                            ! docmax
    real :: dnitr2n2o                         ! Fraction of denitrification lost as N2O. Range of possible values: 0.002 - 0.047 (Xu-Ri and Prentice, 2008)
  end type params_ntransform_type

  type( params_ntransform_type ) :: params_ntransform

  !----------------------------------------------------------------
  ! Module-internal (private) variables
  !----------------------------------------------------------------
  real, dimension(nlu)  :: dn2o             ! soil N2O emissions [gN/m2/d]
  real, dimension(nlu)  :: dn2              ! soil N2 emissions [gN/m2/d]
  real, dimension(nlu)  :: dno              ! soil NO emissions [gN/m2/d]
  real, dimension(nlu)  :: dnloss           ! total N loss (gaseous+leaching) [gN/m2/d]
  real, dimension(nlu)  :: ddenitr          ! gross denitrification [gN/m2/d]
  real, dimension(nlu)  :: dnitr            ! gross nitrification [gN/m2/d]
  real, dimension(nlu)  :: dnvol            ! N volatilisation[gN/m2/d]
  real, dimension(nlu)  :: dnleach          ! N leaching [gN/m2/d]

  real, dimension(nlu,maxgrid), save :: no_w,  no_d      ! NO in wet and dry microsites (split done in ntransform) [gN/m2]
  real, dimension(nlu,maxgrid), save :: n2o_w, n2o_d     ! N2O in wet and dry microsites (split done in ntransform) [gN/m2]
  real, dimension(nlu,maxgrid), save :: n2_w             ! N2 in wet microsites (split done in ntransform) [gN/m2]
  
  real, dimension(nlu,maxgrid), save :: pno2              ! NO2 [gN/m2]

  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! daily
  real, allocatable, dimension(:,:,:) :: outdno3
  real, allocatable, dimension(:,:,:) :: outdnh4
  real, allocatable, dimension(:,:,:) :: outdnloss     ! daily total N loss (gaseous+leacing) (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outddenitr    ! daily amount of N denitrified (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outdnitr      ! daily amount of N nitrified (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outdnvol      ! daily amount of N volatilised (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outdnleach    ! daily amount of N leached (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outdn2o       ! daily N2O emitted (gN/m2/d)

  ! annual
  real, dimension(nlu,maxgrid) :: outano3
  real, dimension(nlu,maxgrid) :: outanh4
  real, dimension(nlu,maxgrid) :: outanloss            ! annual total N loss (gaseous+leacing) (gN/m2/yr)
  real, dimension(nlu,maxgrid) :: outadenitr           ! annual denitrified N (gN/m2/yr)
  real, dimension(nlu,maxgrid) :: outan2o              ! annual N2O emitted (gaseous+leacing) (gN/m2/yr)

contains

  subroutine ntransform( dm, mo, jpngr, dnhxdep, dnoydep, aprec )
    !////////////////////////////////////////////////////////////////
    !  Litter and SOM decomposition and nitrogen mineralisation.
    !  1st order decay of litter and SOM _pools, governed by temperature
    !  and soil moisture following LPJ (Sitch et al., 2003) and 
    !  Xu-Ri & Prentice (XXX).
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_params_core, only: pft_start, pft_end
    use md_rates
    use md_waterbal, only: soilphys, psoilphys
    use md_soiltemp, only: dtemp_soil
    use md_plant, only: pexud
    use md_interface


    ! XXX try: this is wrong: dw1 is only plant available water. 
    ! should be water-filled pore space = ( (porosity - ice) - (total fluid water volume) ) / dz

    ! arguments
    integer, intent(in) :: mo            ! month
    integer, intent(in) :: dm            ! day of the current month
    integer, intent(in) :: jpngr         ! grid cell number
    real, intent(in)    :: dnhxdep       ! daily N deposition as NHx [gN/d]
    real, intent(in)    :: dnoydep       ! daily N deposition as NOy [gN/d]
    real, intent(in)    :: aprec         ! annual total precipitation [mm/d]
    
    ! local variables
    integer    :: lu                     ! gridcell unit counter variable
    
    real, save :: ph_soil
    real, save :: nh3max
    
    real       :: dnmax                  ! labile carbon availability modifier
    real       :: ftemp_vol              ! temperature rate modifier for ammonia volatilization
    real       :: ftemp_nitr             ! temperature rate modifier for nitrification
    real       :: ftemp_denitr           ! temperature rate modifier for denitrification
    real       :: ftemp_diffus           ! temperature rate modifier for gas difussion from soil
    real       :: fph                    ! soil-pH modifier
    real       :: fwet                   ! fraction of pools in wet microsites (subject to denitrification)
    real       :: fdry                   ! fraction of pools in dry microsites (subject to nitrification)
    
    real       :: no3_inc, n2o_inc, no_inc, no2_inc, n2_inc      ! pool increments, temporary variables
    real       :: tmp                                            ! temporary variable
        
    real       :: nh4_w, no3_w, no2_w    ! anaerobic pools
    real       :: nh4_d, no3_d, no2_d    ! aerobic pools
    real       :: doc_w, no, n2o, n2     ! anaerobic pools
    
    real       :: doc_d                  ! aerobic pools

    ! Variables N balance test
    logical, parameter :: baltest_trans = .false.  ! set to false to do mass conservation test during transient simulation
    logical :: verbose = .false.  ! set to true to activate verbose mode
    logical :: baltest
    real :: nbal_before_1, nbal_after_1, nbal1, nbal_before_2, nbal_after_2, nbal2
    real :: no3bal_0, no3bal_1, nh4bal_0, nh4bal_1
    real, parameter :: eps = 9.999e-8    ! numerical imprecision allowed in mass conservation tests

    if (baltest_trans .and. .not. interface%steering%spinup) then
      baltest = .true.
      verbose = .true.
    else
      baltest = .false.
    end if

    !-------------------------------------------------------------------------
    ! Record for balances
    !-------------------------------------------------------------------------
    ! all pools plus all losses summed up
    if (verbose) print*,'              with state variables:'
    if (verbose) print*,'              ninorg = ', pno3(1,jpngr)%n14 + pnh4(1,jpngr)%n14 + no_w(1,jpngr) + no_d(1,jpngr) + n2o_w(1,jpngr) + n2o_d(1,jpngr) + n2_w(1,jpngr) + pno2(1,jpngr)
    if (verbose) print*,'              nloss  = ', dnloss(1)
    if (verbose) print*,'              dndep  = ', dnoydep + dnhxdep
    if (baltest) nbal_before_1 = pno3(1,jpngr)%n14 + pnh4(1,jpngr)%n14 + dnloss(1) + no_w(1,jpngr) + no_d(1,jpngr) + n2o_w(1,jpngr) + n2o_d(1,jpngr) + n2_w(1,jpngr) + pno2(1,jpngr) + dnoydep + dnhxdep
    if (baltest) nbal_before_2 = pno3(1,jpngr)%n14 + pnh4(1,jpngr)%n14 + ddenitr(1) + dnitr(1) + dnvol(1) + dnleach(1) + dnoydep + dnhxdep
    if (verbose) print*,'executing ntransform() ... '

    !///////////////////////////////////////////////////////////////////////
    ! INITIALIZATION 
    !-----------------------------------------------------------------------
    dnloss(:)  = 0.0
    dnvol(:)   = 0.0
    ddenitr(:) = 0.0
    dnitr(:)   = 0.0
    dnleach(:) = 0.0
    
    if ( dm==1 .and. mo==1 ) then
      !///////////////////////////////////////////////////////////////////////
      ! ANNUAL INITIALIZATION 
      !-----------------------------------------------------------------------
      ! Calculate soil PH using empirical relationship with annual precip
      ! Eq.5, Tab.5, XP08 (ntransform.cpp:65) (c++:aprec in mm/yr; F: mm/yr)
      !------------------------------------------------------------------
      ph_soil = 3810.0/(762.0+aprec)+3.8

      ! Deprotonation of NH4 to NH3 depends on soil pH
      !------------------------------------------------------------------
      if (ph_soil>6.0) then
        nh3max = 1.0
      else
        nh3max = 0.00001
      endif
      
    endif
          
    ! LOOP OVER GRIDCELL LAND UNITS
    do lu=1,nlu

      !-------------------------------------------------------------------------
      ! Add N deposition to inorganic pools
      !-------------------------------------------------------------------------
      ! pno3(lu,jpngr)%n14 = pno3(lu,jpngr)%n14 + dnoydep
      ! pnh4(lu,jpngr)%n14 = pnh4(lu,jpngr)%n14 + dnhxdep

      pno3(lu,jpngr)%n14 = pno3(lu,jpngr)%n14 + dnoydep
      pnh4(lu,jpngr)%n14 = pnh4(lu,jpngr)%n14 + dnhxdep

      ! ! xxx try: 
      ! pno3(lu,jpngr)%n14 = pno3(lu,jpngr)%n14 + 10.0 / 365.0
      ! pnh4(lu,jpngr)%n14 = pnh4(lu,jpngr)%n14 + 10.0 / 365.0


      !-------------------------------------------------------------------------
      ! Record for balances
      !-------------------------------------------------------------------------
      ! all pools plus all losses summed up
      if (verbose) print*,'              before:'
      if (verbose) print*,'              no3 = ', pno3(lu,jpngr)%n14
      if (verbose) print*,'              no4 = ', pnh4(lu,jpngr)%n14
      if (baltest) no3bal_0 = pno3(lu,jpngr)%n14
      if (baltest) nh4bal_0 = pnh4(lu,jpngr)%n14
 

      ! must rather be wtot_up which includes water below permanent wilting point (see waterbalance.F).
      !------------------------------------------------------------------
      ! reference temperature: 25°C
      ftemp_vol = min( 1.0, ftemp( dtemp_soil(lu,jpngr), "lloyd_and_taylor", ref_temp=25.0 ) )   


      !///////////////////////////////////////////////////////////////////////
      ! AMMONIUM VOLATILIZATION (ntransform.cpp:41)
      !-----------------------------------------------------------------------
      ! use mw1 for monthly timestep and wpool for daily, because this is updated daily
      ! XXX nh3max is not considered in the equations presented in the paper! XXX
      fph                = exp( 2.0 * ( ph_soil - 10.0 ) )
      dnvol(lu)          = nh3max * ftemp_vol**2 * fph * soilphys(lu)%wscal * ( 1.0 - soilphys(lu)%wscal ) * pnh4(lu,jpngr)%n14
      pnh4(lu,jpngr)%n14 = pnh4(lu,jpngr)%n14 - dnvol(lu)
      dnloss(lu)         = dnloss(lu) + dnvol(lu)

      ! if (nh4>0.0) print*,'fvol ', dnvol(lu) / nh4 

      !///////////////////////////////////////////////////////////////////////
      ! NITRATE LEACHING
      !-----------------------------------------------------------------------
      ! Reduce NO3 by fraction dnleach(lu)
      !------------------------------------------------------------------      
      dnleach(lu)        = pno3(lu,jpngr)%n14 * soilphys(lu)%fleach
      pno3(lu,jpngr)%n14 = pno3(lu,jpngr)%n14 - dnleach(lu)
      dnloss(lu)         = dnloss(lu) + dnleach(lu)


      !///////////////////////////////////////////////////////////////////////
      ! SUBSTRATE PARTITIONING (ntransform.cpp:95)
      !------------------------------------------------------------------
      ! Nitrification (aerobic) and denitrification (anaerobic) can occur
      ! simulataneously in different microsites. Substrate is thus parti-
      ! tioned according to the water content.
      
      ! wet (anaerobic) fraction
      !------------------------------------------------------------------
      ! print*,'ntransform wscal ', soilphys(lu)%wscal

      fwet  = soilphys(lu)%wscal / 3.3
      nh4_w = fwet * pnh4(lu,jpngr)%n14
      no3_w = fwet * pno3(lu,jpngr)%n14
      no2_w = fwet * pno2(lu,jpngr)

      doc_w = pexud(lu,jpngr)%c12 * fwet

      ! dry (aerobic) fraction
      !------------------------------------------------------------------
      fdry  = 1.0 - fwet
      nh4_d = fdry * pnh4(lu,jpngr)%n14
      no3_d = fdry * pno3(lu,jpngr)%n14
      no2_d = fdry * pno2(lu,jpngr)

      ! doc_d = sum( pexud(pft_start(lu):pft_end(lu),jpngr)%c12 ) * fdry
      doc_w = pexud(lu,jpngr)%c12 * fdry


      !///////////////////////////////////////////////////////////////////////
      ! NITRIFICATION in aerobic microsites (ntransform.cpp:123)
      !------------------------------------------------------------------
      ftemp_nitr = max( min( (((70.0-dtemp_soil(lu,jpngr))/(70.0-38.0))**12.0) * exp(12.0*(dtemp_soil(lu,jpngr)-38.0)/(70.0-38.0)), 1.0), 0.0)

      
      ! gross nitrification rate (Eq.1, Tab.8, XP08)
      !------------------------------------------------------------------
      no3_inc    = params_ntransform%maxnitr * ftemp_nitr * nh4_d
      dnitr(lu)  = no3_inc
      nh4_d      = nh4_d - no3_inc   
  
      ! NO from nitrification (Eq.3, Tab.8, XP08)
      !------------------------------------------------------------------
      no_inc         = params_ntransform%non * no3_inc
      no3_inc        = no3_inc - no_inc      
      no_d(lu,jpngr) = no_d(lu,jpngr) + no_inc
      
      ! N2O from nitrification (Eq.4, Tab.8, XP08)
      !------------------------------------------------------------------
      n2o_inc         = params_ntransform%n2on * no3_inc
      no3_inc         = no3_inc - n2o_inc
      n2o_d(lu,jpngr) = n2o_d(lu,jpngr) + n2o_inc
      no3_d           = no3_d + no3_inc
            
      ! if N loss is defined w.r.t. reduction in NH4 and NO3 pools, then this is the correct formulation:
      dnloss(lu) = dnloss(lu) + n2o_inc + no_inc

      ! xxx debug
      if (baltest) no3bal_1 = no3_w + no3_d - no3_inc
      if (baltest) nh4bal_1 = nh4_w + nh4_d + dnitr(lu)

      if (baltest) nbal1 = no3bal_1 - no3bal_0
      if (baltest) nbal2 = nh4bal_1 - nh4bal_0
      if (verbose) print*,'              --- preliminary balance after nitrification '
      if (verbose) print*,'              ', nbal1
      if (verbose) print*,'              ', nbal2
      if (baltest .and. abs(nbal1)>eps) stop 'balance 1 not satisfied'
      if (baltest .and. abs(nbal2)>eps) stop 'balance 2 not satisfied'


      !///////////////////////////////////////////////////////////////////////
      ! DENITRIFICATION (ntransform.cpp:177) in anaerobic microsites
      !------------------------------------------------------------------
      ! reference temperature: 22°C
      ftemp_denitr = ftemp( dtemp_soil(lu,jpngr), "lloyd_and_taylor", ref_temp=22.0 )

      
      ! Effect of labile carbon availability on denitrification (Eq.2, Tab.9, XP08)
      ! doc is last year's doc because it is only available at the end of the month
      ! while this SR is calculated daily, even when _dailymode==0.
      !------------------------------------------------------------------
      dnmax = params_ntransform%docmax * doc_w / ( params_ntransform%kdoc + doc_w )                     ! dnmax < 1 for all doc_w 

      ! ! xxx try:
      ! dnmax = 0.5

      ! print*,'fMM DOC ',  params_ntransform%docmax * doc_w / ( params_ntransform%kdoc + doc_w )
      
      ! Denitrification ratio, NO3->NO2 (Eq.3, Tab.9, XP08)
      !------------------------------------------------------------------
      no2_inc     = min( dnmax * ftemp_denitr * no3_w / ( params_ntransform%kn + no3_w ) * 1000.0, no3_w )
      if (no2_inc>no3_w) stop 'no2_inc > no3_w'
      
      no3_w       = no3_w - no2_inc
      no2_w       = no2_w + no2_inc
      ddenitr(lu) = no2_inc
      
      ! if N loss is defined w.r.t. reduction in NH4 and NO3 pools, then this is the correct formulation:
      dnloss(lu) = dnloss(lu) + no2_inc


      ! Transformation NO2->N2 (Eq.4., Tab.9, XP08)
      !------------------------------------------------------------------
      n2_inc = min( dnmax * ftemp_denitr * no2_w / ( params_ntransform%kn + no2_w ) * 1000.0, no2_w )
      if (n2_inc>no2_w) stop 'n2_inc > no2_w'

      no2_w = no2_w - n2_inc
      
      ! N2O from denitrification (Eq.6, Tab.9, XP08)
      !------------------------------------------------------------------
      ! n2o_inc = 0.018d0*ftemp_denitr*(1.01d0-0.21d0*soilphys(lu)%wscal)*n2_inc  !Colin says 0.018 was used here. Code I got had 0.015
      ! Factor reduced from 1.8% to 1.2% to get ~6.5 TgN/yr N2O emissions
      ! n2o_inc = dnitr2n2o*ftemp_denitr*(1.01-0.21*soilphys(lu)%wscal)*n2_inc
      ! XXX try: Changed this to from 0.21 to 1.0 in order to get plausible results
      n2o_inc = params_ntransform%dnitr2n2o * ftemp_denitr * ( 1.01 - 0.8 * soilphys(lu)%wscal ) * n2_inc
      n2_inc  = n2_inc - n2o_inc
      n2o_w(lu,jpngr) = n2o_w(lu,jpngr) + n2o_inc

      ! NO from denitrification (Eq.5, Tab.9, XP08)
      !------------------------------------------------------------------
      ! no_inc = 0.0001*ftemp_denitr*(1.01-0.21*soilphys(lu)%wscal)*n2_inc
      ! XXX try: Changed this to from 0.21 to 1.0 in order to get plausible results
      no_inc = 0.0001 * ftemp_denitr * ( 1.01 - 0.8 * soilphys(lu)%wscal ) * n2_inc
      n2_inc = n2_inc - no_inc
      no_w(lu,jpngr) = no_w(lu,jpngr) + no_inc

      ! N2 from denitrification
      !------------------------------------------------------------------
      n2_w(lu,jpngr) = n2_w(lu,jpngr) + n2_inc


      !///////////////////////////////////////////////////////////////////////
      ! UPDATE POOLS (ntransform.cpp:389)
      !------------------------------------------------------------------
      ! nh4, no3 and no2 (sum of sub-_pools in wet and dry microsites) are
      ! defined as global variables. They are split into sub-_pools at the
      ! beginning of ntransform (see "substrate partitioning"), processed
      ! through "denitrification" and "nitrification" and added up again
      ! here. In contrast, each sub-pool for wet and dry conditions of no,
      ! n2o and n2 is defined as a global (common) variable, while the sum
      ! of the sup-_pools (no, n2o and n2) is defined locally and only used
      ! for the diffusion/emission (see below).
      !------------------------------------------------------------------
      pnh4(lu,jpngr)%n14 = nh4_w + nh4_d
      pno3(lu,jpngr)%n14 = no3_w + no3_d
      pno2(lu,jpngr)     = no2_w + no2_d

      no  = no_w(lu,jpngr) + no_d(lu,jpngr)
      n2o = n2o_w(lu,jpngr) + n2o_d(lu,jpngr)
      n2  = n2_w(lu,jpngr)

      !///////////////////////////////////////////////////////////////////////
      ! Diffusion of NO, N2O and N2 from the soil (ntransform.cpp:281)
      !------------------------------------------------------------------
      ! reference temperature: 25°C. Corresponds to Eq.1, Tab.10, Xu-Ri & Prentice, 2008
      ftemp_diffus = min( 1.0, ftemp( dtemp_soil(lu,jpngr), "lloyd_and_taylor", ref_temp=25.0 ))

      ! Total gaseous escape
      !------------------------------------------------------------------
      dno(lu) = ftemp_diffus*(1.0-soilphys(lu)%wscal)*no
      dn2o(lu)= ftemp_diffus*(1.0-soilphys(lu)%wscal)*n2o
      dn2(lu) = ftemp_diffus*(1.0-soilphys(lu)%wscal)*n2

      ! if N loss is defined w.r.t. gaseous escape, then this is the correct formulation:
      ! dnloss(lu) = dnloss(lu) + dno(lu) + dn2o(lu) + dn2(lu)

      ! Gaseous escape of pools at dry microsites
      !------------------------------------------------------------------
      tmp = ftemp_diffus*(1.0-soilphys(lu)%wscal)*no_d(lu,jpngr)
      no_d(lu,jpngr) = no_d(lu,jpngr)-tmp
      
      tmp = ftemp_diffus*(1.0-soilphys(lu)%wscal)*n2o_d(lu,jpngr)
      n2o_d(lu,jpngr) = n2o_d(lu,jpngr)-tmp

      ! Gaseous escape of pools at wet microsites
      !------------------------------------------------------------------
      tmp = ftemp_diffus*(1.0-soilphys(lu)%wscal)*no_w(lu,jpngr)
      no_w(lu,jpngr) = no_w(lu,jpngr)-tmp
      
      tmp = ftemp_diffus*(1.0-soilphys(lu)%wscal)*n2o_w(lu,jpngr)
      n2o_w(lu,jpngr) = n2o_w(lu,jpngr)-tmp
                 
      n2_w(lu,jpngr) = n2_w(lu,jpngr) - dn2(lu)
      
    enddo                                                 ! lu

    !-------------------------------------------------------------------------
    ! Test mass conservation
    !-------------------------------------------------------------------------
    ! all pools plus all losses summed up
    if (baltest) nbal_after_1 = pno3(1,jpngr)%n14 + pnh4(1,jpngr)%n14 + dnloss(1) + no_w(1,jpngr) + no_d(1,jpngr) + n2o_w(1,jpngr) + n2o_d(1,jpngr) + n2_w(1,jpngr) + pno2(1,jpngr)
    if (baltest) nbal_after_2 = pno3(1,jpngr)%n14 + pnh4(1,jpngr)%n14 + ddenitr(1) + dnitr(1) + dnvol(1) + dnleach(1) - no3_inc
    if (baltest) nbal1 = nbal_after_1 - nbal_before_1
    if (baltest) nbal2 = nbal_after_2 - nbal_before_2
    if (verbose) print*,'              ==> returned:'
    if (verbose) print*,'              ninorg = ', pno3(1,jpngr)%n14 + pnh4(1,jpngr)%n14 + no_w(1,jpngr) + no_d(1,jpngr) + n2o_w(1,jpngr) + n2o_d(1,jpngr) + n2_w(1,jpngr) + pno2(1,jpngr)
    if (verbose) print*,'              nloss  = ', dnloss(1)
    if (verbose) print*,'   --- balance: '
    if (verbose) print*,'       d( ninorg + loss )', nbal1
    if (verbose) print*,'       d( ninorg + loss )', nbal2
    if (baltest .and. abs(nbal1)>eps) stop 'balance 1 not satisfied'
    if (baltest .and. abs(nbal2)>eps) stop 'balance 2 not satisfied'

  end subroutine ntransform

  
  subroutine getpar_modl_ntransform()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads waterbalance module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! maximum nitrification rate
    params_ntransform%maxnitr = getparreal( 'params/params_ntransform_xuri.dat', 'maxnitr' )

    ! maximum NO from nitrification (day-1)
    params_ntransform%non = getparreal( 'params/params_ntransform_xuri.dat', 'non' )

    ! maximum N2O from nitrification (day-1)
    params_ntransform%n2on = getparreal( 'params/params_ntransform_xuri.dat', 'n2on' )

    ! Michaelis-Menten coefficient [gN/m2]. Use this value if soil represents top 100 cm 
    params_ntransform%kn = getparreal( 'params/params_ntransform_xuri.dat', 'kn' )

    ! Michaelis-Menten coefficient [gC/m2]. Use this value if soil represents top 100 cm 
    params_ntransform%kdoc = getparreal( 'params/params_ntransform_xuri.dat', 'kdoc' )

    ! docmax
    params_ntransform%docmax = getparreal( 'params/params_ntransform_xuri.dat', 'docmax' )

    ! Fraction of denitrification lost as N2O. Range of possible values: 0.002 - 0.047 (Xu-Ri and Prentice, 2008)
    params_ntransform%dnitr2n2o = getparreal( 'params/params_ntransform_xuri.dat', 'dnitr2n2o' )

  end subroutine getpar_modl_ntransform

  
  subroutine initglobal_ntransform()
    !////////////////////////////////////////////////////////////////
    ! Subroutine initialises pool variables
    !----------------------------------------------------------------
    ! public variables
    pnh4(:,:)  = nitrogen(10.0)  ! start from non-zero to allow growth
    pno3(:,:)  = nitrogen(10.0)  ! start from non-zero to allow growth

    ! module-specific variables
    pno2(:,:)      = 0.0
    no_w(:,:)     = 0.0
    no_d(:,:)     = 0.0
    n2o_w(:,:)    = 0.0
    n2o_d(:,:)    = 0.0
    n2_w(:,:)     = 0.0

  end subroutine initglobal_ntransform


  subroutine initdaily_ntransform()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    !----------------------------------------------------------------
    use md_interface

    dn2o   (:) = 0.0
    dn2    (:) = 0.0
    dno    (:) = 0.0
    dnloss (:) = 0.0
    ddenitr(:) = 0.0
    dnitr  (:) = 0.0
    dnvol  (:) = 0.0
    dnleach(:) = 0.0

  end subroutine initdaily_ntransform


  subroutine initio_ntransform()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use md_interface

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(interface%params_siml%runname)

    !----------------------------------------------------------------
    ! DAILY OUTPUT
    !----------------------------------------------------------------
    if (interface%params_siml%loutntransform) then

      ! SOIL NO3
      filnam=trim(prefix)//'.d.no3.out'
      open(107,file=filnam,err=888,status='unknown')

      ! SOIL NH4
      filnam=trim(prefix)//'.d.nh4.out'
      open(506,file=filnam,err=888,status='unknown')

      ! DAILY TOTAL N LOSS (gN/m2/d)
      filnam=trim(prefix)//'.d.nloss.out'
      open(500,file=filnam,err=888,status='unknown')

      ! DAILY TOTAL N LOSS (gN/m2/d)
      filnam=trim(prefix)//'.d.nvol.out'
      open(501,file=filnam,err=888,status='unknown')

      ! DAILY TOTAL N LOSS (gN/m2/d)
      filnam=trim(prefix)//'.d.denitr.out'
      open(502,file=filnam,err=888,status='unknown')

      ! DAILY TOTAL N LOSS (gN/m2/d)
      filnam=trim(prefix)//'.d.nitr.out'
      open(503,file=filnam,err=888,status='unknown')

      ! DAILY TOTAL N LOSS (gN/m2/d)
      filnam=trim(prefix)//'.d.nleach.out'
      open(504,file=filnam,err=888,status='unknown')

      ! DAILY N2O EMISSIONS (gN/m2/d)
      filnam=trim(prefix)//'.d.n2o.out'
      open(505,file=filnam,err=888,status='unknown')

      !----------------------------------------------------------------
      ! ANNUAL OUTPUT
      !----------------------------------------------------------------
      ! ANNUAL TOTAL N LOSS (gN/m2/yr)
      filnam=trim(prefix)//'.a.nloss.out'
      open(550,file=filnam,err=888,status='unknown')

      ! ANNUAL TOTAL DENITRIFIED N (gN/m2/yr)
      filnam=trim(prefix)//'.a.denitr.out'
      open(553,file=filnam,err=888,status='unknown')

      ! ANNUAL N2O EMISSIONS (gN/m2/yr)
      filnam=trim(prefix)//'.a.n2o.out'
      open(551,file=filnam,err=888,status='unknown')

      ! SOIL NO3 (mean over days)
      filnam=trim(prefix)//'.a.no3.out'
      open(316,file=filnam,err=888,status='unknown')

      ! SOIL NH4 (mean over days)
      filnam=trim(prefix)//'.a.nh4.out'
      open(552,file=filnam,err=888,status='unknown')

    end if

    return

    888  stop 'INITIO_NTRANSFORM: error opening output files'

  end subroutine initio_ntransform


  subroutine initoutput_ntransform()
    !////////////////////////////////////////////////////////////////
    !  Initialises waterbalance-specific output variables
    !----------------------------------------------------------------
    use md_interface

    if (interface%params_siml%loutntransform) then

      if (interface%steering%init) allocate( outdnloss ( nlu,ndayyear,maxgrid ) ) ! daily total N loss (gaseous+leacing) (gN/m2/d)
      if (interface%steering%init) allocate( outddenitr( nlu,ndayyear,maxgrid ) ) ! daily amount of N denitrified (gN/m2/d)
      if (interface%steering%init) allocate( outdnitr  ( nlu,ndayyear,maxgrid ) ) ! daily amount of N nitrified (gN/m2/d)
      if (interface%steering%init) allocate( outdnvol  ( nlu,ndayyear,maxgrid ) ) ! daily amount of N volatilised (gN/m2/d)
      if (interface%steering%init) allocate( outdnleach( nlu,ndayyear,maxgrid ) ) ! daily amount of N leached (gN/m2/d)
      if (interface%steering%init) allocate( outdn2o   ( nlu,ndayyear,maxgrid ) ) ! daily N2O emitted (gN/m2/d)
      if (interface%steering%init) allocate( outdno3   ( nlu,ndayyear,maxgrid ) ) ! daily total inorganic N (gN/m2)
      if (interface%steering%init) allocate( outdnh4   ( nlu,ndayyear,maxgrid ) ) ! daily total inorganic N (gN/m2)

      outdnloss(:,:,:)  = 0.0
      outddenitr(:,:,:) = 0.0
      outdnitr(:,:,:)   = 0.0
      outdnvol(:,:,:)   = 0.0
      outdnleach(:,:,:) = 0.0
      outdn2o(:,:,:)    = 0.0
      outdno3(:,:,:)    = 0.0
      outdnh4(:,:,:)    = 0.0
      
      outanloss(:,:)  = 0.0
      outadenitr(:,:) = 0.0
      outan2o(:,:)    = 0.0
      outano3(:,:)    = 0.0
      outanh4(:,:)    = 0.0

    end if

  end subroutine initoutput_ntransform


  subroutine getout_daily_ntransform( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy    
    integer, intent(in) :: doy    

    if (interface%params_siml%loutntransform) then
      !----------------------------------------------------------------
      ! DAILY
      ! Collect daily output variables
      !----------------------------------------------------------------
      outdnloss(:,doy,jpngr)  = dnloss(:)
      outddenitr(:,doy,jpngr) = ddenitr(:)
      outdnitr(:,doy,jpngr)   = dnitr(:)
      outdnvol(:,doy,jpngr)   = dnvol(:)
      outdnleach(:,doy,jpngr) = dnleach(:)
      outdn2o(:,doy,jpngr)    = dn2o(:)
      outdno3(:,doy,jpngr)    = pno3(:,jpngr)%n14
      outdnh4(:,doy,jpngr)    = pnh4(:,jpngr)%n14

      !----------------------------------------------------------------
      ! ANNUAL SUM OVER DAILY VALUES
      ! Collect annual output variables as sum of daily values
      !----------------------------------------------------------------
      outano3(:,jpngr)    = outano3(:,jpngr) + pno3(:,jpngr)%n14 / ndayyear
      outanh4(:,jpngr)    = outanh4(:,jpngr) + pnh4(:,jpngr)%n14 / ndayyear
      outanloss(:,jpngr)  = outanloss(:,jpngr) + dnloss(:)
      outadenitr(:,jpngr) = outadenitr(:,jpngr) + ddenitr(:)
      outan2o(:,jpngr)    = outan2o(:,jpngr) + dn2o(:)

    end if

  end subroutine getout_daily_ntransform


  subroutine writeout_ascii_ntransform( year )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use md_params_core, only: npft
    use md_interface

    ! arguments
    integer, intent(in) :: year       ! simulation year

    ! local variables
    real :: itime
    integer :: day, jpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii_ntransform: think of something ...'
    jpngr = 1

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    !-------------------------------------------------------------------------
    if (interface%params_siml%loutntransform) then

      if ( .not. interface%steering%spinup &
        .and. interface%steering%outyear>=interface%params_siml%daily_out_startyr &
        .and. interface%steering%outyear<=interface%params_siml%daily_out_endyr ) then

        ! Write daily output only during transient simulation
        do day=1,ndayyear

          ! Define 'itime' as a decimal number corresponding to day in the year + year
          itime = real(interface%steering%outyear) + real(day-1)/real(ndayyear)

          if (nlu>1) stop 'writeout_ascii_ntransform: write out lu-area weighted sum'
          if (npft>1) stop 'writeout_ascii_ntransform: think of something for ccost output'

          write(506,999) itime, sum(outdnh4(:,day,jpngr))
          write(107,999) itime, sum(outdno3(:,day,jpngr))
          write(500,999) itime, sum(outdnloss(:,day,jpngr))
          write(501,999) itime, sum(outdnvol(:,day,jpngr))
          write(502,999) itime, sum(outddenitr(:,day,jpngr))
          write(503,999) itime, sum(outdnitr(:,day,jpngr))
          write(504,999) itime, sum(outdnleach(:,day,jpngr))
          write(505,999) itime, sum(outdn2o(:,day,jpngr))

        end do
      end if

      !-------------------------------------------------------------------------
      ! ANNUAL OUTPUT
      ! Write annual value, summed over all PFTs / LUs
      ! xxx implement taking sum over PFTs (and gridcells) in this land use category
      !-------------------------------------------------------------------------
      itime = real(interface%steering%outyear)

      write(316,999) itime, sum(outano3(:,jpngr))
      write(552,999) itime, sum(outanh4(:,jpngr))
      write(550,999) itime, sum(outanloss(:,jpngr))
      write(553,999) itime, sum(outadenitr(:,jpngr))
      write(551,999) itime, sum(outan2o(:,jpngr))

    end if

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_ntransform


end module md_ntransform
