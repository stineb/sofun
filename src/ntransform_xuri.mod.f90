module _ntransform
  !////////////////////////////////////////////////////////////////
  ! INORGANIC NITROGEN DYNAMICS MODULE AFTER XURI & PRENTICE 2008
  ! Contains the "main" subroutine 'ntransform' and all necessary 
  ! subroutines for handling input/output. 
  ! Every module that implements 'ntransform' must contain this list 
  ! of subroutines (names that way).
  !   - getpar_modl_ntransform
  !   - initio_ntransform
  !   - initoutput_ntransform
  !   - getout_daily_ntransform
  !   - getout_monthly_ntransform
  !   - writeout_ascii_ntransform
  !   - ntransform
  ! Required module-independent model state variables (necessarily 
  ! updated by 'ntransform') are:
  !   - inorganic N pool, sum of NO3 and NH4 ('pninorg')
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _params_core, only: nlu, maxgrid, ndayyear
  use _params_siml, only: loutntransform

  implicit none

  ! NTRANSFORM PARAMETERS
  ! read from file, therefore not declared as parameter
  real :: MAXNITR                           ! maximum nitrification rate
  real :: NON                               ! maximum NO from nitrification (day-1)
  real :: N2ON                              ! maximum N2O from nitrification (day-1)
  real :: KN                                ! Michaelis-Menten coefficient [gN/m2]. Use this value if soil represents top 100 cm 
  real :: KDOC                              ! Michaelis-Menten coefficient [gC/m2]. Use this value if soil represents top 100 cm 
  real :: DOCMAX                            ! DOCMAX
  real :: DNITR2N2O                         ! Fraction of denitrification lost as N2O. Range of possible values: 0.002 - 0.047 (Xu-Ri and Prentice, 2008)

  ! MODULE-SPECIFIC VARIABLES
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
  
  real, dimension(nlu,maxgrid), save :: no2              ! NO2 [gN/m2]
  real, dimension(nlu,maxgrid), save :: fno3             ! fraction: no3/(no3+nh4)

  ! OUTPUT VARIABLES
  ! daily
  real, allocatable, dimension(:,:,:) :: outdnloss     ! daily total N loss (gaseous+leacing) (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outddenitr    ! daily amount of N denitrified (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outdnitr      ! daily amount of N nitrified (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outdnvol      ! daily amount of N volatilised (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outdnleach    ! daily amount of N leached (gN/m2/d)
  real, allocatable, dimension(:,:,:) :: outdn2o       ! daily N2O emitted (gN/m2/d)

  ! annual
  real, dimension(nlu,maxgrid) :: outanloss              ! annual total N loss (gaseous+leacing) (gN/m2/yr)
  real, dimension(nlu,maxgrid) :: outan2o                ! annual N2O emitted (gaseous+leacing) (gN/m2/yr)

contains

  subroutine ntransform( dm, mo, jpngr, dndep, aprec )
    !////////////////////////////////////////////////////////////////
    !  Litter and SOM decomposition and nitrogen mineralisation.
    !  1st order decay of litter and SOM _pools, governed by temperature
    !  and soil moisture following LPJ (Sitch et al., 2003) and 
    !  Xu-Ri & Prentice (XXX).
    !  June 2014
    !  b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use _rates
    use _vars_core, only: pninorg, psoilphys, soilphys, ddoc

    ! XXX beni: this is wrong: dw1 is only plant available water. 
    ! should be water-filled pore space = ( (porosity - ice) - (total fluid water volume) ) / dz

    ! arguments
    integer, intent(in) :: mo                              ! month
    integer, intent(in) :: dm                              ! day of the current month
    integer, intent(in) :: jpngr                           ! grid cell number
    real, intent(in)    :: dndep                           ! daily N deposition [gN/d]
    real, intent(in)    :: aprec                           ! annual total precipitation [mm/d]
    
    ! local variables
    integer       :: lu                                             ! gridcell unit counter variable
    
    real, save    :: ph_soil
    real, save    :: nh3max
    
    real          :: dnmax                                          ! labile carbon availability modifier
    real          :: ftemp_ninorg                                   ! temperature modifier
    
    real          :: no3_inc, n2o_inc, no_inc, no2_inc, n2_inc      ! temporary variables
    real          :: tmp                                            ! temporary variable
    
    real          :: nh4      ! ammonium [gN/m2]
    real          :: no3      ! nitrate [gN/m2]
    
    real          :: nh4_w, no3_w, no2_w                            ! anaerobic _pools
    real          :: nh4_d, no3_d, no2_d                            ! aerobic _pools
    real          :: doc_w, no, n2o, n2                             ! anaerobic _pools
    
    real          :: doc_d                                          ! aerobic _pools

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
      ! Add N deposition to inorganic pool
      !-------------------------------------------------------------------------
      ! write(0,*) 'adding N deposition ', dndep(doy)
      pninorg(lu,jpngr)%n14 = pninorg(lu,jpngr)%n14 + dndep
              
      !-------------------------------------------------------------------------
      ! Define NO3 and NH4 from total inorganic N and previous day's shares
      !-------------------------------------------------------------------------
      no3 = pninorg(lu,jpngr)%n14 * fno3(lu,jpngr)
      nh4 = pninorg(lu,jpngr)%n14 * (1.0 - fno3(lu,jpngr))

      ! Daily updated soil moisture and soil temperature are required for
      ! ntransform.

      ! must rather be wtot_up which includes water below permanent wilting point (see waterbalance.F).
      !------------------------------------------------------------------
      ! reference temperature: 25°C
      ftemp_ninorg = min( 1.0, ftemp( psoilphys(lu,jpngr)%temp, "lloyd_and_taylor", ref_temp=25.0 ) )   


      !///////////////////////////////////////////////////////////////////////
      ! AMMONIUM VOLATILIZATION (ntransform.cpp:41)
      !-----------------------------------------------------------------------
      ! use mw1 for monthly timestep and wpool for daily, because this is updated daily
      ! XXX nh3max is not considered in the equations presented in the paper! XXX
      dnvol(lu)  = ( nh3max * ( soilphys(lu)%wscal * ( 1.0 - soilphys(lu)%wscal ) ) * ftemp_ninorg * ftemp_ninorg * exp( 2.0 * ( ph_soil - 10.0 ) ) ) * nh4
      nh4        = nh4 - dnvol(lu)
      dnloss(lu) = dnloss(lu) + dnvol(lu)

      
      !///////////////////////////////////////////////////////////////////////
      ! NITRATE LEACHING
      !-----------------------------------------------------------------------
      ! Reduce NO3 by fraction dnleach(lu)
      !------------------------------------------------------------------
      dnleach(lu) = no3 * soilphys(lu)%ro / psoilphys(lu,jpngr)%wcont
      no3         = no3 - dnleach(lu)
      dnloss(lu)  = dnloss(lu) + dnleach(lu)

      !///////////////////////////////////////////////////////////////////////
      ! SUBSTRATE PARTITIONING (ntransform.cpp:95)
      !------------------------------------------------------------------
      ! Nitrification (aerobic) and denitrification (anaerobic) can occur
      ! simulataneously in different microsites. Substrate is thus parti-
      ! tioned according to the water content.
      
      ! wet (anaerobic) fraction
      !------------------------------------------------------------------
      nh4_w = soilphys(lu)%wscal / 3.3 * nh4
      no3_w = soilphys(lu)%wscal / 3.3 * no3
      no2_w = soilphys(lu)%wscal / 3.3 * no2(lu,jpngr)

      doc_w = ddoc(lu) * soilphys(lu)%wscal / 3.3
      
      ! write(0,*) 'mo, dm, ddoc(lu) ', mo, dm, ddoc(lu)

      ! dry (aerobic) fraction
      !------------------------------------------------------------------
      nh4_d = ( 1.0 - soilphys(lu)%wscal / 3.3 ) * nh4
      no3_d = ( 1.0 - soilphys(lu)%wscal / 3.3 ) * no3
      no2_d = ( 1.0 - soilphys(lu)%wscal / 3.3 ) * no2(lu,jpngr)

      doc_d = ddoc(lu) * ( 1.0 - soilphys(lu)%wscal / 3.3 )

    
      !///////////////////////////////////////////////////////////////////////
      ! NITRIFICATION in aerobic microsites (ntransform.cpp:123)
      !------------------------------------------------------------------
      ftemp_ninorg = max( min( (((70.0-psoilphys(lu,jpngr)%temp)/(70.0-38.0))**12.0) * exp(12.0*(psoilphys(lu,jpngr)%temp-38.0)/(70.0-38.0)), 1.0), 0.0)

      
      ! gross nitrification rate (Eq.1, Tab.8, XP08)
      !------------------------------------------------------------------
      no3_inc    = MAXNITR * ftemp_ninorg * nh4_d
      dnitr(lu)  = no3_inc
      nh4_d      = nh4_d - no3_inc   
      dnloss(lu) = dnloss(lu) + no3_inc

      
      ! NO from nitrification (Eq.3, Tab.8, XP08)
      !------------------------------------------------------------------
      no_inc         = NON * no3_inc
      no3_inc        = no3_inc - no_inc
      
      no_d(lu,jpngr) = no_d(lu,jpngr) + no_inc

      
      ! N2O from nitrification (Eq.4, Tab.8, XP08)
      !------------------------------------------------------------------
      n2o_inc         = N2ON * no3_inc
      no3_inc         = no3_inc - n2o_inc
      
      n2o_d(lu,jpngr) = n2o_d(lu,jpngr) + n2o_inc
      no3_d           = no3_d + no3_inc
      
      dn2o(lu) = n2o_inc
      
      !///////////////////////////////////////////////////////////////////////
      ! DENITRIFICATION (ntransform.cpp:177) in anaerobic microsites
      !------------------------------------------------------------------
      ! reference temperature: 22°C
      ftemp_ninorg = ftemp( psoilphys(lu,jpngr)%temp, "lloyd_and_taylor", ref_temp=22.0 )

      
      ! Effect of labile carbon availability on denitrification (Eq.2, Tab.9, XP08)
      ! doc is last year's doc because it is only available at the end of the month
      ! while this SR is calculated daily, even when _dailymode==0.
      !------------------------------------------------------------------
      dnmax = DOCMAX * doc_w / ( KDOC + doc_w )                     ! dnmax < 1 for all doc_w 

      
      ! Denitrification ratio, NO3->NO2 (Eq.3, Tab.9, XP08)
      !------------------------------------------------------------------
      no2_inc     = min( dnmax * ftemp_ninorg * no3_w / ( KN + no3_w ) * 1000.0, no3_w )
      if (no2_inc>no3_w) stop 'no2_inc > no3_w'
      
      no3_w       = no3_w - no2_inc
      no2_w       = no2_w + no2_inc
      ddenitr(lu) = no2_inc
      dnloss(lu)  = dnloss(lu) + no2_inc

      
      ! Transformation NO2->N2 (Eq.4., Tab.9, XP08)
      !------------------------------------------------------------------
      n2_inc = min( dnmax * ftemp_ninorg * no2_w / ( KN + no2_w ) * 1000.0, no2_w )
      if (n2_inc>no2_w) stop 'n2_inc > no2_w'

      no2_w = no2_w - n2_inc

      
      ! N2O from denitrification (Eq.6, Tab.9, XP08)
      !------------------------------------------------------------------
      ! n2o_inc = 0.018d0*ftemp_ninorg*(1.01d0-0.21d0*soilphys(lu)%wscal)*n2_inc  !Colin says 0.018 was used here. Code I got had 0.015
      ! Factor reduced from 1.8% to 1.2% to get ~6.5 TgN/yr N2O emissions
      ! n2o_inc = DNITR2N2O*ftemp_ninorg*(1.01-0.21*soilphys(lu)%wscal)*n2_inc
      ! XXX try: Changed this to from 0.21 to 1.0 in order to get plausible results
      n2o_inc = DNITR2N2O * ftemp_ninorg * ( 1.01 - 0.8 * soilphys(lu)%wscal ) * n2_inc
      n2_inc  = n2_inc - n2o_inc
      
      n2o_w(lu,jpngr) = n2o_w(lu,jpngr) + n2o_inc
      dn2o(lu) = dn2o(lu) + n2o_inc


      ! NO from denitrification (Eq.5, Tab.9, XP08)
      !------------------------------------------------------------------
      ! no_inc = 0.0001*ftemp_ninorg*(1.01-0.21*soilphys(lu)%wscal)*n2_inc
      ! XXX try: Changed this to from 0.21 to 1.0 in order to get plausible results
      no_inc = 0.0001 * ftemp_ninorg * ( 1.01 - 0.8 * soilphys(lu)%wscal ) * n2_inc
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
      nh4 = nh4_w + nh4_d
      no3 = no3_w + no3_d
      no2 = no2_w + no2_d

      pninorg(lu,jpngr)%n14 = nh4 + no3
      if ( (nh4+no3)>0.0 ) then
        fno3(lu,jpngr) = no3 / (nh4 + no3)
      else
        fno3(lu,jpngr) = 0.0
      end if

      no  = no_w(lu,jpngr) + no_d(lu,jpngr)
      n2o = n2o_w(lu,jpngr) + n2o_d(lu,jpngr)
      n2  = n2_w(lu,jpngr)

              
      !///////////////////////////////////////////////////////////////////////
      ! Diffusion of NO, N2O and N2 from the soil (ntransform.cpp:281)
      !------------------------------------------------------------------
      ! reference temperature: 25°C. Corresponds to Eq.1, Tab.10, Xu-Ri & Prentice, 2008
      ftemp_ninorg = min( 1.0, ftemp( psoilphys(lu,jpngr)%temp, "lloyd_and_taylor", ref_temp=25.0 ))

      ! Nitrification _fluxes
      !------------------------------------------------------------------
      tmp = ftemp_ninorg*(1.0-soilphys(lu)%wscal)*no_d(lu,jpngr)
      no_d(lu,jpngr) = no_d(lu,jpngr)-tmp
      
      tmp = ftemp_ninorg*(1.0-soilphys(lu)%wscal)*n2o_d(lu,jpngr)
      n2o_d(lu,jpngr) = n2o_d(lu,jpngr)-tmp

      
      ! Denitrification _fluxes
      !------------------------------------------------------------------
      tmp = ftemp_ninorg*(1.0-soilphys(lu)%wscal)*no_w(lu,jpngr)
      no_w(lu,jpngr) = no_w(lu,jpngr)-tmp
      
      tmp = ftemp_ninorg*(1.0-soilphys(lu)%wscal)*n2o_w(lu,jpngr)
      n2o_w(lu,jpngr) = n2o_w(lu,jpngr)-tmp
                 

      ! xxx try: soilphys(lu)%wscal was always too high therefore no n2o escaped
      ! ! Total _fluxes (XXX should be equal to sum of nitrification and
      ! ! denitrification _fluxes. XXX)
      ! !------------------------------------------------------------------
      ! dno(lu) = ftemp_ninorg*(1.0-soilphys(lu)%wscal)*no
      ! dn2o(lu)= ftemp_ninorg*(1.0-soilphys(lu)%wscal)*n2o
      ! dn2(lu) = ftemp_ninorg*(1.0-soilphys(lu)%wscal)*n2

      ! write(0,*) 'n2o ', n2o

      n2_w(lu,jpngr) = n2_w(lu,jpngr)-dn2(lu)
      
    enddo                                                 ! lu

  end subroutine ntransform

  
  subroutine getpar_modl_ntransform()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads waterbalance module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use _sofunutils, only: getparreal

    ! maximum nitrification rate
    MAXNITR = getparreal( 'params/params_ntransform_xuri.dat', 'MAXNITR' )

    ! maximum NO from nitrification (day-1)
    NON = getparreal( 'params/params_ntransform_xuri.dat', 'NON' )

    ! maximum N2O from nitrification (day-1)
    N2ON = getparreal( 'params/params_ntransform_xuri.dat', 'N2ON' )

    ! Michaelis-Menten coefficient [gN/m2]. Use this value if soil represents top 100 cm 
    KN = getparreal( 'params/params_ntransform_xuri.dat', 'KN' )

    ! Michaelis-Menten coefficient [gC/m2]. Use this value if soil represents top 100 cm 
    KDOC = getparreal( 'params/params_ntransform_xuri.dat', 'KDOC' )

    ! DOCMAX
    DOCMAX = getparreal( 'params/params_ntransform_xuri.dat', 'DOCMAX' )

    ! Fraction of denitrification lost as N2O. Range of possible values: 0.002 - 0.047 (Xu-Ri and Prentice, 2008)
    DNITR2N2O = getparreal( 'params/params_ntransform_xuri.dat', 'DNITR2N2O' )

  end subroutine getpar_modl_ntransform

  
  subroutine init_global_ntransform()
    !////////////////////////////////////////////////////////////////
    ! Subroutine initialises pool variables, specific for this module
    !----------------------------------------------------------------
    no2(:,:)      = 1.0  ! xxx try
    fno3(:,:)     = 0.0
    no_w(:,:)     = 0.0
    no_d(:,:)     = 0.0
    n2o_w(:,:)    = 0.0
    n2o_d(:,:)    = 0.0
    n2_w(:,:)     = 0.0

  end subroutine init_global_ntransform


  subroutine initio_ntransform()
    !////////////////////////////////////////////////////////////////
    ! OPEN ASCII OUTPUT FILES FOR OUTPUT
    !----------------------------------------------------------------
    use _params_siml, only: runname

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(runname)

    !----------------------------------------------------------------
    ! DAILY OUTPUT
    !----------------------------------------------------------------
    if (loutntransform) then

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

    end if

    !----------------------------------------------------------------
    ! ANNUAL OUTPUT
    !----------------------------------------------------------------
    ! ANNUAL TOTAL N LOSS (gN/m2/yr)
    filnam=trim(prefix)//'.a.nloss.out'
    open(550,file=filnam,err=888,status='unknown')

    ! ANNUAL N2O EMISSIONS (gN/m2/yr)
    filnam=trim(prefix)//'.a.n2o.out'
    open(551,file=filnam,err=888,status='unknown')

    return

    888  stop 'INITIO_NTRANSFORM: error opening output files'

  end subroutine initio_ntransform


  subroutine initoutput_ntransform()
    !////////////////////////////////////////////////////////////////
    !  Initialises waterbalance-specific output variables
    !----------------------------------------------------------------

    if (loutntransform) then

      allocate( outdnloss ( nlu,ndayyear,maxgrid ) ) ! daily total N loss (gaseous+leacing) (gN/m2/d)
      allocate( outddenitr( nlu,ndayyear,maxgrid ) ) ! daily amount of N denitrified (gN/m2/d)
      allocate( outdnitr  ( nlu,ndayyear,maxgrid ) ) ! daily amount of N nitrified (gN/m2/d)
      allocate( outdnvol  ( nlu,ndayyear,maxgrid ) ) ! daily amount of N volatilised (gN/m2/d)
      allocate( outdnleach( nlu,ndayyear,maxgrid ) ) ! daily amount of N leached (gN/m2/d)
      allocate( outdn2o   ( nlu,ndayyear,maxgrid ) ) ! daily N2O emitted (gN/m2/d)

      outdnloss(:,:,:)  = 0.0
      outddenitr(:,:,:) = 0.0
      outdnitr(:,:,:)   = 0.0
      outdnvol(:,:,:)   = 0.0
      outdnleach(:,:,:) = 0.0
      outdn2o(:,:,:)    = 0.0

    end if
      
    outanloss(:,:)    = 0.0
    outan2o(:,:)      = 0.0

  end subroutine initoutput_ntransform


  subroutine getout_daily_ntransform( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    !  SR called daily to sum up output variables.
    !----------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy    
    integer, intent(in) :: doy    

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    !----------------------------------------------------------------
    if (loutntransform) then
      outdnloss(:,doy,jpngr)  = dnloss(:)
      outddenitr(:,doy,jpngr) = ddenitr(:)
      outdnitr(:,doy,jpngr)   = dnitr(:)
      outdnvol(:,doy,jpngr)   = dnvol(:)
      outdnleach(:,doy,jpngr) = dnleach(:)
      outdn2o(:,doy,jpngr)    = dn2o(:)
    end if

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    outanloss(:,jpngr)    = outanloss(:,jpngr) + dnloss(:)
    outan2o(:,jpngr)      = outan2o(:,jpngr) + dn2o(:)

  end subroutine getout_daily_ntransform


  subroutine writeout_ascii_ntransform( year, spinup )
    !/////////////////////////////////////////////////////////////////////////
    ! WRITE WATERBALANCE-SPECIFIC VARIABLES TO OUTPUT
    !-------------------------------------------------------------------------
    use _params_core, only: npft
    use _params_siml, only: firstyeartrend, spinupyears, daily_out_startyr, &
      daily_out_endyr, outyear

    ! arguments
    integer, intent(in) :: year       ! simulation year
    logical, intent(in) :: spinup     ! true during spinup years

    ! local variables
    real :: itime
    integer :: myday, myjpngr

    ! xxx implement this: sum over gridcells? single output per gridcell?
    if (maxgrid>1) stop 'writeout_ascii_ntransform: think of something ...'
    myjpngr = 1

    !-------------------------------------------------------------------------
    ! DAILY OUTPUT
    !-------------------------------------------------------------------------
    if (loutntransform) then
      if ( .not. spinup .and. outyear>=daily_out_startyr .and. outyear<=daily_out_endyr ) then
        ! Write daily output only during transient simulation
        do myday=1,ndayyear

          ! Define 'itime' as a decimal number corresponding to day in the year + year
          itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(myday-1)/real(ndayyear)

          if (nlu>1) stop 'writeout_ascii_ntransform: write out lu-area weighted sum'
          if (npft>1) stop 'writeout_ascii_ntransform: think of something for ccost output'

          write(500,999) itime, sum(outdnloss(:,myday,myjpngr))
          write(501,999) itime, sum(outdnvol(:,myday,myjpngr))
          write(502,999) itime, sum(outddenitr(:,myday,myjpngr))
          write(503,999) itime, sum(outdnitr(:,myday,myjpngr))
          write(504,999) itime, sum(outdnleach(:,myday,myjpngr))
          write(505,999) itime, sum(outdn2o(:,myday,myjpngr))

        end do
      end if
    end if

    !-------------------------------------------------------------------------
    ! ANNUAL OUTPUT
    ! Write annual value, summed over all PFTs / LUs
    ! xxx implement taking sum over PFTs (and gridcells) in this land use category
    !-------------------------------------------------------------------------
    itime = real(year) + real(firstyeartrend) - real(spinupyears)

    write(550,999) itime, sum(outanloss(:,myjpngr))
    write(551,999) itime, sum(outan2o(:,myjpngr))

    ! write(0,*) 'outan2o written to output ', sum(outan2o(:,myjpngr))

    return
    
    999 format (F20.8,F20.8)

  end subroutine writeout_ascii_ntransform


end module _ntransform
