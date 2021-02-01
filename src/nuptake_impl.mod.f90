module md_nuptake_impl
  !////////////////////////////////////////////////////////////////
  ! Contains functions for using empirical relationships to model
  ! N uptake, driven by GPP 
  !--------------------------------------------------------------
  ! load modules
  use md_params_core
  use md_tile, only: tile_type, tile_fluxes_type
  use md_interface, only: interface
  use md_grid, only: gridtype, domaininfo_type

  implicit none
  private
  public nuptake_impl, get_preds_nimpl, getpar_nimpl, initio_nc_nimpl, initoutput_nimpl, &
    getout_annual_nimpl, writeout_nc_nimpl

  !----------------------------------------------------------------
  ! Variables without memory (not necessarily just fluxes; just define the type) 
  !----------------------------------------------------------------
  type plant_nimpl_fluxes_type
    !Forest model
    real :: anpp        ! annual total net primary production (gC m-2 yr-1)
    real :: aanpp       ! annual aboveground net primary production (gC m-2 yr-1)
    real :: abnpp       ! annual belowground net primary production (gC m-2 yr-1)
    real :: alnpp       ! annual net primary production for leaf production (gC m-2 yr-1)
    real :: awnpp       ! annual net primary production for wood production (gC m-2 yr-1)
    real :: leafcn      ! the ratio of leaf nitrogen per mass to leaf carbon per mass (unitness)
    real :: nre         ! nitrogen resorption efficiency of leaf
    real :: alnf        ! annual leaf nitrogen flux (gN m-2 yr-2)
    real :: awnf        ! annual wood nitrogen flux (gN m-2 yr-2)
    real :: abnf        ! annual belowground(root) nitrogen flux (gN m-2 yr-2)
    real :: nuptake     ! annual nitrogen uptake in ecosystem (gN m-2 yr-1)
    real :: avcmax25    ! annual maximum carboxylation capacity (umol m-2 s-1)
    !Grassland model    
    real :: anpp_grass  ! annual total net primary production in grassland (gC m-2 yr-1)
    real :: aanpp_grass  ! annual aboveground net primary production in grassland (gC m-2 yr-1)
    real :: abnpp_grass  ! annual belowground net primary production in grassland (gC m-2 yr-1)
    real :: alnf_grass        ! annual leaf nitrogen flux in grassland (gN m-2 yr-2)
    real :: abnf_grass        ! annual belowground(root) nitrogen flux in grassland  (gN m-2 yr-2)
    real :: nuptake_grass     ! annual nitrogen uptake in ecosystem in grassland (gN m-2 yr-1)
  end type plant_nimpl_fluxes_type

  type tile_nimpl_fluxes_type
    type(plant_nimpl_fluxes_type), dimension(npft) :: plant
  end type tile_nimpl_fluxes_type

  !----------------------------------------------------------------
  ! object containing module-specific fluxes (create this new object)
  !----------------------------------------------------------------
  type(tile_nimpl_fluxes_type), dimension(nlu) :: tile_nimpl_fluxes

  !----------------------------------------------------------------
  ! Coefficients of statistical models
  !----------------------------------------------------------------
  type coef_nimpl_type

    !!!Forest model below
    ! BP:GPP model
    real :: cnsoil_bp
    real :: age_bp
    real :: fapar_bp
    real :: alpha_bp
    real :: intersect_bp

    ! ANPP:GPP model
    real :: cnsoil_anpp
    real :: age_anpp
    real :: fapar_anpp
    real :: alpha_anpp
    real :: intersect_anpp

    ! ALNPP:NPP model
    real :: ppfd_alnpp
    real :: tg_alnpp
    real :: vpd_alnpp
    real :: intersect_alnpp

    ! Leaf C:N model
    real :: vcmax25_leafcn
    real :: lma_leafcn
    real :: cmass_leafcn

    ! nre model
    real :: tg_nre
    real :: vpd_nre
    real :: intersect_nre    

    ! Constant ratio
    real :: root_cn
    real :: wood_cn

    !!!Grassland model below
    ! NPP:GPP model
    real :: cue_grass

    ! ANPP:NPP model
    real :: tg_anpp_grass
    real :: alpha_anpp_grass
    real :: intersect_anpp_grass

    ! Constant ratio of root in grassland
    real :: root_cn_grass


  end type coef_nimpl_type  

  type(coef_nimpl_type) :: coef_nimpl

  !----------------------------------------------------------------
  ! Predictor fields
  !----------------------------------------------------------------
  type preds_nimpl_type
    real :: cnsoil
    real :: age
    real :: fapar
    real :: alpha
    real :: ppfd
    real :: tg
    real :: vpd
    real :: lma
  end type preds_nimpl_type

  type(preds_nimpl_type), dimension(:), allocatable :: preds_nimpl

  !----------------------------------------------------------------
  ! Specify file and variable names for NetCDF reading (Vcmax25 may need additional input from SOFUN)
  !----------------------------------------------------------------
  character(len=100), parameter :: filnam_cnsoil = "./input/global/nimpl/CNrt.nc"
  character(len=100), parameter :: varnam_cnsoil = "CNrt"
  character(len=100), parameter :: filnam_age    = "./input/global/nimpl/age.nc"
  character(len=100), parameter :: varnam_age    = "age"
  character(len=100), parameter :: filnam_fapar  = "./input/global/nimpl/fAPAR.nc"
  character(len=100), parameter :: varnam_fapar  = "fAPAR"
  character(len=100), parameter :: filnam_alpha  = "./input/global/nimpl/alpha.nc"
  character(len=100), parameter :: varnam_alpha  = "alpha"
  character(len=100), parameter :: filnam_ppfd   = "./input/global/nimpl/PPFD.nc"
  character(len=100), parameter :: varnam_ppfd   = "PPFD"
  character(len=100), parameter :: filnam_Tg     = "./input/global/nimpl/Tg.nc"
  character(len=100), parameter :: varnam_Tg     = "Tg"
  character(len=100), parameter :: filnam_vpd    = "./input/global/nimpl/vpd.nc"
  character(len=100), parameter :: varnam_vpd    = "vpd"
  character(len=100), parameter :: filnam_lma    = "./input/global/nimpl/LMA.nc"
  character(len=100), parameter :: varnam_lma    = "LMA"
  
  !----------------------------------------------------------------
  ! Module-specific output variables
  !----------------------------------------------------------------
  ! annual, forest
  real, dimension(:), allocatable :: outanpp
  real, dimension(:), allocatable :: outaanpp
  real, dimension(:), allocatable :: outabnpp
  real, dimension(:), allocatable :: outalnpp
  real, dimension(:), allocatable :: outawnpp
  real, dimension(:), allocatable :: outleafcn
  real, dimension(:), allocatable :: outnre
  real, dimension(:), allocatable :: outalnf
  real, dimension(:), allocatable :: outawnf
  real, dimension(:), allocatable :: outabnf
  real, dimension(:), allocatable :: outnuptake
  real, dimension(:), allocatable :: outavcmax25
  ! annual, grassland
  real, dimension(:), allocatable :: outanpp_grass
  real, dimension(:), allocatable :: outaanpp_grass
  real, dimension(:), allocatable :: outabnpp_grass
  real, dimension(:), allocatable :: outalnf_grass
  real, dimension(:), allocatable :: outabnf_grass
  real, dimension(:), allocatable :: outnuptake_grass
  !----------------------------------------------------------------
  ! Module-specific NetCDF output file and variable names
  !----------------------------------------------------------------
  !Forest
  character(len=256) :: ncoutfilnam_anpp
  character(len=256) :: ncoutfilnam_aanpp
  character(len=256) :: ncoutfilnam_abnpp
  character(len=256) :: ncoutfilnam_alnpp
  character(len=256) :: ncoutfilnam_awnpp
  character(len=256) :: ncoutfilnam_leafcn
  character(len=256) :: ncoutfilnam_nre
  character(len=256) :: ncoutfilnam_lnf
  character(len=256) :: ncoutfilnam_wnf
  character(len=256) :: ncoutfilnam_bnf
  character(len=256) :: ncoutfilnam_nuptake
  character(len=256) :: ncoutfilnam_annualvcmax25
  !Grassland
  character(len=256) :: ncoutfilnam_anpp_grass
  character(len=256) :: ncoutfilnam_aanpp_grass
  character(len=256) :: ncoutfilnam_abnpp_grass
  character(len=256) :: ncoutfilnam_lnf_grass
  character(len=256) :: ncoutfilnam_bnf_grass
  character(len=256) :: ncoutfilnam_nuptake_grass

  !Forest
  character(len=*), parameter :: NPP_NAME = "npp"
  character(len=*), parameter :: ANPP_NAME = "anpp"
  character(len=*), parameter :: BNPP_NAME = "bnpp"
  character(len=*), parameter :: LNPP_NAME = "lnpp"
  character(len=*), parameter :: WNPP_NAME = "wnpp"
  character(len=*), parameter :: LEAFCN_NAME = "leafcn"
  character(len=*), parameter :: NRE_NAME = "nre"
  character(len=*), parameter :: LNF_NAME = "lnf"
  character(len=*), parameter :: WNF_NAME = "wnf"
  character(len=*), parameter :: BNF_NAME = "bnf"
  character(len=*), parameter :: nuptake_NAME = "nuptake"
  character(len=*), parameter :: ANNUALVCMAX25_NAME = "annualvcmax25"
  !Grassland
  character(len=*), parameter :: NPP_NAME_GRASS = "npp_grass"
  character(len=*), parameter :: ANPP_NAME_GRASS = "anpp_grass"
  character(len=*), parameter :: BNPP_NAME_GRASS = "bnpp_grass"
  character(len=*), parameter :: LNF_NAME_GRASS = "lnf_grass"
  character(len=*), parameter :: BNF_NAME_GRASS = "bnf_grass"
  character(len=*), parameter :: nuptake_NAME_GRASS = "nuptake_grass"


contains

  subroutine nuptake_impl( jpngr, dogridcell, tile, tile_fluxes, init )
    !////////////////////////////////////////////////////////////////
    ! Determines all the downstream fluxes as a function of GPP
    ! using pre-defined statistical relationships
    !--------------------------------------------------------------
    ! arguments
    integer, intent(in) :: jpngr
    logical, intent(in) :: dogridcell
    type(tile_type), dimension(nlu), intent(in) :: tile
    type(tile_fluxes_type), dimension(nlu), intent(in) :: tile_fluxes
    logical, intent(in) :: init

    ! local variable
    integer :: lu, pft

    !--------------------------------------------------------------
    ! Predict using statistical models
    ! Note that some data was log-transfromed (cnsoil, age, PPFD, vpd, Vcmax25, lma) while some are not (alpha, fAPAR, Tg)
    ! All ratios were using logit function (logit(y)=x), therefore we should convert them into y = 1/(1+ exp(-x))
    ! leaf c/n model were using log function, so it should be exp in advance
    !--------------------------------------------------------------
    ! Make tile_nimpl_fluxes a field
    if (dogridcell) then
      lu = 1
      do pft = 1,npft
        if ((preds_nimpl(jpngr)%cnsoil > 0.0).and.(preds_nimpl(jpngr)%age > 0.0).and.(preds_nimpl(jpngr)%fapar > 0.0).and.(preds_nimpl(jpngr)%alpha > 0.0)) then
          tile_nimpl_fluxes(lu)%plant(pft)%anpp  = tile_fluxes(lu)%plant(pft)%agpp * (1/(1 + EXP(-(coef_nimpl%cnsoil_bp * LOG(preds_nimpl(jpngr)%cnsoil) + coef_nimpl%age_bp * LOG(preds_nimpl(jpngr)%age) + coef_nimpl%fapar_bp * preds_nimpl(jpngr)%fapar +coef_nimpl%alpha_bp * preds_nimpl(jpngr)%alpha + coef_nimpl%intersect_bp))))
          tile_nimpl_fluxes(lu)%plant(pft)%aanpp = tile_fluxes(lu)%plant(pft)%agpp * (1/(1 + EXP(-(coef_nimpl%cnsoil_anpp * LOG(preds_nimpl(jpngr)%cnsoil) + coef_nimpl%age_anpp * LOG(preds_nimpl(jpngr)%age) + coef_nimpl%fapar_anpp * preds_nimpl(jpngr)%fapar +coef_nimpl%alpha_anpp * preds_nimpl(jpngr)%alpha + coef_nimpl%intersect_anpp))))
          tile_nimpl_fluxes(lu)%plant(pft)%abnpp = tile_nimpl_fluxes(lu)%plant(pft)%anpp - tile_nimpl_fluxes(lu)%plant(pft)%aanpp
          else
            tile_nimpl_fluxes(lu)%plant(pft)%anpp  = 0
            tile_nimpl_fluxes(lu)%plant(pft)%aanpp  = 0
            tile_nimpl_fluxes(lu)%plant(pft)%abnpp = 0
            !print*,'tile_fluxes(lu)%plant(pft)%agpp',tile_fluxes(lu)%plant(pft)%agpp 
            !print*,'jpngr',jpngr
            !print*,'tile_nimpl_fluxes(lu)%plant(pft)%anpp',tile_nimpl_fluxes(lu)%plant(pft)%anpp
            !print*,'jpngr',jpngr            
            !print*,'preds_nimpl(jpngr)%cnsoil ',preds_nimpl(jpngr)%cnsoil 
            !print*,'jpngr',jpngr
            !print*,'preds_nimpl(jpngr)%age ',preds_nimpl(jpngr)%age
            !print*,'jpngr',jpngr 
            !print*,'preds_nimpl(jpngr)%fapar',preds_nimpl(jpngr)%fapar
            !print*,'jpngr',jpngr 
            !print*,'preds_nimpl(jpngr)%alpha ',preds_nimpl(jpngr)%alpha
            !print*,'jpngr',jpngr
        end if

        !tile_nimpl_fluxes(lu)%plant(pft)%abnpp = tile_nimpl_fluxes(lu)%plant(pft)%anpp - tile_nimpl_fluxes(lu)%plant(pft)%aanpp
        !print*,'4'
        if ((preds_nimpl(jpngr)%ppfd > 0.0).and.(preds_nimpl(jpngr)%vpd > 0.0)) then
          tile_nimpl_fluxes(lu)%plant(pft)%alnpp = tile_nimpl_fluxes(lu)%plant(pft)%aanpp * (1/(1+EXP(-(coef_nimpl%ppfd_alnpp * LOG(preds_nimpl(jpngr)%ppfd) + coef_nimpl%tg_alnpp * preds_nimpl(jpngr)%tg + coef_nimpl%vpd_alnpp * LOG(preds_nimpl(jpngr)%vpd) + coef_nimpl%intersect_alnpp))))
          else
            tile_nimpl_fluxes(lu)%plant(pft)%alnpp = 0
        end if
        !print*,'5'
        if ((preds_nimpl(jpngr)%cnsoil > 0.0).and.(preds_nimpl(jpngr)%age > 0.0).and.(preds_nimpl(jpngr)%fapar > 0.0).and.(preds_nimpl(jpngr)%alpha > 0.0).and.(preds_nimpl(jpngr)%ppfd > 0.0).and.(preds_nimpl(jpngr)%vpd > 0.0)) then
          tile_nimpl_fluxes(lu)%plant(pft)%awnpp = tile_nimpl_fluxes(lu)%plant(pft)%aanpp - tile_nimpl_fluxes(lu)%plant(pft)%alnpp
          else
            tile_nimpl_fluxes(lu)%plant(pft)%awnpp = 0
        end if
        !model 1: Nmass ~ vcmax25/lma - r2 = 0.066
        if ((preds_nimpl(jpngr)%lma > 0.0)) then
          tile_nimpl_fluxes(lu)%plant(pft)%leafcn = (coef_nimpl%lma_leafcn/coef_nimpl%cmass_leafcn) + (coef_nimpl%vcmax25_leafcn/coef_nimpl%cmass_leafcn) * (tile_fluxes(lu)%canopy%avcmax25_max)/(preds_nimpl(jpngr)%lma)!it is actually leaf n/c here...    
          !print*,'tile_nimpl_fluxes(lu)%plant(pft)%leafcn',tile_nimpl_fluxes(lu)%plant(pft)%leafcn
          !print*,'jpngr',jpngr
          else
            tile_nimpl_fluxes(lu)%plant(pft)%leafcn = 0
        end if      
          !model 1b: Nmass ~ vcmax25/lma + soil C/N.
          !tile_nimpl_fluxes(lu)%plant(pft)%leafcn = (0.017 + (0.005)*((tile_fluxes(lu)%plant(pft)%avcmax25_max)/(preds_nimpl(jpngr)%lma)) + (-0.0003)*(preds_nimpl(jpngr)%cnsoil))/0.4638 !it is actually leaf n/c here...
          !model 1c: Nmass ~ vcmax25 + lma + soil C/N. 
          !tile_nimpl_fluxes(lu)%plant(pft)%leafcn = (0.02922 + (-0.0004462)*(preds_nimpl(jpngr)%cnsoil) + 0.00002715*(tile_fluxes(lu)%plant(pft)%avcmax25_max) + (-0.00006713)*(preds_nimpl(jpngr)%lma))/0.4638 !it is actually leaf n/c here...
        !model 2: Nmass ~ Tg + PPFD + LMA
        !tile_nimpl_fluxes(lu)%plant(pft)%leafcn = (0.02035+ (-0.0001095)*preds_nimpl(jpngr)%tg + (0.00001875)*(preds_nimpl(jpngr)%ppfd) + (-0.00006863)*(preds_nimpl(jpngr)%lma))/0.4638  !it is actually leaf n/c here...the alternative model for PPFD + Tg + LMA     
        
        !nre
        if (preds_nimpl(jpngr)%vpd > 0.0) then
          tile_nimpl_fluxes(lu)%plant(pft)%nre = (1/(1+EXP(-(coef_nimpl%tg_nre * preds_nimpl(jpngr)%tg + coef_nimpl%vpd_nre * LOG(preds_nimpl(jpngr)%vpd) + coef_nimpl%intersect_nre))))
          else
            tile_nimpl_fluxes(lu)%plant(pft)%nre = 0
        end if

        !exclude nre in lnf????  
        tile_nimpl_fluxes(lu)%plant(pft)%alnf = tile_nimpl_fluxes(lu)%plant(pft)%alnpp * tile_nimpl_fluxes(lu)%plant(pft)%leafcn * (1.0 - tile_nimpl_fluxes(lu)%plant(pft)%nre)!it is actually leaf n/c here. 
        !tile_nimpl_fluxes(lu)%plant(pft)%alnf = tile_nimpl_fluxes(lu)%plant(pft)%alnpp * tile_nimpl_fluxes(lu)%plant(pft)%leafcn!it is actually leaf n/c here. 
        
        tile_nimpl_fluxes(lu)%plant(pft)%awnf = tile_nimpl_fluxes(lu)%plant(pft)%awnpp / coef_nimpl%wood_cn
        tile_nimpl_fluxes(lu)%plant(pft)%abnf = tile_nimpl_fluxes(lu)%plant(pft)%abnpp / coef_nimpl%root_cn

        !exclude nre in nuptake?????????          
        tile_nimpl_fluxes(lu)%plant(pft)%nuptake = tile_nimpl_fluxes(lu)%plant(pft)%alnpp * tile_nimpl_fluxes(lu)%plant(pft)%leafcn * (1.0 - tile_nimpl_fluxes(lu)%plant(pft)%nre) + (tile_nimpl_fluxes(lu)%plant(pft)%awnpp / coef_nimpl%wood_cn) + (tile_nimpl_fluxes(lu)%plant(pft)%abnpp / coef_nimpl%root_cn)
        !tile_nimpl_fluxes(lu)%plant(pft)%nuptake = tile_nimpl_fluxes(lu)%plant(pft)%alnpp * tile_nimpl_fluxes(lu)%plant(pft)%leafcn + (tile_nimpl_fluxes(lu)%plant(pft)%awnpp / coef_nimpl%wood_cn) + (tile_nimpl_fluxes(lu)%plant(pft)%abnpp / coef_nimpl%root_cn)

        !!Grassland models
        tile_nimpl_fluxes(lu)%plant(pft)%anpp_grass  = tile_fluxes(lu)%plant(pft)%agpp * coef_nimpl%cue_grass

        if ((preds_nimpl(jpngr)%alpha > 0.0).and.(preds_nimpl(jpngr)%tg > 0.0)) then
          tile_nimpl_fluxes(lu)%plant(pft)%aanpp_grass = tile_nimpl_fluxes(lu)%plant(pft)%anpp_grass * (1/(1+EXP(-(coef_nimpl%alpha_anpp_grass * preds_nimpl(jpngr)%alpha + coef_nimpl%tg_anpp_grass * preds_nimpl(jpngr)%tg + coef_nimpl%intersect_anpp_grass))))
          tile_nimpl_fluxes(lu)%plant(pft)%abnpp_grass = tile_nimpl_fluxes(lu)%plant(pft)%anpp_grass - tile_nimpl_fluxes(lu)%plant(pft)%aanpp_grass         
          else
            tile_nimpl_fluxes(lu)%plant(pft)%aanpp_grass = 0
            tile_nimpl_fluxes(lu)%plant(pft)%abnpp_grass = 0
        end if        

        tile_nimpl_fluxes(lu)%plant(pft)%alnf_grass = tile_nimpl_fluxes(lu)%plant(pft)%aanpp_grass * tile_nimpl_fluxes(lu)%plant(pft)%leafcn * (1.0 - tile_nimpl_fluxes(lu)%plant(pft)%nre)!it is actually leaf n/c here. 
        tile_nimpl_fluxes(lu)%plant(pft)%abnf_grass = tile_nimpl_fluxes(lu)%plant(pft)%abnpp_grass / coef_nimpl%root_cn_grass
        tile_nimpl_fluxes(lu)%plant(pft)%nuptake_grass = (tile_nimpl_fluxes(lu)%plant(pft)%abnpp_grass / coef_nimpl%root_cn_grass)  + tile_nimpl_fluxes(lu)%plant(pft)%aanpp_grass * tile_nimpl_fluxes(lu)%plant(pft)%leafcn * (1.0 - tile_nimpl_fluxes(lu)%plant(pft)%nre)

      end do
    else

      lu = 1
      !print*,'1 b'
      tile_nimpl_fluxes(lu)%plant(:)%anpp  = dummy
      !print*,'2 b'
      tile_nimpl_fluxes(lu)%plant(:)%aanpp = dummy
      !print*,'3 b'
      tile_nimpl_fluxes(lu)%plant(:)%abnpp = dummy
      !print*,'4 b'
      tile_nimpl_fluxes(lu)%plant(:)%alnpp = dummy
      !print*,'5 b'
      tile_nimpl_fluxes(lu)%plant(:)%awnpp = dummy
      !print*,'6 b'
      tile_nimpl_fluxes(lu)%plant(:)%leafcn = dummy
      !print*,'7 b'
      tile_nimpl_fluxes(lu)%plant(:)%nre = dummy
      tile_nimpl_fluxes(lu)%plant(:)%alnf = dummy
      tile_nimpl_fluxes(lu)%plant(:)%awnf = dummy
      tile_nimpl_fluxes(lu)%plant(:)%abnf = dummy

      !!Grassland models
      tile_nimpl_fluxes(lu)%plant(:)%anpp_grass = dummy
      tile_nimpl_fluxes(lu)%plant(:)%aanpp_grass = dummy
      tile_nimpl_fluxes(lu)%plant(:)%abnpp_grass = dummy
      tile_nimpl_fluxes(lu)%plant(:)%alnf_grass = dummy
      tile_nimpl_fluxes(lu)%plant(:)%abnf_grass = dummy
      tile_nimpl_fluxes(lu)%plant(:)%nuptake_grass = dummy       

    end if

  end subroutine nuptake_impl


  subroutine get_preds_nimpl( domaininfo, grid )
    !////////////////////////////////////////////////////////////////
    ! Some explanations XXX 
    !----------------------------------------------------------------
    ! arguments
    use md_params_core, only: npft

    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid

    ! allocate memory
    allocate( preds_nimpl(domaininfo%maxgrid) )

    !--------------------------------------------------------------
    ! Read predictor fields from files, populates 'preds_nimpl'
    !--------------------------------------------------------------
    ! NPP:GPP model
    call get_preds_nc_byvar( trim(filnam_cnsoil), trim(varnam_cnsoil), domaininfo, grid, preds_nimpl(:)%cnsoil )
    call get_preds_nc_byvar( trim(filnam_age),    trim(varnam_age),    domaininfo, grid, preds_nimpl(:)%age )
    call get_preds_nc_byvar( trim(filnam_fapar),  trim(varnam_fapar),  domaininfo, grid, preds_nimpl(:)%fapar )
    call get_preds_nc_byvar( trim(filnam_alpha),  trim(varnam_alpha),  domaininfo, grid, preds_nimpl(:)%alpha )

    ! ! ANPP:GPP model
    ! call get_preds_nc_byvar( trim(filnam_cnsoil), trim(varnam_cnsoil), domaininfo, grid, preds_nimpl(:)%cnsoil )
    ! call get_preds_nc_byvar( trim(filnam_age),    trim(varnam_age),    domaininfo, grid, preds_nimpl(:)%age )
    ! call get_preds_nc_byvar( trim(filnam_fapar),  trim(varnam_fapar),  domaininfo, grid, preds_nimpl(:)%fapar )
    ! call get_preds_nc_byvar( trim(filnam_alpha),  trim(varnam_alpha),  domaininfo, grid, preds_nimpl(:)%alpha )

    ! ALNPP:NPP model
    call get_preds_nc_byvar( trim(filnam_ppfd), trim(varnam_ppfd), domaininfo, grid, preds_nimpl(:)%ppfd )
    call get_preds_nc_byvar( trim(filnam_tg),   trim(varnam_tg),   domaininfo, grid, preds_nimpl(:)%tg )
    call get_preds_nc_byvar( trim(filnam_vpd),  trim(varnam_vpd),  domaininfo, grid, preds_nimpl(:)%vpd )

    ! Leaf C:N model
    ! call get_preds_nc_byvar( trim(filnam_vcmax25), trim(varnam_vcmax25), domaininfo, grid, preds_nimpl(:)%vcmax25 )
    call get_preds_nc_byvar( trim(filnam_lma),     trim(varnam_lma),     domaininfo, grid, preds_nimpl(:)%lma )

  end subroutine get_preds_nimpl


  subroutine getpar_nimpl()
    !////////////////////////////////////////////////////////////////
    !Extract all coeffients and constant in a dat file
    !--------------------------------------------------------------
    use md_sofunutils, only: getparreal
    !!Forest models
    ! BP/GPP model
    coef_nimpl%cnsoil_bp  = getparreal( 'params/params_nimpl.dat', 'cnsoil_bp' )
    coef_nimpl%age_bp     = getparreal( 'params/params_nimpl.dat', 'age_bp' )
    coef_nimpl%fapar_bp  = getparreal( 'params/params_nimpl.dat', 'fapar_bp' )
    coef_nimpl%alpha_bp     = getparreal( 'params/params_nimpl.dat', 'alpha_bp' )
    coef_nimpl%intersect_bp     = getparreal( 'params/params_nimpl.dat', 'intersect_bp' )

    ! ANPP/GPP model
    coef_nimpl%cnsoil_anpp  = getparreal( 'params/params_nimpl.dat', 'cnsoil_anpp' )
    coef_nimpl%age_anpp     = getparreal( 'params/params_nimpl.dat', 'age_anpp' )
    coef_nimpl%fapar_anpp  = getparreal( 'params/params_nimpl.dat', 'fapar_anpp' )
    coef_nimpl%alpha_anpp     = getparreal( 'params/params_nimpl.dat', 'alpha_anpp' )
    coef_nimpl%intersect_anpp     = getparreal( 'params/params_nimpl.dat', 'intersect_anpp' )

    ! ALNPP/NPP model
    coef_nimpl%ppfd_alnpp = getparreal( 'params/params_nimpl.dat', 'ppfd_alnpp' )
    coef_nimpl%tg_alnpp = getparreal( 'params/params_nimpl.dat', 'tg_alnpp' )
    coef_nimpl%vpd_alnpp = getparreal( 'params/params_nimpl.dat', 'vpd_alnpp' )
    coef_nimpl%intersect_alnpp     = getparreal( 'params/params_nimpl.dat', 'intersect_alnpp' )
    
    ! Leaf C:N model
    coef_nimpl%vcmax25_leafcn = getparreal( 'params/params_nimpl.dat', 'vcmax25_leafcn' )
    coef_nimpl%lma_leafcn = getparreal( 'params/params_nimpl.dat', 'lma_leafcn' )
    coef_nimpl%cmass_leafcn     = getparreal( 'params/params_nimpl.dat', 'cmass_leafcn' )

    ! leaf nre model
    coef_nimpl%tg_nre = getparreal( 'params/params_nimpl.dat', 'tg_nre' )
    coef_nimpl%vpd_nre = getparreal( 'params/params_nimpl.dat', 'vpd_nre' )
    coef_nimpl%intersect_nre     = getparreal( 'params/params_nimpl.dat', 'intersect_nre' )

    ! Constant
    coef_nimpl%root_cn = getparreal( 'params/params_nimpl.dat', 'root_cn' )
    coef_nimpl%wood_cn = getparreal( 'params/params_nimpl.dat', 'wood_cn' )

    !!Grassland models
    ! BP/GPP model
    coef_nimpl%cue_grass = getparreal( 'params/params_nimpl.dat', 'cue_grass' )
    ! ANPP/NPP model
    coef_nimpl%tg_anpp_grass = getparreal( 'params/params_nimpl.dat', 'tg_anpp_grass' )
    coef_nimpl%alpha_anpp_grass = getparreal( 'params/params_nimpl.dat', 'alpha_anpp_grass' )
    coef_nimpl%intersect_anpp_grass = getparreal( 'params/params_nimpl.dat', 'intersect_anpp_grass' )
    ! Constant of root C/N
    coef_nimpl%root_cn_grass = getparreal( 'params/params_nimpl.dat', 'root_cn_grass' )

  end subroutine getpar_nimpl


  subroutine get_preds_nc_byvar( filnam, varname, domaininfo, grid, pred )
    !////////////////////////////////////////////////////////////////
    ! xxx add explanation here
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: check

    ! arguments
    character(len=*) :: filnam
    character(len=*) :: varname
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(inout) :: grid
    real, dimension(domaininfo%maxgrid), intent(out) :: pred

    ! local variables
    integer :: ncid, varid
    integer :: latdimid, londimid
    integer :: nlat_arr, nlon_arr
    real, allocatable, dimension(:)   :: lon_arr
    real, allocatable, dimension(:)   :: lat_arr
    real, allocatable, dimension(:,:) :: pred_arr

    integer :: i, pft, jpngr, ilon_arr, ilat_arr, n_noinfo
    integer, dimension(domaininfo%maxgrid) :: ilon
    integer, dimension(domaininfo%maxgrid) :: ilat
    real, allocatable, dimension(:) :: tmp
    real :: ncfillvalue
    real :: dlat, dlon
    character(len=3), parameter :: lonname = "lon"
    character(len=3), parameter :: latname = "lat"
    ! character(len=100), parameter :: dimname_pft = "z"
    ! character(len=100), parameter :: varname = "pftcover"
    ! character(len=100), parameter :: filnam = "./input/global/landcover/modis_landcover_halfdeg_2010_FILLED.nc"

    !----------------------------------------------------------------  
    ! Get vegetation cover information from file
    !----------------------------------------------------------------
    print*,'getting predictor from ', trim(filnam), ' ...'

    ! Read arrays of all months of current year from file  
    call check( nf90_open( trim(filnam), NF90_NOWRITE, ncid ) )

    ! get dimension ID for latitude
    call check( nf90_inq_dimid( ncid, trim(latname), latdimid ) )

    ! Get latitude information: nlat
    call check( nf90_inquire_dimension( ncid, latdimid, len = nlat_arr ) )

    ! get dimension ID for longitude
    call check( nf90_inq_dimid( ncid, trim(lonname), londimid ) )

    ! Get latitude information: nlon
    call check( nf90_inquire_dimension( ncid, londimid, len = nlon_arr ) )

    ! for index association, get ilon and ilat vectors
    ! Allocate array sizes now knowing nlon and nlat 
    allocate( lon_arr(nlon_arr) )
    allocate( lat_arr(nlat_arr) )

    !! Get longitude and latitude values
    !call check( nf90_get_var( ncid, londimid, lon_arr ) )
    !call check( nf90_get_var( ncid, latdimid, lat_arr ) )

    ! xxx try:
    lon_arr = (/ (i, i = 1,nlon_arr) /)
    lon_arr = (lon_arr - 1) * 0.5 - 180.0 + 0.25 
    lat_arr = (/ (i, i = 1,nlat_arr) /)
    lat_arr = (lat_arr - 1) * 0.5 - 90.0 + 0.25 

    ! Check if the resolution of the climate input files is identical to the model grid resolution
    dlon = lon_arr(2) - lon_arr(1)
    dlat = lat_arr(2) - lat_arr(1)

    !print*,'dlon', dlon
    !print*,'dlat', dlat

    if (dlon/=domaininfo%dlon) stop 'Longitude resolution of nimpl predictor file is not identical with model grid.'
    if (dlat/=domaininfo%dlat) stop 'latitude resolution of nimpl predictor file is not identical with model grid.'

    ! get index associations
    do jpngr=1,domaininfo%maxgrid
      ilon_arr = 1
      do while (grid(jpngr)%lon/=lon_arr(ilon_arr))
        ilon_arr = ilon_arr + 1
      end do
      ilon(jpngr) = ilon_arr

      ilat_arr = 1
      do while (grid(jpngr)%lat/=lat_arr(ilat_arr))
        ilat_arr = ilat_arr + 1
      end do
      ilat(jpngr) = ilat_arr
    end do

    ! allocate size of output array
    allocate( pred_arr(nlon_arr,nlat_arr) )
    ! allocate( tmp(npft_in))

    ! Get the varid of the data variable, based on its name
    call check( nf90_inq_varid( ncid, trim(varname), varid ) )

    ! Read the array
    call check( nf90_get_var( ncid, varid, pred_arr, start=(/1, 1/), count=(/nlon_arr, nlat_arr/) ) )

    ! Get _FillValue from file (assuming that all are the same for WATCH-WFDEI)
    call check( nf90_get_att( ncid, varid, "_FillValue", ncfillvalue ) )

    ! close NetCDF files
    call check( nf90_close( ncid ) )

    ! read from array to define field    
    do jpngr=1,domaininfo%maxgrid
      
      pred(jpngr) = pred_arr(ilon(jpngr),ilat(jpngr))

      ! is data actually available?
      if (pred(jpngr)==ncfillvalue) then
        pred(jpngr) = dummy
        !grid(jpngr)%dogridcell = .false.
      end if

    end do

    ! deallocate memory again (the problem is that climate input files are of unequal length in the record dimension)
    deallocate( pred_arr )
    deallocate( lon_arr )
    deallocate( lat_arr )

    return

  end subroutine get_preds_nc_byvar


  subroutine initio_nc_nimpl()
    !////////////////////////////////////////////////////////////////
    ! Initialises module-specific NetCDF output files.
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: init_nc_3D_time, check
    
    ! local variables
    character(len=256) :: prefix

    character(len=*), parameter :: TITLE = "SOFUN GP-model output, module md_nimpl"
    !character(len=12) :: coef_age_bp_char

    integer :: jpngr, doy

    character(len=4)  :: year_char
    write(year_char,999) interface%steering%outyear

    ! convert parameter values to charaters
    !write(coef_age_bp_char,888) coef_nimpl%age_bp
    ! xxx complement. this is like as.character()

    prefix = "./output_nc/"//trim(interface%params_siml%runname)

    if ( .not. interface%steering%spinup ) then
      !----------------------------------------------------------------
      ! Annual NPP output file (anpp)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_anpp = trim(prefix)//'.'//year_char//".a.npp.nc"
        print*,'initialising ', trim(ncoutfilnam_anpp), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_anpp), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = NPP_NAME, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual net primary productivivty", &
                          title    = TITLE & 
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if
      !----------------------------------------------------------------
      ! Annual ANPP output file (aanpp)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_aanpp = trim(prefix)//'.'//year_char//".a.anpp.nc"
        print*,'initialising ', trim(ncoutfilnam_aanpp), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_aanpp), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = ANPP_NAME, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual aboveground net primary productivivty", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )      
      end if
      !----------------------------------------------------------------
      ! Annual BNPP output file (abnpp)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_abnpp = trim(prefix)//'.'//year_char//".a.bnpp.nc"
        print*,'initialising ', trim(ncoutfilnam_abnpp), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_abnpp), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = BNPP_NAME, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual belowground net primary productivivty", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )      
      end if
      !----------------------------------------------------------------
      ! Annual leaf NPP output file (alnpp)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_alnpp = trim(prefix)//'.'//year_char//".a.lnpp.nc"
        print*,'initialising ', trim(ncoutfilnam_alnpp), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_alnpp), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = LNPP_NAME, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual leaf net primary productivivty", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )   
      end if
      !----------------------------------------------------------------
      ! Annual wood NPP output file (awnpp)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_awnpp = trim(prefix)//'.'//year_char//".a.wnpp.nc"
        print*,'initialising ', trim(ncoutfilnam_awnpp), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_awnpp), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = WNPP_NAME, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual wood net primary productivivty", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if
      !----------------------------------------------------------------
      ! Annual leaf C/N output file (leafcn)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_leafcn = trim(prefix)//'.'//year_char//".a.leafcn.nc"
        print*,'initialising ', trim(ncoutfilnam_leafcn), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_leafcn), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = LEAFCN_NAME, &
                          varunits = "unitness", &
                          longnam  = "leaf carbon to nitrogen ratio", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )                              
      end if
      !----------------------------------------------------------------
      ! Annual leaf nitrogen efficiency (nre)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_nre = trim(prefix)//'.'//year_char//".a.nre.nc"
        print*,'initialising ', trim(ncoutfilnam_nre), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_nre), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = NRE_NAME, &
                          varunits = "unitness", &
                          longnam  = "leaf nitrogen resorption efficiency", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )                              
      end if
      !----------------------------------------------------------------
      ! Annual leaf nitrogen flux output file (lnf)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_lnf = trim(prefix)//'.'//year_char//".a.lnf.nc"
        print*,'initialising ', trim(ncoutfilnam_lnf), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_lnf), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = LNF_NAME, &
                          varunits = "gN m-2 yr-1", &
                          longnam  = "Nitrogen flux at leaf compartment", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )   
      end if
      !----------------------------------------------------------------
      ! Annual wood nitrogen flux output file (wnf)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_wnf = trim(prefix)//'.'//year_char//".a.wnf.nc"
        print*,'initialising ', trim(ncoutfilnam_wnf), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_wnf), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = WNF_NAME, &
                          varunits = "gN m-2 yr-1", &
                          longnam  = "Nitrogen flux at wood compartment", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )   
      end if
      !----------------------------------------------------------------
      ! Annual belowground nitrogen flux output file (bnf)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_bnf = trim(prefix)//'.'//year_char//".a.bnf.nc"
        print*,'initialising ', trim(ncoutfilnam_bnf), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_bnf), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = BNF_NAME, &
                          varunits = "gN m-2 yr-1", &
                          longnam  = "Nitrogen flux at belowground compartment", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if
      !----------------------------------------------------------------
      ! Annual N uptake in ecosystem (nuptake)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_nuptake = trim(prefix)//'.'//year_char//".a.nuptake.nc"
        print*,'initialising ', trim(ncoutfilnam_nuptake), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_nuptake), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = nuptake_NAME, &
                          varunits = "gN m-2 yr-1", &
                          longnam  = "Annual nitrogen uptake in ecosystem", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if
      !----------------------------------------------------------------
      ! Annual vcmax25 (annualvcmax25)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_annualvcmax25 = trim(prefix)//'.'//year_char//".a.annualvcmax25.nc"
        print*,'initialising ', trim(ncoutfilnam_annualvcmax25), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_annualvcmax25), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = ANNUALVCMAX25_NAME, &
                          varunits = "umol m-2 s-1", &
                          longnam  = "Annual maximum carboxylation capacity", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if
      !Grassland models
      !----------------------------------------------------------------
      ! Annual NPP output file in grassland (anpp_grass)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_anpp_grass = trim(prefix)//'.'//year_char//".a.npp_grass.nc"
        print*,'initialising ', trim(ncoutfilnam_anpp_grass), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_anpp_grass), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = NPP_NAME_GRASS, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual net primary productivivty in grassland", &
                          title    = TITLE & 
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if
      !----------------------------------------------------------------
      ! Annual ANPP output file in grassland (aanpp_grass)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_aanpp_grass = trim(prefix)//'.'//year_char//".a.anpp_grass.nc"
        print*,'initialising ', trim(ncoutfilnam_aanpp_grass), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_aanpp_grass), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = ANPP_NAME_GRASS, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual aboveground net primary productivivty in grassland", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )      
      end if
      !----------------------------------------------------------------
      ! Annual BNPP output file in grassland (abnpp_grass)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_abnpp_grass = trim(prefix)//'.'//year_char//".a.bnpp_grass.nc"
        print*,'initialising ', trim(ncoutfilnam_abnpp_grass), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_abnpp_grass), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = BNPP_NAME_GRASS, &
                          varunits = "gC m-2 yr-1", &
                          longnam  = "annual belowground net primary productivivty in grassland", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )      
      end if
      !----------------------------------------------------------------
      ! Annual leaf nitrogen flux output file in grassland  (lnf_grass)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_lnf_grass = trim(prefix)//'.'//year_char//".a.lnf_grass.nc"
        print*,'initialising ', trim(ncoutfilnam_lnf_grass), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_lnf_grass), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = LNF_NAME_GRASS, &
                          varunits = "gN m-2 yr-1", &
                          longnam  = "Nitrogen flux at leaf compartment in grassland", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )   
      end if
      !----------------------------------------------------------------
      ! Annual belowground nitrogen flux output file in grassland (bnf_grass)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_bnf_grass = trim(prefix)//'.'//year_char//".a.bnf_grass.nc"
        print*,'initialising ', trim(ncoutfilnam_bnf_grass), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_bnf_grass), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = BNF_NAME_GRASS, &
                          varunits = "gN m-2 yr-1", &
                          longnam  = "Nitrogen flux at belowground compartment in grassland", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if
      !----------------------------------------------------------------
      ! Annual N uptake in ecosystem in grassland (nuptake_grass)
      !----------------------------------------------------------------
      if (interface%params_siml%loutnimpl) then
        ncoutfilnam_nuptake_grass = trim(prefix)//'.'//year_char//".a.nuptake_grass.nc"
        print*,'initialising ', trim(ncoutfilnam_nuptake_grass), '...'
        call init_nc_3D_time(  filnam  = trim(ncoutfilnam_nuptake_grass), &
                          nlon     = interface%domaininfo%nlon, &
                          nlat     = interface%domaininfo%nlat, &
                          lon      = interface%domaininfo%lon, &
                          lat      = interface%domaininfo%lat, &
                          outyear  = interface%steering%outyear, &
                          outdt    = 365, &
                          outnt    = 1, &
                          varnam   = nuptake_NAME_GRASS, &
                          varunits = "gN m-2 yr-1", &
                          longnam  = "Annual nitrogen uptake in ecosystem in grassland", &
                          title    = TITLE &
                          !globatt1_nam = "coef_age_bp", globatt1_val = coef_age_bp_char &
                          ! XXX add more attributes XXX
                          )
      end if      
    end if

    888  format (F12.6)  ! for numeric (decimals)
    999  format (I4.4)   ! for integers
    
  end subroutine initio_nc_nimpl


  subroutine initoutput_nimpl( ngridcells )
    !////////////////////////////////////////////////////////////////
    ! Initialises module-specific output variables
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: ngridcells

    ! annual
    if (interface%params_siml%loutnimpl) then

      if (interface%steering%init) then
        !Forest
        allocate( outanpp(ngridcells) )
        allocate( outaanpp(ngridcells) )
        allocate( outabnpp(ngridcells) )
        allocate( outalnpp(ngridcells) )
        allocate( outawnpp(ngridcells) )
        allocate( outleafcn(ngridcells) )
        allocate( outnre(ngridcells) )
        allocate( outalnf(ngridcells) )
        allocate( outawnf(ngridcells) )
        allocate( outabnf(ngridcells) )
        allocate( outnuptake(ngridcells) )
        allocate( outavcmax25(ngridcells) )
        !Grassland
        allocate( outanpp_grass(ngridcells) )  
        allocate( outaanpp_grass(ngridcells) )        
        allocate( outabnpp_grass(ngridcells) )
        allocate( outalnf_grass(ngridcells) )
        allocate( outabnf_grass(ngridcells) )
        allocate( outnuptake_grass(ngridcells) )
      end if
      !Forest
      outanpp(:) = 0.0
      outaanpp(:) = 0.0
      outabnpp(:) = 0.0
      outalnpp(:) = 0.0
      outawnpp(:) = 0.0
      outleafcn(:) = 0.0
      outnre(:) = 0.0
      outalnf(:) = 0.0
      outawnf(:) = 0.0
      outabnf(:) = 0.0
      outnuptake(:) = 0.0
      outavcmax25(:) = 0.0
      !Grassland
      outanpp_grass(:) = 0.0
      outaanpp_grass(:) = 0.0
      outabnpp_grass(:) = 0.0
      outalnf_grass(:) = 0.0
      outabnf_grass(:) = 0.0
      outnuptake_grass(:) = 0.0    
    end if

  end subroutine initoutput_nimpl


  subroutine getout_annual_nimpl( jpngr, tile, tile_fluxes )
    !////////////////////////////////////////////////////////////////
    ! Called once a year to gather annual output variables.
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !----------------------------------------------------------------
    use md_interface

    ! arguments
    integer, intent(in) :: jpngr
    type(tile_type), dimension(nlu), intent(in) :: tile
    !type(tile_fluxes_type), dimension(nlu), intent(in) :: tile_fluxes
    type(tile_fluxes_type), dimension(nlu), intent(in) :: tile_fluxes

    ! local variables
    integer :: pft, lu

    ! outanrlarea(jpngr) = anrlarea
    if (interface%params_siml%loutnimpl) then

      lu = 1
      if ( abs(sum(tile(lu)%plant(:)%fpc_grid) - 1.0) > eps ) stop 'getout_annual_nimpl(): fpc_grid does not sum up to 1.0'
      !print*,'jpngr',jpngr
      !print*,'sum(tile(lu)%plant(:)%fpc_grid)',sum(tile(lu)%plant(:)%fpc_grid)
      ! grid-cell average across PFTs, weighted by their FPC grid 
      !Forest
      outanpp(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%anpp * 1.0 ) ! because npp was derived from gpp, while gpp has already been accounted by fpc_grid, there no need to multiply again!
      outaanpp(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%aanpp * 1.0 )
      outabnpp(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%abnpp * 1.0)
      outalnpp(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%alnpp * 1.0)
      outawnpp(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%awnpp * 1.0)
      outleafcn(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%leafcn * 0.5) !Because here leafcn was derived from canopy vcmax25 (has been applied by weighted-sum method already). So c3 and c4 leaf c/n is the same here, so just halve them here.
      outnre(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%nre * 0.5)! There is no difference (classfication) between c3 and c4! They are the same value of c3 and c4, so let's just half-weight them.
      outalnf(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%alnf * 1.0)
      outawnf(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%awnf * 1.0)
      outabnf(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%abnf * 1.0)
      outnuptake(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%nuptake * 1.0 )
      outavcmax25(jpngr) = sum( tile_fluxes(lu)%plant(:)%avcmax25_max * tile(lu)%plant(:)%fpc_grid) !same above, in the end this result should be = canopy vcmax25
      !Grassland
      outanpp_grass(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%anpp_grass * 1.0)
      outaanpp_grass(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%aanpp_grass * 1.0)
      outabnpp_grass(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%abnpp_grass * 1.0)
      outalnf_grass(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%alnf_grass * 1.0)
      outabnf_grass(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%abnf_grass * 1.0)
      outnuptake_grass(jpngr) = sum( tile_nimpl_fluxes(lu)%plant(:)%nuptake_grass * 1.0)
    end if

  end subroutine getout_annual_nimpl


  subroutine writeout_nc_nimpl()
    !/////////////////////////////////////////////////////////////////////////
    ! Writes module-specific NetCDF output
    !
    ! This is designed for use within SOFUN and requires arguments as
    ! derived-types, defined elsewhere. For other applications, implement 
    ! the function calls (e.g., calc_dgpp()) differently and 
    ! comment/delete this subroutine.
    !-------------------------------------------------------------------------
    use netcdf
    use md_io_netcdf, only: write_nc_2D, write_nc_3D_time, check
    
    if ( .not. interface%steering%spinup ) then
      !Forest model
      !-------------------------------------------------------------------------
      ! Annual NPP
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_anpp), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_anpp), &
                                                              NPP_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outanpp(:) &
                                                              )
    
      !-------------------------------------------------------------------------
      ! Annual ANPP
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_aanpp), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_aanpp), &
                                                              ANPP_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outaanpp(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual BNPP
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_abnpp), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_abnpp), &
                                                              BNPP_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outabnpp(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual LNPP
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_alnpp), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_alnpp), &
                                                              LNPP_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outalnpp(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual WNPP
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_awnpp), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_awnpp), &
                                                              WNPP_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outawnpp(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual leafcn
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_leafcn), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_leafcn), &
                                                              LEAFCN_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outleafcn(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual nre
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_nre), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_nre), &
                                                              NRE_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outnre(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual lnf
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_lnf), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_lnf), &
                                                              LNF_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outalnf(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual wnf
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_wnf), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_wnf), &
                                                              WNF_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outawnf(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual bnf
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_bnf), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_bnf), &
                                                              BNF_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outabnf(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual nuptake
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_nuptake), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_nuptake), &
                                                              nuptake_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outnuptake(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual annualvcmax25 --> not used anymore
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_annualvcmax25), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_annualvcmax25), &
                                                              ANNUALVCMAX25_NAME, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outavcmax25(:) &
                                                              )
      !Grassland model
      !-------------------------------------------------------------------------
      ! Annual NPP in grassland
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_anpp_grass), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_anpp_grass), &
                                                              NPP_NAME_GRASS, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outanpp_grass(:) &
                                                              )
    
      !-------------------------------------------------------------------------
      ! Annual ANPP in grassland
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_aanpp_grass), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_aanpp_grass), &
                                                              ANPP_NAME_GRASS, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outaanpp_grass(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual BNPP in grassland
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_abnpp_grass), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_abnpp_grass), &
                                                              BNPP_NAME_GRASS, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outabnpp_grass(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual lnf in grassland
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_lnf_grass), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_lnf_grass), &
                                                              LNF_NAME_GRASS, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outalnf_grass(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual bnf in grassland
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_bnf_grass), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_bnf_grass), &
                                                              BNF_NAME_GRASS, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outabnf_grass(:) &
                                                              )
      !-------------------------------------------------------------------------
      ! Annual nuptake in grassland
      !-------------------------------------------------------------------------
      if (interface%params_siml%loutnimpl) print*,'writing ', trim(ncoutfilnam_nuptake_grass), '...'
      if (interface%params_siml%loutnimpl) call write_nc_2D( trim(ncoutfilnam_nuptake_grass), &
                                                              nuptake_NAME_GRASS, &
                                                              interface%domaininfo%maxgrid, &
                                                              interface%domaininfo%nlon, &
                                                              interface%domaininfo%nlat, &
                                                              interface%grid(:)%ilon, &
                                                              interface%grid(:)%ilat, &
                                                              interface%grid(:)%dogridcell, &
                                                              outnuptake_grass(:) &
                                                              )
    end if

  end subroutine writeout_nc_nimpl


end module md_nuptake_impl
