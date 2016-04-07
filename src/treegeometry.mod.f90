module md_vegdynamics
  !////////////////////////////////////////////////////////////////
  ! Module contains tree geometry variables and functions.
  ! Variables have 'jpngr' dimension, thus have "memory" from one
  ! year to the next like pool variables.
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: npft, maxgrid
  implicit none

  logical, dimension(npft,maxgrid) :: ispresent   ! boolean whether PFT is present
  real, dimension(npft,maxgrid)    :: height      ! tree height (m)
  real, dimension(npft,maxgrid)    :: crownarea   ! individual's tree crown area
  real, dimension(npft,maxgrid)    :: nind        ! number of individuals [1/m2]
  real, dimension(npft,maxgrid)    :: lai_ind
  real, dimension(npft,maxgrid)    :: fpc_grid
  real, dimension(npft,maxgrid)    :: fapar_ind 

contains

  subroutine update_fpc( pft, jpngr )
    !////////////////////////////////////////////////////////////////
    ! Updates fpc_grid, fapar_ind, and lai_ind
    !----------------------------------------------------------------
    use md_pools, only: pleaf
    use md_params_modl, only: tree, islu, lu_category, grass
    
    implicit none

    ! ARGUMENTS
    integer, intent(in) :: pft, jpngr

    ! LOCAL VARIABLES
    ! logical :: distribute
    ! real :: share
    real :: cleaf_tot

    !-------------------------------------------------------------------------
    ! Non-linearity of Beer-Law causes very high FPC values when 2 Grasses are present.
    ! (Beer Law does NOT make sense for grasses, anyway.)
    ! Thus, use sum of all grass/moss-leaf masses and calculate FPC based on the sum.
    ! Then compute each PFT's FPC as the product of total-grass FPC times each PFT's leaf mass.
    !-------------------------------------------------------------------------
    lai_ind(pft,jpngr)  = get_lai_ind( pleaf(pft,jpngr)%c%c12, crownarea(pft,jpngr), pft )
    fapar_ind(pft,jpngr)  = get_fapar( lai_ind(pft,jpngr) )
    fpc_grid(pft,jpngr) = get_fpc_grid( crownarea(pft,jpngr), nind(pft,jpngr), fapar_ind(pft,jpngr) )    


    !    distribute = .false.
    !    if (grass(pft)) then
    !      cleaf_tot = sum( pleaf(:,jpngr)%c%c12, mask=(grass(:).and.islu(:,lu_category(pft))) )
    !      distribute = .true.
    !    else
    !      cleaf_tot = pleaf(pft,jpngr)%c%c12
    !    end if    

    !    xxx lai_ind calculation here : this does not work : pft is passed on from this pft, whereas sum of all pfts leaf mass is used xxx!    

    !    lai_ind(pft,jpngr)  = get_lai_ind( cleaf_tot, crownarea(pft,jpngr), pft )
    !    fapar_ind(pft,jpngr)  = get_fapar( lai_ind(pft,jpngr) )
    !    fpc_grid(pft,jpngr) = get_fpc_grid( crownarea(pft,jpngr), nind(pft,jpngr), fapar_ind(pft,jpngr) )    

    !    if ( distribute ) then
    !      if (cleaf_tot>0.0) then
    !        share = pleaf(pft,jpngr)%c%c12 / cleaf_tot
    !        lai_ind(pft,jpngr)  = lai_ind(pft,jpngr) * share
    !        fapar_ind(pft,jpngr)  = fapar_ind(pft,jpngr) * share
    !        fpc_grid(pft,jpngr) = fpc_grid(pft,jpngr) * share
    !      else
    !        lai_ind(pft,jpngr)  = 0.0
    !        fapar_ind(pft,jpngr)  = 0.0
    !        fpc_grid(pft,jpngr) = 0.0
    !      end if
    !    end if

    return

  end subroutine update_fpc


  function get_nind( tree, fpc_grid, crownarea, fapar_ind ) result( nind )
    !////////////////////////////////////////////////////////////////
    !----------------------------------------------------------------
    use md_pools, only: pleaf

    implicit none

    ! ARGUMENTS
    logical, intent(in) :: tree
    real,intent(in) :: fpc_grid
    real,intent(in) :: crownarea
    real,intent(in) :: fapar_ind

    ! FUNCTION RETURN VARIABLE
    real, intent(out) :: nind

    if (tree) then
      nind = fpc_grid / ( crownarea * fapar_ind )
    else
      nind = 1.0
    end if

  end function get_nind


  function get_fpc_grid( crownarea, nind, fapar_ind ) result( fpc_grid )
    !////////////////////////////////////////////////////////////////
    ! Function returns total fractional plant cover of a PFT
    ! Eq. 8 in Sitch et al., 2003
    !----------------------------------------------------------------

    ! arguments
    real, intent(in) :: crownarea
    real, intent(in) :: nind
    real, intent(in) :: fapar_ind

    ! function return variable
    real, intent(out) :: fpc_grid

    fpc_grid = crownarea * nind * fapar_ind

  end function get_fpc_grid


  function get_fapar( lai_ind ) result( fapar_ind )
    !////////////////////////////////////////////////////////////////
    ! Function returns fractional plant cover an individual
    ! Eq. 7 in Sitch et al., 2003
    !----------------------------------------------------------------
    use md_params_modl, only: kbeer

    ! arguments
    real, intent(in) :: lai_ind

    ! function return variable
    real, intent(out) :: fapar_ind

    fapar_ind = ( 1.0 - exp( -1.0 * kbeer * lai_ind) )

  end function get_fapar


  function get_lai_ind( pleaf, crownarea, pft ) result( lai_ind )
    !////////////////////////////////////////////////////////////////
    ! Function returns leaf area index of an individual
    ! Eq. 5 in Sitch et al., 2003
    !----------------------------------------------------------------
    use md_params_modl, only: sla
    implicit none

    ! arguments
    real, intent(in)    :: pleaf
    real, intent(in)    :: crownarea
    integer, intent(in) :: pft

    ! function return variable
    real :: lai_ind

#if _check_sanity
    if (crownarea<0.0) stop 'GET_LAI_IND: crownarea < 0'
#endif

    if (crownarea==0.0) then
      lai_ind = 0.0
    else
      lai_ind = pleaf * sla(pft) / crownarea
    endif

  end function get_lai_ind


  subroutine update_tree_geometry( pft, jpngr, psapw_temp, ntoc_sm )
    !////////////////////////////////////////////////////////////////
    ! Subroutine updating tree geometry.
    ! Requires global variables:
    !   - psapw, pwood
    ! Requires global model parameters:
    !   - wooddens, allom, crownarea_max, reinickerp, latosa, sla
    ! Updates global variables (geometry-related):
    !   - height, crownarea, psapw
    ! Summary of calculation:
    ! ( Cwood+Csap', Cleaf, Croot ) --> ( LA, SA, H, D, CA, Csap)
    ! Allometry is independent of sapwood/heartwood ratio,
    ! only sum of the two is used (VI)
    !----------------------------------------------------------------
    use md_classdefs
    use md_pools, only: pleaf, psapw, pwood
    use md_params_modl, only: wooddens, allom1, allom2, allom3, &
      crownarea_max, reinickerp, latosa, sla, pi, tree

    ! arguments
    integer, intent(in)          :: pft, jpngr
    type(orgpool), intent(in)    :: psapw_temp
    real, intent(in)             :: ntoc_sm

    ! local variables
    real    :: stemdiam

    if (tree(pft)) then

      ! (I)   LA = latosa * SA
      ! (II)  Cleaf = lmtorm * Croot
      ! (III) H = allom2 * D**allom3
      ! (IV)  CA = min (allom1 * D**reinickerp, crownarea_max)
      ! (V)   LA = Cleaf * sla
      ! (VI)  wooddens = ( Cwood + Csap' ) / ( (D/2)**2 * pi * H )

      ! 1. solve (VI) for D using (III)
      stemdiam = (4.0*(psapw_temp%c%c12+pwood(pft,jpngr)%c%c12)/ &
        wooddens(pft)/pi/allom2(pft))**(1.0/(2.0+allom3))

      ! 2. solve (III) for H
      height(pft,jpngr) = allom2(pft)*stemdiam**allom3

      ! 3. (IV) for CA
      crownarea(pft,jpngr) = min(crownarea_max(pft), allom1(pft)*stemdiam**reinickerp)

      ! 4. Recalculate sapwood mass, transferring excess sapwood to heartwood
      ! compartment, if necessary to satisfy (I)
      psapw(pft,jpngr)%c%c12 = pleaf(pft,jpngr)%c%c12*height(pft,jpngr)*wooddens(pft)*sla(pft)/latosa(pft)
      psapw(pft,jpngr)%n%n14 = psapw(pft,jpngr)%c%c12*ntoc_sm
      call orgcp( orgminus( psapw_temp, psapw(pft,jpngr)), pwood(pft,jpngr) )

    else

      height(pft,jpngr) = 0.0
      crownarea(pft,jpngr) = 1.0

    end if

  end subroutine update_tree_geometry
  

  function get_geometry_sapl( grass, lai0, sla, latosa, reinickerp, allom1, &
    allom2, allom3, wooddens, woodtosapw, lmtorm0, cton_pro )
    !////////////////////////////////////////////////////////////////////
    ! Function to determine sapling tree geometry. Returns an array 
    ! containing (lm_sapl, sm_sapl, hm_sapl), given geometry parameters.
    !--------------------------------------------------------------------
    use md_classdefs, only: orgpool
    use md_params_core, only: pi
    implicit none

    ! arguments
    logical, intent(in) :: grass
    real, intent(in) :: lai0                        ! sapling (or grass on ((interface%steering%init))ialisation) LAI (=lai_sapl)
    real, intent(in) :: sla                         ! specific leaf area [m2/gC]
    real, intent(in) :: latosa                      ! ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b) (m2 * m-2)
    real, intent(in) :: reinickerp                  ! Reinicker-p for geometry
    real, intent(in) :: allom1, allom2, allom3  ! allometry parameters
    real, intent(in) :: wooddens                    ! wood density (gC * m-3)
    real, intent(in) :: woodtosapw                  ! sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    real, intent(in) :: lmtorm0                     ! leaf to root ratio under non-water stressed conditionss
    real, intent(in) :: cton_pro                    ! C:N ratio of new production

    ! local variables
    real :: x
    type(orgpool) :: lm_sapl
    type(orgpool) :: sm_sapl
    type(orgpool) :: hm_sapl
    type(orgpool) :: rm_sapl
    real :: stemdiam
    real :: height_sapl

    ! function return variable
    type(orgpool), dimension(4) :: get_geometry_sapl


    if (grass) then
      !///////////////////////////////////////////////////////
      ! GRASS "SAPLING" GEOMETRY
      !------------------------------------------------------- 
      lm_sapl%c%c12 = lai0 / sla
      sm_sapl%c%c12 = 0.0
      hm_sapl%c%c12 = 0.0

    else
      !///////////////////////////////////////////////////////
      ! TREE SAPLING GEOMETRY
      !------------------------------------------------------- 

      ! Calculate leafmass for a sapling individual
      !  (1) lai = leafmass * sla / (crown area)
      !  (2) (leaf area) = latosa * (sapwood xs area)
      !         (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
      !  (3) (crown area) = allom1 * (stem diameter) ** reinickerp
      !         (Reinickes theory)
      ! From (1),
      !  (4) leafmass = lai * (crown area) / sla
      ! From (1) & (3),
      !  (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
      ! From (2),
      !  (6) leafmass = latosa * (sapwood xs area) / sla
      !  (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
      ! From (6) and (7),
      !  (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
      ! From (8),
      !  (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
      ! (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
      ! Define x,
      ! (11) x = [ (sapwood diameter)+(heartwood diameter) ] / 
      !          (sapwood diameter)
      ! From (10) & (11),
      ! (12) (stem diameter) = x * (sapwood diameter)
      ! From (5), (9) & (12),
      ! (13) leafmass = lai * allom1 * x**reinickerp * 
      !               (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
      ! From (13),
      ! (14) leafmass = [ lai * allom1 * x**reinickerp *
      !      (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(2/(2-reinickerp))
      x = woodtosapw
       
      lm_sapl%c%c12 = (lai0*allom1*x**reinickerp*  &
        (4.0*sla/pi/latosa)**(reinickerp*0.5) / &
        sla)**(2.0/(2.0-reinickerp)) !eqn 14

      ! Calculate sapling stem diameter
      ! From (9) & (12),
      ! (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5
      stemdiam = x*(4.0*lm_sapl%c%c12*sla/pi/latosa)**0.5 !Eqn 15
    
      !     Calculate sapling height
      ! (16) height = allom2 * (stem diameter)**allom3 (source?)
      height_sapl = allom2*stemdiam**allom3 !Eqn 16
    
      ! Calculate sapling sapwood mass
      ! (17) (sapwood volume) = height * (sapwood xs area)
      ! (18) (sapwood xs area) = leafmass * sla / latosa
      ! From (17) & (18),
      ! (19) (sapwood volume) = height * leafmass * sla / latosa
      ! (20) (sapwood mass) = (wood density) * (sapwood volume)
      ! From (19) & (20),
      ! (21) (sapwood mass) = (wood density) * height * leafmass * sla /
      ! latosa
      sm_sapl%c%c12 = wooddens*height_sapl*lm_sapl%c%c12*sla/latosa !Eqn 21

      ! Calculate sapling heartwood mass
      !     From (11),
      ! (22) (heartwood mass) = (x-1) * (sapwood mass)
      hm_sapl%c%c12 = (x-1.0)*sm_sapl%c%c12 !Eqn 22

    end if

    ! root mass of sapling is a function of leaf mass            
    rm_sapl%c%c12 = ( 1.0 / lmtorm0 ) * lm_sapl%c%c12

    ! N content of sapling
    lm_sapl%n%n14 = lm_sapl%c%c12 / cton_pro
    rm_sapl%n%n14 = rm_sapl%c%c12 / cton_pro 
    sm_sapl%n%n14 = sm_sapl%c%c12 / cton_pro
    hm_sapl%n%n14 = hm_sapl%c%c12 / cton_pro

    get_geometry_sapl = (/lm_sapl, sm_sapl, hm_sapl, rm_sapl/)

  end function get_geometry_sapl

end module md_vegdynamics

