!//////////////////////////////////////////////////////////////////////////////
!
!     SUBROUTINE WATERBALANCE
!     Calculations of soil hydrology, plant water balance, daily canopy
!     conductance and PFT daily water stress factor
!
!     instead of applying this model:
!     use latent heat flux data from FLUXNET station
!     or simple 1-bucket model, different implementations in Tyler's, 
!       Wang Han's, etc. code     
!
!//////////////////////////////////////////////////////////////////////////////

      subroutine waterbalance(init,jpngr,d,m,icount,rootprop,dgp,dgc,
     $     kperc,awc,drunoff_drain,drunoff_surf,dwscal,daet,intercp_tot,
     $     daep,daept,dmoss_cap)

      IMPLICIT NONE

      real, dimension(nbucket,nlu) :: aettotal  !  active total transpiration from layers 1, 2


# 191
!     -------------------------------------------------------------------------
!     Initialisations
!     -------------------------------------------------------------------------
      aettotal(1,:)=0.d0
      aettotal(2,:)=0.d0
      drunoff_surf(:)=0.d0
      daet(:)=0.d0
      perc=0.d0
      qdrain=0.0d0
      intercp_tot(:) = 0.0d0
      cover(:) = 0.0d0 

# 319
!     We assume that fractional water holding capacity exhibits only minor
!     differences between lower and upper bucket

!     -------------------------------------------------------------------------
!     START PFT LOOP
!     -------------------------------------------------------------------------

      do pft=1,npft

!     -------------------------------------------------------------------------
!     CALCULATE ACTUAL CANOPY CONDUCTANCE, POTENTIAL WATER SCALAR,
!     INTERCEPTION LOSS AND ACTUAL EVAPOTRANSPIRATION (AET) FOR EACH PFT
!     -------------------------------------------------------------------------
        if (present(pft,jpngr)) then
          lu=lu_category(pft)

!     -------------------------------------------------------------------------
!     INTRODUCE the calculation of the available water here
!     -------------------------------------------------------------------------
!     Calculate effective supply function in root zone today
!     Eqn 24, Haxeltine & Prentice 1996
!     rootprop(2,pft): fraction of active fraction of roots uptaking water 
!     from lower soil layer
!     rootprop(1,pft): the remainder in the top layer

          wr=wpool(1,lu,jpngr)*rootprop(1,pft)
     $         +wpool(2,lu,jpngr)*rootprop(2,pft)
          

# 359
!     -------------------------------------------------------------------------
!     SUPPLY
!     -------------------------------------------------------------------------
          
!     The calculation of 'supply' has changed from Gerten et al. (2004)
!     Eq. 6 to what was in the webversion (F77), but neither Sitch et
!     al. 2003 nor Gerten et al. 2004 mention it. LPJ-GUESS uses a
!     similar Eqn.  
          
!     R. Wania: Supply should be scaled by FPC, because Emax is in mm/d
!     which can be seen as L/m2/day. Emax is the maximum transpiration
!     per grid cell. If a PFT has an FPC of only 0.5, then the supply
!     rate will be only half.

# 374
          supply=emax(pft)*wr
# 378
!          print*,'supply',supply
!          print*,'emax',emax(pft)
!          print*,'wr',wr
!          print*,'dphen',dphen(d,pft)
!          print*,'fpc_grid',fpc_grid(pft,jpngr) 

!     -------------------------------------------------------------------------
!     Calculate PFT-specific INTERCEPTION LOSS (Kergoat 1998)
!     -------------------------------------------------------------------------
# 388
!     no interception
          int_store(pft) = 0.0d0
          wet = 0.0d0
          intercep(pft) = 0.0d0

# 403
!     -------------------------------------------------------------------------
!     DEMAND
!     -------------------------------------------------------------------------

!     Calculate AET demand function and potential demand 
!     assuming full leaf cover (Eqn 23, Haxeltine & Prentice 1996)
!     Note fraction of daylength
!     The demand is calculated differently by Gerten than by Sitch.
!     Sitch uses *(1-exp(-dgp/gm)),
!     whereas Gerten uses /(1+gm/dgp), which leads to a slower
!     increase in the potential demand. 
!     This is Gerten's version:
!               demand=(1.0-wet)*dpet(d)*alpha/
!     *              (1+gm/(dgp(d,pft)*dphen(d,pft)))

!     The problem with both, Gerten and Sitch's version, is that the
!     scaling of demand by fpc_grid is implicit in dgp. In this
!     equation, dgp influences demand in a non-linear way however. This
!     causes the following problem:
!     Assume we have 2 PFTs with the same physiological parameters.
!     We allocate some numbers to the variables in the eqn.
!     wet = 0.0
!     dpet(d) = 10
!     alpham = 1.391
!     gm = 3.26
!     non-fpc-scaled dgp = 5
!     dphen = 1
c
!     Situation 1:
!     PFT1: fpc_grid = 0.9
!     PFT2: fpc_grid = 0.1
!     demand(pft1) = (1-0)*10*1.391/(1+3.26/5*0.9*1) = 8.1
!     demand(pft2) = (1-0)*10*1.391/(1+3.26/5*0.1*1) = 1.8
!     total demand = 9.9

!     Situation 2:
!     PFT1: fpc_grid = 0.6
!     PFT2: fpc_grid = 0.4
!     demand(pft1) = (1-0)*10*1.391/(1+3.26/5*0.6*1) = 6.7
!     demand(pft2) = (1-0)*10*1.391/(1+3.26/5*0.4*1) = 5.3
!     total demand = 12.0
!     
!     The total demand differs because of a different balance of FPCs.
!     A suggestion to solve this by Rita Wania (19/12/2008)
!     1. Take out the scaling by fpc_grid from dgp
!     2. Scale the entire demand equation by fpc_grid.
!     -------------------------------------------------------------------------
# 451
          demand=dpet(d)*alpham*(1.0d0-
     *         exp(-dgp(d,pft)*dphen(d,pft)/gm))
          demandpot=dpet(d)*alpham*(1.0d0-exp(-dgp(d,pft)/gm))
          
!     Calculate daily potential water scalar (= water stress factor)
          if (demand.gt.0.0d0) then !note that always demandpot >= demand
            dwscal(pft)=min(supply/demandpot,1.d0)
          else
            dwscal(pft)=1.0d0
          endif

# 493
!     Calculate actual canopy conductance (gc), according to balance
!     between supply and demand Eqn 25, Haxeltine & Prentice 1996
!     R. Wania: Note that aet can exceed dpet because dpet is
!     equilibrium evaporation, not potential evaporation!
!     -------------------------------------------------------------------------

            if(supply.ge.demand) then
              dgc(d,pft)=dgp(d,pft)*dphen(d,pft) 
            else 
              if(dpet(d).gt.0.0d0) then
!     version before bug correction r.2650 by S. Keel
c                dgc(d,pft)=MAX(0.0d0,(gm*supply/((1.0d0-wet)
c     $               * dpet(d)*alpham-supply)))
!     version after bug correction r.2650 by S. Keel
                dgc(d,pft)=MAX(0.0d0,(gm*fpc_grid(pft,jpngr)
     $                *dphen(d,pft)*supply/((1.0d0-wet)
     $               * dpet(d)*alpham*fpc_grid(pft,jpngr)*dphen(d,pft)
     $                -supply)))
              else
                dgc(d,pft)=0.0d0
              endif
            endif
            aet(pft)=MIN(supply,demand)  
            
!            print*,'supply,demand',supply,demand

!     -------------------------------------------------------------------------
!     AET
!     -------------------------------------------------------------------------

# 528
            
!     Accumulate total AET and interception loss
!     -------------------------------------------------------------------------
            if (wr.eq.0.0d0) then
              beta_w(1)=0.0d0
              beta_w(2)=0.0d0
            else
              beta_w(1)=rootprop(1,pft)*wpool(1,lu,jpngr)/wr
              beta_w(2)=rootprop(2,pft)*wpool(2,lu,jpngr)/wr
            endif
            
            aettotal(1,lu)=aettotal(1,lu)+beta_w(1)*aet(pft)
            aettotal(2,lu)=aettotal(2,lu)+beta_w(2)*aet(pft)
            intercp_tot(lu)=intercp_tot(lu)+intercep(pft)
     $           * fpc_grid(pft,jpngr)
            

!     Calculate area covered by PFTs
!     -------------------------------------------------------------------------
            cover(lu)=cover(lu)+(fpc_grid(pft,jpngr)*dphen(d,pft))


!     Compute total transpiration, output variable for gpp
!     -------------------------------------------------------------------------
!     identical to daet=aettotal(1,lu)+aettotal(2,lu) sum over lu classes

            daet(lu)=daet(lu)+aet(pft)   
          
          endif                 !present
      enddo                     !pft
      
!     ---------------
!     END OF PFT LOOP
!     ---------------
      

!     Compute total transpiration and effective rainfall
!     and meltwater from snow
!     ---------------------------------------------------
      do lu=1,nlu
# 571
          prec(lu)=dprec(d)-intercp_tot(lu) + dmelt(d) !dmelt regarded as lu independent
# 577
      enddo                     ! nlu




!     =========================== Simple hydrology ============================
!     kuno 2005-02-24:	split water balance for landuse

!     calculate daily percolation from soil layer 1 to 2
!     Eqn 31, Haxeltine & Prentice 1996

      do lu=1,nlu

          perc=kperc(1)*wpool(1,lu,jpngr)**kperc(2)
          daep(lu)=0.0d0               ! no evaporation from soil
          qdrain = 0.0d0               ! no drainage into sub-soil

!     Update soil water content
!     Eqns 26,27, Haxeltine & Prentice 1996
          wpool(1,lu,jpngr) = wpool(1,lu,jpngr)
     $         +(
     $            prec(lu)             ! precipitation + snow melt - interception
     $           -perc                 ! percolation from soil layers 1 -> 2
     $           -aettotal(1,lu)       ! active total transpiration
     $           -daep(lu)             ! daily actual evaporation from soil
     $           )/awc(1)

          wpool(2,lu,jpngr) = wpool(2,lu,jpngr)
     $         +(
     $           perc                 ! percolation from soil layers 1 -> 2
     $          -qdrain               ! drainage into sub-soil
     $          -aettotal(2,lu)       ! active total transpiration
     $          )/awc(2)
          
          if (wpool(1,lu,jpngr).le.1.0d-8) wpool(1,lu,jpngr)=0.d0
          if (wpool(2,lu,jpngr).le.1.0d-8) wpool(2,lu,jpngr)=0.d0
          

!     Transfer excess soil water to runoff
          if (wpool(2,lu,jpngr).ge.1.d0) then
            drunoff_drain(lu)=((wpool(2,lu,jpngr)-1.d0)*awc(2))
            wpool(2,lu,jpngr)=1.0d0
          endif
          
          if (wpool(1,lu,jpngr).ge.1.0d0) then
            drunoff_surf(lu)=(wpool(1,lu,jpngr)-1.0d0)*awc(1)
            wpool(1,lu,jpngr)=1.0d0
          endif

      enddo                     !lu



!     ------------------------------------------------------------------
!     N-LEACHING (changed from soilwater.cpp:226)
!     N-leached is now the drainage plus surface runoff as a fraction of
!     total soil water (both given in mm)
!     ------------------------------------------------------------------
      do lu=1,nlu
        f_leach(d,lu) = min(
     $       1.0d0,
     $       max(
     $       0.0d0,
     $       1.d0-((awc(1)+awc(2))/(awc(1)+awc(2)
     $       +drunoff_drain(lu)+drunoff_surf(lu)))
     $       )
     $       )
      enddo

# 1479
      
      return
      end




