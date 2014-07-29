c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE ALLOCATION
c     Allocation of annual C increment to leaf, stem and fine root
c     compartments, update of individual structure

      subroutine allocation(jpngr,year)

      implicit none
#include "para.inc"
#include "modelpara.inc"
#include "fluxLU.inc"
#include "climateLU.inc"
#include "poolsLU.inc"
#include "pftparLU.inc"

c     ARGUMENTS
      INTEGER jpngr, year
      
c     PARAMETERS
      INTEGER nseg
      PARAMETER(nseg=20)
      INTEGER jmax
      PARAMETER(jmax=40)        !max number of iterations in search for alloc soln

c     LOCAL VARIABLES
      INTEGER pft
      INTEGER j                 !counter in bisection loop

      REAL*8 lminc_ind          !individual leafmass increment this year
      REAL*8 sminc_ind          !individual sapmass increment this year
      REAL*8 rminc_ind          !individual fineroot mass increment this year
      REAL*8 bm_inc_ind         !individual total biomass increment this year
      REAL*8 sap_xsa            !cross sectional area of sapwood  
      REAL*8 stemdiam           !stem diameter 
      REAL*8 x1,x2,rtbis,dx,xmid,sign !working vars in bisection
      REAL*8 fx1,fmid
      REAL*8 lminc_ind_min      !min leafmass increment to maintain current sapwood
      REAL*8 rminc_ind_min      !min rootmass increment to support new leafmass
      REAL*8 lmtorm             !ratio of leafmass to fine rootmass
      REAL*8 fpc_ind            !individual FPC
      REAL*8 crownarea_max      !maximum crown area (m2)
      REAL*8 allom_fact         !pft dependent factor in allom eq.
      REAL*8 allom_exp          !pft dependent exponent in allom eq.
      REAL*8 tree_fact          !pft dependent tree properties in allom eq.
#if DyN
      REAL*8 plantC             !sum of C in all plant compartments
      REAL*8 plantN             !sum of N in all plant compartments
      REAL*8 f_cl               !Partitioning coefficient of C allocation to leafs
      REAL*8 f_cr               !Partitioning coefficient of C allocation to roots
      REAL*8 f_cs               !Partitioning coefficient of C allocation to sapwood
      REAL*8 f_nl               !Partitioning coefficient of N allocation to leafs
      REAL*8 f_nr               !Partitioning coefficient of N allocation to roots
      REAL*8 f_ns               !Partitioning coefficient of N allocation to sapwood
      LOGICAL test
      REAL*8 temp
#endif

c     Loop over PFTs
c     -------------------------------------------------------------------------
      do pft=1,npft
         
        if (present(pft,jpngr)) then
          
          if(nind(pft,jpngr).eq.0) stop
          bm_inc_ind=bm_inc(pft,jpngr,1)/nind(pft,jpngr)
          
c     calculate this year's leaf to fine root mass ratio from mean annual 
c     water scalar and pft specific parameter
          
          lmtorm=pftpar(pft,18)*wscal(pft,jpngr)
                    
          if (tree(pft)) then
c     -------------------------------------------------------------------------
c     TREE ALLOCATION
c     -------------------------------------------------------------------------
            
            crownarea_max=pftpar(pft,20)
            
c     Work out minimum leaf production to maintain current sapmass
            
            lminc_ind_min=latosa(pft)*sm_ind(pft,jpngr,1)/
     $           (wooddens(pft)*height(pft,jpngr)*sla(pft))-
     $           lm_ind(pft,jpngr,1) !eqn (27)
            
            rminc_ind_min=latosa(pft)*sm_ind(pft,jpngr,1)/
     $           (wooddens(pft)*height(pft,jpngr)*
     $           sla(pft)*lmtorm)-rm_ind(pft,jpngr,1) !eqn (30)
            
            if (rminc_ind_min.gt.0.0d0
     $           .and.lminc_ind_min.gt.0.0d0
     $           .and.rminc_ind_min+lminc_ind_min.le.bm_inc_ind) then !(see growth.cpp:535)
c     Normal allocation (positive increment to all living C
c     compartments)
              x1=0.0d0
              x2=(bm_inc_ind-(lm_ind(pft,jpngr,1)/
     $             lmtorm-rm_ind(pft,jpngr,1)))/
     $             (1.0d0+1.0d0/lmtorm)
              dx=(x2-x1)/dble(nseg)

c     allocate constants in allom eq. for simplification and numerical stability
              allom_fact = allom2(pft)**(2.0d0/allom3)           
              allom_exp = dnint(1.0d0+2.0d0/allom3) !must be a whole number
              tree_fact = sla(pft)*wooddens(pft)/latosa(pft)

              if (lm_ind(pft,jpngr,1).eq.0.0d0) x1=x1+dx !to avoid division by zero
               
c     evaluate f(x1)=LHS of eqn (22) at x1
               
              fx1=allom_fact*((sm_ind(pft,jpngr,1)+
     $             bm_inc_ind-x1-((lm_ind(pft,jpngr,1)+x1)/
     $             lmtorm)+rm_ind(pft,jpngr,1)+hm_ind(pft,jpngr,1))
     $             /wooddens(pft))/(pi/4.0d0)-((sm_ind(pft,jpngr,1)+
     $             bm_inc_ind-x1-((lm_ind(pft,jpngr,1)+x1)/lmtorm)+
     $             rm_ind(pft,jpngr,1))/((lm_ind(pft,jpngr,1)+x1)*
     $             tree_fact))**allom_exp
              
c     Find approximate location of leftmost root on the interval
c     (x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
c     change in sign of f(xmid) relative to f(x1).
              
              fmid=fx1
              xmid=x1

              do while (fmid*fx1.gt.0.0d0.and.xmid.lt.x2)
                
                xmid=xmid+dx
                fmid=allom_fact*((sm_ind(pft,jpngr,1)+
     $               bm_inc_ind-xmid-((lm_ind(pft,jpngr,1)+xmid)/
     $               lmtorm)+rm_ind(pft,jpngr,1)+
     *               hm_ind(pft,jpngr,1))/wooddens(pft))/(pi/4.0d0)-
     *               ((sm_ind(pft,jpngr,1)+bm_inc_ind-xmid-
     *               ((lm_ind(pft,jpngr,1)+xmid)/lmtorm)+
     *               rm_ind(pft,jpngr,1))/((lm_ind(pft,jpngr,1)
     $               +xmid)*tree_fact))**allom_exp
                
              enddo
              
              x1=xmid-dx
              x2=xmid
              
c     Apply bisection method to find root on the new interval (x1,x2)
              
              fx1=allom_fact*((sm_ind(pft,jpngr,1)+
     $             bm_inc_ind-x1-((lm_ind(pft,jpngr,1)+x1)/lmtorm)+
     $             rm_ind(pft,jpngr,1)+hm_ind(pft,jpngr,1))
     *             /wooddens(pft))/(pi/4.0d0)-((sm_ind(pft,jpngr,1)+
     $             bm_inc_ind-x1-((lm_ind(pft,jpngr,1)+x1)/lmtorm)+
     $             rm_ind(pft,jpngr,1))/((lm_ind(pft,jpngr,1)+x1)*
     $             tree_fact))**allom_exp
              
              if (fx1.ge.0.0d0) then
                sign=-1.0d0
              else
                sign=1.0d0
              endif
              
              rtbis=x1
              dx=x2-x1
              
c     Bisection loop
c     Search iterates on value of xmid until xmid lies within
c     xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0
              
              fmid=1.0d0        !dummy value to guarantee entry to loop
              j=0.0d0           !number of iterations so far (maximum tries=jmax)
              
              do while (dx.ge.xacc.and.abs(fmid).gt.yacc)
                
                dx=dx*0.5d0
                xmid=rtbis+dx
                
c     calculate fmid=f(xmid) [eqn (22)]
                
                fmid=allom_fact*((sm_ind(pft,jpngr,1)+
     $               bm_inc_ind-xmid-((lm_ind(pft,jpngr,1)+xmid)/
     $               lmtorm)+rm_ind(pft,jpngr,1)+
     *               hm_ind(pft,jpngr,1))/wooddens(pft))/(pi/4.0d0)-
     *               ((sm_ind(pft,jpngr,1)+bm_inc_ind-xmid-
     *               ((lm_ind(pft,jpngr,1)+xmid)/lmtorm)+
     *               rm_ind(pft,jpngr,1))/((lm_ind(pft,jpngr,1)+
     *               xmid)*tree_fact))**allom_exp
                
                if (fmid*sign.le.0.0d0) rtbis=xmid
                j=j+1
                
              enddo
              
c     Now rtbis contains numerical solution for lminc_ind given
c     eqn (22)
              lminc_ind=rtbis                               ! (growth.cpp:608)
              
c     Calculate increments in other compartments using allometry
c     relationships

              rminc_ind=(lm_ind(pft,jpngr,1)+lminc_ind)/
     $             lmtorm-rm_ind(pft,jpngr,1) !eqn (9)
              sminc_ind=bm_inc_ind-lminc_ind-rminc_ind !eqn (1)
              
            else                !normal
              
              
              lminc_ind=(bm_inc_ind-lm_ind(pft,jpngr,1)/
     $             lmtorm+rm_ind(pft,jpngr,1))/
     *             (1.0d0+1.0d0/lmtorm) !eqn (33)
              
              if (lminc_ind.ge.0.0d0) then
                
c     Positive allocation to leafmass
                
                rminc_ind=bm_inc_ind-lminc_ind !eqn (31)
                
c     Add killed roots (if any) to below-ground litter
                 
                if (rminc_ind.lt.0.0d0) then

C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     THE FOLLOWING TWO STATEMENTS ARE MISSING IN C++ CODE (growth.cpp:633)
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                  lminc_ind=bm_inc_ind
                  rminc_ind=((lm_ind(pft,jpngr,1) + lminc_ind) 
     $                 / lmtorm)-rm_ind(pft,jpngr,1)
                  
                  call addc(
     $                 litter_bg(pft,jpngr,:),
     $                 (-rminc_ind)*nind(pft,jpngr),
     $                 rm_ind(pft,jpngr,2:ncvar)
     $                 )                                    ! (growth.cpp:633)
#if DyN
c     Update N pool and litter following carbon allocation (growth.cpp:1100)
                  temp = max(rm_ind(pft,jpngr,n_id)
     $                 +rminc_ind*afntoc_rm(pft,jpngr),0.d0)
                  call addn(
     $                 litter_bg(pft,jpngr,:),
     $                 (rm_ind(pft,jpngr,n_id)-temp)*nind(pft,jpngr)
     $                 )
                  rm_ind(pft,jpngr,n_id) = temp
#endif
                   
c     kuno note 2005-01-19: for rminc_ind to be negative in (31), lmtorm has
c     to increase from previous year ( lmtorm/lmtorm(alt) > [lminc_ind+lm_ind]/lm_ind ). 
c     This is sufficient for rminc_ind to always be negative in the above code.
c     Thus, here rm_ind can only die off, which is ok.
                   
                endif
                 
              else
                 
c     Negative allocation to leaf mass
                 
                rminc_ind=bm_inc_ind
                lminc_ind=(rm_ind(pft,jpngr,1)+
     $               rminc_ind)*lmtorm-lm_ind(pft,jpngr,1)
c     from eqn (9)
                 
c     Add killed leaves to litter
                 
                call addc(
     $               litter_ag_fast(pft,jpngr,:), !target pool
     $               (-lminc_ind*nind(pft,jpngr)), !12C to add
     $               lm_ind(pft,jpngr,2:ncvar) !13C signature
     $               )                                      ! (growth.cpp:645)

#if DyN
c     Update N pool and litter following carbon allocation (growth.cpp:1099)
                temp = max(lm_ind(pft,jpngr,n_id)
     $               +lminc_ind*afntoc_lm(pft,jpngr),0.d0)
                call addn(
     $               litter_ag_fast(pft,jpngr,:),
     $               (lm_ind(pft,jpngr,n_id)-temp)*nind(pft,jpngr)
     $               )
                lm_ind(pft,jpngr,n_id) = temp
#endif

                 
              endif
               
c     Calculate sminc_ind (must be negative)
               
c     from (25),
c     (34) lm_ind + lminc_ind = latosa * (sm_ind + sminc_ind)
c     / (wooddens * height * sla)
c     from (34),
c     (35) sminc_ind = (lm_ind + lminc_ind) * wooddens * height * sla
c     / latosa - sm_ind
               
              sminc_ind=(lm_ind(pft,jpngr,1)+lminc_ind)*
     $             wooddens(pft)*height(pft,jpngr)*
     *             sla(pft)/latosa(pft)-sm_ind(pft,jpngr,1) !eqn (35)
               
c     kuno note 2005-01-19:	if the allometric relations were satisfied at this
c     point, we would have from (35)
c     sminc_ind=lminc_ind*sm_ind/lm_ind
c     but - they don't seem to be (??)
c     Thus, I don't see on what condition sminc_ind gets negative,
c     so the following is probably ok (avoids unsensical mass fluxes, but
c     violates the condition (A) (not enough sapwood to support leaves).
               
#if version_sg
c     prevent negative allocation to heartwood sg 22.10.01
              if (sminc_ind.gt.0.0d0) sminc_ind=0.0d0
c     this can lead to "catastrophe" with fj version, maybe with sg version as well?
c     Is the budget significantly violated when sm growth is not compensated with hm??
c     Can this lead to escalation??? No, budget is not violated, it is only a redistribution
c     of mass between sapwood and heartwood.
#else
c     do not prevent negative allocation to heartwood (version_fj)
#endif
               
c     Convert killed sapwood to heartwood
              call addc(
     $             hm_ind(pft,jpngr,:),
     $             (-sminc_ind),
     $             sm_ind(pft,jpngr,2:ncvar)
     $             )                                        ! (grwoth.cpp:656)
               
            endif               !normal
             
c     Increment C compartments
            call addc(
     $           lm_ind(pft,jpngr,:),
     $           lminc_ind,
     $           bm_inc(pft,jpngr,2:ncvar)
     $           )
            call addc(
     $           rm_ind(pft,jpngr,:),
     $           rminc_ind,
     $           bm_inc(pft,jpngr,2:ncvar)
     $           )
            call addc(
     $           sm_ind(pft,jpngr,:),
     $           sminc_ind,
     $           bm_inc(pft,jpngr,2:ncvar)
     $           )                                          ! (grwoth.cpp:1091 ff)

#if DyN
c     --------------------------------------------------------------------
c     N ALLOCATION for trees
c     Tab.3, Xu-Ri & Prentice (2008) (XP08)
c     (growth.cpp:717 ff)
c     XXX This code is consistent with Tab.3, but not with c++ code!
c     --------------------------------------------------------------------
            
            if (lm_ind(pft,jpngr,1).gt.0.d0) then
               
c     Update fraction of C mass in different compartments
               plantC = lm_ind(pft,jpngr,1)+sm_ind(pft,jpngr,1)
     $              +hm_ind(pft,jpngr,1)+rm_ind(pft,jpngr,1)
               f_cl = lm_ind(pft,jpngr,1)/plantC
               f_cr = rm_ind(pft,jpngr,1)/plantC
               f_cs = (sm_ind(pft,jpngr,1)+hm_ind(pft,jpngr,1))/plantC
               
c     Calculate partitioning coefficients (fraction of total N allocated to different compartments)
c     (Eq.2, Tab.3, XP08)
               f_nl = 1/(1+(rcton_LS*f_cs+rcton_LR*f_cr)/f_cl)
               f_nr=rcton_LR*f_nl*f_cr/f_cl
               f_ns=rcton_LS*f_nl*f_cs/f_cl
               
            else
               
               f_nl=0.d0
               f_nr=0.d0
               f_ns=1.d0
               
            endif

            if (abs(f_nl+f_nr+f_ns-1.d0).gt.1.d-7) then
               print*,'sum of fractions to be allocated is not = 0 !'
               stop
            endif

c     XXX c++ code is inconsistent with equations of Tab.3, XP08 in terms of sm. Because sm_ind became
c     negative, when implemented as in c++, the equations of Tab.3, XP08 are considered here.
            plantN = lm_ind(pft,jpngr,n_id)+sm_ind(pft,jpngr,n_id)
     $           +rm_ind(pft,jpngr,n_id)+bm_inc(pft,jpngr,n_id) !Eq.1, Tab.3, XP08 (growth.cpp:1111)
            
c$$$            plantN = lm_ind(pft,jpngr,n_id)+sm_ind(pft,jpngr,n_id)
c$$$     $           +hm_ind(pft,jpngr,n_id)+rm_ind(pft,jpngr,n_id)
c$$$     $           +bm_inc(pft,jpngr,n_id)                    !Eq.1, Tab.3, XP08 (growth.cpp:1111) XXX Is hm_ind
!     really allocatable? XXX
                        
            lm_ind(pft,jpngr,n_id) = plantN*f_nl
            rm_ind(pft,jpngr,n_id) = plantN*f_nr
            sm_ind(pft,jpngr,n_id) = plantN*f_ns            !Eq.3, Tab.3, XP08
            
c     Define N/C ratio (instead of C/N) for efficiency (most operations are of type Nmass=Cmass*[N/C ratio]).
            if (lm_ind(pft,jpngr,1).gt.0.d0) then 
               afntoc_lm(pft,jpngr) = lm_ind(pft,jpngr,n_id)/lm_ind(pft,jpngr,1) !Eq.4, Tab.3, XP08 (growth.cpp:1120)
            else
               afntoc_lm(pft,jpngr) = 1/cton_pro(pft)
            endif

            if (rm_ind(pft,jpngr,1).gt.0.d0) then 
               afntoc_rm(pft,jpngr) = rm_ind(pft,jpngr,n_id)/rm_ind(pft,jpngr,1)
            else
               afntoc_rm(pft,jpngr) = 1/cton_pro(pft)
            endif

            if (sm_ind(pft,jpngr,1).gt.0.d0) then
               afntoc_sm(pft,jpngr) =sm_ind(pft,jpngr,n_id)/sm_ind(pft,jpngr,1)
            else
               afntoc_sm(pft,jpngr) = 1/cton_pro(pft)
            endif            
#endif
            
c     spahni 2010-02-10: check if compartments are still positive
            if (lm_ind(pft,jpngr,1).lt.1.0d-15.or.
     $           rm_ind(pft,jpngr,1).lt.1.0d-15.or.
     $           sm_ind(pft,jpngr,1).lt.1.0d-15) then
              present(pft,jpngr) = .false. !remove pft
              call kill(pft,jpngr) !move any carbon (and N) to litter
            endif
            
             
c     Calculate new height, diameter and crown area
            if (present(pft,jpngr)) then
              sap_xsa=lm_ind(pft,jpngr,1)*sla(pft)/latosa(pft) !eqn (5)
              height(pft,jpngr)=sm_ind(pft,jpngr,1)/sap_xsa
     $             /wooddens(pft)
              stemdiam=(height(pft,jpngr)/allom2(pft))**(1.0d0/allom3) !eqn (C)
              crownarea(pft,jpngr)=min(allom1(pft)*stemdiam**reinickerp,
     $             crownarea_max) !eqn (D)
            else
              sap_xsa = 0.0d0
              height(pft,jpngr) = 0.0d0
              stemdiam = 0.0d0
              crownarea(pft,jpngr) = 0.0d0
           endif


          elseif (grass(pft)) then            
c     -------------------------------------------------------------------------
c     GRASS ALLOCATION
c     Distribute this year's production among leaves and fine roots
c     according to leaf to rootmass ratio [eqn (33)]
c     Relocation of C from one compartment to the other not allowed:
c     negative increment in either compartment transferred to litter
c     -------------------------------------------------------------------------

             lminc_ind=(bm_inc_ind-lm_ind(pft,jpngr,1)/
     $            lmtorm+rm_ind(pft,jpngr,1))/
     *            (1.0d0+1.0d0/lmtorm)
             rminc_ind=bm_inc_ind-lminc_ind
             
             if (lminc_ind.ge.0.0d0) then
C     comments to comply with sitch/FJ version, kuno as of 2005-01-18
c     kuno 2005-05-31:	changed to preprocessor flag
#if version_sg
c     Add killed roots (if any) to below-ground litter
              if (rminc_ind.lt.0.0d0) then
c     change sg 17.1.00 
                 lminc_ind = bm_inc_ind 
                  
                 call addc(
     $                  litter_bg(pft,jpngr,:),
     $                  (-rminc_ind)*nind(pft,jpngr),
     $                  rm_ind(pft,jpngr,2:ncvar)
     $                  )
#if DyN
                temp = max(rm_ind(pft,jpngr,n_id)
     $               +rminc_ind*afntoc_rm(pft,jpngr),0.d0)
                call addn(
     $               litter_bg(pft,jpngr,:),
     $               (rm_ind(pft,jpngr,n_id)-temp)*nind(pft,jpngr)
     $               )
                rm_ind(pft,jpngr,n_id) = temp
#endif
                
              endif
#endif
              
            else                !lminc_ind < 0.0d0
                   
c     Negative allocation to leaf mass
               
              rminc_ind=bm_inc_ind
              lminc_ind=(rm_ind(pft,jpngr,1)+
     $             rminc_ind)*lmtorm-lm_ind(pft,jpngr,1)
c     from eqn (9)
              
c     Add killed leaves to litter
c     spahni 2010-02-11: removed *2 in equation below - typo ?
              call addc(
     $             litter_ag_fast(pft,jpngr,:),
     $             (-lminc_ind)*nind(pft,jpngr),
     $             lm_ind(pft,jpngr,2:ncvar)
     $             )

#if DyN
              temp = max(lm_ind(pft,jpngr,n_id)
     $             +lminc_ind*afntoc_lm(pft,jpngr),0.d0)
              call addn(
     $             litter_ag_fast(pft,jpngr,:),
     $             (lm_ind(pft,jpngr,n_id)-temp)*nind(pft,jpngr)
     $             )
              lm_ind(pft,jpngr,n_id) = temp
#endif              
              
            endif
            
c     Increment C compartments (growth.cpp:1162)
            
            call addc(
     $           lm_ind(pft,jpngr,:),
     $           lminc_ind,
     $           bm_inc(pft,jpngr,2:ncvar)
     $           )
            
            call addc(
     $           rm_ind(pft,jpngr,:),
     $           rminc_ind,
     $           bm_inc(pft,jpngr,2:ncvar)
     $           )


#if DyN
c     --------------------------------------------------------------------
c     N ALLOCATION for Grasses
c     Tab.3, Xu-Ri & Prentice (2008) (XP08)
c     (growth.cpp:717 ff)
c     --------------------------------------------------------------------

            if (lm_ind(pft,jpngr,1).gt.0.d0) then
               
c     Update fraction of C mass in different compartments
               plantC = lm_ind(pft,jpngr,1)+rm_ind(pft,jpngr,1)
               f_cl = lm_ind(pft,jpngr,1)/plantC
               f_cr = rm_ind(pft,jpngr,1)/plantC
               
c     Calculate partitioning coefficients (fraction of total N allocated to different compartments)
c     (Eq.2, Tab.3, XP08)
               f_nl=1/(1+(rcton_LR*f_cr)/f_cl);
               f_nr=rcton_LR*f_nl*f_cr/f_cl;
               
            else
               
               f_nl=0.0
               f_nr=1.0
               
            endif
            
            plantN = lm_ind(pft,jpngr,n_id)+rm_ind(pft,jpngr,n_id)
     $           +bm_inc(pft,jpngr,n_id)                    !Eq.1, Tab.3, XP08 (growth.cpp:1181)
            
            lm_ind(pft,jpngr,n_id) = plantN*f_nl
            rm_ind(pft,jpngr,n_id) = plantN*f_nr            !Eq.3, Tab.3, XP08 (growth.cpp:1184)
            
c     Define N/C ratio (instead of C/N) for efficiency (most operations are of type Nmass=Cmass*[N/C ratio]).
c     condition introduced by bs to avoid fpe.
            if (lm_ind(pft,jpngr,1).gt.0.d0) then
               afntoc_lm(pft,jpngr) = lm_ind(pft,jpngr,n_id)
     $              /lm_ind(pft,jpngr,1)                    !Eq.4, Tab.3, XP08 (growth.cpp:1189)
            else
               afntoc_lm(pft,jpngr) = 1/cton_pro(pft)
            endif

            if (rm_ind(pft,jpngr,1).gt.0.d0) then
               afntoc_rm(pft,jpngr) = rm_ind(pft,jpngr,n_id)
     $              /rm_ind(pft,jpngr,1)
            else
               afntoc_rm(pft,jpngr) = 1/cton_pro(pft)
            endif
            
#endif          


c     spahni 2010-02-10: check if compartments are still positive
            if (lm_ind(pft,jpngr,1).lt.1.0d-15.or.
     $           rm_ind(pft,jpngr,1).lt.1.0d-15) then
              present(pft,jpngr) = .false. !remove pft
              call kill(pft,jpngr) !move any carbon (and nitrogen) to litter
            endif

            
          else                  !moss
c     -------------------------------------------------------------------------
c     MOSS ALLOCATION

             lminc_ind = bm_inc_ind                         !always > 0 in this loop
            
             call addc(
     $            lm_ind(pft,jpngr,:),
     $            lminc_ind,
     $            bm_inc(pft,jpngr,2:ncvar)
     $            )

c     spahni 2010-02-10: check if compartments are still positive
            if (lm_ind(pft,jpngr,1).lt.1.0d-15) then
              present(pft,jpngr) = .false. !remove pft
              call kill(pft,jpngr) !move any carbon to litter
            endif

c     spahni 2010-02-10: check if compartments are still positive
            if (lm_ind(pft,jpngr,1).lt.1.0d-15) then
              present(pft,jpngr) = .false. !remove pft
              call kill(pft,jpngr) !move any carbon to litter
            endif

          endif                 ! end: if tree
          
        endif                   !present
        

        if (present(pft,jpngr)) then        
c     Update LAI and FPC
#if (competition_scheme_fj || competition_scheme_sg)          
          if (crownarea(pft,jpngr).gt.0.0d0) then
            lai_ind(pft,jpngr)=(lm_ind(pft,jpngr,1)*sla(pft))/
     $           crownarea(pft,jpngr)
          else
            lai_ind(pft,jpngr)=0.0d0
          endif
          
          fpc_ind=1.0d0-dexp(-0.5d0*lai_ind(pft,jpngr))
          fpc_grid_old(pft,jpngr)=fpc_grid(pft,jpngr)
          fpc_grid(pft,jpngr)=crownarea(pft,jpngr)*
     $         nind(pft,jpngr)*fpc_ind
          fpc_inc(pft,jpngr)=max(0.d0,fpc_grid(pft,jpngr)-
     $         fpc_grid_old(pft,jpngr))
#else                           !take all fpc change into account
c     kuno 2005-01-31: set fpc_grid old after, and fpc_inc before light
c     in order to include all primary production into light competition
          call update_fpc(pft,jpngr)
#endif
          
          
        endif                   !present
         
      enddo                     !pft
       
      return
      end
      
