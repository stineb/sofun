subroutine allocation_annual( jpngr )
  !////////////////////////////////////////////////////////////////
  ! Annual allocation scheme adopted from LPX in case 
  ! _allocation_pipe == 0 or adpted from Makela et al., 2008 else.
  !----------------------------------------------------------------
#include "sofun_module_control.inc"
  use params_core
  use classdefs
  use pools, only: pleaf, proot, psapw, plitt_af, plitt_as, plabl
  use treegeometry, only: ispresent, nind
  use params_modl, only: lmtorm0, tree, grass, rcton_LS, rcton_LR, &
    cton_pro 
  use outvars, only: outaCalloc, outaNalloc, outaCveg2lit, outaNveg2lit
#if _allocation_makela
  use params_modl, only: c_H, alpha_w
#endif

  implicit none

  ! ARGUMENTS
  integer, intent(in) :: jpngr          ! gridcell number


#if _allocation_makela
  !////////////////////////////////////////////////////////////////
  ! Daily allocation scheme. Pipe model (steady state). Adopted
  ! from Makela et al., 2008
  !----------------------------------------------------------------

  ! LOCAL VARIABLES
  integer :: pft 
  real    :: L_ss                          ! Pipe length at quasi-steady state

  do pft=1,npft
    if (ispresent(pft,jpngr)) then

      ! Pipe length at quasi-steady state (Eq.15 in Makela et al., 2008)
      L_ss = c_H * ntoc( pleaf(pft,jpngr) )

      ! Mass of live wood
      psapw(pft,jpngr)%c%c12 = alpha_w * pleaf(pft,jpngr)%c%c12  * L_ss

    end if
  end do
  return

#else
  !////////////////////////////////////////////////////////////////
  ! Annual allocation scheme adopted from LPX.
  !----------------------------------------------------------------

  ! LOCAL VARIABLES
  integer :: pft                        ! PFT counter
  real    :: req                        ! N required if all labile C was to be used for allocation
  real    :: avl                        ! N stored in labile pool
  real    :: nalloc                     ! N actually allocated this year
  real    :: calloc                     ! C actually allocated this year
  real    :: lmtorm                     ! adjusted leaf-to-root allocation ratio
  real    :: calloc_ind                 ! total C increment per individual (=bm_inc_ind)
  real    :: nalloc_ind                 ! total N increment per individual
  real    :: lminc_ind                  ! leaf C mass increment per individual
  real    :: rminc_ind                  ! root C mass increment per individual
  real    :: plantC                     ! total plant C after C allocation
  real    :: plantN                     ! total plant N after N allocation
  real    :: f_cl                       ! fraction of total plant C in leaves
  real    :: f_cr                       ! fraction of total plant C in roots
  real    :: f_nl                       ! fraction of total plant N in leaves
  real    :: f_nr                       ! fraction of total plant N in roots
  real    :: tmp
  type(orgpool) :: orgtmp

  integer, parameter :: nseg = 20
  integer, parameter :: jmax = 40       ! max number of iterations in search for alloc soln

#if _check_sanity  
  real    :: cbal_before
  real    :: cbal_after
  real    :: nbal_before
  real    :: nbal_after

  cbal_before = plabl(1,1)%c%c12 + pleaf(1,1)%c%c12 + proot(1,1)%c%c12
  nbal_before = plabl(1,1)%n%n14 + pleaf(1,1)%n%n14 + proot(1,1)%n%n14
  !print*,'A: ALLOCATION'
  !print*,'pleaf(1,1)   ',pleaf(1,1)
  !print*,'proot(1,1)   ',proot(1,1) 
  !print*,'plabl(1,1)   ',plabl(1,1) 
  !print*,'C:N in plabl ',cton( plabl(1,1) )
#endif

  ! PFT LOOP
  do pft=1,npft

    if (tree(pft)) stop &
      'ALLOCATION: not implemented for tree PFTs. abort.'
     
    if (ispresent(pft,jpngr)) then
      
      ! sanity check`
      if (nind(pft,jpngr)==0) stop

      ! xxx debug
      !print*,'plabl',plabl

      ! determine how much of labile pool to allocate to growth this year
      req = plabl(pft,jpngr)%c%c12/cton_pro(pft)
      avl = plabl(pft,jpngr)%n%n14
      !req = dbminc%c%c12/cton_pro(pft)
      !avl = dbminc%n%n14

      !print*,'req, avl ',req,avl
      nalloc = min(req,avl)
      calloc = cton_pro(pft)*nalloc

      ! OUTPUT COLLECTION
      outaCalloc(pft,jpngr) = calloc
      outaNalloc(pft,jpngr) = nalloc

      ! unit conversion to per-individual
      calloc_ind = calloc / nind(pft,jpngr)
      nalloc_ind = nalloc / nind(pft,jpngr)

      !print*,'C:N of alloc     ',calloc / nalloc
      !print*,'C:N of alloc_ind ',calloc_ind / nalloc_ind

      !print*,'C:N of whole plant ',cton( orgplus(pleaf(pft,jpngr),proot(pft,jpngr)), default=0.0 )

      !print*,'anpp         ',anpp(pft)
      ! remove C and N to be allocated from labile pool
      call csub( carbon(calloc), plabl(pft,jpngr)%c )
      call nsub( nitrogen(nalloc), plabl(pft,jpngr)%n )

      ! calculate this year's leaf to fine root mass ratio from mean annual 
      ! water scalar and pft specific parameter
      ! xxx try
      !lmtorm = lmtorm0 * wscal(pft,jpngr)
      lmtorm = lmtorm0(pft) * 1.0  

      ! xxx debug
      !lmtorm = 1.0
      !print*,'ALLOCATION: lmtorm set to ',lmtorm 

      if (grass(pft)) then
        !-------------------------------------------------------------------------
        ! GRASS ALLOCATION
        ! Distribute this year's production among leaves and fine roots
        ! according to leaf to rootmass ratio [eqn (33)]
        ! Relocation of C from one compartment to the other not allowed:
        ! negative increment in either compartment transferred to litter
        !-------------------------------------------------------------------------
        ! allocation to leaves
        lminc_ind = ( proot(pft,jpngr)%c%c12 + calloc_ind - pleaf(pft,jpngr)%c%c12/lmtorm ) / ( 1.0 + 1.0/lmtorm )

        ! allocation to roots
        rminc_ind = calloc_ind - lminc_ind
         
        if (lminc_ind<0.0) then   
          ! Negative allocation to leaf mass
          !print*,'Negative allocation to leaf mass'
          rminc_ind = calloc_ind
          lminc_ind = lmtorm * ( proot(pft,jpngr)%c%c12 + rminc_ind ) - pleaf(pft,jpngr)%c%c12
            
          ! Add killed leaves C and N to litter (lminc_ind is neg.!)
          ! 'tmp' is fraction of leaves killed
          tmp = (-1.0)*nind(pft,jpngr)*lminc_ind/pleaf(pft,jpngr)%c%c12
          
          call orgmvRec( orgfrac( tmp, pleaf(pft,jpngr) ), pleaf(pft,jpngr), plitt_af(pft,jpngr), outaCveg2lit(pft,jpngr), outaNveg2lit(pft,jpngr) )

          ! not used : call ncp( (-lminc_ind)*ntoc(pleaf(pft,jpngr))*nind(pft,jpngr), plitt_af(pft,jpngr)%n ) 

          ! ! 'tmp' is the amount of N necessary for (reduced) leaf C, assuming const. leaf-C:N 
          ! tmp = max( pleaf(pft,jpngr,n_id) + lminc_ind*afntoc_lm(pft,jpngr), 0.0 )
          ! call ncp( (pleaf(pft,jpngr)%n-tmp)*nind(pft,jpngr), plitt_af(pft,jpngr)%n )
          ! pleaf(pft,jpngr)%n = tmp
          
        endif

        ! Increment C compartments (growth.cpp:1162). When lminc_ind is neg., this is a subtraction;
        ! note that orgcp is used above and not orgmv.
        call ccp( carbon(lminc_ind), pleaf(pft,jpngr)%c )
        call ccp( carbon(rminc_ind), proot(pft,jpngr)%c )

        !print*,'lminc_ind',lminc_ind
        !print*,'rminc_ind',rminc_ind
        
        !--------------------------------------------------------------------
        ! N ALLOCATION for Grasses
        ! Tab.3, Xu-Ri & Prentice (2008) (XP08)
        ! (growth.cpp:717 ff)
        !--------------------------------------------------------------------
        if (pleaf(pft,jpngr)%c%c12>0.0 .and. nalloc_ind>0.0) then

          ! Update fraction of C mass in different compartments
          !plantC = pleaf(pft,jpngr)%c%c12 + proot(pft,jpngr)%c%c12
          !f_cl   = pleaf(pft,jpngr)%c%c12 / plantC
          !f_cr   = proot(pft,jpngr)%c%c12 / plantC
           
          ! Calculate partitioning coefficients (fraction of total N allocated to different compartments)
          ! (Eq.2, Tab.3, XP08)
          
          ! xxx debug
          !rcton_LR = 1.0
          !print*,'ALLOCATION: rcton_LR set to ',rcton_LR

          tmp = rcton_LR * proot(pft,jpngr)%c%c12 / pleaf(pft,jpngr)%c%c12
          f_nl = (proot(pft,jpngr)%n%n14 + nalloc_ind - tmp * pleaf(pft,jpngr)%n%n14 ) / (nalloc_ind * (1.0+tmp) )
          f_nr = (1.0 - f_nl)

          !f_nl = 1.0 / ( 1.0 + ( rcton_LR * f_cr ) / f_cl )
          !f_nr = rcton_LR * f_nl * f_cr / f_cl

          !print*,'f_nl, f_nr',f_nl, f_nr
           
        else
           
          f_nl = 0.0
          f_nr = 1.0
           
        endif
        
        !plantN = pleaf(pft,jpngr)%n%n14 + proot(pft,jpngr)%n%n14 + nalloc_ind                   !Eq.1, Tab.3, XP08 (growth.cpp:1181)

        !if (plantN>0.0) print*,'plantC/plantN',plantc/plantN
        
        !pleaf(pft,jpngr)%n = nitrogen(plantN*f_nl)
        !proot(pft,jpngr)%n = nitrogen(plantN*f_nr)            !Eq.3, Tab.3, XP08 (growth.cpp:1184)
                    
        call ncp( nitrogen(nalloc_ind*f_nl), pleaf(pft,jpngr)%n )
        call ncp( nitrogen(nalloc_ind*f_nr), proot(pft,jpngr)%n )

        !print*,'leaf-C  ',pleaf(pft,jpngr)%c%c12
        !print*,'root-C  ',proot(pft,jpngr)%c%c12
        !print*,'leaf-N  ',pleaf(pft,jpngr)%n%n14
        !print*,'root-N  ',proot(pft,jpngr)%n%n14
        !print*,'leaf-C:N',cton(pleaf(pft,jpngr))
        !print*,'root-C:N',cton(proot(pft,jpngr))
        !print*,'ALLCOATION: C:N of whole plant ',cton( orgplus(pleaf(pft,jpngr),proot(pft,jpngr)), default=0.0 )
        !print*,'prescC:N',cton_pro(pft)
        !print*,'CN: leaf:root',cton(pleaf(pft,jpngr))/cton(proot(pft,jpngr))
        !print*,'rcton_LR',rcton_LR

        ! xxx try
        ! spahni 2010-02-10: check if compartments are still positive
        !if (pleaf(pft,jpngr)%c%c12<1.0d-15 .or.proot(pft,jpngr)%c%c12<1.0d-15) then
          !remove PFT
          !ispresent(pft,jpngr) = .false. 
          ! xxx try: uncomment and define 'kill' SR
          !call kill(pft,jpngr) !move any carbon (and nitrogen) to litter
        !endif

      ! else
      !   !-------------------------------------------------------------------------
      !   ! TREE ALLOCATION
      !   ! Using the algorithm implemented in LPX.
      !   !-------------------------------------------------------------------------            
            
      !   ! Work out minimum leaf production to maintain current sapmass
            
      !   lminc_ind_min = latosa(pft)*psapw(pft,jpngr)%c%c12/(wooddens(pft)*height(pft,jpngr)*sla(pft))-pleaf(pft,jpngr)%c%c12      !eqn (27)
        
      !   rminc_ind_min = latosa(pft)*psapw(pft,jpngr)%c%c12/(wooddens(pft)*height(pft,jpngr)*sla(pft)*lmtorm)-proot(pft,jpngr)%c%c12 !eqn (30)
        
      !   if (rminc_ind_min>0.0 .and. lminc_ind_min>0.0 .and. rminc_ind_min+lminc_ind_min<=calloc_ind) then !(see growth.cpp:535)
      !     ! Normal allocation (positive increment to all living C compartments)
      !     x1=0.0
      !     x2=(calloc_ind-(pleaf(pft,jpngr)%c%c12/lmtorm-proot(pft,jpngr)%c%c12))/(1.0+1.0/lmtorm)
      !     dx=(x2-x1)/real(nseg)

      !     ! allocate constants in allom eq. for simplification and numerical stability
      !     allom_fact = allom2(pft)**(2.0d0/allom3)           
      !     allom_exp = dnint(1.0+2.0d0/allom3)         !must be a whole number
      !     tree_fact = sla(pft)*wooddens(pft)/latosa(pft)

      !     if (pleaf(pft,jpngr)%c%c12==0.0) x1=x1+dx    !to avoid division by zero
          
      !     ! evaluate f(x1)=LHS of eqn (22) at x1
      !     fx1=allom_fact*((psapw(pft,jpngr)%c%c12+ &
      !       calloc_ind-x1-((pleaf(pft,jpngr)%c%c12+x1)/ &
      !       lmtorm)+proot(pft,jpngr)%c%c12+pwood(pft,jpngr)%c%c12) &
      !       /wooddens(pft))/(pi/4.0d0)-((psapw(pft,jpngr)%c%c12+ &
      !       calloc_ind-x1-((pleaf(pft,jpngr)%c%c12+x1)/lmtorm)+ &
      !       proot(pft,jpngr)%c%c12)/((pleaf(pft,jpngr)%c%c12+x1)* &
      !       tree_fact))**allom_exp
          
      !     ! Find approximate location of leftmost root on the interval
      !     ! (x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
      !     ! change in sign of f(xmid) relative to f(x1).
      !     fmid=fx1
      !     xmid=x1

      !     do while (fmid*fx1>0.0.and.xmid<x2)
            
      !       xmid=xmid+dx
      !       fmid=allom_fact*((psapw(pft,jpngr)%c%c12+ &
      !         calloc_ind-xmid-((pleaf(pft,jpngr)%c%c12+xmid)/ &
      !         lmtorm)+proot(pft,jpngr)%c%c12+ &
      !         pwood(pft,jpngr)%c%c12)/wooddens(pft))/(pi/4.0d0)- &
      !         ((psapw(pft,jpngr)%c%c12+calloc_ind-xmid- &
      !         ((pleaf(pft,jpngr)%c%c12+xmid)/lmtorm)+ &
      !         proot(pft,jpngr)%c%c12)/((pleaf(pft,jpngr)%c%c12 &
      !         +xmid)*tree_fact))**allom_exp
            
      !     enddo
          
      !     x1=xmid-dx
      !     x2=xmid
          
      !     ! Apply bisection method to find root on the new interval (x1,x2)
      !     fx1=allom_fact*((psapw(pft,jpngr)%c%c12+ &
      !       calloc_ind-x1-((pleaf(pft,jpngr)%c%c12+x1)/lmtorm)+ &
      !       proot(pft,jpngr)%c%c12+pwood(pft,jpngr)%c%c12) &
      !       /wooddens(pft))/(pi/4.0d0)-((psapw(pft,jpngr)%c%c12+ &
      !       calloc_ind-x1-((pleaf(pft,jpngr)%c%c12+x1)/lmtorm)+ &
      !       proot(pft,jpngr)%c%c12)/((pleaf(pft,jpngr)%c%c12+x1)* &
      !       tree_fact))**allom_exp
          
      !     if (fx1>=0.0) then
      !       fsign=-1.0
      !     else
      !       fsign=1.0
      !     endif
          
      !     rtbis=x1
      !     dx=x2-x1
          
      !     ! Bisection loop
      !     ! Search iterates on value of xmid until xmid lies within
      !     ! xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0
          
      !     fmid=1.0                                    !dummy value to guarantee entry to loop
      !     j=0.0                                       !number of iterations so far (maximum tries=jmax)
          
      !     do while (dx>=xacc.and.abs(fmid)>yacc)
            
      !       dx=dx*0.5d0
      !       xmid=rtbis+dx
            
      !       ! calculate fmid=f(xmid) [eqn (22)]
      !       fmid=allom_fact*((psapw(pft,jpngr)%c%c12+ &
      !         calloc_ind-xmid-((pleaf(pft,jpngr)%c%c12+xmid)/ &
      !         lmtorm)+proot(pft,jpngr)%c%c12+ &
      !         pwood(pft,jpngr)%c%c12)/wooddens(pft))/(pi/4.0d0)- &
      !         ((psapw(pft,jpngr)%c%c12+calloc_ind-xmid- &
      !         ((pleaf(pft,jpngr)%c%c12+xmid)/lmtorm)+ &
      !         proot(pft,jpngr)%c%c12)/((pleaf(pft,jpngr)%c%c12+ &
      !         xmid)*tree_fact))**allom_exp
            
      !       if (fmid*fsign<=0.0) rtbis=xmid
      !       j=j+1
            
      !     enddo
          
      !     ! Now rtbis contains numerical solution for lminc_ind given
      !     ! eqn (22)
      !     lminc_ind=rtbis                               ! (growth.cpp:608)
          
      !     ! Calculate increments in other compartments using allometry
      !     ! relationships
      !     rminc_ind=(pleaf(pft,jpngr)%c%c12+lminc_ind)/ &lmtorm-proot(pft,jpngr)%c%c12  !eqn (9)
      !     sminc_ind=calloc_ind-lminc_ind-rminc_ind      !eqn (1)
          
      !   else                                            !normal
          
      !     lminc_ind=(calloc_ind-pleaf(pft,jpngr)%c%c12/lmtorm+proot(pft,jpngr)%c%c12)/(1.0+1.0/lmtorm)  !eqn (33)
          

      !     if (lminc_ind>=0.0) then
            
      !       ! Positive allocation to leafmass
      !       rminc_ind=calloc_ind-lminc_ind              !eqn (31)
            
      !       ! Add killed roots (if any) to below-ground litter
      !       if (rminc_ind<=0.0) then

      !         lminc_ind=calloc_ind

      !         ! rminc_ind is calculated so that ratio 'lmtorm' is satisfied, given lminc_inc
      !         rminc_ind=((pleaf(pft,jpngr)%c%c12 + lminc_ind) / lmtorm)-proot(pft,jpngr)%c%c12
              
      !         ! from eqn (9)
      !         ! rminc_ind < 0, => add "superfluous" rminc_ind to litter, N along with it                  
      !         ntoc_tmp = ntoc(proot(pft,jpngr))
      !         call orgcp( orgpool( carbon(-rminc_ind), nitrogen(-rminc_ind*ntoc_tmp) ), &
      !                     plitt_bg(pft,jpngr), &
      !                     nind(pft,jpngr) &
      !                   )
      !         proot(pft,jpngr)%n%n14 = proot(pft,jpngr)%n%n14 + rminc_ind*ntoc_tmp


      !         specialcase_rm(pft) = .true.
              
      !         ! kuno note 2005-01-19: for rminc_ind to be negative in (31), lmtorm has
      !         ! to increase from previous year ( lmtorm/lmtorm(alt) > [lminc_ind+lm_ind]/lm_ind ). 
      !         ! This is sufficient for rminc_ind to always be negative in the above code.
      !         ! Thus, here rm_ind can only die off, which is ok.
              
      !       endif
            
      !     else                                          !lminc_ind < 0
            
      !       ! Negative allocation to leaf mass
      !       rminc_ind=calloc_ind
      !       lminc_ind=(proot(pft,jpngr)%c%c12+rminc_ind)*lmtorm-pleaf(pft,jpngr)%c%c12
    
      !       ! from eqn (9)
      !       ! lminc_ind < 0 => add "superfluous" lminc_ind to litter
      !       ntoc_tmp = ntoc(pleaf(pft,jpngr))
      !       call orgcp( orgpool( carbon(-lminc_ind), nitrogen(-lminc_ind*ntoc_tmp) ), &
      !                   plitt_af(pft,jpngr), &
      !                   nind(pft,jpngr) &
      !                 )
      !       pleaf(pft,jpngr)%n%n14 = pleaf(pft,jpngr)%n%n14 + lminc_ind*ntoc_tmp

      !       specialcase_lm(pft) = .true.
            
      !     endif
          
      !     ! Calculate sminc_ind (must be negative)
                
      !     ! from (25),
      !     ! (34) lm_ind + lminc_ind = latosa * (sm_ind + sminc_ind)
      !     ! / (wooddens * height * sla)
      !     ! from (34),
      !     ! (35) sminc_ind = (lm_ind + lminc_ind) * wooddens * height * sla
      !     ! / latosa - sm_ind
          
      !     sminc_ind=(pleaf(pft,jpngr)%c%c12+lminc_ind)*wooddens(pft)*height(pft,jpngr)*sla(pft)/latosa(pft)-psapw(pft,jpngr)%c%c12 !eqn (35)
          
      !     ! kuno note 2005-01-19: if the allometric relations were satisfied at this
      !     ! point, we would have from (35)
      !     ! sminc_ind=lminc_ind*sm_ind/lm_ind
      !     ! but - they don't seem to be (??)
      !     ! Thus, I don't see on what condition sminc_ind gets negative,
      !     ! so the following is probably ok (avoids unsensical mass fluxes, but
      !     ! violates the condition (A) (not enough sapwood to support leaves).
                              
      !     ! Convert killed sapwood to heartwood (sminc_ind .le. 0 in all cases)
      !     ! Don'd do anything with N as it is re-distributed between plant pools later anyway
      !     ! print*,'converted killed sapwood to heartwood.'
      !     ccp( carbon(-sminc_ind), pwood(pft,jpngr) )

      !   endif                                           !normal
        
      !   ! Increment C compartments
      !   ccp( carbon(lminc_ind), pleaf(pft,jpngr) )
      !   ccp( carbon(rminc_ind), proot(pft,jpngr) )
      !   ccp( carbon(sminc_ind), psapw(pft,jpngr) )

      !   ! For this PFT, bm_inc is now allocated.
      !   bm_inc(pft,jpngr,1) = 0.0

      !   ! This is not necessary when using the ks_hydrology flag. ==> Why does this happen?
      !   if (pleaf(pft,jpngr)%c%c12<0.0.or.&
      !       proot(pft,jpngr)%c%c12<0.0.or.&
      !       psapw(pft,jpngr)%c%c12<0.0) then
      !     print*,'REMOVING PFT DUE TO NEGATIVE PLANT POOL'
      !     ispresent(pft,jpngr) = .false.
      !     ! In this case (pft), N has not yet been allocated. I.o. to conserve
      !     ! N balance, put N to be allocated back to Ninorg.
      !     call kill(pft,jpngr)
      !   endif

      !   !--------------------------------------------------------------------
      !   ! N ALLOCATION for trees
      !   ! Tab.3, Xu-Ri & Prentice (2008) (XP08)
      !   ! (growth.cpp:717 ff)
      !   !--------------------------------------------------------------------
      !   if (pleaf(pft,jpngr)%c%c12>0.0) then
          
      !   ! Calculate partitioning coefficients (fraction of total N allocated to different compartments)
      !   ! (Eq.2, Tab.3, XP08). Hold the ratio of CN ratios of leaf vs. sapwood constant. Same for
      !   ! leaf vs. root.
      !     f_nl = 1/(1+(rcton_LS*psapw(pft,jpngr)%c%c12+rcton_LR*proot(pft,jpngr)%c%c12)/pleaf(pft,jpngr)%c%c12)
      !     f_nr = rcton_LR*f_nl*proot(pft,jpngr)%c%c12/pleaf(pft,jpngr)%c%c12
      !     f_ns = rcton_LS*f_nl*psapw(pft,jpngr)%c%c12/pleaf(pft,jpngr)%c%c12
          
      !   else
          
      !     f_nl=0.0
      !     f_nr=0.0
      !     f_ns=1.0
          
      !   endif

      !   if (abs(f_nl+f_nr+f_ns-1.0)>1.d-7) then
      !     print*,'sum of fractions to be allocated is not = 0 !'
      !     stop
      !   endif
        
      !   if (nind(pft,jpngr).gt.1.d-17) then
      !     plantN = lm_ind(pft,jpngr,n_id)+sm_ind(pft,jpngr,n_id)+rm_ind(pft,jpngr,n_id)+bm_inc(pft,jpngr,n_id)/nind(pft,jpngr) 

      !   ! Eq.1, Tab.3, XP08 (growth.cpp:1111). bm_ind is in units of kg/m2, whereas lm_ind, etc. are in units of kg/m2/ind. (tested without division and CN ratio of biomass change (before and after allocation SR) is not correct!)
      !   elseif (bm_inc(pft,jpngr,n_id).gt.0.0) then
      !     stop 'problem in allocation.F'
      !   endif

      !   lm_ind(pft,jpngr,n_id) = plantN*f_nl
      !   rm_ind(pft,jpngr,n_id) = plantN*f_nr
      !   sm_ind(pft,jpngr,n_id) = plantN*f_ns            !Eq.3, Tab.3, XP08

      !   ! N is now allocated. Thus bm_inc can be set to 0.
      !   bm_inc(pft,jpngr,n_id) = 0.0

      !   call update_afntoc(pft,jpngr)
        
      !   ! spahni 2010-02-10: check if compartments are still positive
      !   if (pleaf(pft,jpngr)%c%c12<1.0d-15 .or. &
      !       proot(pft,jpngr)%c%c12<1.0d-15 .or. &
      !       psapw(pft,jpngr)%c%c12<1.0d-15) then
      !     present(pft,jpngr) = .false.                  !remove pft
      !     call kill(pft,jpngr)                          !move any carbon (and N) to litter
      !     print*,'REMOVING PFT DUE TO VERY SMALL PLANT POOL'
      !   endif
        
        
      !   ! Calculate new height, diameter and crown area
      !   if (present(pft,jpngr)) then
      !     sap_xsa = pleaf(pft,jpngr)%c%c12*sla(pft)/latosa(pft) !eqn (5)
      !     height(pft,jpngr)=psapw(pft,jpngr)%c%c12/sap_xsa/wooddens(pft)
      !     stemdiam=(height(pft,jpngr)/allom2(pft))**(1.0/allom3) !eqn (C)
      !     crownarea(pft,jpngr)=min(allom1(pft)*stemdiam**reinickerp,crownarea_max) !eqn (D)
      !   else
      !     sap_xsa = 0.0
      !     height(pft,jpngr) = 0.0
      !     stemdiam = 0.0
      !     crownarea(pft,jpngr) = 0.0
      !   endif

      endif                 ! if grass

    endif                   !ispresent


    ! xxx try
    !if (ispresent(pft,jpngr)) then        
      ! Update LAI and FPC
      !call update_fpc(pft,jpngr)
    !endif                   !ispresent
     
  enddo                     !pft
 
#if _check_sanity  
  ! xxx test
  cbal_after = plabl(1,1)%c%c12 + pleaf(1,1)%c%c12 + proot(1,1)%c%c12
  nbal_after = plabl(1,1)%n%n14 + pleaf(1,1)%n%n14 + proot(1,1)%n%n14
  if (abs(cbal_after-cbal_before)>1.d9) then
    print*,'cbal_before, cbal_after ', cbal_before, cbal_after
    stop 'C balance violated in ALLOCATION'
  endif
  if (abs(nbal_after-nbal_before)>1.d9) then
    print*,'nbal_before, nbal_after ', nbal_before, nbal_after
    stop 'N balance violated in ALLOCATION'
  endif
  !print*,'B: ALLOCATION'
  !print*,'pleaf(1,1)   ',pleaf(1,1)
  !print*,'proot(1,1)   ',proot(1,1) 
  !print*,'plabl(1,1)   ',plabl(1,1) 
  !print*,'C:N in plabl ',cton( plabl(1,1), default=0.0 )
#endif

  return

#endif

end subroutine allocation_annual


