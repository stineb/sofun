      module classdefs
       !////////////////////////////////////////////////////////////////
       !  Module contains Fortran 90 derived-type declarations to define
       !  objects
       !----------------------------------------------------------------

       implicit none

       ! Carbon, so far contains only c12 (to be extended for c13)
       type carbon
         real :: c12
       end type carbon

       type nitrogen
         real :: n14
       end type nitrogen

       ! Organic pools, contain carbon (c12) and nitrogen (n14)
       type orgpool
         type (carbon)   :: c
         type (nitrogen) :: n
       end type orgpool

       !! Plants, contain leaves and roots
       !type plantclass
         !type (orgpool) :: lm ! leafmass
         !type (orgpool) :: rm ! rootmass
       !end type plantclass

       ! Soil, contains a pool with fast and slow turnover
       !type soilclass
         !type (orgpool) :: fs ! fast
         !type (orgpool) :: sl ! slow
       !end type soilclass

       ! Litter, contains ...
       !type litterclass
         !type (orgpool) :: af ! above-ground fast
         !type (orgpool) :: as ! above-ground slow
         !type (orgpool) :: bg ! below-ground
       !end type litterclass

       contains
       !=========================LOW-LEVEL================================


       !--------------------------ORGANIC---------------------------------

       subroutine orgcp(amount,to,d13C)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "copy" organic mass to pool (e.g. for output).
         !  Does NOT substract amount moved ('amount') from source
         !----------------------------------------------------------------

         implicit none
         type (orgpool), intent(in) :: amount
         type (orgpool), intent(inout) :: to
         real, optional, intent(in) :: d13C

         call ccp(amount%c,to%c)
         call ncp(amount%n,to%n)

!         if (present(d13C)) then
!           to%c%c12 = amount%c%c12 + to%c%c12
!           to%n%n14 = amount%n%n14 + to%n%n14
!         else
!           to%c%c12 = amount%c%c12 + to%c%c12
!           to%n%n14 = amount%n%n14 + to%n%n14
!         end if

       end subroutine orgcp


       subroutine orgsub(amount,from)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "substract" organic mass ('amount') from source 
         !  pool ('from'). ONLY substracts, but does NOT add.
         !----------------------------------------------------------------

         implicit none
         type (orgpool), intent(in) :: amount
         type (orgpool), intent(inout) :: from

         call csub(amount%c,from%c)
         call nsub(amount%n,from%n)

       end subroutine orgsub


       subroutine orgmv(amount,from,to,scale)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "move" organic mass ('amount') from source pool 
         !  ('from') to destination pool ('to'). Substracts amount moved 
         !  ('amount') from source. 'orgmv' is the combination of 'orgcp' 
         !  and 'orgsub'
         !----------------------------------------------------------------

         implicit none
         type (orgpool), intent(in) :: amount
         type (orgpool), intent(inout) :: from
         type (orgpool), intent(inout) :: to
         real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

         if (present(scale)) then
           call orgcp(orgfrac(scale,amount),to) 
           call orgsub(amount,from)       
         else
           call orgcp(amount,to)
           call orgsub(amount,from)
         endif  

       end subroutine orgmv


       !--------------------------CARBON----------------------------------

       subroutine cmv(amount,from,to,scale)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "move" only C from organic mass ('amount') from 
         !  source pool ('from') to destination pool ('to'). Substracts 
         !  amount moved ('amount') from source. 'cmv' is the combination 
         !  of 'ccp' and 'csub'. 
         !----------------------------------------------------------------

         implicit none
         type (carbon), intent(in) :: amount
         type (carbon), intent(inout) :: from
         type (carbon), intent(inout) :: to
         real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

         if (present(scale)) then
           call ccp(cfrac(scale,amount),to)
           call csub(amount,from)
         else
           call ccp(amount,to)
           call csub(amount,from)
         endif

       end subroutine cmv


       subroutine ccp(amount,to,d13C)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "copy" carbon to pool (e.g. for output).
         !  Does NOT substract amount moved ('amount') from source
         !----------------------------------------------------------------

         implicit none
         type (carbon), intent(in) :: amount
         type (carbon), intent(inout) :: to
         real, optional, intent(in) :: d13C

         to%c12 = amount%c12 + to%c12

!         if (present(d13C)) then
!           to%c%c12 = amount%c%c12 + to%c%c12
!           to%n%n14 = amount%n%n14 + to%n%n14
!         else
!           to%c%c12 = amount%c%c12 + to%c%c12
!           to%n%n14 = amount%n%n14 + to%n%n14
!         end if

       end subroutine ccp


       subroutine csub(amount,from)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "substract" organic mass ('amount') from source 
         !  pool ('from'). ONLY substracts, but does NOT add.
         !----------------------------------------------------------------

         implicit none
         type (carbon), intent(in) :: amount
         type (carbon), intent(inout) :: from

         if (from%c12>0) then
           from%c12 = from%c12 - amount%c12
         else if (from%c12>0) then
           stop 'in SR csub: trying to remove c12 from empty pool'
         endif
           
       end subroutine csub

       !--------------------------NITROGEN--------------------------------

       subroutine nmv(amount,from,to,scale)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "move" only C from organic mass ('amount') from 
         !  source pool ('from') to destination pool ('to'). Substracts 
         !  amount moved ('amount') from source. 'nmv' is the combination 
         !  of 'ccp' and 'csub'. 
         !----------------------------------------------------------------

         implicit none
         type (nitrogen), intent(in) :: amount
         type (nitrogen), intent(inout) :: from
         type (nitrogen), intent(inout) :: to
         real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

         if (present(scale)) then
           call ncp(nfrac(scale,amount),to)
           call nsub(amount,from)
         else
           call ncp(amount,to)
           call nsub(amount,from)
         endif

       end subroutine nmv


       subroutine ncp(amount,to)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "copy" nitrogen to pool (e.g. for output).
         !  Does NOT substract amount moved ('amount') from source
         !----------------------------------------------------------------

         implicit none
         type (nitrogen), intent(in) :: amount
         type (nitrogen), intent(inout) :: to

         to%n14 = amount%n14 + to%n14

       end subroutine ncp


       subroutine nsub(amount,from)
         !////////////////////////////////////////////////////////////////
         !  Generic SR to "substract" nitrogen ('amount') from source 
         !  pool ('from'). ONLY substracts, but does NOT add.
         !----------------------------------------------------------------

         implicit none
         type (nitrogen), intent(in) :: amount
         type (nitrogen), intent(inout) :: from

         if (from%n14>0) then
           from%n14 = from%n14 - amount%n14
         else if (from%n14>0) then
           stop 'in SR csub: trying to remove n14 from empty pool'
         endif

       end subroutine nsub


       !--------------------------FUNCTIONS--------------------------------

       function orgfrac(frac,from)
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable of type 'orgpool' and size
         !  of a fraction 'frac' of source pool ('from')
         !----------------------------------------------------------------
         implicit none

         ! arguments
         real :: frac
         type (orgpool), intent(in) :: from
         
         ! function return variable
         type (orgpool) :: orgfrac

         orgfrac%c%c12 = frac * from%c%c12
         orgfrac%n%n14 = frac * from%n%n14

       end function orgfrac

       function cfrac(frac,from)
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable of type 'carbon' and size
         !  of a fraction 'frac' of source pool ('from')
         !----------------------------------------------------------------
         implicit none

         ! arguments
         real :: frac
         type (carbon), intent(in) :: from
         
         ! function return variable
         type (carbon) :: cfrac

         cfrac%c12 = frac * from%c12

       end function cfrac


       function nfrac(frac,from)
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable of type 'nitrogen' and size
         !  of a fraction 'frac' of source pool ('from')
         !----------------------------------------------------------------
         implicit none

         ! arguments
         real :: frac
         type (nitrogen), intent(in) :: from
         
         ! function return variable
         type (nitrogen) :: nfrac

         nfrac%n14 = frac * from%n14

       end function nfrac


       function orgplus( pool1, pool2 )
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable sum of two pools of type 
         !  'orgpool'. Sum is of type 'orgpool' as well.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (orgpool), intent(in) :: pool1, pool2
         
         ! function return variable
         type (orgpool) :: orgplus

         orgplus%c = cplus(pool1%c,pool2%c)
         orgplus%n = nplus(pool1%n,pool2%n)

       end function orgplus


       function cplus( pool1, pool2, pool3, pool4, pool5, pool6, pool7, pool8, pool9, pool10 )
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable sum of two pools of type 
         !  'carbon'. Sum is of type 'carbon' as well.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (carbon), intent(in) :: pool1
         type (carbon), intent(in) :: pool2
         type (carbon), intent(in), optional :: pool3
         type (carbon), intent(in), optional :: pool4
         type (carbon), intent(in), optional :: pool5
         type (carbon), intent(in), optional :: pool6
         type (carbon), intent(in), optional :: pool7
         type (carbon), intent(in), optional :: pool8
         type (carbon), intent(in), optional :: pool9
         type (carbon), intent(in), optional :: pool10
         
         ! function return variable
         type (carbon) :: cplus

         cplus%c12 = pool1%c12 + pool2%c12

         if (present(pool3)) then
           cplus%c12 = cplus%c12 + pool3%c12
           if (present(pool4)) then
             cplus%c12 = cplus%c12 + pool4%c12
             if (present(pool5)) then
               cplus%c12 = cplus%c12 + pool5%c12
               if (present(pool6)) then
                 cplus%c12 = cplus%c12 + pool6%c12
                 if (present(pool7)) then
                   cplus%c12 = cplus%c12 + pool7%c12
                   if (present(pool8)) then
                     cplus%c12 = cplus%c12 + pool8%c12
                     if (present(pool9)) then
                       cplus%c12 = cplus%c12 + pool9%c12
                       if (present(pool10)) then
                         cplus%c12 = cplus%c12 + pool10%c12
                       end if
                     end if
                   end if
                 end if
               end if
             end if
           end if
         end if

       end function cplus


       function nplus( pool1, pool2, pool3, pool4, pool5, pool6, pool7, pool8, pool9, pool10 )
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable sum of two pools of type 
         !  'nitrogen'. Sum is of type 'nitrogen' as well.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (nitrogen), intent(in) :: pool1
         type (nitrogen), intent(in) :: pool2
         type (nitrogen), intent(in), optional :: pool3
         type (nitrogen), intent(in), optional :: pool4
         type (nitrogen), intent(in), optional :: pool5
         type (nitrogen), intent(in), optional :: pool6
         type (nitrogen), intent(in), optional :: pool7
         type (nitrogen), intent(in), optional :: pool8
         type (nitrogen), intent(in), optional :: pool9
         type (nitrogen), intent(in), optional :: pool10
         
         ! function return variable
         type (nitrogen) :: nplus

         nplus%n14 = pool1%n14 + pool2%n14

         if (present(pool3)) then
           nplus%n14 = nplus%n14 + pool3%n14
           if (present(pool4)) then
             nplus%n14 = nplus%n14 + pool4%n14
             if (present(pool5)) then
               nplus%n14 = nplus%n14 + pool5%n14
               if (present(pool6)) then
                 nplus%n14 = nplus%n14 + pool6%n14
                 if (present(pool7)) then
                   nplus%n14 = nplus%n14 + pool7%n14
                   if (present(pool8)) then
                     nplus%n14 = nplus%n14 + pool8%n14
                     if (present(pool9)) then
                       nplus%n14 = nplus%n14 + pool9%n14
                       if (present(pool10)) then
                         nplus%n14 = nplus%n14 + pool10%n14
                       end if
                     end if
                   end if
                 end if
               end if
             end if
           end if
         end if

       end function nplus


       function cminus( pool1, pool2 )
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable difference of two pools of 
         !  type 'carbon'. Sum is of type 'carbon' as well.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (carbon), intent(in) :: pool1
         type (carbon), intent(in) :: pool2
         
         ! function return variable
         type (carbon) :: cminus

         cminus%c12 = pool1%c12 - pool2%c12

       end function cminus


       function cton( pool )
         !////////////////////////////////////////////////////////////////
         !  Generic function to return the C:N ratio of an organic pool.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (orgpool), intent(in) :: pool
         
         ! function return variable
         real :: cton

         cton = pool%c%c12 / pool%n%n14

       end function cton


       function ntoc( pool )
         !////////////////////////////////////////////////////////////////
         !  Generic function to return the N:C ratio of an organic pool.
         !  This is equal to the inverse of the 'cton' function.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (orgpool), intent(in) :: pool
         
         ! function return variable
         real :: ntoc

         ntoc = pool%n%n14 / pool%c%c12 

       end function ntoc


       end module classdefs
