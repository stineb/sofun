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

       ! Plants, contain leaves and roots
       type plantclass
         type (orgpool) :: lm ! leafmass
         type (orgpool) :: rm ! rootmass
       end type plantclass

       ! Soil, contains a pool with fast and slow turnover
       type soilclass
         type (orgpool) :: fs ! fast
         type (orgpool) :: sl ! slow
       end type soilclass

       ! Litter, contains ...
       type litterclass
         type (orgpool) :: af ! above-ground fast
         type (orgpool) :: as ! above-ground slow
         type (orgpool) :: bg ! below-ground
       end type litterclass

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

       subroutine cmv(amount,from,to)
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

         call ccp(amount,to)
         call csub(amount,from)

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

         from%c12 = from%c12 - amount%c12

       end subroutine csub

       !--------------------------NITROGEN--------------------------------

       subroutine nmv(amount,from,to)
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

         call ncp(amount,to)
         call nsub(amount,from)

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

         from%n14 = from%n14 - amount%n14

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


       function orgadd(pool1,pool2)
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable sum of two pools of type 
         !  'orgpool'. Sum is of type 'orgpool' as well.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (orgpool), intent(in) :: pool1, pool2
         
         ! function return variable
         type (orgpool) :: orgadd

         orgadd%c = cadd(pool1%c,pool2%c)
         orgadd%n = nadd(pool1%n,pool2%n)

       end function orgadd


       function cadd(pool1,pool2)
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable sum of two pools of type 
         !  'carbon'. Sum is of type 'carbon' as well.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (carbon), intent(in) :: pool1, pool2
         
         ! function return variable
         type (carbon) :: cadd

         cadd%c12 = pool1%c12 + pool2%c12

       end function cadd


       function nadd(pool1,pool2)
         !////////////////////////////////////////////////////////////////
         !  Generic function to return variable sum of two pools of type 
         !  'nitrogen'. Sum is of type 'nitrogen' as well.
         !----------------------------------------------------------------
         implicit none

         ! arguments
         type (nitrogen), intent(in) :: pool1, pool2
         
         ! function return variable
         type (nitrogen) :: nadd

         nadd%n14 = pool1%n14 + pool2%n14

       end function nadd


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
