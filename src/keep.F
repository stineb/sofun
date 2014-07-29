       !===============================================

       module definitions

       implicit none

       type orgpool
         real :: c12
         real :: n14
       end type orgpool

       type plantclass
         type (orgpool) :: leafmass
         type (orgpool) :: rootmass
       end type plantclass

       type soilclass
         type (orgpool) :: fast
         type (orgpool) :: slow
       end type soilclass

       type litterclass
         type (orgpool) :: coarse
         type (orgpool) :: fine
         type (orgpool) :: below
       end type litterclass

       end module definitions

       !===============================================

       program main

       use definitions
       implicit none

       integer, parameter :: nlu = 2
       integer, parameter :: npft = 1
       integer :: lu

       type (orgpool), dimension(nlu) :: dlitter
       type (soilclass), dimension(nlu) :: soil
       type (litterclass), dimension(nlu) :: litter
       type (plantclass), dimension(npft) :: plant

       do lu=1,nlu

         soil(lu)%fast%c12 = 180
         soil(lu)%fast%n14 = 10

         soil(lu)%slow%c12 = 1800
         soil(lu)%slow%n14 = 100

         litter(lu)%fine%c12 = 90
         litter(lu)%fine%n14 = 5

         litter(lu)%coarse%c12 = 90
         litter(lu)%coarse%n14 = 5

         litter(lu)%below%c12 = 90
         litter(lu)%below%n14 = 5

         dlitter(lu)%c12 = 18
         dlitter(lu)%n14 = 1

         print*,'before'
         print*,'soil ', soil(lu)
         print*,'litter ', litter(lu)

         ! finally, adding mass to pool is invoked by
         !soil = dlitter + soil  ! dlitter .moveto. soil
         call orgmv(dlitter(lu),litter(lu)%coarse,soil(lu)%slow)

         call orgmv(dlitter(lu),litter(lu)%fine,soil(lu)%fast)

         print*,'after'
         print*,'soil ', soil(lu)
         print*,'litter ', litter(lu)

       end do

       contains

       subroutine orgcp(amount,to,d13C)

         implicit none
         type (orgpool), intent(in) :: amount
         type (orgpool), intent(inout) :: to
         real, optional, intent(in) :: d13C

         to%c12 = amount%c12 + to%c12
         to%n14 = amount%n14 + to%n14

!         if (present(d13C)) then
!           to%c12 = amount%c12 + to%c12
!           to%n14 = amount%n14 + to%n14
!         else
!           to%c12 = amount%c12 + to%c12
!           to%n14 = amount%n14 + to%n14
!         end if

       end subroutine orgcp


       subroutine orgmv(amount,from,to)

         implicit none
         type (orgpool), intent(in) :: amount
         type (orgpool), intent(inout) :: from
         type (orgpool), intent(inout) :: to

         to%c12 = to%c12 + amount%c12
         to%n14 = to%n14 + amount%n14

         from%c12 = from%c12 - amount%c12
         from%n14 = from%n14 - amount%n14

       end subroutine orgmv

       end program

       !===============================================
