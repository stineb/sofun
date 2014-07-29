       module definitions
       type vector
          real :: v(3)
       end type vector

       TYPE point
              REAL :: X, Y
       END TYPE point

       TYPE circle
              TYPE (point) :: Center
              REAL :: Radius
       END TYPE circle

       end module definitions

       !==================================

       program main

       use definitions
       implicit none

       call define

       end program main

       !==================================

       subroutine define
       use definitions
       implicit none
       TYPE (circle) :: C
       C = circle(point(0., 0.), 1.)
       print*,'C',C
       end subroutine define


!       function vadd(v1,v2) result (v3)

!              implicit none

!              type vector
!                 real*8 :: v(3)
!              end type vector
              
!              type (vector), intent(in) :: v1,v2
!              type (vector)  :: v3
              
!              v3%v = v1%v + v2%v ! Note the use of array operations
              
!       end function vadd
