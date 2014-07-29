!      Main program
       PROGRAM try

!      DERIVED.F90
!      Define a derived-type structure,
!      type variables, and assign values
       
       TYPE member
         INTEGER age
         CHARACTER (LEN = 20) name
       END TYPE member

       TYPE (member) :: george
       TYPE (member) :: ernie
       
       george     = member( 33, 'George Brown' )
       ernie%age  = 56
       ernie%name = 'Ernie Brown'
       
       WRITE (*,*) george
       WRITE (*,*) ernie
 
!         CALL outer

       END PROGRAM try

       SUBROUTINE outer
          REAL x,y                  ! declaration of 'y' is not actually necessary
          
          x = 5

          CALL inner
          
          PRINT*, 'outer: x =',x
          PRINT*, 'outer: y =',y    ! value y is not defined inside outer 

       CONTAINS
          SUBROUTINE inner
             REAL y
             
             y = x + 1.

             PRINT*, 'inner: x =',x  ! inner has access to 'x' by host association
             PRINT*, 'inner: y =',y
             
          END SUBROUTINE inner       ! SUBROUTINE mandatory
       END SUBROUTINE outer
