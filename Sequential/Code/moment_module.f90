!::::::::::::::::::::::::
! MADE BY GENIS LLEOPART
!::::::::::::::::::::::::

! Contains the subroutine momentum, which calculate
! the total intertial momentum of our system

! Variables IN:
!   N 	--> Number of particles | INTEGER
!   vel --> Velocities of the particles | DIMENSION(N,3), REAL(8)

! Variables OUT:
!   P_total --> Total momentum | DIMENSION(3), REAL(8)
module moment_module
implicit none
contains

subroutine momentum(N,vel,P_total)
implicit none
INTEGER, INTENT(IN) :: N
REAL(8), DIMENSION(N,3), INTENT(IN) :: vel
REAL(8), DIMENSION(3), INTENT(OUT) :: P_total
INTEGER :: i

P_total = 0.0D0
Do i = 1,N,1
  P_total(:) = P_total(:) + vel(i,:)
end do

end subroutine momentum


end module moment_module
