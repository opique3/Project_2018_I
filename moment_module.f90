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
