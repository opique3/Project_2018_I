module velocities_module
implicit none
contains

subroutine IN_velocities(N,T,seed,vel)
implicit none
INTEGER, INTENT(IN) :: N, seed
REAL(8), INTENT(IN) :: T
REAL(8), DIMENSION(N,3), INTENT(out) :: vel
INTEGER :: i
REAL(8) :: x1, x2, w, sigma

call srand(seed)
sigma = sqrt(T)

Do i = 1,N,1
  Do while(w .ge. 1.0 .or. w .eq. 0.0)
    x1 = 2.*rand() - 1.
    x2 = 2.*rand() - 1.
    w = x1*x1 + x2*x2
  end do

  w = sigma*dsqrt((-2.*log(w))/w)
  vel(i,1) = x1*w
  vel(i,2) = x2*w

  Do while(w .ge. 1.0 .or. w .eq. 0.0)
    x1 = 2.*rand() - 1.
    x2 = 2.*rand() - 1.
    w = x1*x1 + x2*x2
  end do

  w = sigma*dsqrt((-2.*log(w))/w)
  vel(i,3) = x1*w

end do

end subroutine IN_velocities

end module velocities_module
