module genis_module
implicit none
contains

subroutine IN_positions(M,N,density,pos)
implicit none
INTEGER, INTENT(IN) :: M, N
REAL(8), INTENT(IN) :: density
REAL(8), DIMENSION(N,3), INTENT(OUT) :: pos
INTEGER :: NN, ix, iy, iz, num1, num2
REAL(8) :: a, a2, L
REAL(8), DIMENSION(3) :: e1, e2, e3

L = (N/density)**(1./3.)
a = L/M
a2 = a/2.
NN = N/2

e1 = (/ a, 0.0D0, 0.0D0 /)
e2 = (/ 0.0D0, a, 0.0D0 /)
e3 = (/ 0.0D0, 0.0D0, a /)

Do iz = 1,M,1
  Do iy = 1,M,1
    Do ix = 1,M,1
        num1 = ix + (iy-1)*M + (iz-1)*M*M
        num2 = num1 + NN
        pos(num1,:) = ix*e1(:) + iy*e2(:) + iz*e3(:)
        pos(num2,:) = pos(num1,:) + a2*e1(:) + a2*e2(:) + a2*e3(:)
    end do
  end do
end do

end subroutine IN_positions

subroutine IN_velocities(N,T,seed,vel)
implicit none
INTEGER, INTENT(IN) :: N, seed
REAL(8), DIMENSION(N,3), INTENT(IN) :: vel
INTEGER :: i
REAL(8) :: x1, x2, w, sigma

call srand(seed)
sigma = dsqrt(T)

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


end module genis_module
