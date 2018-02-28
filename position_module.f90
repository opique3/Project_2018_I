module positions
implicit none
contains

subroutine IN_positions(M,N,density,pos)
implicit none
INTEGER, INTENT(IN) :: M, N
REAL(8), INTENT(IN) :: density
REAL(8), DIMENSION(N,3), INTENT(OUT) :: pos
INTEGER :: ix, iy, iz, num
REAL(8) :: a, a2, L
REAL(8), DIMENSION(3) :: e1, e2, e3

L = (N/density)**(1./3.)
a = L/M
e1 = (/ a, 0.0D0, 0.0D0 /)
e2 = (/ 0.0D0, a, 0.0D0 /)
e3 = (/ 0.0D0, 0.0D0, a /)

Do iz = 1,M,1
  Do iy = 1,M,1
    Do ix = 1,M,1
        num = ix + (iy-1)*M + (iz-1)*M*M
        pos(num1,:) = ix*e1(:) + iy*e2(:) + iz*e3(:)
    end do
  end do
end do

end subroutine IN_positions

end module positions
