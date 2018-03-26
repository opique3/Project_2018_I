!:::::::::::::::::::::::::
! MADE BY GENIS LLEOPART
!:::::::::::::::::::::::::

! This module contain the subroutine positions, which
! generate the initial positions of the particles 
! following a Simple Cubic Geometry

! Variables IN:
!   nPart --> Number of particles | INTEGER
!   boxSize --> Length of our box simulation | REAL(8)

! Variables OUT:
!   pos --> Particle positions | DIMENSION(nPart,3), REAL(8)

module positions_module
implicit none
contains

subroutine SC_init_conditions(nPart, pos, boxSize)
implicit none
integer, intent(in)                             :: nPart
real(8), intent(in)                             :: boxSize
real(8), dimension(nPart,3), intent(out)        :: pos
real(8)                                         :: dist
real(8), dimension(3)                           :: e1, e2, e3, r
integer                                         :: i, j, k

! Generate the length of the basis box and the basis vectors
dist = boxSize/dfloat(nPart)**(1./3.)
r  = (/dist, dist, dist/) 
e1 = (/dist, 0.0D0, 0.0D0/)
e2 = (/0.0D0, dist, 0.D0/)
e3 = (/0.0D0, 0.0D0, dist/)

! Generate structures of 8 particles, following a SC
do i = 0, nPart - 8, 8
        pos(i+1,:) = r
        pos(i+2,:) = e1 + r
        pos(i+3,:) = e2 + r
        pos(i+4,:) = e3 + r
        pos(i+5,:) = e1 + e2 + r
        pos(i+6,:) = e1 + e3 + r
        pos(i+7,:) = e2 + e3 + r
        pos(i+8,:) = e1 + e2 + e3 + r

       r(1) = r(1) + 2*dist
       if (r(1) > boxSize) then
               r(1) = dist
               r(2) = r(2) + 2*dist
       end if
       if (r(2) > boxSize)  then
               r(2) = dist
               r(3) = r(3) + 2*dist
       end if
end do

end subroutine SC_init_conditions

end module positions_module
