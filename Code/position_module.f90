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

! La idea aqui es construïr els vectors directors que generen la cela SC
! on les seves components coincideixin amb les posicions dels atoms a la cel·la
! unitat
dist = boxSize/dfloat(nPart)**(1./3.)
r  = (/dist, dist, dist/) 
e1 = (/dist, 0.0D0, 0.0D0/)
e2 = (/0.0D0, dist, 0.D0/)
e3 = (/0.0D0, 0.0D0, dist/)

! Aquesta forma de col·locar les particules en una cel·la SC
! fa que el numero total de particules hagi de ser multiple de 8
! ja que son el les que hi ha en la cel·la unitat.
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
