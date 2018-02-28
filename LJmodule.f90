module lj_module
use pbc_module
implicit none
contains

!Subrutina pel càlcul de les forces conseqüència de la interacció LJ entre les nostres partícules.
subroutine LJ_pot(nPart, pos, eps, sig, boxSize, cutOff, F, V)
implicit none
integer, intent(in)                             :: nPart
real(8), dimension(nPart,3), intent(in)         :: pos
real(8), intent(in)                             :: eps, sig, boxSize, cutOff
real(8), dimension(nPart,3), intent(out)        :: F
real(8)                                         :: V
real(8), dimension(3)                           :: dist
real(8)                                         :: rij, dV
integer                                         :: i, j, k

V = 0.
F(:,:) = 0.
do i = 1, nPart, 1; do j = i + 1, nPart, 1
        dist(:) = pos(i,:) - pos(j,:)
        call pbc(dist, boxSize)
        rij = dsqrt(dot_product(dist,dist))
        if (rij < cutOff) then
                dist(:) = dist(:)/rij
                V  = V + 4.*eps*((sig/rij)**12. - (sig/rij)**6.)
                dV = 4*eps*(12.*sig**12./rij**13. - 6.*sig**6./rij**7)
                F(i,:) = F(i,:) + dV*dist(:)
                F(j,:) = F(j,:) - dV*dist(:)
        end if
end do; end do
end subroutine LJ_pot
end module lj_module
