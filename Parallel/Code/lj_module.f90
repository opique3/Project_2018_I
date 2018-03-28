!------------------------
!MADE BY ORIOL PIQUÉ
!-------------------------

module lj_module
use pbc_module
use mpi
implicit none
contains

!This module contains the subroutine that calculates the resulting forces that appear because of the LJ potential interaction between the particles.


!Variables(in):

!Number of particles (nPart)
!Positions array (pos)
!Epsilon (eps)
!Sigma (sig)
!Size of the simulation box (boxSize)
!Cutoff distance (cutOff)


!Variables(out):

!Forces array (F)
!Potential (V)


subroutine LJ_pot(nPart, pos, eps, sig, boxSize, cutOff, F, V_partial, myFirstPart, myLastPart)
implicit none
integer, intent(in)                             :: nPart
real(8), dimension(:,:), intent(in)             :: pos
real(8), intent(in)                             :: eps, sig, boxSize, cutOff
real(8), dimension(:,:), intent(out)            :: F
real(8)                                         :: V_partial
real(8), dimension(3)                           :: dist
real(8)                                         :: rij, dV
integer                                         :: i, j, k
integer						:: ierror
integer, parameter				:: rMaster = 0

V_partial = 0.
F(:,:) = 0.
do i = myFirstPart, myLastPart, 1; do j = 1, nPart, 1
        dist(:) = pos(i,:) - pos(j,:)
        call pbc(dist, boxSize)
        rij = dsqrt(dot_product(dist,dist))
        if (rij < cutOff) then
                dist(:) = dist(:)/rij
                V_partial  = V_partial + 4.*eps*((sig/rij)**12. - (sig/rij)**6.)
                dV = 4*eps*(12.*sig**12./rij**13. - 6.*sig**6./rij**7)
                F(i,:) = F(i,:) + dV*dist(:)
        end if
end do; end do

call mpi_reduce(V_partial, V, 1, mpi_real8, mpi_sum, rMaster, mpi_comm_world, ierror) 

end subroutine LJ_pot
end module lj_module
