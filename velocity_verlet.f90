module Verlet_algorithm
contains
subroutine velocity_verlet(time, dt, pos, vel, nPart, eps, sig, boxSize, cutOff, V, F)
implicit none
integer, intent(in)                             :: nPart
real(8), intent(in)                             :: dt, eps, sig, boxSize, cutOff
real(8), intent(inout)                          :: time
real(8), intent(out)                            :: V
real(8), dimension(nPart,3), intent(inout)      :: pos, vel, F
real(8), dimension(nPart,3)                     :: F_aux
integer                                         :: i, j, k
real(8), dimension(3)                           :: vec

if (time == 0) call LJ_pot(nPart, pos, eps, sig, boxSize, cutOff, F_aux, V)
if (time /= 0) F_aux(:,:) = F(:,:)
time = time + dt

do i = 1, nPart, 1
        vec(:) = pos(i,:) + vel(i,:)*dt + F_aux(i,:)*dt**3./2.
        call pbc(vec, boxSize)
        pos(i,:) = vec(:)
end do

call LJ_pot(nPart, pos, eps, sig, boxSize, cutOff, F, V)

do i = 1, nPart, 1
        vel(i,:) = vel(i,:) + (F_aux(i,:) + F(i,:))*dt/2.
end do
end subroutine velocity_verlet

endmodule
