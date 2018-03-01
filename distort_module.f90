module distort_module
use pbc_module
implicit none
contains

subroutine distort_geometry(nPart, pos, boxSize, seed)
implicit none
integer, intent(in)                             :: nPart, seed
real(8), dimension(nPart,3), intent(inout)      :: pos
real(8), intent(in)                             :: boxSize
real(8)                                         :: randN, maxDisp
real(8), dimension(3)                           :: vec
integer                                         :: i, j, k

call srand(seed)
! La distortion es, com a maxim, un 1/3 de la distancia maxima
! que hi ha entre les particules veines.
! S'aplica a cada dimensió
maxDisp = boxSize/(5*dfloat(nPart)**(1./3.))
do i = 1, nPart, 1
        do j = 1, 3, 1
                randN = 2*rand() - 1. ! Numero entre -1 i 1
                randN = maxDisp*randN
                pos(i,j) = pos(i,j) + randN
                if (pos(i,j) > boxSize) pos(i,j) = pos(i,j) - boxSize
                if (pos(i,j) < 0.)      pos(i,j) = pos(i,j) + boxSize
        end do
        vec(:) = pos(i,:)
        call pbc(vec, boxSize)
        pos(i,:) = vec(:)
end do
end subroutine distort_geometry

end module distort_module

