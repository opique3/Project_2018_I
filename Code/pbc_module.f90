!------------------------
!MADE BY ORIOL PIQUÉ
!-------------------------

module pbc_module
implicit none
contains

!This module contains the subroutine that implements periodic boundary conditions when calculating the potential.


!Variables(in):

!Size of the simulation box (boxSize)
!Potential (V)


!Variables(out):

!Potential (V)


subroutine pbc(V, boxSize)
implicit none
real(8), intent(in)                             :: boxSize
real(8), dimension(3), intent(inout)            :: V

V(:) = V(:) - nint(V(:)/boxSize)*boxSize
end subroutine pbc
end module pbc_module
