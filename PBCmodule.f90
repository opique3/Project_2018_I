!--------------------
!MADE BY ORIOL PIQUÉ
!--------------------

module PBCmodule
implicit none
contains

!Subrutina per les condicions periòdiques de contorn per una sola partícula.
subroutine pbc(v, boxSize)
implicit none
real(8), intent(in)                             :: boxSize
real(8), dimension(3), intent(inout)            :: v

v(:) = v(:) - nint(v(:)/boxSize)*boxSize
end subroutine pbc
