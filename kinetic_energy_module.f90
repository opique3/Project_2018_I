module kintetic_energy_module
contains
subroutine kinetic_energy(vel, KE, Tinst, nPart)
implicit none
integer, intent(in)                             :: nPart
real(8), dimension(nPart,3), intent(in)         :: vel
real(8), intent(out)                            :: KE, Tinst
real(8), dimension(3)                           :: vec
real(8)                                         :: modV
integer                                         :: i

KE = 0.0D0
do i = 1, nPart, 1
        vec(:) = vel(i,:)
        modV = dsqrt(dot_product(vec, vec))
        KE = KE + modV**2.
end do
KE = KE/2.0
Tinst = 2.0*KE/(3.0*float(nPart))
end subroutine kinetic_energy

end module kintetic_energy_module
