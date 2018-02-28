program test_pos
use positions
implicit none
integer :: num_part, num_nodes, i
real(8) :: density
real(8), allocatable, dimension (:,:) :: pos

open(unit=99, file='initial_positions.xyz')


num_part = 9**3
num_nodes = (num_part*0.5)**(1./3.)
print *, num_nodes
num_nodes = int(num_nodes)
density = 0.88

allocate(pos(num_part,3))
pos = 0.0

call IN_positions(num_nodes, num_part, density, pos)

write(99,*) num_part
write(99,*) ''
do i = 1, num_part
    write(99, *) 'C ' ,pos(i,:)
end do

end program test_pos
