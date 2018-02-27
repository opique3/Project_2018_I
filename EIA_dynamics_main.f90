program eia_dynamics
use moment
use positions
use velocities
use andersen_thermostate
use Verlet_algorithm
use 
use
use

implicit none
!timing variables
real :: start_time, end_time


call cpu_time(start_time)



!First read the parameters of the simulation
call read_input()

!Initialize the system
call init_positions()
call init_velocity

!Melt the crystal to obtain a liquid
do i = 1, nint(0.1*num_iterations) 


!

call cpu_time(end_time)
print *, 'Execution time = ', end_time-start_time
end program eia_dynamics