program dynamics
!use the necessary modules
use moment_module
use positions_module
use velocities_module
use andersen_therm_module
use vel_verlet_module
use lj_module
use pbc_module
use read_data_module
use print_positions_module
use kinetic_energy_module
use print_data_module

implicit none
character(25)                           :: fName, fff
integer                                 :: i, j, k, l, m, n
integer                                 :: fStat, un, num_nodes
real(4)                                 :: start, finish
real(8)                                 :: dt, boxSize, cutOff, T, density, dtAux
integer                                 :: nPartDim, nPart, nSteps
real(8), allocatable, dimension(:,:)    :: pos, F, vel
real(8)                                 :: V, eps, sig, time, KE, Tinst
real(8), allocatable, dimension(:)      :: total_momentum
integer                                 :: seed, trjCount, thermCount
integer                                 :: initUn, finUn, trajUn, dataUn, velUn
    
call cpu_time(start)
call get_command_argument(1, fName, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program'
        call exit()
end if 
un = 100
open(unit=un, file=trim(fName), status='old') 
    
! Llegir les dades del fitxer de entrada i alocatar la variable que contindrà
! la posició de totes les partícules.
call readData(un, dt, boxSize, cutOff, nPartDim, T, eps, sig, nSteps, density, seed)
nPart = nPartDim**3
allocate(pos(nPart,3), F(nPart,3), vel(nPart,3), total_momentum(3))
if (mod(float(nPart),8.) /= 0.) then
        print*, 'Total number of particles must be'
        print*, 'divisible by 8 to fit a SC cell'
        print*, 'Exitting Program'
        call exit()
end if
    
    
! Generar la posició inicial de totes les partícules conforme a un cristall 
! d'una cel·la bcc. I fondre el cristall amb una temperatura molt alta
num_nodes = int(nPart**(1/3))
call IN_positions(num_nodes, nPart, density, pos)
call IN_velocities(nPart, T, seed, vel)
dtAux = 0.001*dt
time = 0.0D0
do i = 1, int(0.001*nSteps)
    call velocity_verlet(time, dtAux, pos, vel, nPart, eps, sig, boxSize, cutOff, V, F)
    call andersen_thermo(dt, 100.0d0, nPart, i, vel)
end do

!call distort_geometry(nPart, pos, boxSize, seed)
! Genera una distribució de velocitats gaussiana amb l'algoritme de Box_Muller
call IN_velocities(nPart, T, seed, vel)
    
! Obre els arxius on s'imprimiran els resultats de la simulació
initUn = 101; finUn = 102; trajUn = 103; dataUn = 104; velUn = 105
open(unit=initUn, file='initial.out')   ! Coord. Inicials
open(unit=finUn , file='final.out')     ! Coord. Finals
open(unit=trajUn, file='traj.xyz')      ! Trajectoria
open(unit=dataUn, file='data.out')     ! T, Ken, V, temps...
open(unit=velUn,  file='velocity.out')     ! T, Ken, V, temps...

call print_positions(initUn, nPart, pos, time)
! Començem la MD inicialitzant el temps a 0 i fent el numero de pasos
! especificats
time = 0.0D0; trjCount = 0; thermCount = 0
write(velUn, *) nSteps/100 - 1, nPart
write(trajUn,*) nSteps/100 - 1, nPart
call print_positions(trajUn, nPart, pos, time)

call print_positions(velUn,  nPart, vel, time)

do i = 1, nSteps, 1
        call velocity_verlet(time, dt, pos, vel, nPart, eps, sig, boxSize, cutOff, V, F)
        if (trjCount == 100) then
                call print_positions(trajUn, nPart, pos, time)
                call print_positions(velUn,  nPart, vel, time)
                call kinetic_energy(vel, KE, Tinst, nPart)
                call momentum(nPart, vel, total_momentum)
                call print_data(time, V, KE, Tinst, total_momentum, dataUn)
                trjCount = 0
        end if
        if (thermCount == 10) then
                call andersen_thermo(dt, T, nPart, i, vel)
                thermCount = 0
        end if
        trjCount = trjCount + 1
        thermCount = thermCount + 1
end do
call print_positions(finUn, nPart, pos, time)


call cpu_time(finish)
write(dataUn,*) "# CPU TIME: ", finish - start
close(un); close(initUn); close(finUn); close(trajUn); close(dataUn)
contains

end program dynamics
    
