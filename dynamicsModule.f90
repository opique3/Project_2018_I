module dynamicsModule
implicit none
contains

subroutine readData(un, dt, boxSize, cutOff, nPartDim, T, eps, sig, nSteps, seed)
implicit none
real(8), intent(out)                            :: dt, boxSize, cutOff, T
real(8), intent(out)                            :: eps, sig
integer, intent(out)                            :: nPartDim, un, nSteps, seed
real(8)                                         :: density
read(un,*) dt                   ! Increment de temps
!read(un,*) boxSize              ! Tamany de la caixa de simulació
read(un,*) density              ! Densitat en unitats adimensionals 
read(un,*) cutOff               ! Distància màxima a la que es calcula LJ
read(un,*) nPartDim             ! Numero total de partícules en cada dimensio
read(un,*) T                    ! Temperatura de la simulació
read(un,*) eps                  ! Epsilon per al potencial LJ
read(un,*) sig                  ! Sigma per al potencial LJ
read(un,*) nSteps               ! Numero de passos de temps que fa la simulació
read(un,*) seed                 ! llavor per la generació de numeros aleatoris.
close(un)

! D = N/V --> V = N/D --> L = (N/D)**(1/3)
boxSize = nPartDim/(density)**(1./3.)
print*, 'BS = ', boxSize

end subroutine readData

subroutine print_positions(un, nPart, pos, time)
implicit none
integer, intent(in)                             :: un, nPart
real(8), dimension(nPart,3), intent(in)         :: pos
real(8), intent(in)                             :: time
integer                                         :: i

write(un,*) nPart
write(un,*) 'Simulation Time = ', time,'seconds'
do i = 1, nPart, 1
        write(un,*) 'C', pos(i,:)
end do

end subroutine print_positions

subroutine SC_init_conditions(nPart, pos, boxSize)
implicit none
integer, intent(in)                             :: nPart
real(8), intent(in)                             :: boxSize
real(8), dimension(nPart,3), intent(out)        :: pos
real(8)                                         :: dist
real(8), dimension(3)                           :: e1, e2, e3, r
integer                                         :: i, j, k

! La idea aqui es construïr els vectors directors que generen la cela SC
! on les seves components coincideixin amb les posicions dels atoms a la cel·la
! unitat
dist = boxSize/dfloat(nPart)**(1./3.)
r  = (/dist, dist, dist/) 
e1 = (/dist, 0.0D0, 0.0D0/)
e2 = (/0.0D0, dist, 0.D0/)
e3 = (/0.0D0, 0.0D0, dist/)

! Aquesta forma de col·locar les particules en una cel·la SC
! fa que el numero total de particules hagi de ser multiple de 8
! ja que son el les que hi ha en la cel·la unitat.
do i = 0, nPart - 8, 8
        pos(i+1,:) = r
        pos(i+2,:) = e1 + r
        pos(i+3,:) = e2 + r
        pos(i+4,:) = e3 + r
        pos(i+5,:) = e1 + e2 + r
        pos(i+6,:) = e1 + e3 + r
        pos(i+7,:) = e2 + e3 + r
        pos(i+8,:) = e1 + e2 + e3 + r

       r(1) = r(1) + 2*dist
       if (r(1) > boxSize) then
               r(1) = dist
               r(2) = r(2) + 2*dist
       end if
       if (r(2) > boxSize)  then
               r(2) = dist
               r(3) = r(3) + 2*dist
       end if
end do

end subroutine SC_init_conditions

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

! Subrutina per les condicions periodiques de contorn per una sola particula.
! S'utilitza en el calcul del potencial LJ de la subrutina posterior.
subroutine pbc(v, boxSize)
implicit none
real(8), intent(in)                             :: boxSize
real(8), dimension(3), intent(inout)            :: v

v(:) = v(:) - nint(v(:)/boxSize)*boxSize
end subroutine pbc

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

subroutine boxMuller_polar(seed, vel, nPart, T)
implicit none
integer, intent(in)                             :: seed, nPart
real(8), dimension(nPart,3), intent(out)        :: vel
real(8), intent(in)                             :: T
real(8)                                         :: x1, x2, y1, y2, w
integer                                         :: i, j, k

call srand(seed)
do i = 1, nPart/2, 1; do j = 1, 3, 1
        w = 1.0D0
        do while (w >= 1.0D0)
                x1 = 2.0*rand() - 1.0
                x2 = 2.0*rand() - 1.0
                w = x1**2. + x2**2.
        end do
        w  = dsqrt(-2.0*dlog(w)/w)
        y1 = x1*w
        y2 = x2*w
        vel(i,j) = y1*dsqrt(T)
        vel(i + nPart/2,j) = y2*dsqrt(T)
end do; end do

end subroutine boxMuller_polar

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

subroutine andersen_thermo(dt, T, nPart, seed, vel, eps)
implicit none
integer, intent(in)                             :: nPart, seed
real(8), intent(in)                             :: dt, T, eps
real(8), dimension(nPart,3), intent(inout)      :: vel
real(8)                                         :: x1, x2, y1, y2, w
integer                                         :: i, j, k

call srand(seed)
do i = 1, nPart, 1
        if (rand() < 0.1) then
                ! canvia la velocitat d'aquesta particula segons una distribució
                ! gaussiana com maxwell boltzman.
                ! Aixó ho fem utilitzant el algoritme Box_muller, pero com
                ! aquest genera parells de numeros aleatoris i tenim tres
                ! components en el vector de velocitats, podem fer dues vegades
                ! el BM i descartar un dels numeros aleatoris.
                w = 1.0D0
                do while (w >= 1.0D0)
                        x1 = 2.0*rand() - 1.0
                        x2 = 2.0*rand() - 1.0
                        w = x1**2. + x2**2.
                end do
                w  = dsqrt(-2.0*dlog(w)/w)
                y1 = x1*w
                y2 = x2*w
                vel(i,1) = y1*dsqrt(T)
                vel(i,2) = y2*dsqrt(T)
                w = 1.0D0
                do while (w >= 1.0D0)
                        x1 = 2.0*rand() - 1.0
                        x2 = 2.0*rand() - 1.0
                        w = x1**2. + x2**2.
                end do
                w  = dsqrt(-2.0*dlog(w)/w)
                y1 = x1*w
                vel(i,3) = y1*dsqrt(T)
        end if
end do
end subroutine andersen_thermo

end module dynamicsModule



















