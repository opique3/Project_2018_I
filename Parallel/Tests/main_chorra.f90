program main_chorra
use mpi
use send_rec_module
implicit none
real(8), dimension(:,:), allocatable            :: A, B
integer                                         :: nPart, numProcs, ierror
integer                                         :: myFirstPart, myLastPart
integer                                         :: status, i, rank
integer                                         :: numParts

call mpi_init(ierror)
call mpi_comm_rank(mpi_comm_world, rank, ierror)
call mpi_comm_size(mpi_comm_world, numProcs, ierror)

nPart = 100
numParts = nPart/numProcs



myFirstPart = rank*numParts + 1
myLastPart = myFirstPart + numParts - 1
if (rank == numProcs - 1) myLastPart = nPart
print*, 'nPARTS', rank, numParts, myLastPart - myFirstPart, myFirstPart, myLastPart

allocate(A(nPart,3))
allocate(B(myLastPart - myFirstPart + 1,3))

B(:,:) = 100.0D0
B(:,:) = B(:,:) + rank


if (rank == 0) then
        call send_recv_array(B, myFirstPart, myLastPart, rank, nPart, status, full_dats=A)
else
        call send_recv_array(B, myFirstPart, myLastPart, rank, nPart, status)
end if

if (rank == 0) then
        do i = 1, nPart, 1
                print*, A(i,:)
        end do
end if

call mpi_finalize(ierror)

end program main_chorra
