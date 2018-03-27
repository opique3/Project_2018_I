module send_rec_module
use mpi
implicit none
contains

subroutine send_recv_array(dats, myFirstPart, myLastPart, ierror, nPart)
implicit none
integer, intent(in)                                     :: ierror, nPart
integer, intent(in)                                     :: myFirstPart, myLastPart
real(8), dimension(:,:), intent(inout)                  :: dats
real(8), optional, dimension(:,:), intent(out)          :: full_dats
integer, parameter                                      :: rMaster = 0
integer                                                 :: i, nDats, numProcs

call mpi_comm_size(mpi_comm_world, numProcs, ierror)
numParts = nPart/numProcs
if (rank == rMaster) then
        if (.not.present(full_dats)) then
                print*, "ERROR WHILE CALLING send_recv_array()"
                print*, "'full_dats' not given for the MASTER THREAD"
                print*, "EXITTING PROGRAM"
                call exit()
        end if

        do i = 0, numProcs - 1, 1
                procFirstPart = i*numParts + 1
                if (i /= numProcs - 1) procLastPart = procFirstPart + numParts - 1
                if (i == numProcs - 1) procLastPart = nPart
                if (i == rMaster) then
                        full_dats(procFirstPart:procLastPart,:) = dats(:,:)
                        cycle
                end if
                nDats = myLastPart - myFirstPart
                call mpi_recv(&
                & dats,                 &
                & 3*nDats,              &
                & mpi_real,             &
                & i,                    &
                & mpi_any_tag,          &
                & mpi_comm_world,       &
                & status,               &
                & ierror
                &)
                full_dats(procFirstPart:procLastPart,:) = dats(:,:)
                
        end do
else if (rank /= rMaster) then
        nDats = size(dats)
        call mpi_send(&
        & dats,                 &
        & 3*nDats,              &
        & mpi_real,             &
        & rMaster,              &
        & 2001,                 &
        & mpi_comm_world,       &
        & ierror                &
        &)

end if
end subroutine send_recv
end module send_rec_module

