module send_rec_module
use mpi
implicit none
contains

subroutine send_recv_array(dats, myFirstPart, myLastPart, ierror)
use mpi
implicit none
integer, intent(in)                                     :: ierror
integer, intent(in)                                     :: myFirstPart, myLastPart
real(8), dimension(:,:), intent(inout)                  :: dats
real(8), optional, dimension(:,:), intent(out)          :: full_dats
integer, parameter                                      :: rMaster = 0
integer                                                 :: i, nDats

if (rank == rMaster) then
        if (.not.present(full_dats)) then
                print*, "ERROR WHILE CALLING send_recv_array()"
                print*, "'full_dats' not given for the MASTER THREAD"
                print*, "EXITTING PROGRAM"
                call exit()
        end if

        if (.not.present()) then
                print*, "ERROR WHILE CALLING send_recv_array()"
                print*, "'nTotDats' not given for the MASTER THREAD"
                print*, "EXITTING PROGRAM"
                call exit()
        end if

        do i = 0, nProcs - 1, 1
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
                & nDats,                &
                & mpi_real,             &
                & i,                    &
                & mpi_any_tag,          &
                & mpi_comm_world,       &
                & status,               &
                & ierror
                &)
                full_dats(procFirstPart:procLastPart,:) = dats(:,:)
                
        end if
else if (rank /= rMaster) then
        nDats = myLastPart - myFirstPart
        call mpi_send(&
        & dats,                 &
        & nDats,                &
        & mpi_real,             &
        & rMaster,              &
        & 2001,                 &
        & mpi_comm_world,       &
        & ierror                &
        &)

end if
end subroutine send_recv
end module send_rec_module




