module test_namelist_mod

implicit none

contains

subroutine test_namelist()

    use mpi
    use namelist_read_mod, only : read_namelist_as_str

    character(*), parameter :: namelist_text = "$TST1"        // new_line("A") // &
                                              "ipar1 = 57,"   // new_line("A") // &
                                              "fpar1 = 5.7"   // new_line("A") // &
                                              "$END"
    character(*), parameter :: namelist_fname = ".test.namelist"
    character(:), allocatable :: namstr_as_read
    integer(kind=4) ipar1, ipar1_master
    real(kind=8)    fpar1, fpar1_master
    namelist /tst1/ ipar1, fpar1

    integer myid, Np, ierr

    integer unit_nam_wr

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    if (myid==0) then
        open(newunit=unit_nam_wr, file=namelist_fname, action="write", &
             form="unformatted",access="stream")
        write(unit_nam_wr) namelist_text
        close(unit_nam_wr)
    end if

    call read_namelist_as_str(namstr_as_read,namelist_fname, myid, master_id = 0)

    read(namstr_as_read, tst1)

    if(myid == 0) then
        print *, "namelist params:", ipar1, fpar1
        ipar1_master = ipar1
        fpar1_master = fpar1
    end if

    call MPI_BCAST(ipar1_master,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fpar1_master,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

    ipar1 = abs(ipar1 - ipar1_master)
    fpar1 = abs(fpar1 - fpar1_master)

    call MPI_REDUCE(ipar1,ipar1_master,1,MPI_INT,   MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(fpar1,fpar1_master,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then
        if(ipar1_master == 0 .and. fpar1_master < 1e-16 .and. ierr==MPI_SUCCESS) then
            print *, "namelist read test passed"
        else
            print *, "namelist read_test_failed"
        end if
        !delete test namelist file
        open(newunit=unit_nam_wr, file=namelist_fname, action="write")
        close(unit_nam_wr,status="delete")
    end if

end subroutine test_namelist

end module test_namelist_mod
