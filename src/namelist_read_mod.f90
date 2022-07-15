module namelist_read_mod

implicit none

contains

subroutine read_namelist_as_str(namelist_str, namelist_fname, mpi_id, master_id)
    use mpi

    character(:),    intent(out), allocatable :: namelist_str
    character(*),    intent(in) :: namelist_fname
    integer(kind=4), intent(in) :: mpi_id
    integer(kind=4), intent(in), optional :: master_id

    integer(kind=4) :: namelist_len, unit_nam
    integer :: master_id_loc, ierr

    master_id_loc = 0
    if(present(master_id)) master_id_loc = master_id

    if (master_id_loc==mpi_id) then
        print *, "----reading namelist: ", namelist_fname," ----------"
        inquire(file=namelist_fname,size=namelist_len)
        open(newunit=unit_nam, file=namelist_fname, action="read", &
             form="unformatted",access="stream")

        allocate(character(namelist_len) :: namelist_str)
        read(unit_nam) namelist_str

        close(unit_nam)
        print '(A)', namelist_str !I think it will be usefull for future model use
                               !as it will allow to get used namelists from log file
        print *, "-----namelist end---------------"
    end if

    !broadcast namelist length
    call MPI_BCAST(namelist_len, 1, MPI_INT, master_id_loc, MPI_COMM_WORLD, ierr)

    if (master_id_loc /= mpi_id) allocate(character(namelist_len) :: namelist_str)

    call MPI_BCAST(namelist_str, namelist_len, MPI_CHARACTER, master_id_loc, MPI_COMM_WORLD, ierr)

end subroutine read_namelist_as_str

end module namelist_read_mod
