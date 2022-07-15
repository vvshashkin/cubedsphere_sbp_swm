module mpi_paneled_outputer_mod

! use outputer_abstract_mod, only : outputer_t
! use grid_field_mod,        only : grid_field_t
! use partition_mod,         only : partition_t
! use mpi
!
!     implicit none
!
! type, public, extends(outputer_t) :: mpi_paneled_outputer_t
!
!     integer(kind=4)               :: rec_num = 1
!     integer(kind=mpi_offset_kind) :: record_disp ! displacement in bytes between records
!     integer(kind=4), allocatable  :: type_local(:), type_global(:) ! mpi subarrays
!
! contains
!
!     procedure, public :: write => mpi_io_write
!
! end type mpi_paneled_outputer_t
!
! contains
!
! subroutine mpi_io_write(this, f, partition, file_name, rec_num)
!
!     class(mpi_paneled_outputer_t), intent(inout) :: this
!     type(grid_field_t),            intent(inout) :: f
!     type(partition_t),             intent(in)    :: partition
!     character(*),                  intent(in)    :: file_name
!     integer(kind=4),               intent(in), &
!                                    optional      :: rec_num
!
!     real(kind=4), allocatable :: buf(:,:,:,:)
!     integer(kind=mpi_offset_kind) :: disp
!     integer(kind=4) :: file_info, file_handle
!     integer(kind=4) :: t, ierr, irec
!     integer(kind=4) :: is, ie, js, je, ks, ke
!
!     irec = 1
!
!     if (present(rec_num)) irec=rec_num
!
!     call MPI_Info_create(file_info, ierr)
!     call mpi_file_open(mpi_comm_world, file_name, &
!         MPI_MODE_WRONLY + MPI_MODE_CREATE, file_info, file_handle, ierr)
!
!     do t = partition%ts, partition%te
!
!         disp = int((irec-1)*this%record_disp,kind=mpi_offset_kind)
!
!         call MPI_File_set_view(file_handle, disp, mpi_real4, this%type_global(t), 'native',file_info, ierr)
!
!         is = partition%tile(t)%is; ie = partition%tile(t)%ie
!         js = partition%tile(t)%js; je = partition%tile(t)%je
!         ks = partition%tile(t)%ks; ke = partition%tile(t)%ke
!
!         allocate(buf(is:ie,js:je,1,ks:ke))
!         buf(is:ie,js:je,1,ks:ke) = real(f%block(t)%p(is:ie,js:je,ks:ke),4)
!
!         call MPI_File_write(file_handle, buf, 1, this%type_local(t), mpi_status_ignore, ierr)
!
!         deallocate(buf)
!
!     end do
!
!     call MPI_File_close(file_handle, ierr)
!
!
! end subroutine mpi_io_write

end module mpi_paneled_outputer_mod
