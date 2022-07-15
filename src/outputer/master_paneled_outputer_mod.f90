module master_paneled_outputer_mod

use outputer_abstract_mod, only : outputer_t
use grid_field_mod,        only : grid_field_t
use exchange_abstract_mod, only : exchange_t
use domain_mod,            only : domain_t
use tiles_mod,             only : tiles_t

implicit none

type, public, extends(outputer_t) :: master_paneled_outputer_t
    class(exchange_t), allocatable :: gather_exch
    type(tiles_t)                  :: tiles
    type(grid_field_t)             :: buf
    integer(kind=4)                :: master_id
    integer(kind=4)                :: rec_num = 1

contains
    procedure, public :: write => master_paneled_write
end type master_paneled_outputer_t

contains

subroutine master_paneled_write(this, f, domain, file_name, rec_num)

    class(master_paneled_outputer_t), intent(inout) :: this
    type(grid_field_t),               intent(inout) :: f
    character(*),                     intent(in)    :: file_name
    type(domain_t),                   intent(in)    :: domain
    integer(kind=4),  optional,       intent(in)    :: rec_num

    real(kind=4), allocatable :: buffer(:,:,:)
    integer(kind=4) irec, reclen
    integer(kind=4) :: t, i, j, k

    if (domain%parcomm%myid == this%master_id) then
        do t = domain%partition%ts, domain%partition%te
            do k = this%tiles%tile(t)%ks, this%tiles%tile(t)%ke
                do j = this%tiles%tile(t)%js, this%tiles%tile(t)%je
                    do i = this%tiles%tile(t)%is, this%tiles%tile(t)%ie
                        this%buf%tile(t)%p(i,j,k) = f%tile(t)%p(i,j,k)
                    end do
                end do
            end do
        end do
        call this%gather_exch%do(this%buf, domain%parcomm)
    else
        call this%gather_exch%do(f, domain%parcomm)
    end if

    if (domain%parcomm%myid == this%master_id) then
            irec = 1
            if (present(rec_num)) irec = rec_num

            reclen  = domain%partition%num_panels*this%tiles%Nx*this%tiles%Ny
            allocate(buffer(this%tiles%Nx,this%tiles%Ny,1:domain%partition%num_panels))

            open(newunit=this%out_stream, file = trim(file_name), &
                 access="direct", recl = reclen)

            do k = 1, this%tiles%Nz
                do t = 1, this%tiles%Nt
                    do j = this%tiles%tile(t)%js, this%tiles%tile(t)%je
                        do i = this%tiles%tile(t)%is, this%tiles%tile(t)%ie
                            buffer(i,j,domain%partition%panel_map(t)) = &
                            real(this%buf%tile(t)%p(i,j,k),4)
                        end do
                    end do
                end do

                write(this%out_stream, rec = (irec-1)*(this%tiles%Nz-1+1)+(k-1+1)) buffer
            end do
            deallocate(buffer)
            close(this%out_stream)
    end if

end subroutine master_paneled_write

end module master_paneled_outputer_mod
