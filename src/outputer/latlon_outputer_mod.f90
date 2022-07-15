module latlon_outputer_mod

use outputer_abstract_mod, only : outputer_t, outputer_vector_t
use grid_field_mod,        only : grid_field_t
use exchange_abstract_mod, only : exchange_t
use domain_mod,            only : domain_t
use tiles_mod,             only : tiles_t
use abstract_regrid_mod,   only : regrid_t, regrid_vec_t

implicit none

type, public, extends(outputer_t) :: latlon_outputer_t
    class(exchange_t), allocatable :: gather_exch
    type(grid_field_t)             :: exchange_buf
    integer(kind=4)                :: master_id
    type(tiles_t)                  :: tiles
    class(regrid_t),   allocatable :: regrid
    type(domain_t)                 :: regrid_domain
    type(grid_field_t)             :: regrid_work
    integer(kind=4)                :: Nlon, Nlat
    real(kind=8),      allocatable :: buffer(:,:,:)

contains
    procedure, public :: write => latlon_write
end type latlon_outputer_t

type, public, extends(outputer_vector_t) :: latlon_vec_outputer_t
    class(exchange_t), allocatable     :: gather_exch_u, gather_exch_v
    type(grid_field_t)                 :: exchange_buf_u, exchange_buf_v
    integer(kind=4)                    :: master_id
    type(tiles_t)                      :: tiles_u, tiles_v
    class(regrid_vec_t),   allocatable :: regrid
    type(domain_t)                     :: regrid_domain
    type(grid_field_t)                 :: regrid_work_u, regrid_work_v
    integer(kind=4)                    :: Nlon, Nlat
    real(kind=8),      allocatable     :: ulatlon(:,:,:), vlatlon(:,:,:)

contains
    procedure, public :: write => latlon_vec_write
end type latlon_vec_outputer_t

contains

subroutine latlon_write(this, f, domain, file_name, rec_num)

    class(latlon_outputer_t),         intent(inout) :: this
    type(grid_field_t),               intent(inout) :: f
    character(*),                     intent(in)    :: file_name
    type(domain_t),                   intent(in)    :: domain
    integer(kind=4),  optional,       intent(in)    :: rec_num

    integer(kind=4) :: i, j, k, t, ks, ke, panel_ind
    integer(kind=4) :: ts, te, js, je, is, ie

    if (domain%parcomm%myid == this%master_id) then
        do t = domain%partition%ts, domain%partition%te
            call this%tiles%tile(t)%getind(is, ie, js, je, ks, ke)
            do k = ks, ke
                do j = js, je
                    do i = is, ie
                        this%exchange_buf%tile(t)%p(i,j,k) = f%tile(t)%p(i,j,k)
                    end do
                end do
            end do
        end do
        call this%gather_exch%do(this%exchange_buf, domain%parcomm)
    else
        call this%gather_exch%do(f, domain%parcomm)
    end if

    if(domain%parcomm%myid == this%master_id) then

        open(newunit=this%out_stream, file = trim(file_name), &
             access="direct", recl = this%Nlat*this%Nlon)

        ts = this%exchange_buf%ts
        te = this%exchange_buf%te
        ks = this%exchange_buf%tile(ts)%ks
        ke = this%exchange_buf%tile(ts)%ke
        !map exchanged grid function to work buffer before regrid
        do k = ks,ke

            do t=ts,te
                is = this%exchange_buf%tile(t)%is
                ie = this%exchange_buf%tile(t)%ie
                js = this%exchange_buf%tile(t)%js
                je = this%exchange_buf%tile(t)%je
                panel_ind = domain%partition%panel_map(t)
                do j=js, je; do i=is, ie
                    this%regrid_work%tile(panel_ind)%p(i,j,1) = this%exchange_buf%tile(t)%p(i,j,k)
                end do; end do
            end do

            call this%regrid%do_regrid(this%buffer,this%regrid_work,this%regrid_domain)
            write(this%out_stream,rec=(rec_num-1)*(ke-ks+1)+(k-ks+1)) real(this%buffer,4)
            !print *, k, size(this%buffer), this%Nlat*this%Nlon
        end do
        close(this%out_stream)
    end if

end subroutine latlon_write

subroutine latlon_vec_write(this, u, v, domain, file_name_u, file_name_v, rec_num)

    class(latlon_vec_outputer_t), intent(inout) :: this
    type(grid_field_t),           intent(inout) :: u, v
    character(*),                 intent(in)    :: file_name_u, file_name_v
    type(domain_t),               intent(in)    :: domain
    integer(kind=4),              intent(in)    :: rec_num

    integer(kind=4) :: i, j, k, t, ks, ke, panel_ind
    integer(kind=4) :: ts, te, js, je, is, ie

     if (domain%parcomm%myid == this%master_id) then
         do t = domain%partition%ts, domain%partition%te
             call this%tiles_u%tile(t)%getind(is, ie, js, je, ks, ke)
             do k = ks, ke
                 do j = js, je
                     do i = is, ie
                         this%exchange_buf_u%tile(t)%p(i,j,k) = u%tile(t)%p(i,j,k)
                     end do
                 end do
             end do
         end do
         call this%gather_exch_u%do(this%exchange_buf_u, domain%parcomm)
     else
         call this%gather_exch_u%do(u, domain%parcomm)
     end if

     if (domain%parcomm%myid == this%master_id) then
         do t = domain%partition%ts, domain%partition%te
             call this%tiles_v%tile(t)%getind(is, ie, js, je, ks, ke)
             do k = ks, ke
                 do j = js, je
                     do i = is, ie
                         this%exchange_buf_v%tile(t)%p(i,j,k) = v%tile(t)%p(i,j,k)
                     end do
                 end do
             end do
         end do
         call this%gather_exch_v%do(this%exchange_buf_v, domain%parcomm)
     else
         call this%gather_exch_v%do(v, domain%parcomm)
     end if

    if(domain%parcomm%myid == this%master_id) then

        open(newunit=this%out_stream_u, file = trim(file_name_u), &
             access="direct", recl = this%Nlat*this%Nlon)
        open(newunit=this%out_stream_v, file = trim(file_name_v), &
             access="direct", recl = this%Nlat*this%Nlon)
        ts = this%exchange_buf_u%ts
        te = this%exchange_buf_v%te
        ks = this%exchange_buf_u%tile(ts)%ks
        ke = this%exchange_buf_v%tile(ts)%ke
        !map exchanged grid function to work buffer before regrid
        do k = ks,ke

            do t=ts,te
                is = this%exchange_buf_u%tile(t)%is
                ie = this%exchange_buf_u%tile(t)%ie
                js = this%exchange_buf_u%tile(t)%js
                je = this%exchange_buf_u%tile(t)%je
                panel_ind = domain%partition%panel_map(t)
                do j=js, je; do i=is, ie
                    this%regrid_work_u%tile(panel_ind)%p(i,j,1) = this%exchange_buf_u%tile(t)%p(i,j,k)
                end do; end do
            end do
            do t=ts,te
                is = this%exchange_buf_v%tile(t)%is
                ie = this%exchange_buf_v%tile(t)%ie
                js = this%exchange_buf_v%tile(t)%js
                je = this%exchange_buf_v%tile(t)%je
                panel_ind = domain%partition%panel_map(t)
                do j=js, je; do i=is, ie
                    this%regrid_work_v%tile(panel_ind)%p(i,j,1) = this%exchange_buf_v%tile(t)%p(i,j,k)
                end do; end do
            end do

            call this%regrid%do_regrid(this%ulatlon,this%vlatlon, &
                                       this%regrid_work_u,this%regrid_work_v, &
                                       this%regrid_domain)
            write(this%out_stream_u,rec=(rec_num-1)*(ke-ks+1)+(k-ks)+1) real(this%ulatlon,4)
            write(this%out_stream_v,rec=(rec_num-1)*(ke-ks+1)+(k-ks)+1) real(this%vlatlon,4)

        end do
        close(this%out_stream_u)
        close(this%out_stream_v)
    end if

end subroutine latlon_vec_write

end module latlon_outputer_mod
