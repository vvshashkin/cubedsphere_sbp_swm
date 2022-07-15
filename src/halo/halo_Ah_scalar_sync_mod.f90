module halo_Ah_scalar_sync_mod

use halo_mod,          only : halo_t
use exchange_halo_mod, only : exchange_t

implicit none

type, extends(halo_t) :: halo_Ah_scalar_sync_t

    class(exchange_t), allocatable  :: exch_halo

    contains

    procedure :: get_halo_scalar => make_Ah_scalar_sync

end type

contains

subroutine make_Ah_scalar_sync(this,f,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(halo_Ah_scalar_sync_t),  intent(inout) :: this
    class(grid_field_t),           intent(inout) :: f
    type(domain_t),                intent(in)    :: domain
    integer(kind=4),               intent(in)    :: halo_width

    integer(kind=4)  :: t

    call this%exch_halo%do(f, domain%parcomm)

    do t = f%ts, f%te
        call sync_Ah_scalar_tile(f%tile(t),domain%mesh_xy%tile(t))
    end do
end subroutine make_Ah_scalar_sync

subroutine sync_Ah_scalar_tile(f, mesh)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: f
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    if(is == 1) then
        do k=ks,ke
            if(js == 1) then
                f%p(1,1,k) = (f%p(1,1,k)+f%p(0,1,k)+f%p(1,0,k))/3.0_8
            end if
            do j = max(2,js),min(mesh%ny-1,je)
                f%p(1,j,k) = 0.5_8*(f%p(0,j,k)+f%p(1,j,k))
            end do
            if(je == mesh%ny) then
                f%p(1,je,k) = (f%p(1,je,k)+f%p(0,je,k)+f%p(1,je+1,k))/3.0_8
            end if
        end do
    end if
    if(ie == mesh%nx) then
        do k=ks,ke
            if(js == 1) then
                f%p(ie,1,k) = (f%p(ie,1,k)+f%p(ie+1,1,k)+f%p(ie,0,k))/3.0_8
            end if
            do j = max(2,js),min(mesh%ny-1,je)
                f%p(ie,j,k) = 0.5_8*(f%p(ie,j,k)+f%p(ie+1,j,k))
            end do
            if(je == mesh%ny) then
                f%p(ie,je,k) = (f%p(ie,je,k)+f%p(ie+1,je,k)+f%p(ie,je+1,k))/3.0_8
            end if
        end do
    end if
    if(js == 1) then
        do k=ks,ke
            do i = max(2,is),min(mesh%nx-1,ie)
                f%p(i,1,k) = 0.5_8*(f%p(i,1,k)+f%p(i,0,k))
            end do
        end do
    end if
    if(je == mesh%ny) then
        do k=ks,ke
            do i = max(2,is),min(mesh%nx-1,ie)
                f%p(i,je,k) = 0.5_8*(f%p(i,je,k)+f%p(i,je+1,k))
            end do
        end do
    end if

end subroutine sync_Ah_scalar_tile

end module halo_Ah_scalar_sync_mod
