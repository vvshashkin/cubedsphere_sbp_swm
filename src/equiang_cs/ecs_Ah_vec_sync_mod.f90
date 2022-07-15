!Syncronize Ah-grid vector components across edges of cubed sphere
!1)get vector components from heighbouring face
!2)transform them to the coordinate system of current face
!3) make half sums with components from current face at the same point
!Transformation weights are precomputed in ecs_halo_Ah_vec_sync_factory_mod

module ecs_halo_Ah_vec_sync_mod

use halo_mod,          only : halo_vec_t
use exchange_halo_mod, only : exchange_t
use parcomm_mod,       only : parcomm_global

implicit none

type, extends(halo_vec_t) :: ecs_halo_Ah_vec_sync_t

    integer(kind=4)                        :: ts, te
    character(len=:),          allocatable :: components_type
    class(exchange_t),         allocatable :: exch_edges
    type(tile_sync_t),         allocatable :: tile(:)

    contains

    procedure :: get_halo_vector => sync_ecs_Ah_vector

end type ecs_halo_Ah_vec_sync_t

type tile_sync_t
    real(kind=8), dimension(:), allocatable :: qb, qr, qt, ql
end type tile_sync_t

contains

subroutine sync_ecs_Ah_vector(this,u,v,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_Ah_vec_sync_t),    intent(inout) :: this
    class(grid_field_t),              intent(inout) :: u,v
    type(domain_t),                   intent(in)    :: domain
    integer(kind=4),                  intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_edges%do_vec(u,v, domain%parcomm)

    do t = this%ts, this%te
        if(this%components_type == "contravariant") then
            call syncronize_grad_on_edges(u%tile(t), v%tile(t), this%tile(t), &
                                          domain%mesh_xy%tile(t))
        else if(this%components_type == "covariant") then
            call syncronize_grad_on_edges(v%tile(t), u%tile(t), this%tile(t), &
                                          domain%mesh_xy%tile(t))
        else
            call parcomm_global%abort("ecs_Ah vector sync"// &
                                       " -unknown type of vector components"// &
                                       this%components_type)
        end if
    end do
end subroutine sync_ecs_Ah_vector

subroutine syncronize_grad_on_edges(u,v,q,mesh)

    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t),          intent(inout) :: u, v
    type(tile_sync_t),           intent(in)    :: q
    type(tile_mesh_t),           intent(in)    :: mesh

    integer(kind=4) :: ks, ke, js, je, is, ie, i, j, k

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    if(is==1 .and. js==1) then
        do k=ks,ke
            u%p(1,1,k) = (u%p(1,1,k)+u%p(0,1,k)+u%p(1,0,k)+q%qb(1)*v%p(1,0,k))/3.0_8
            v%p(1,1,k) = (v%p(1,1,k)+v%p(0,1,k)+v%p(1,0,k)+q%ql(1)*u%p(0,1,k))/3.0_8
        end do
    end if
    if(ie==mesh%nx .and. js==1) then
        do k=ks,ke
            u%p(ie,1,k) = (u%p(ie,1,k)+u%p(ie+1,1,k)+u%p(ie,0,k)+q%qb(ie)*v%p(ie,0,k))/3.0_8
            v%p(ie,1,k) = (v%p(ie,1,k)+v%p(ie+1,1,k)+v%p(ie,0,k)+q%qr(1)*u%p(ie+1,1,k))/3.0_8
        end do
    end if
    if(is==1 .and. je==mesh%ny) then
        do k=ks,ke
            u%p(1,je,k) = (u%p(1,je,k)+u%p(0,je,k)+u%p(1,je+1,k)+q%qt(1)*v%p(1,je+1,k))/3.0_8
            v%p(1,je,k) = (v%p(1,je,k)+v%p(0,je,k)+v%p(1,je+1,k)+q%ql(je)*u%p(0,je,k))/3.0_8
        end do
    end if
    if(ie==mesh%nx .and. je==mesh%ny) then
        do k=ks,ke
            u%p(ie,je,k) = (u%p(ie,je,k)+u%p(ie+1,je,k)+u%p(ie,je+1,k)+q%qt(ie)*v%p(ie,je+1,k))/3.0_8
            v%p(ie,je,k) = (v%p(ie,je,k)+v%p(ie+1,je,k)+v%p(ie,je+1,k)+q%qr(je)*u%p(ie+1,je,k))/3.0_8
        end do
    end if
    if(is == 1) then
        do k=ks,ke; do j=max(js,2),min(je,mesh%ny-1)
            u%p(1,j,k) = 0.5_8*(u%p(0,j,k)+u%p(1,j,k))
            v%p(1,j,k) = 0.5_8*(v%p(1,j,k)+v%p(0,j,k)+q%ql(j)*u%p(0,j,k))
        end do; end do
    end if
    if(ie == mesh%nx) then
        do k=ks,ke; do j=max(js,2),min(je,mesh%ny-1)
            u%p(ie,j,k) = 0.5_8*(u%p(ie,j,k)+u%p(ie+1,j,k))
            v%p(ie,j,k) = 0.5_8*(v%p(ie,j,k)+v%p(ie+1,j,k)+q%qr(j)*u%p(ie+1,j,k))
        end do; end do
    end if
    if(js==1) then
        do k=ks,ke; do i=max(is,2),min(ie,mesh%nx-1)
            v%p(i,1,k) = 0.5_8*(v%p(i,0,k)+v%p(i,1,k))
            u%p(i,1,k) = 0.5_8*(u%p(i,0,k)+u%p(i,1,k)+q%qb(i)*v%p(i,0,k))
        end do; end do
    end if
    if(je == mesh%ny) then
        do k=ks,ke; do i=max(is,2),min(ie,mesh%nx-1)
            v%p(i,je,k) = 0.5_8*(v%p(i,je,k)+v%p(i,je+1,k))
            u%p(i,je,k) = 0.5_8*(u%p(i,je,k)+u%p(i,je+1,k)+q%qt(i)*v%p(i,je+1,k))
        end do; end do
    end if
end subroutine syncronize_grad_on_edges

end module  ecs_halo_Ah_vec_sync_mod
