module vector_advection_Ah_mod

use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use halo_mod,                      only : halo_vec_t
use exchange_abstract_mod,         only : exchange_t
use abstract_vector_advection_mod, only : vector_advection_operator_t
use abstract_v_nabla_mod,          only : v_nabla_operator_t

implicit none

type, public, extends(vector_advection_operator_t) :: vector_advection_Ah_t
    class(halo_vec_t),         allocatable :: sync_edges_cov
    class(halo_vec_t),         allocatable :: sync_edges_contra
    class(exchange_t),         allocatable :: exch_uv_interior
    class(v_nabla_operator_t), allocatable :: v_nabla_op
contains
    procedure :: calc_vec_advection
    procedure :: calc_vec_advection_contra
end type vector_advection_Ah_t

contains

subroutine calc_vec_advection(this, u_tend, v_tend, u, v, ut, vt, domain)
    class(vector_advection_Ah_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u,  v!covariant components
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    call this%exch_uv_interior%do_vec(u,v,domain%parcomm)

    call this%v_nabla_op%calc_v_nabla(u_tend, u, ut, vt, domain%mesh_u)
    call this%v_nabla_op%calc_v_nabla(v_tend, v, ut, vt, domain%mesh_v)

    do t = domain%mesh_xy%ts, domain%mesh_xy%te
        call add_metric_terms_tile(u_tend%tile(t), v_tend%tile(t),               &
                                   u%tile(t), v%tile(t), ut%tile(t), vt%tile(t), &
                                   domain%mesh_xy%tile(t), domain%mesh_xy%scale)
    end do

    call this%sync_edges_cov%get_halo_vector(u_tend,v_tend,domain,0)
end subroutine calc_vec_advection

subroutine add_metric_terms_tile(u_tend, v_tend, u, v, ut, vt, mesh, scale)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(in)    :: u, v, ut, vt
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: scale

    type(tile_field_t),     intent(inout) :: u_tend, v_tend

    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke


    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k = ks, ke

        do j = js, je
            do i = is, ie
                u_tend%p(i,j,k) =  u_tend%p(i,j,k)+&
                                  (u%p(i,j,k)*ut%p(i,j,k)*mesh%G(1,1,1,i,j,k)+ &
                                   v%p(i,j,k)*ut%p(i,j,k)*mesh%G(1,1,2,i,j,k)+ &
                                   u%p(i,j,k)*vt%p(i,j,k)*mesh%G(1,2,1,i,j,k)+ &
                                   v%p(i,j,k)*vt%p(i,j,k)*mesh%G(1,2,2,i,j,k)) / scale
            end do
        end do

        do j = js, je
            do i = is, ie
                v_tend%p(i,j,k) =  v_tend%p(i,j,k)+&
                                  (u%p(i,j,k)*ut%p(i,j,k)*mesh%G(2,1,1,i,j,k)+ &
                                   v%p(i,j,k)*ut%p(i,j,k)*mesh%G(2,1,2,i,j,k)+ &
                                   u%p(i,j,k)*vt%p(i,j,k)*mesh%G(2,2,1,i,j,k)+ &
                                   v%p(i,j,k)*vt%p(i,j,k)*mesh%G(2,2,2,i,j,k)) / scale
            end do
        end do

    end do

end subroutine add_metric_terms_tile

subroutine calc_vec_advection_contra(this, u_tend, v_tend, ut, vt, domain)
    class(vector_advection_Ah_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    call this%exch_uv_interior%do_vec(ut,vt,domain%parcomm)

    call this%v_nabla_op%calc_v_nabla(u_tend, ut, ut, vt, domain%mesh_u)
    call this%v_nabla_op%calc_v_nabla(v_tend, vt, ut, vt, domain%mesh_v)

    do t = domain%mesh_xy%ts, domain%mesh_xy%te
        call add_metric_terms_contra_tile(u_tend%tile(t), v_tend%tile(t),               &
                                          ut%tile(t), vt%tile(t),                       &
                                          domain%mesh_xy%tile(t), domain%mesh_xy%scale)
    end do

    call this%sync_edges_contra%get_halo_vector(u_tend,v_tend,domain,0)
end subroutine calc_vec_advection_contra

subroutine add_metric_terms_contra_tile(u_tend, v_tend, ut, vt, mesh, scale)

    use mesh_mod, only : tile_mesh_t

    type(tile_field_t),     intent(in)    :: ut, vt
    type(tile_mesh_t),      intent(in)    :: mesh
    real(kind=8),           intent(in)    :: scale

    type(tile_field_t),     intent(inout) :: u_tend, v_tend

    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8)    :: hx

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k = ks, ke

        do j = js, je
            do i = is, ie
                u_tend%p(i,j,k) = u_tend%p(i,j,k)-&
                                  (ut%p(i,j,k)*ut%p(i,j,k)*mesh%G(1,1,1,i,j,k)+ &
                                   2.0_8*ut%p(i,j,k)*vt%p(i,j,k)*mesh%G(1,2,1,i,j,k)+ &
                                   vt%p(i,j,k)*vt%p(i,j,k)*mesh%G(2,2,1,i,j,k)) / scale
            end do
        end do

        do j = js, je
            do i = is, ie
                v_tend%p(i,j,k) = v_tend%p(i,j,k)-&
                                  (ut%p(i,j,k)*ut%p(i,j,k)*mesh%G(1,1,2,i,j,k)+ &
                                   2.0_8*vt%p(i,j,k)*ut%p(i,j,k)*mesh%G(1,2,2,i,j,k)+ &
                                   vt%p(i,j,k)*vt%p(i,j,k)*mesh%G(2,2,2,i,j,k)) / scale
            end do
        end do

    end do

end subroutine add_metric_terms_contra_tile

end module vector_advection_Ah_mod
