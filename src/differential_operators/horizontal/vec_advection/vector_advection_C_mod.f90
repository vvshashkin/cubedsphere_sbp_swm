module vector_advection_C_mod

use abstract_v_nabla_mod,          only : v_nabla_operator_t
use abstract_vector_advection_mod, only : vector_advection_operator_t
use grid_field_mod,                only : grid_field_t, tile_field_t
use domain_mod,                    only : domain_t
use halo_mod,                      only : halo_vec_t
use abstract_interpolators2d_mod,  only : interpolator2d_vec2vec_t
use parcomm_mod,                   only : parcomm_global

implicit none

type, public, extends(vector_advection_operator_t) :: vector_advection_C_t
    class(v_nabla_operator_t), allocatable       :: v_nabla_op
    type(grid_field_t)                           :: u_at_v, v_at_u, uh, vh
    class(interpolator2d_vec2vec_t), allocatable :: interp_v2h_op
    class(interpolator2d_vec2vec_t), allocatable :: interp_h2v_op
    class(halo_vec_t), allocatable               :: halo_uv, tendency_edge_sync
contains
    procedure :: calc_vec_advection
    procedure :: calc_vec_advection_contra
end type vector_advection_C_t

contains

subroutine calc_vec_advection(this, u_tend, v_tend, u, v, ut, vt, domain)
    class(vector_advection_C_t),  intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u,  v!covariant components
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    call parcomm_global%abort("vector_advection_c_t, calc_vec_advection not implemented")

end subroutine calc_vec_advection

subroutine calc_vec_advection_contra(this, u_tend, v_tend, ut, vt, domain)
    class(vector_advection_C_t),  intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: ut, vt!contravariant components
    type(grid_field_t),           intent(inout) :: u_tend, v_tend!advective tendencies

    integer(kind=4) :: t

    call this%interp_v2h_op%interp2d_vec2vec(this%uh, this%vh, ut, vt, domain)
    call this%interp_h2v_op%interp2d_vec2vec(this%v_at_u, this%u_at_v, this%vh, this%uh, domain)

    call this%halo_uv%get_halo_vector(ut, vt, domain, 3)

    call this%v_nabla_op%calc_v_nabla(u_tend, ut, ut, this%v_at_u, domain%mesh_u)
    call this%v_nabla_op%calc_v_nabla(v_tend, vt, this%u_at_v, vt, domain%mesh_v)

    do t = domain%mesh_o%ts, domain%mesh_o%te
        call add_metric_terms_1comp_contra_tile(u_tend%tile(t), ut%tile(t), ut%tile(t), this%v_at_u%tile(t), &
                                                domain%mesh_u%scale, domain%mesh_u%tile(t),"u")
        call add_metric_terms_1comp_contra_tile(v_tend%tile(t), vt%tile(t), this%u_at_v%tile(t), vt%tile(t), &
                                                domain%mesh_v%scale, domain%mesh_v%tile(t),"v")
    end do

    call this%tendency_edge_sync%get_halo_vector(u_tend, v_tend, domain, 1)
end subroutine calc_vec_advection_contra

subroutine add_metric_terms_1comp_contra_tile(uvt_tend, uvt, ut, vt, scale, mesh, component)

    use mesh_mod, only : tile_mesh_t
    use tile_mod, only : tile_t

    type(tile_field_t),     intent(in)    :: uvt, ut, vt
    real(kind=8),           intent(in)    :: scale
    type(tile_mesh_t),      intent(in)    :: mesh
    character(len=*),       intent(in)    :: component

    type(tile_field_t),     intent(inout) :: uvt_tend

    integer(kind=4) :: i, j, k
    integer(kind=4) :: component_num
    integer(kind=4) :: is, ie, js, je, ks, ke
    real(kind=8)    :: hx
    real(kind=8)    :: dx, dy, zl, zr

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    hx = mesh%hx

    if(component=="u") then
        component_num = 1
    else if(component=="v") then
        component_num = 2
    else
        call parcomm_global%abort("vector_advection_C_mod, unknown vector component type:"//component)
    end if

    do k = ks, ke
        do j = js, je
            do i = is, ie
                uvt_tend%p(i,j,k) = uvt_tend%p(i,j,k) -&
                                    (ut%p(i,j,k)*ut%p(i,j,k)*mesh%G(1,1,component_num,i,j,k)+ &
                                     2.0_8*ut%p(i,j,k)*vt%p(i,j,k)*mesh%G(1,2,component_num,i,j,k)+ &
                                     vt%p(i,j,k)*vt%p(i,j,k)*mesh%G(2,2,component_num,i,j,k)) / scale
            end do
        end do
    end do

end subroutine add_metric_terms_1comp_contra_tile

end module vector_advection_C_mod
